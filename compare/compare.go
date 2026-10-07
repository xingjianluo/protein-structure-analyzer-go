// Package compare computes cross-structure interface differences for a batch
// of related structures. Results are deterministic: independent of submission
// order and stable across repeated runs.
package compare

import (
	"fmt"
	"hash/fnv"
	"math"
	"sort"
	"strconv"
	"strings"
	"sync"

	"github.com/xingjianluo/protein-structure-analyzer-go/contact"
	"github.com/xingjianluo/protein-structure-analyzer-go/structure"
)

// DistanceChangeTolerance is the minimum max-min spread (Angstrom) of one
// contact's distance across structures for it to be reported as changed.
const DistanceChangeTolerance = 0.5

// InterfacePresence records in which structures one chain-pair interface
// appears, and how strong it is in each.
type InterfacePresence struct {
	ChainA string
	ChainB string
	// PresentIn lists source IDs (sorted) where the interface exists.
	PresentIn []string
	// ContactCounts maps source ID to the contact count there.
	ContactCounts map[string]int
	// MinDists maps source ID to the minimum contact distance there.
	MinDists map[string]float64
}

// ContactChange describes how one residue-level contact differs across the
// batch: either it is absent from some structures, or its distance varies by
// more than DistanceChangeTolerance.
type ContactChange struct {
	ChainA string
	ChainB string
	SeqA   int
	ICodeA byte
	SeqB   int
	ICodeB byte
	// NameA/NameB are the residue names in the lexicographically smallest
	// source where the contact is present.
	NameA string
	NameB string
	// PresentIn lists source IDs (sorted) where this contact exists.
	PresentIn []string
	// MinDistBySource maps source ID to the contact distance there.
	MinDistBySource map[string]float64
	// MaxSpread is max-min of the distance across present structures.
	MaxSpread float64
}

// Result is the deterministic outcome of a batch comparison.
type Result struct {
	// Fingerprint identifies options+inputs; an identical fingerprint means
	// an identical conclusion.
	Fingerprint string
	Options     contact.Options
	Sources     []string
	// Common interfaces are present in every compared structure.
	Common []InterfacePresence
	// Partial interfaces appear in some but not all structures.
	Partial []InterfacePresence
	// Changes lists residue contacts whose presence set or distance varies.
	Changes []ContactChange
}

// BatchError reports which inputs of a comparison batch failed. It is
// distinguishable by type; each failure is a *contact.InputError.
type BatchError struct {
	Failures []error
}

func (e *BatchError) Error() string {
	parts := make([]string, len(e.Failures))
	for i, f := range e.Failures {
		parts[i] = f.Error()
	}
	return "batch comparison failed: " + strings.Join(parts, "; ")
}

// Comparer computes batch comparisons and caches the last result keyed by an
// options+inputs fingerprint, so a repeated identical submission is served
// from cache while any change of options or data invalidates it.
type Comparer struct {
	analyzer *contact.Analyzer

	mu         sync.RWMutex
	lastKey    string
	lastResult *Result
}

// New creates a Comparer.
func New(a *contact.Analyzer) *Comparer {
	if a == nil {
		a = contact.NewAnalyzer(contact.DefaultCutoff)
	}
	return &Comparer{analyzer: a}
}

// fingerprintStructure hashes the effective-model content that feeds contact
// calculation: chain IDs, residue keys, names, hetero flags and atom
// coordinates. O(atoms), deterministic for identical content.
func fingerprintStructure(h interface{ Write([]byte) (int, error) }, s *structure.Structure) {
	if s.Effective == nil {
		fmt.Fprintf(h, "src=%s no-effective-model\n", s.Source)
		return
	}
	fmt.Fprintf(h, "src=%s model=%d\n", s.Source, s.Effective.Serial)
	for _, c := range s.Effective.Chains {
		fmt.Fprintf(h, "chain %s\n", c.ID)
		for _, r := range c.Residues {
			fmt.Fprintf(h, "res %d %c %s het=%v std=%v\n", r.Seq, r.ICode, r.Name, r.Hetero, r.IsStandardAA)
			for _, a := range r.Atoms {
				fmt.Fprintf(h, "atom %s %x %x %x\n", a.Name,
					math.Float64bits(a.Coord.X), math.Float64bits(a.Coord.Y), math.Float64bits(a.Coord.Z))
			}
		}
	}
}

// batchFingerprint builds the cache/identity key from options and the sorted
// per-structure content hashes. Submission order does not affect it.
func batchFingerprint(opts contact.Options, structs []*structure.Structure) string {
	type unit struct {
		source string
		hash   uint64
	}
	units := make([]unit, 0, len(structs))
	for _, s := range structs {
		h := fnv.New64a()
		fingerprintStructure(h, s)
		units = append(units, unit{source: s.Source, hash: h.Sum64()})
	}
	sort.Slice(units, func(i, j int) bool { return units[i].source < units[j].source })
	var b strings.Builder
	b.WriteString("opts{")
	b.WriteString(opts.Fingerprint())
	b.WriteString("}")
	for _, u := range units {
		b.WriteString("|")
		b.WriteString(u.source)
		b.WriteString("#")
		b.WriteString(strconv.FormatUint(u.hash, 16))
	}
	return b.String()
}

// contactKey identifies one residue-level contact by endpoint positions.
type contactKey struct {
	chainA string
	seqA   int
	icodeA byte
	chainB string
	seqB   int
	icodeB byte
}

// Compare analyzes the batch under opts and returns the cross-structure
// comparison. The result is independent of the submission order of structs.
// If any input cannot be analyzed, a *BatchError is returned and no result is
// produced (the cached previous result is untouched).
func (c *Comparer) Compare(structs []*structure.Structure, opts contact.Options) (*Result, error) {
	if len(structs) == 0 {
		return nil, &BatchError{Failures: []error{&contact.AnalysisError{Reason: "empty batch"}}}
	}
	key := batchFingerprint(opts, structs)

	c.mu.RLock()
	if c.lastResult != nil && c.lastKey == key {
		res := c.lastResult
		c.mu.RUnlock()
		return res, nil
	}
	c.mu.RUnlock()

	// Analyze every input; any failure aborts the whole batch.
	reports := map[string]*contact.InterfaceReport{}
	var failures []error
	for i, s := range structs {
		if s == nil {
			failures = append(failures, &contact.InputError{Index: i, Source: "", Err: contact.ErrEmptyStructure})
			continue
		}
		rep, err := c.analyzer.AnalyzeWithOptions(s, opts)
		if err != nil {
			failures = append(failures, &contact.InputError{Index: i, Source: s.Source, Err: err})
			continue
		}
		reports[s.Source] = rep
	}
	if len(failures) > 0 {
		return nil, &BatchError{Failures: failures}
	}

	res := FromReports(opts, reports)

	c.mu.Lock()
	c.lastKey = key
	c.lastResult = res
	c.mu.Unlock()
	return res, nil
}

// fingerprintReports builds the public result identity from the options and
// the canonical content of every per-structure report. Any change of options
// or of any input's analysis output yields a different fingerprint.
func fingerprintReports(opts contact.Options, reports map[string]*contact.InterfaceReport, sources []string) string {
	h := fnv.New64a()
	fmt.Fprintf(h, "opts{%s}\n", opts.Fingerprint())
	for _, src := range sources {
		rep := reports[src]
		fmt.Fprintf(h, "src=%s model=%d cutoff=%g\n", src, rep.Model, rep.Cutoff)
		for _, it := range rep.Interfaces {
			fmt.Fprintf(h, "iface %s %s %d %x\n", it.ChainA, it.ChainB, it.ContactCount, math.Float64bits(it.MinDist))
		}
		for _, ct := range rep.Contacts {
			fmt.Fprintf(h, "contact %s %d %c %s %s %d %c %s %x\n",
				ct.ChainA, ct.SeqA, ct.ICodeA, ct.NameA,
				ct.ChainB, ct.SeqB, ct.ICodeB, ct.NameB, math.Float64bits(ct.MinDist))
		}
	}
	return "cmp-" + strconv.FormatUint(h.Sum64(), 16)
}

// FromReports aggregates per-structure interface reports into the
// deterministic batch comparison. It is exported so that callers which
// already analyzed their inputs (e.g. the store layer) reuse the exact same
// aggregation and ordering rules as Compare.
func FromReports(opts contact.Options, reports map[string]*contact.InterfaceReport) *Result {
	sources := make([]string, 0, len(reports))
	for src := range reports {
		sources = append(sources, src)
	}
	sort.Strings(sources)
	return buildResult(fingerprintReports(opts, reports, sources), opts, reports, sources)
}

// buildResult aggregates per-structure reports into the deterministic
// comparison. sources must be sorted.
func buildResult(key string, opts contact.Options, reports map[string]*contact.InterfaceReport, sources []string) *Result {
	n := len(sources)

	// --- interface-level aggregation ---
	type ifaceAgg struct {
		present []string
		counts  map[string]int
		mins    map[string]float64
	}
	ifaces := map[string]*ifaceAgg{}
	// --- contact-level aggregation ---
	type contactAgg struct {
		nameA, nameB string
		present      []string
		dists        map[string]float64
	}
	contacts := map[contactKey]*contactAgg{}

	for _, src := range sources {
		rep := reports[src]
		for _, it := range rep.Interfaces {
			k := it.ChainA + "\x00" + it.ChainB
			g := ifaces[k]
			if g == nil {
				g = &ifaceAgg{counts: map[string]int{}, mins: map[string]float64{}}
				ifaces[k] = g
			}
			g.present = append(g.present, src)
			g.counts[src] = it.ContactCount
			g.mins[src] = it.MinDist
		}
		for _, ct := range rep.Contacts {
			k := contactKey{ct.ChainA, ct.SeqA, ct.ICodeA, ct.ChainB, ct.SeqB, ct.ICodeB}
			g := contacts[k]
			if g == nil {
				g = &contactAgg{nameA: ct.NameA, nameB: ct.NameB, dists: map[string]float64{}}
				contacts[k] = g
			}
			g.present = append(g.present, src)
			g.dists[src] = ct.MinDist
		}
	}

	res := &Result{
		Fingerprint: key,
		Options:     opts,
		Sources:     sources,
	}
	for k, g := range ifaces {
		i := strings.IndexByte(k, 0)
		p := InterfacePresence{
			ChainA: k[:i], ChainB: k[i+1:],
			PresentIn: g.present, ContactCounts: g.counts, MinDists: g.mins,
		}
		if len(g.present) == n {
			res.Common = append(res.Common, p)
		} else {
			res.Partial = append(res.Partial, p)
		}
	}
	sort.Slice(res.Common, func(i, j int) bool {
		a, b := res.Common[i], res.Common[j]
		if a.ChainA != b.ChainA {
			return a.ChainA < b.ChainA
		}
		return a.ChainB < b.ChainB
	})
	sort.Slice(res.Partial, func(i, j int) bool {
		a, b := res.Partial[i], res.Partial[j]
		if a.ChainA != b.ChainA {
			return a.ChainA < b.ChainA
		}
		return a.ChainB < b.ChainB
	})

	for k, g := range contacts {
		lo, hi := math.Inf(1), math.Inf(-1)
		for _, d := range g.dists {
			if d < lo {
				lo = d
			}
			if d > hi {
				hi = d
			}
		}
		spread := hi - lo
		if len(g.present) == n && spread <= DistanceChangeTolerance {
			continue // unchanged everywhere: not reported
		}
		res.Changes = append(res.Changes, ContactChange{
			ChainA: k.chainA, SeqA: k.seqA, ICodeA: k.icodeA,
			ChainB: k.chainB, SeqB: k.seqB, ICodeB: k.icodeB,
			NameA: g.nameA, NameB: g.nameB,
			PresentIn: g.present, MinDistBySource: g.dists, MaxSpread: spread,
		})
	}
	sort.Slice(res.Changes, func(i, j int) bool {
		a, b := res.Changes[i], res.Changes[j]
		if a.ChainA != b.ChainA {
			return a.ChainA < b.ChainA
		}
		if a.ChainB != b.ChainB {
			return a.ChainB < b.ChainB
		}
		if a.SeqA != b.SeqA {
			return a.SeqA < b.SeqA
		}
		if a.ICodeA != b.ICodeA {
			return a.ICodeA < b.ICodeA
		}
		if a.SeqB != b.SeqB {
			return a.SeqB < b.SeqB
		}
		return a.ICodeB < b.ICodeB
	})
	return res
}
