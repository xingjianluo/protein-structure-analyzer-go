package contact

import (
	"math"
	"sort"

	"github.com/xingjianluo/protein-structure-analyzer-go/structure"
)

// DefaultCutoff is the heavy-atom distance cutoff in Angstrom.
const DefaultCutoff = 4.5

// ResiduePair identifies two contacting residues from different chains.
type ResiduePair struct {
	ChainA  string
	SeqA    int
	ICodeA  byte
	NameA   string
	ChainB  string
	SeqB    int
	ICodeB  byte
	NameB   string
	MinDist float64
}

// ChainInterface aggregates contacts between two chains.
type ChainChainInterface struct {
	ChainA       string
	ChainB       string
	ContactCount int
	MinDist      float64
}

// InterfaceReport is the deterministic result of an interface analysis.
type InterfaceReport struct {
	Source     string
	Model      int
	Cutoff     float64
	Contacts   []ResiduePair
	Interfaces []ChainChainInterface
}

// Analyzer computes residue contacts and chain interfaces.
type Analyzer struct {
	cutoff float64
}

// NewAnalyzer creates an Analyzer with the given distance cutoff. A non-positive
// value selects DefaultCutoff.
func NewAnalyzer(cutoff float64) *Analyzer {
	if cutoff <= 0 {
		cutoff = DefaultCutoff
	}
	return &Analyzer{cutoff: cutoff}
}

// Cutoff returns the configured cutoff.
func (a *Analyzer) Cutoff() float64 { return a.cutoff }

func canonicalChainOrder(a, b string) (string, string) {
	if a <= b {
		return a, b
	}
	return b, a
}

// Analyze computes the interface report with DefaultOptions: non-standard
// records (HETATM, modified residues, metal ions) are excluded and every
// chain participates. Use AnalyzeWithOptions to change that.
func (a *Analyzer) Analyze(s *structure.Structure) (*InterfaceReport, error) {
	return a.AnalyzeWithOptions(s, DefaultOptions())
}

// AnalyzeWithOptions computes the interface report restricted by opts:
// non-standard residues are excluded unless opts.IncludeHetero is set, and
// only the selected chains participate when opts.Chains is non-empty. The
// result depends solely on the structure content and opts: chains are
// processed in model order, and the output lists are sorted into a canonical
// order.
func (a *Analyzer) AnalyzeWithOptions(s *structure.Structure, opts Options) (*InterfaceReport, error) {
	if s == nil || s.Effective == nil {
		return nil, ErrEmptyStructure
	}
	norm := opts.normalized()
	if norm.Cutoff != a.cutoff {
		// Per-call cutoff overrides the analyzer default without mutation.
		a = &Analyzer{cutoff: norm.Cutoff}
	}
	allow := norm.chainSet()
	rep := &InterfaceReport{Source: s.Source, Model: s.Effective.Serial, Cutoff: a.cutoff}
	var chains []*structure.Chain
	for _, c := range s.Effective.Chains {
		if allow == nil || allow[c.ID] {
			chains = append(chains, c)
		}
	}
	usable := func(r *structure.Residue) bool {
		if !r.UsableForContact() {
			return false
		}
		if !norm.IncludeHetero && (r.Hetero || !r.IsStandardAA) {
			return false
		}
		return true
	}
	type agg struct {
		count int
		min   float64
	}
	ifaceAgg := map[string]*agg{}

	for i := 0; i < len(chains); i++ {
		for j := i + 1; j < len(chains); j++ {
			ci, cj := chains[i], chains[j]
			for _, ri := range ci.Residues {
				if !usable(ri) {
					continue
				}
				for _, rj := range cj.Residues {
					if !usable(rj) {
						continue
					}
					d := ResidueMinDistance(ri, rj)
					if d <= a.cutoff {
						ca, cb := canonicalChainOrder(ci.ID, cj.ID)
						pair := ResiduePair{
							ChainA: ca, SeqA: ri.Seq, ICodeA: ri.ICode, NameA: ri.Name,
							ChainB: cb, SeqB: rj.Seq, ICodeB: rj.ICode, NameB: rj.Name,
							MinDist: d,
						}
						if ci.ID != ca {
							pair = ResiduePair{
								ChainA: ca, SeqA: rj.Seq, ICodeA: rj.ICode, NameA: rj.Name,
								ChainB: cb, SeqB: ri.Seq, ICodeB: ri.ICode, NameB: ri.Name,
								MinDist: d,
							}
						}
						rep.Contacts = append(rep.Contacts, pair)
						key := ca + "\x00" + cb
						g := ifaceAgg[key]
						if g == nil {
							g = &agg{min: math.Inf(1)}
							ifaceAgg[key] = g
						}
						g.count++
						if d < g.min {
							g.min = d
						}
					}
				}
			}
		}
	}

	for key, g := range ifaceAgg {
		// split key without importing strings for the separator
		idA, idB := key, ""
		for k := 0; k < len(key); k++ {
			if key[k] == 0 {
				idA, idB = key[:k], key[k+1:]
				break
			}
		}
		rep.Interfaces = append(rep.Interfaces, ChainChainInterface{
			ChainA: idA, ChainB: idB, ContactCount: g.count, MinDist: g.min,
		})
	}

	sort.Slice(rep.Contacts, func(x, y int) bool {
		p, q := rep.Contacts[x], rep.Contacts[y]
		if p.ChainA != q.ChainA {
			return p.ChainA < q.ChainA
		}
		if p.ChainB != q.ChainB {
			return p.ChainB < q.ChainB
		}
		if p.SeqA != q.SeqA {
			return p.SeqA < q.SeqA
		}
		if p.ICodeA != q.ICodeA {
			return p.ICodeA < q.ICodeA
		}
		if p.SeqB != q.SeqB {
			return p.SeqB < q.SeqB
		}
		if p.ICodeB != q.ICodeB {
			return p.ICodeB < q.ICodeB
		}
		return p.NameA < q.NameA
	})
	sort.Slice(rep.Interfaces, func(x, y int) bool {
		p, q := rep.Interfaces[x], rep.Interfaces[y]
		if p.ChainA != q.ChainA {
			return p.ChainA < q.ChainA
		}
		return p.ChainB < q.ChainB
	})
	return rep, nil
}

// AnalyzeMany analyzes every structure and returns a source-keyed map plus
// per-input errors. Results never depend on slice order.
func AnalyzeMany(a *Analyzer, structs []*structure.Structure) (map[string]*InterfaceReport, []error) {
	out := map[string]*InterfaceReport{}
	var errs []error
	for i, s := range structs {
		rep, err := a.Analyze(s)
		if err != nil {
			errs = append(errs, &InputError{Index: i, Source: sourceOf(s), Err: err})
			continue
		}
		out[rep.Source] = rep
	}
	return out, errs
}

func sourceOf(s *structure.Structure) string {
	if s == nil {
		return ""
	}
	return s.Source
}
