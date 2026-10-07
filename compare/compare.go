// Package compare computes deterministic cross-structure interface
// comparisons for a batch of related structures.
package compare

import (
	"fmt"
	"math"
	"sort"

	"github.com/xingjianluo/protein-structure-analyzer-go/contact"
	"github.com/xingjianluo/protein-structure-analyzer-go/structure"
)

// SignificantDistDelta is the minimum difference in interface minimum
// distance (Angstrom) treated as a significant change.
const SignificantDistDelta = 0.5

// Options selects which residues and chains participate in the comparison.
// The zero value compares all chains using standard amino acids only.
type Options struct {
	// IncludeHetero also admits non-standard records (HETATM: metal ions,
	// modified residues, ligands). Default false: they are excluded.
	IncludeHetero bool
	// Chains restricts the comparison to these chain IDs. Empty means all.
	Chains []string
}

// Input is one structure submitted for batch comparison.
type Input struct {
	// ID identifies the input in the report (typically the file path).
	ID     string
	Struct *structure.Structure
}

// DuplicateIDError is returned when two inputs share the same ID, which
// would make the canonical ordering ambiguous.
type DuplicateIDError struct{ ID string }

func (e *DuplicateIDError) Error() string {
	return fmt.Sprintf("compare: duplicate input ID %q", e.ID)
}

// FileClass distinguishes how an input was treated, by type rather than by
// matching message text.
type FileClass int

const (
	ClassUnknown FileClass = iota
	// ClassCompared: the structure took part in the comparison.
	ClassCompared
	// ClassDegraded: parsed with data-quality problems; excluded from the
	// comparison so incomplete data is never treated as normal.
	ClassDegraded
	// ClassExcluded: unusable input (nil or no effective model).
	ClassExcluded
)

func (c FileClass) String() string {
	switch c {
	case ClassCompared:
		return "compared"
	case ClassDegraded:
		return "degraded"
	case ClassExcluded:
		return "excluded"
	default:
		return "unknown"
	}
}

// FileStatus records the classification of one input.
type FileStatus struct {
	ID      string
	Class   FileClass
	Reasons []structure.Issue
}

// StructureContactStat aggregates one chain-pair interface in one structure.
type StructureContactStat struct {
	ID       string
	Present  bool
	Contacts int
	MinDist  float64
}

// InterfacePresence shows in which structures a chain-pair interface exists.
type InterfacePresence struct {
	ChainA     string
	ChainB     string
	PresentIn  []string
	AbsentFrom []string
	// Shared is true when the interface exists in every compared structure.
	Shared bool
}

// ContactChange aggregates per-structure statistics for one chain pair.
type ContactChange struct {
	ChainA string
	ChainB string
	Stats  []StructureContactStat
	// Changed is true when presence differs, contact counts differ, or the
	// minimum distance moves by more than SignificantDistDelta across the
	// compared structures.
	Changed bool
}

// Report is the deterministic result of a batch comparison. Every list is
// sorted into a canonical order, so the same inputs always yield an
// identical report regardless of submission order or repetition.
type Report struct {
	// Options echoes the normalized switches actually used.
	Options Options
	// Files lists every submitted input, sorted by ID, with its class.
	Files []FileStatus
	// Compared is how many inputs actually took part.
	Compared int
	// Interfaces lists every observed chain pair, shared and partial,
	// sorted by (ChainA, ChainB).
	Interfaces []InterfacePresence
	// Changes lists per-structure contact statistics per chain pair,
	// sorted by (ChainA, ChainB).
	Changes []ContactChange
}

// Comparer computes batch comparisons with a fixed contact cutoff.
type Comparer struct {
	analyzer *contact.Analyzer
}

// NewComparer creates a Comparer. A non-positive cutoff selects
// contact.DefaultCutoff.
func NewComparer(cutoff float64) *Comparer {
	return &Comparer{analyzer: contact.NewAnalyzer(cutoff)}
}

// classify decides how one input is treated.
func classify(s *structure.Structure) FileClass {
	if s == nil || s.Effective == nil {
		return ClassExcluded
	}
	if s.IsDegraded() {
		return ClassDegraded
	}
	return ClassCompared
}

// Compare analyzes the batch and returns the deterministic report. The
// result depends only on the set of inputs and the options: inputs are
// sorted by ID before any computation, and every output list is sorted.
func (c *Comparer) Compare(inputs []Input, opts Options) (*Report, error) {
	nopts := normalize(opts)

	// Canonical input order: by ID, stable for ties (rejected below).
	sorted := append([]Input{}, inputs...)
	sort.SliceStable(sorted, func(i, j int) bool { return sorted[i].ID < sorted[j].ID })
	for i := 1; i < len(sorted); i++ {
		if sorted[i].ID == sorted[i-1].ID {
			return nil, &DuplicateIDError{ID: sorted[i].ID}
		}
	}

	rep := &Report{Options: nopts}
	type analyzed struct {
		id  string
		rep *contact.InterfaceReport
	}
	var compared []analyzed

	for _, in := range sorted {
		st := FileStatus{ID: in.ID, Class: classify(in.Struct)}
		if in.Struct != nil && st.Class != ClassExcluded {
			st.Reasons = append(st.Reasons, in.Struct.Warnings...)
		}
		rep.Files = append(rep.Files, st)
		if st.Class != ClassCompared {
			continue
		}
		filtered := filterStructure(in.Struct, nopts)
		ir, err := c.analyzer.Analyze(filtered)
		if err != nil {
			// No effective model after filtering (e.g. chain selection
			// matched nothing): treat as excluded, never silently compared.
			rep.Files[len(rep.Files)-1].Class = ClassExcluded
			continue
		}
		compared = append(compared, analyzed{id: in.ID, rep: ir})
	}
	rep.Compared = len(compared)

	// Union of chain pairs across compared structures.
	type pairKey struct{ a, b string }
	pairs := map[pairKey]map[string]contact.ChainChainInterface{}
	var order []pairKey
	for _, an := range compared {
		for _, it := range an.rep.Interfaces {
			k := pairKey{it.ChainA, it.ChainB}
			m := pairs[k]
			if m == nil {
				m = map[string]contact.ChainChainInterface{}
				pairs[k] = m
				order = append(order, k)
			}
			m[an.id] = it
		}
	}
	sort.Slice(order, func(i, j int) bool {
		if order[i].a != order[j].a {
			return order[i].a < order[j].a
		}
		return order[i].b < order[j].b
	})

	for _, k := range order {
		per := pairs[k]
		pres := InterfacePresence{ChainA: k.a, ChainB: k.b}
		ch := ContactChange{ChainA: k.a, ChainB: k.b}
		var firstCount, firstDist = -1, 0.0
		uniform := true
		for _, an := range compared {
			stat := StructureContactStat{ID: an.id}
			if it, ok := per[an.id]; ok {
				stat.Present = true
				stat.Contacts = it.ContactCount
				stat.MinDist = it.MinDist
				pres.PresentIn = append(pres.PresentIn, an.id)
				if firstCount < 0 {
					firstCount, firstDist = it.ContactCount, it.MinDist
				} else if it.ContactCount != firstCount ||
					math.Abs(it.MinDist-firstDist) > SignificantDistDelta {
					uniform = false
				}
			} else {
				pres.AbsentFrom = append(pres.AbsentFrom, an.id)
				uniform = false
			}
			ch.Stats = append(ch.Stats, stat)
		}
		pres.Shared = len(pres.AbsentFrom) == 0
		ch.Changed = !uniform
		rep.Interfaces = append(rep.Interfaces, pres)
		rep.Changes = append(rep.Changes, ch)
	}
	return rep, nil
}
