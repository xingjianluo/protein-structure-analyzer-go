package compare

import (
	"sort"

	"github.com/xingjianluo/protein-structure-analyzer-go/structure"
)

// normalize returns a canonical copy of the options: the chain selection is
// sorted and de-duplicated so that logically identical selections always
// produce identical reports.
func normalize(opts Options) Options {
	out := Options{IncludeHetero: opts.IncludeHetero}
	if len(opts.Chains) > 0 {
		seen := map[string]bool{}
		for _, c := range opts.Chains {
			if !seen[c] {
				seen[c] = true
				out.Chains = append(out.Chains, c)
			}
		}
		sort.Strings(out.Chains)
	}
	return out
}

// chainAllowed reports whether the chain ID passes the selection.
func (o Options) chainAllowed(id string) bool {
	if len(o.Chains) == 0 {
		return true
	}
	i := sort.SearchStrings(o.Chains, id)
	return i < len(o.Chains) && o.Chains[i] == id
}

// filterStructure builds a filtered view of s honoring the options: only
// selected chains are kept, and unless IncludeHetero is set only standard
// amino-acid residues participate. Residue and atom data are shared with the
// original (analysis never mutates them).
func filterStructure(s *structure.Structure, opts Options) *structure.Structure {
	out := &structure.Structure{Source: s.Source}
	if s.Effective == nil {
		return out
	}
	m := &structure.Model{Serial: s.Effective.Serial}
	for _, ch := range s.Effective.Chains {
		if !opts.chainAllowed(ch.ID) {
			continue
		}
		nc := &structure.Chain{ID: ch.ID}
		for _, r := range ch.Residues {
			if !opts.IncludeHetero && !r.IsStandardAA {
				continue
			}
			nc.Residues = append(nc.Residues, r)
		}
		m.Chains = append(m.Chains, nc)
	}
	out.Models = []*structure.Model{m}
	out.Effective = m
	return out
}
