package contact

import (
	"sort"
	"strconv"
	"strings"
)

// Options controls which residues and chains participate in interface
// determination. The zero value is NOT the default; use DefaultOptions.
type Options struct {
	// Cutoff is the heavy-atom distance cutoff in Angstrom (<=0 selects DefaultCutoff).
	Cutoff float64
	// IncludeHetero, when false (the default), excludes non-standard records
	// (HETATM, modified residues, metal ions) from interface determination.
	IncludeHetero bool
	// Chains restricts the analysis to the listed chain IDs. Nil or empty
	// means every chain participates.
	Chains []string
}

// DefaultOptions returns the default analysis options: standard amino acids
// only, all chains, default cutoff.
func DefaultOptions() Options {
	return Options{Cutoff: DefaultCutoff}
}

// normalized returns a canonical copy: cutoff resolved and chain IDs sorted
// and deduplicated so that logically identical option sets compare equal.
func (o Options) normalized() Options {
	out := o
	if out.Cutoff <= 0 {
		out.Cutoff = DefaultCutoff
	}
	if len(out.Chains) > 0 {
		cp := append([]string{}, out.Chains...)
		sort.Strings(cp)
		uniq := cp[:0]
		for i, c := range cp {
			if i == 0 || cp[i-1] != c {
				uniq = append(uniq, c)
			}
		}
		out.Chains = append([]string{}, uniq...)
	}
	return out
}

// chainSet returns the normalized chain allow-list, or nil for "all chains".
func (o Options) chainSet() map[string]bool {
	n := o.normalized()
	if len(n.Chains) == 0 {
		return nil
	}
	m := make(map[string]bool, len(n.Chains))
	for _, c := range n.Chains {
		m[c] = true
	}
	return m
}

// Fingerprint returns a stable string identity of the normalized options.
// Any change to cutoff, hetero switch, or chain selection yields a different
// fingerprint, so cached results keyed by it are invalidated automatically.
func (o Options) Fingerprint() string {
	n := o.normalized()
	var b strings.Builder
	b.WriteString("cutoff=")
	b.WriteString(strconv.FormatFloat(n.Cutoff, 'g', -1, 64))
	b.WriteString(";hetero=")
	b.WriteString(strconv.FormatBool(n.IncludeHetero))
	b.WriteString(";chains=")
	b.WriteString(strings.Join(n.Chains, ","))
	return b.String()
}
