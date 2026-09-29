package structure

// Heavy atoms exclude hydrogen (element "H" or "D").
func isHeavy(el string) bool {
	return el != "H" && el != "D" && el != ""
}

// CountIssues returns how many issues of each kind were collected.
func (s *Structure) CountIssues() map[IssueKind]int {
	m := map[IssueKind]int{}
	for _, w := range s.Warnings {
		m[w.Kind]++
	}
	return m
}

// IsDegraded reports whether the structure parsed but contains explicit data
// quality problems (missing atoms, alt-loc collapse, malformed records, ...).
func (s *Structure) IsDegraded() bool {
	for _, w := range s.Warnings {
		if w.Kind == IssueMissingAtom || w.Kind == IssueMalformedRecord ||
			w.Kind == IssueTruncatedRecord || w.Kind == IssueAlternateConformation {
			return true
		}
	}
	return false
}

// HeavyAtomCount returns the number of heavy atoms in the effective model.
func (s *Structure) HeavyAtomCount() int {
	if s.Effective == nil {
		return 0
	}
	n := 0
	for _, c := range s.Effective.Chains {
		for _, r := range c.Residues {
			for _, a := range r.Atoms {
				if isHeavy(a.Element) {
					n++
				}
			}
		}
	}
	return n
}

// ContactUsableResidueCount counts residues that may participate in contacts.
func (s *Structure) ContactUsableResidueCount() int {
	if s.Effective == nil {
		return 0
	}
	n := 0
	for _, c := range s.Effective.Chains {
		for _, r := range c.Residues {
			if r.UsableForContact() {
				n++
			}
		}
	}
	return n
}
