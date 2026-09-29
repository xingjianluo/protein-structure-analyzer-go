package structure

// standardAminoAcids is the set of the 20 canonical amino-acid residue names.
var standardAminoAcids = map[string]bool{
	"GLY": true, "ALA": true, "VAL": true, "LEU": true, "ILE": true,
	"PRO": true, "PHE": true, "TYR": true, "TRP": true, "SER": true,
	"THR": true, "CYS": true, "MET": true, "ASN": true, "GLN": true,
	"ASP": true, "GLU": true, "LYS": true, "ARG": true, "HIS": true,
}

// IsStandardAminoAcid reports whether name is one of the 20 standard residues.
func IsStandardAminoAcid(name string) bool {
	return standardAminoAcids[name]
}

// backboneAtoms are the backbone heavy atoms expected for every amino acid.
var backboneAtoms = []string{"N", "CA", "C", "O"}

// sidechainAtoms maps each standard residue to its expected side-chain heavy
// atoms (PDB atom names, including alternate names for ambiguous protons are
// not relevant: heavy atoms only).
var sidechainAtoms = map[string][]string{
	"GLY": {},
	"ALA": {"CB"},
	"VAL": {"CB", "CG1", "CG2"},
	"LEU": {"CB", "CG", "CD1", "CD2"},
	"ILE": {"CB", "CG1", "CG2", "CD1"},
	"PRO": {"CB", "CG", "CD"},
	"PHE": {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
	"TYR": {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"},
	"TRP": {"CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
	"SER": {"CB", "OG"},
	"THR": {"CB", "OG1", "CG2"},
	"CYS": {"CB", "SG"},
	"MET": {"CB", "CG", "SD", "CE"},
	"ASN": {"CB", "CG", "OD1", "ND2"},
	"GLN": {"CB", "CG", "CD", "OE1", "NE2"},
	"ASP": {"CB", "CG", "OD1", "OD2"},
	"GLU": {"CB", "CG", "CD", "OE1", "OE2"},
	"LYS": {"CB", "CG", "CD", "CE", "NZ"},
	"ARG": {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"},
	"HIS": {"CB", "CG", "ND1", "CD2", "CE1", "NE2"},
}

// ExpectedHeavyAtoms returns the complete set of expected heavy-atom names for
// a standard residue. Returns nil for non-standard residues.
func ExpectedHeavyAtoms(name string) []string {
	if !IsStandardAminoAcid(name) {
		return nil
	}
	out := make([]string, 0, len(backboneAtoms)+len(sidechainAtoms[name]))
	out = append(out, backboneAtoms...)
	out = append(out, sidechainAtoms[name]...)
	return out
}
