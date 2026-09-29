package structure

import (
	"errors"
	"testing"
)

const mmCIFTwoChains = `data_test
_entry.id test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . GLY A 1 ? 0.000 0.000 0.000 1.00 1
ATOM 2 C CA . GLY A 1 ? 1.500 0.000 0.000 1.00 1
ATOM 3 C C . GLY A 1 ? 3.000 0.000 0.000 1.00 1
ATOM 4 O O . GLY A 1 ? 3.500 1.000 0.000 1.00 1
ATOM 5 N N . GLY B 1 ? 0.000 2.000 0.000 1.00 1
ATOM 6 C CA . GLY B 1 ? 1.500 2.000 0.000 1.00 1
ATOM 7 C C . GLY B 1 ? 3.000 2.000 0.000 1.00 1
ATOM 8 O O . GLY B 1 ? 3.500 3.000 0.000 1.00 1
#
`

func TestParseMMCIFFixture(t *testing.T) {
	p := writeTemp(t, "two.cif", mmCIFTwoChains)
	s, err := ParseMMCIFFile(p)
	if err != nil {
		t.Fatalf("parse mmcif: %v", err)
	}
	if len(s.Effective.Chains) != 2 {
		t.Fatalf("want 2 chains, got %d", len(s.Effective.Chains))
	}
	if s.HeavyAtomCount() != 8 {
		t.Fatalf("want 8 heavy atoms, got %d", s.HeavyAtomCount())
	}
	t.Logf("mmCIF parsed: chains=%v heavy=%d",
		[]string{s.Effective.Chains[0].ID, s.Effective.Chains[1].ID}, s.HeavyAtomCount())
}

func TestParseMMCIFMultiModel(t *testing.T) {
	c := mmCIFTwoChains +
		"ATOM 9 C CA . GLY A 1 ? 9.000 9.000 9.000 1.00 2\n"
	p := writeTemp(t, "multi.cif", c)
	s, err := ParseMMCIFFile(p)
	if err != nil {
		t.Fatalf("parse: %v", err)
	}
	if !s.HasIssue(IssueMultiModel) {
		t.Fatalf("want multi-model issue, got %+v", s.Warnings)
	}
	if s.Effective.Serial != 1 {
		t.Fatalf("want model 1 selected, got %d", s.Effective.Serial)
	}
	t.Logf("mmCIF multi-model selects serial %d deterministically", s.Effective.Serial)
}

func TestParseMMCIFNoAtoms(t *testing.T) {
	p := writeTemp(t, "none.cif", "data_empty\n_entry.id x\n")
	_, err := ParseMMCIFFile(p)
	var pe *ParseError
	if !errors.As(err, &pe) {
		t.Fatalf("want *ParseError, got %v", err)
	}
	t.Logf("mmCIF without atom_site rejected: %v", err)
}
