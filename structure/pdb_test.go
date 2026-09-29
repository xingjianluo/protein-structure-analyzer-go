package structure

import (
	"fmt"
	"strings"
	"testing"
)

// atomLine 构造一条符合 PDB 定长列格式的 ATOM/HETATM 记录。
func atomLine(serial int, name, altLoc, resName, chain string, seq int, icode string, x, y, z, occ float64, het bool) string {
	rec := "ATOM  "
	if het {
		rec = "HETATM"
	}
	if altLoc == "" {
		altLoc = " "
	}
	if icode == "" {
		icode = " "
	}
	if chain == "" {
		chain = " "
	}
	return fmt.Sprintf("%s%5d %-4s%s%3s %s%4d%s   %8.3f%8.3f%8.3f%6.2f",
		rec, serial, name, altLoc, resName, chain, seq, icode, x, y, z, occ)
}

func parse(t *testing.T, lines ...string) (*Model, *ParseReport, error) {
	t.Helper()
	m, rep, err := ParsePDB("test.pdb", strings.NewReader(strings.Join(lines, "\n")+"\n"))
	if m != nil {
		t.Logf("parsed: atoms=%d chains=%d issues=%d", rep.Atoms, len(m.Chains), len(rep.Issues))
		for _, is := range rep.Issues {
			t.Logf("  issue [%s] line %d: %s", is.Kind, is.Line, is.Detail)
		}
	}
	return m, rep, err
}

func TestParseBasicSingleModel(t *testing.T) {
	m, rep, err := parse(t,
		atomLine(1, "N", "", "GLY", "A", 1, "", 0, 0, 0, 1.0, false),
		atomLine(2, "CA", "", "GLY", "A", 1, "", 1.5, 0, 0, 1.0, false),
		atomLine(3, "C", "", "GLY", "A", 1, "", 2.5, 1.0, 0, 1.0, false),
	)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if rep.Atoms != 3 || len(m.Chains) != 1 || len(m.Chains[0].Residues) != 1 {
		t.Fatalf("bad model: atoms=%d chains=%+v", rep.Atoms, m.Chains)
	}
	if len(rep.Issues) != 0 {
		t.Fatalf("expected no issues, got %v", rep.Issues)
	}
}

func TestParseMultiModelKeepsFirst(t *testing.T) {
	m, rep, err := parse(t,
		"MODEL        1",
		atomLine(1, "CA", "", "ALA", "A", 1, "", 0, 0, 0, 1.0, false),
		"ENDMDL",
		"MODEL        2",
		atomLine(2, "CA", "", "ALA", "A", 1, "", 9, 9, 9, 1.0, false),
		"ENDMDL",
	)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if !rep.HasKind(IssueMultiModel) {
		t.Fatal("expected MULTI_MODEL issue")
	}
	got := m.Chains[0].Residues[0].Atoms[0]
	t.Logf("decision: kept atom at (%.1f,%.1f,%.1f) from model 1", got.X, got.Y, got.Z)
	if got.X != 0 {
		t.Fatalf("expected first model coordinates, got x=%v", got.X)
	}
}

func TestParseAltLocKeepsHighestOccupancy(t *testing.T) {
	m, rep, err := parse(t,
		atomLine(1, "CA", "A", "LEU", "A", 5, "", 0, 0, 0, 0.4, false),
		atomLine(2, "CA", "B", "LEU", "A", 5, "", 1, 1, 1, 0.6, false),
	)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if !rep.HasKind(IssueAltLoc) {
		t.Fatal("expected ALT_LOC issue")
	}
	res := m.Chains[0].Residues[0]
	if len(res.Atoms) != 1 {
		t.Fatalf("expected 1 atom after altLoc resolution, got %d", len(res.Atoms))
	}
	t.Logf("decision: kept altLoc=%q occ=%.1f", res.Atoms[0].AltLoc, res.Atoms[0].Occupancy)
	if res.Atoms[0].AltLoc != "B" || res.Atoms[0].Occupancy != 0.6 {
		t.Fatalf("expected altLoc B occ 0.6, got %+v", res.Atoms[0])
	}
}

func TestParseInsertionCode(t *testing.T) {
	m, rep, err := parse(t,
		atomLine(1, "CA", "", "GLY", "A", 10, "", 0, 0, 0, 1.0, false),
		atomLine(2, "CA", "", "SER", "A", 10, "A", 5, 0, 0, 1.0, false),
	)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if !rep.HasKind(IssueInsertionCode) {
		t.Fatal("expected INSERTION_CODE issue")
	}
	if len(m.Chains[0].Residues) != 2 {
		t.Fatalf("insertion-code residue must be distinct, got %d residues", len(m.Chains[0].Residues))
	}
	if m.Chains[0].Residues[1].ICode != "A" {
		t.Fatalf("expected ICode A, got %q", m.Chains[0].Residues[1].ICode)
	}
}

func TestParseNonStandardResidue(t *testing.T) {
	m, rep, err := parse(t,
		atomLine(1, "CA", "", "GLY", "A", 1, "", 0, 0, 0, 1.0, false),
		atomLine(2, "ZN", "", "ZN", "A", 2, "", 5, 0, 0, 1.0, true),
	)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if !rep.HasKind(IssueNonStandardRes) {
		t.Fatal("expected NON_STANDARD_RES issue")
	}
	res := m.Chains[0].Residues[1]
	t.Logf("decision: residue %s flagged NonStandard=%v", res.Name, res.NonStandard)
	if !res.NonStandard {
		t.Fatal("HETATM residue must be flagged NonStandard")
	}
}

func TestParseMissingBackboneAtoms(t *testing.T) {
	m, rep, err := parse(t,
		atomLine(1, "CA", "", "GLY", "A", 1, "", 0, 0, 0, 1.0, false),
		atomLine(2, "CB", "", "ALA", "A", 2, "", 5, 0, 0, 1.0, false), // 缺 CA
	)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if !rep.HasKind(IssueMissingAtoms) {
		t.Fatal("expected MISSING_ATOMS issue")
	}
	res := m.Chains[0].Residues[1]
	t.Logf("decision: residue %s:%d Incomplete=%v", res.ChainID, res.Seq, res.Incomplete)
	if !res.Incomplete {
		t.Fatal("residue without CA must be flagged Incomplete")
	}
}

func TestParseMalformedAndTruncated(t *testing.T) {
	_, rep, err := parse(t,
		atomLine(1, "CA", "", "GLY", "A", 1, "", 0, 0, 0, 1.0, false),
		"ATOM      2  CA  GLY A   1      xx.xxx   1.000   1.000  1.00", // 畸形坐标
		"ATOM      3  CA  GLY A   1",                                   // 截断
	)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if !rep.HasKind(IssueMalformedRecord) {
		t.Fatal("expected MALFORMED_RECORD issue")
	}
	if !rep.HasKind(IssueTruncatedRecord) {
		t.Fatal("expected TRUNCATED_RECORD issue")
	}
	if rep.Atoms != 1 {
		t.Fatalf("bad records must be skipped, atoms=%d", rep.Atoms)
	}
}

func TestParseNoAtomsIsError(t *testing.T) {
	_, _, err := parse(t, "HEADER    EMPTY FILE", "END")
	if err == nil {
		t.Fatal("expected ParseError for empty file")
	}
	pe, ok := err.(*ParseError)
	if !ok {
		t.Fatalf("expected *ParseError, got %T", err)
	}
	t.Logf("error reason: %s", pe.Reason)
}

func TestParseDeterministicAcrossAtomOrder(t *testing.T) {
	l1 := atomLine(1, "CA", "", "GLY", "A", 1, "", 0, 0, 0, 1.0, false)
	l2 := atomLine(2, "CA", "", "ALA", "B", 1, "", 3, 0, 0, 1.0, false)
	m1, _, err1 := ParsePDB("a", strings.NewReader(l1+"\n"+l2+"\n"))
	m2, _, err2 := ParsePDB("b", strings.NewReader(l2+"\n"+l1+"\n"))
	if err1 != nil || err2 != nil {
		t.Fatalf("errors: %v %v", err1, err2)
	}
	c1, c2 := Contacts(m1, DefaultContactOptions()), Contacts(m2, DefaultContactOptions())
	t.Logf("contacts forward=%+v reversed=%+v", c1, c2)
	if len(c1) != 1 || len(c2) != 1 || c1[0] != c2[0] {
		t.Fatalf("atom order must not affect result: %+v vs %+v", c1, c2)
	}
	if c1[0].ChainA != "A" || c1[0].ChainB != "B" {
		t.Fatalf("contact direction must be canonical: %+v", c1[0])
	}
}
