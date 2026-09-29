package structure

import (
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// writeTemp writes a PDB fixture into t's temp dir and returns its path.
func writeTemp(t *testing.T, name, content string) string {
	t.Helper()
	dir := t.TempDir()
	p := filepath.Join(dir, name)
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatalf("write fixture %s: %v", name, err)
	}
	t.Logf("fixture %s:\n%s", name, content)
	return p
}

// atomLine builds a fixed-width PDB ATOM/HETATM line.
func atomLine(rec, atom, el, res, chain string, seq int, x, y, z float64) string {
	const f = "%-6s%5d %-4s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  "
	return truncPad(
		fmt.Sprintf(f, rec, 1, atom, " ", res, chain, seq, " ", x, y, z, 1.0, 0.0, el),
		80)
}

func truncPad(s string, n int) string {
	if len(s) >= n {
		return s[:n]
	}
	return s + strings.Repeat(" ", n-len(s))
}

// two tiny chains: chain A GLY 1 at origin area, chain B GLY 1 nearby.
func twoChainPDB() string {
	var b strings.Builder
	b.WriteString(atomLine("ATOM", "N", "N", "GLY", "A", 1, 0, 0, 0) + "\n")
	b.WriteString(atomLine("ATOM", "CA", "C", "GLY", "A", 1, 1.5, 0, 0) + "\n")
	b.WriteString(atomLine("ATOM", "C", "C", "GLY", "A", 1, 3.0, 0, 0) + "\n")
	b.WriteString(atomLine("ATOM", "O", "O", "GLY", "A", 1, 3.5, 1.0, 0) + "\n")
	b.WriteString(atomLine("ATOM", "N", "N", "GLY", "B", 1, 0, 2.0, 0) + "\n")
	b.WriteString(atomLine("ATOM", "CA", "C", "GLY", "B", 1, 1.5, 2.0, 0) + "\n")
	b.WriteString(atomLine("ATOM", "C", "C", "GLY", "B", 1, 3.0, 2.0, 0) + "\n")
	b.WriteString(atomLine("ATOM", "O", "O", "GLY", "B", 1, 3.5, 3.0, 0) + "\n")
	b.WriteString("END\n")
	return b.String()
}

func TestParsePDBTwoChainsOK(t *testing.T) {
	p := writeTemp(t, "two.pdb", twoChainPDB())
	s, err := ParsePDBFile(p)
	if err != nil {
		t.Fatalf("expected parse success, got %v", err)
	}
	if len(s.Effective.Chains) != 2 {
		t.Fatalf("want 2 chains, got %d", len(s.Effective.Chains))
	}
	if s.IsDegraded() {
		t.Fatalf("clean fixture must not be degraded: %+v", s.Warnings)
	}
	if got := s.HeavyAtomCount(); got != 8 {
		t.Fatalf("want 8 heavy atoms, got %d", got)
	}
	t.Logf("parsed %s: %d chains, %d heavy atoms, %d warnings",
		p, len(s.Effective.Chains), s.HeavyAtomCount(), len(s.Warnings))
}

func TestParsePDBMultiModelUsesFirst(t *testing.T) {
	var b strings.Builder
	b.WriteString("MODEL        1\n")
	b.WriteString(atomLine("ATOM", "CA", "C", "GLY", "A", 1, 0, 0, 0) + "\n")
	b.WriteString("ENDMDL\n")
	b.WriteString("MODEL        2\n")
	b.WriteString(atomLine("ATOM", "CA", "C", "GLY", "A", 1, 9, 9, 9) + "\n")
	b.WriteString("ENDMDL\n")
	p := writeTemp(t, "multi.pdb", b.String())

	s, err := ParsePDBFile(p)
	if err != nil {
		t.Fatalf("parse: %v", err)
	}
	if !s.HasIssue(IssueMultiModel) {
		t.Fatalf("expected multi-model issue, got %+v", s.Warnings)
	}
	if s.Effective.Serial != 1 {
		t.Fatalf("want effective model 1, got %d", s.Effective.Serial)
	}
	ca := s.Effective.Chains[0].Residues[0].Atoms[0]
	if ca.Coord.X != 0 || ca.Coord.Y != 0 {
		t.Fatalf("expected first model coords (0,0,0), got %+v", ca.Coord)
	}
	t.Logf("multi-model input: selected model %d with coord %+v (reason: %s)",
		s.Effective.Serial, ca.Coord, s.Warnings[0])
}

func TestParsePDBAlternateConformationDeterministic(t *testing.T) {
	// CA has two alt locs A (occ 0.3) and B (occ 0.7): B must win.
	mk := func() string {
		var b strings.Builder
		l1 := strings.Replace(atomLine("ATOM", "CA", "C", "ALA", "A", 1, 0, 0, 0),
			" 1.00  0.00", " 0.30  0.00", 1)
		l1 = l1[:16] + "A" + l1[17:]
		l2 := strings.Replace(atomLine("ATOM", "CA", "C", "ALA", "A", 1, 5, 0, 0),
			" 1.00  0.00", " 0.70  0.00", 1)
		l2 = l2[:16] + "B" + l2[17:]
		b.WriteString(l1 + "\n" + l2 + "\n")
		return b.String()
	}
	p := writeTemp(t, "alt.pdb", mk())
	s1, err := ParsePDBFile(p)
	if err != nil {
		t.Fatalf("parse: %v", err)
	}
	if !s1.HasIssue(IssueAlternateConformation) {
		t.Fatalf("expected altloc issue, got %+v", s1.Warnings)
	}
	res := s1.Effective.Chains[0].Residues[0]
	if res.AltID != 'B' {
		t.Fatalf("want chosen altloc B (higher occupancy), got %q", res.AltID)
	}
	if len(res.Atoms) != 1 || res.Atoms[0].Coord.X != 5 {
		t.Fatalf("want single CA at x=5, got %+v", res.Atoms)
	}
	// Reproducibility: parse again and compare.
	s2, err := ParsePDBFile(p)
	if err != nil {
		t.Fatalf("second parse: %v", err)
	}
	if s2.Effective.Chains[0].Residues[0].AltID != s1.Effective.Chains[0].Residues[0].AltID {
		t.Fatal("altloc selection not reproducible")
	}
	t.Logf("altloc decision: picked '%c' by highest occupancy, reproducible across parses", res.AltID)
}

func TestParsePDBMissingAtomExcludesResidue(t *testing.T) {
	// ALA with backbone but missing CB -> residue flagged and unusable.
	var b strings.Builder
	for _, a := range [][3]any{
		{"N", "N", []float64{0, 0, 0}},
		{"CA", "C", []float64{1.5, 0, 0}},
		{"C", "C", []float64{3.0, 0, 0}},
		{"O", "O", []float64{3.5, 1.0, 0}},
	} {
		c := a[2].([]float64)
		b.WriteString(atomLine("ATOM", a[0].(string), a[1].(string), "ALA", "A", 1, c[0], c[1], c[2]) + "\n")
	}
	p := writeTemp(t, "missing.pdb", b.String())
	s, err := ParsePDBFile(p)
	if err != nil {
		t.Fatalf("parse: %v", err)
	}
	if !s.HasIssue(IssueMissingAtom) {
		t.Fatalf("expected missing-atom issue, got %+v", s.Warnings)
	}
	r := s.Effective.Chains[0].Residues[0]
	if r.UsableForContact() {
		t.Fatal("residue with missing CB must not participate in contacts")
	}
	if len(r.Missing) != 1 || r.Missing[0] != "CB" {
		t.Fatalf("want Missing=[CB], got %v", r.Missing)
	}
	t.Logf("missing-atom residue excluded from contacts; missing=%v", r.Missing)
}

func TestParsePDBNonStandardAndInsertionCode(t *testing.T) {
	var b strings.Builder
	// MSE as HETATM is non-standard but kept; insertion code A on seq 1.
	b.WriteString(strings.Replace(
		atomLine("HETATM", "SE", "SE", "MSE", "A", 1, 0, 0, 0),
		"%-6s", "%-6s", 1))
	b.WriteString("\n")
	// insertion-code residue 1A
	line := atomLine("ATOM", "CA", "C", "GLY", "A", 2, 2, 0, 0)
	line = line[:26] + "A" + line[27:]
	b.WriteString(line + "\n")
	p := writeTemp(t, "nonstd.pdb", b.String())
	s, err := ParsePDBFile(p)
	if err != nil {
		t.Fatalf("parse: %v", err)
	}
	if !s.HasIssue(IssueNonStandardResidue) {
		t.Fatalf("want non-standard issue, got %+v", s.Warnings)
	}
	if !s.HasIssue(IssueInsertionCode) {
		t.Fatalf("want insertion-code issue, got %+v", s.Warnings)
	}
	t.Logf("non-standard/icode issues distinguishable: %d warnings", len(s.Warnings))
}

func TestParsePDBMalformedAndTruncated(t *testing.T) {
	good := atomLine("ATOM", "CA", "C", "GLY", "A", 1, 0, 0, 0)
	cases := []struct {
		name    string
		content string
		want    IssueKind
	}{
		{"truncated", good[:40] + "\n", IssueTruncatedRecord},
		{"bad-x", good[:30] + "   abc   " + good[39:] + "\n", IssueMalformedRecord},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			p := writeTemp(t, tc.name+".pdb", tc.content+"END\n")
			_, err := ParsePDBFile(p)
			// One dropped record means no usable records at all => ParseError.
			if err == nil {
				t.Fatalf("expected ParseError for %s", tc.name)
			}
			var pe *ParseError
			if !errors.As(err, &pe) {
				t.Fatalf("expected *ParseError, got %T", err)
			}
			t.Logf("%s rejected as ParseError: %v", tc.name, err)
		})
	}

	// A bad record alongside good records: bad one dropped with its own issue
	// kind, good structure still usable.
	mixed := good + "\n" + good[:40] + "\n" +
		atomLine("ATOM", "N", "N", "GLY", "B", 1, 0, 0, 0) + "\nEND\n"
	// Add the rest of GLY backbone atoms for both residues so they are usable.
	mixed = strings.ReplaceAll(mixed, "END\n", "")
	mixed += atomLine("ATOM", "CA", "C", "GLY", "A", 1, 1.5, 0, 0) + "\n"
	mixed += atomLine("ATOM", "C", "C", "GLY", "A", 1, 3.0, 0, 0) + "\n"
	mixed += atomLine("ATOM", "O", "O", "GLY", "A", 1, 3.5, 1, 0) + "\n"
	mixed += atomLine("ATOM", "CA", "C", "GLY", "B", 1, 1.5, 0, 0) + "\n"
	mixed += atomLine("ATOM", "C", "C", "GLY", "B", 1, 3.0, 0, 0) + "\n"
	mixed += atomLine("ATOM", "O", "O", "GLY", "B", 1, 3.5, 1, 0) + "\nEND\n"
	p := writeTemp(t, "mixed.pdb", mixed)
	s, err := ParsePDBFile(p)
	if err != nil {
		t.Fatalf("parse mixed: %v", err)
	}
	if !s.HasIssue(IssueTruncatedRecord) {
		t.Fatalf("want truncated-record warning in %+v", s.Warnings)
	}
	t.Logf("mixed fixture: truncated record dropped distinctly, parse continued")
}

func TestParsePDBEmpty(t *testing.T) {
	p := writeTemp(t, "empty.pdb", "HEADER nothing\nEND\n")
	_, err := ParsePDBFile(p)
	var pe *ParseError
	if !errors.As(err, &pe) {
		t.Fatalf("want *ParseError for atom-less file, got %v", err)
	}
	if !strings.Contains(err.Error(), "no ATOM/HETATM") {
		t.Fatalf("unexpected reason: %v", err)
	}
	t.Logf("empty structure reason distinguishable: %v", err)
}
