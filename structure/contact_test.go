package structure

import (
	"strings"
	"testing"
)

// twoChainFixture 构造两条链：A 链两个残基在原点附近，B 链一个残基，
// 与 A:1 距离 3.0（接触），与 A:2 距离取决于 gap 参数。
func twoChainFixture(gap float64) []string {
	return []string{
		atomLine(1, "CA", "", "GLY", "A", 1, "", 0, 0, 0, 1.0, false),
		atomLine(2, "CB", "", "GLY", "A", 1, "", 1.0, 0, 0, 1.0, false),
		atomLine(3, "CA", "", "ALA", "A", 2, "", 0, gap, 0, 1.0, false),
		atomLine(4, "CA", "", "SER", "B", 1, "", 3.0, 0, 0, 1.0, false),
		atomLine(5, "CB", "", "SER", "B", 1, "", 3.0, 1.0, 0, 1.0, false),
	}
}

func mustParse(t *testing.T, lines []string) *Model {
	t.Helper()
	m, _, err := ParsePDB("fixture", strings.NewReader(strings.Join(lines, "\n")+"\n"))
	if err != nil {
		t.Fatalf("parse fixture: %v", err)
	}
	return m
}

func TestContactWithinCutoff(t *testing.T) {
	m := mustParse(t, twoChainFixture(10.0))
	contacts := Contacts(m, DefaultContactOptions())
	for _, c := range contacts {
		t.Logf("contact %s:%d(%s) <-> %s:%d(%s) min=%.2f pairs=%d",
			c.ChainA, c.SeqA, c.NameA, c.ChainB, c.SeqB, c.NameB, c.MinDistance, c.AtomPairs)
	}
	if len(contacts) != 1 {
		t.Fatalf("expected 1 contact, got %d", len(contacts))
	}
	c := contacts[0]
	if c.ChainA != "A" || c.SeqA != 1 || c.ChainB != "B" || c.SeqB != 1 {
		t.Fatalf("unexpected contact pair: %+v", c)
	}
	if c.MinDistance < 1.9 || c.MinDistance > 2.1 {
		t.Fatalf("expected min distance ~2.0, got %.3f", c.MinDistance)
	}
}

func TestNoContactBeyondCutoff(t *testing.T) {
	lines := []string{
		atomLine(1, "CA", "", "GLY", "A", 1, "", 0, 0, 0, 1.0, false),
		atomLine(2, "CA", "", "SER", "B", 1, "", 100, 0, 0, 1.0, false),
	}
	m := mustParse(t, lines)
	if got := Contacts(m, DefaultContactOptions()); len(got) != 0 {
		t.Fatalf("expected no contacts at 100A separation, got %+v", got)
	}
	if got := Interfaces(m, DefaultContactOptions()); len(got) != 0 {
		t.Fatalf("expected no interfaces, got %+v", got)
	}
}

func TestInterfaceAggregation(t *testing.T) {
	m := mustParse(t, twoChainFixture(3.0)) // A:2 也在阈值内
	ifaces := Interfaces(m, DefaultContactOptions())
	if len(ifaces) != 1 {
		t.Fatalf("expected 1 interface, got %d", len(ifaces))
	}
	ifa := ifaces[0]
	t.Logf("interface %s<->%s contacts=%d hint=%.0f", ifa.ChainA, ifa.ChainB, len(ifa.Contacts), ifa.BuriedAreaHint)
	if ifa.ChainA != "A" || ifa.ChainB != "B" || len(ifa.Contacts) != 2 {
		t.Fatalf("bad interface: %+v", ifa)
	}
}

func TestIncompleteResidueExcluded(t *testing.T) {
	lines := []string{
		atomLine(1, "CA", "", "GLY", "A", 1, "", 0, 0, 0, 1.0, false),
		atomLine(2, "CB", "", "ALA", "A", 2, "", 0, 30.0, 0, 1.0, false), // 缺 CA -> Incomplete
		atomLine(3, "CA", "", "SER", "B", 1, "", 0, 32.0, 0, 1.0, false), // 仅与 A:2 在阈值内
	}
	m, rep, err := ParsePDB("fixture", strings.NewReader(strings.Join(lines, "\n")+"\n"))
	if err != nil {
		t.Fatalf("parse: %v", err)
	}
	if !rep.HasKind(IssueMissingAtoms) {
		t.Fatal("fixture should trigger MISSING_ATOMS")
	}
	contacts := Contacts(m, DefaultContactOptions())
	for _, c := range contacts {
		t.Logf("contact %s:%d <-> %s:%d", c.ChainA, c.SeqA, c.ChainB, c.SeqB)
		if c.SeqA == 2 || c.SeqB == 2 {
			t.Fatalf("incomplete residue must not participate: %+v", c)
		}
	}
	if len(contacts) != 0 {
		t.Fatalf("expected 0 contacts (only A:2 was in range), got %+v", contacts)
	}
}

func TestNonStandardResidueExcludedByDefault(t *testing.T) {
	lines := []string{
		atomLine(1, "CA", "", "GLY", "A", 1, "", 0, 0, 0, 1.0, false),
		atomLine(2, "CA", "", "MSE", "B", 1, "", 2.0, 0, 0, 1.0, true), // 非标准
	}
	m := mustParse(t, lines)
	if got := Contacts(m, DefaultContactOptions()); len(got) != 0 {
		t.Fatalf("non-standard residue must be excluded by default: %+v", got)
	}
	opts := DefaultContactOptions()
	opts.IncludeHetAtoms = true
	got := Contacts(m, opts)
	t.Logf("with IncludeHetAtoms: %+v", got)
	if len(got) != 1 {
		t.Fatalf("expected 1 contact with IncludeHetAtoms, got %+v", got)
	}
}

func TestRepeatedComputationIdentical(t *testing.T) {
	m := mustParse(t, twoChainFixture(3.0))
	first := Interfaces(m, DefaultContactOptions())
	for i := 0; i < 5; i++ {
		got := Interfaces(m, DefaultContactOptions())
		if len(got) != len(first) {
			t.Fatalf("run %d: length differs", i)
		}
		for j := range got {
			if got[j].ChainA != first[j].ChainA || got[j].ChainB != first[j].ChainB ||
				len(got[j].Contacts) != len(first[j].Contacts) {
				t.Fatalf("run %d: result differs at %d", i, j)
			}
			for k := range got[j].Contacts {
				if got[j].Contacts[k] != first[j].Contacts[k] {
					t.Fatalf("run %d: contact %d differs", i, k)
				}
			}
		}
	}
	t.Logf("5 repeated runs identical, %d interfaces", len(first))
}
