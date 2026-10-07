package store

import (
	"errors"
	"fmt"
	"testing"

	"github.com/xingjianluo/protein-structure-analyzer-go/contact"
)

// pdbWithMetal writes chain A (GLY) plus a ZN HETATM on chain B, 3 A away.
func pdbWithMetal() string {
	s := pdbChains("A", 3)
	s = fmt.Sprintf("HETATM%5d %-4s%1s%-3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n",
		9, " ZN", " ", "ZN", "B", 1, 0.0, 3.0, 0.0, 1.0, 0.0, "ZN") + s
	return s
}

func TestSnapshotCarriesComparisonAndOptions(t *testing.T) {
	p1 := writeFix(t, "c1.pdb", pdbChains("AB", 3))
	p2 := writeFix(t, "c2.pdb", pdbChains("AB", 3))

	st := New(contact.NewAnalyzer(4.5))
	if err := st.Load([]string{p2, p1}); err != nil {
		t.Fatal(err)
	}
	snap := st.Current()
	if snap.Comparison == nil {
		t.Fatal("snapshot must carry the batch comparison")
	}
	cmp := snap.Comparison
	if len(cmp.Common) != 1 || cmp.Common[0].ChainA != "A" || cmp.Common[0].ChainB != "B" {
		t.Fatalf("want common A-B interface, got %+v", cmp.Common)
	}
	if len(cmp.Sources) != 2 || cmp.Sources[0] >= cmp.Sources[1] {
		t.Fatalf("sources must be sorted, got %v", cmp.Sources)
	}
	if snap.Options.Fingerprint() != contact.DefaultOptions().Fingerprint() {
		t.Fatalf("snapshot options not recorded: %+v", snap.Options)
	}
	t.Logf("comparison fingerprint=%s sources=%v common=%s-%s presentIn=%v",
		cmp.Fingerprint, cmp.Sources, cmp.Common[0].ChainA, cmp.Common[0].ChainB, cmp.Common[0].PresentIn)
}

func TestOptionChangeInvalidatesSnapshotConclusion(t *testing.T) {
	p := writeFix(t, "metal.pdb", pdbWithMetal())
	st := New(contact.NewAnalyzer(4.5))

	if err := st.Load([]string{p}); err != nil {
		t.Fatal(err)
	}
	before := st.Current()
	if len(before.Comparison.Common) != 0 || len(before.Comparison.Partial) != 0 {
		t.Fatalf("default must exclude ZN from interfaces, got %+v", before.Comparison.Common)
	}
	t.Log("default load: ZN HETATM excluded => no interface")

	err := st.LoadWithOptions([]string{p}, contact.Options{Cutoff: 4.5, IncludeHetero: true})
	if err != nil {
		t.Fatal(err)
	}
	after := st.Current()
	if after == before {
		t.Fatal("option change must commit a new snapshot")
	}
	if after.Comparison.Fingerprint == before.Comparison.Fingerprint {
		t.Fatal("option change must change the comparison fingerprint")
	}
	if len(after.Comparison.Common) != 1 {
		t.Fatalf("with IncludeHetero want 1 interface, got %+v", after.Comparison.Common)
	}
	if !after.Options.IncludeHetero {
		t.Fatal("snapshot must reflect the IncludeHetero switch")
	}
	t.Logf("IncludeHetero=true: interface %s-%s appears; fingerprint %s -> %s",
		after.Comparison.Common[0].ChainA, after.Comparison.Common[0].ChainB,
		before.Comparison.Fingerprint, after.Comparison.Fingerprint)
}

func TestFailedBatchReloadRetainsPreviousComparison(t *testing.T) {
	good := writeFix(t, "g1.pdb", pdbChains("AB", 3))
	good2 := writeFix(t, "g2.pdb", pdbChains("AB", 3))
	// Unopenable (missing) file and an unparsable one both trigger rollback.
	missing := good + ".does-not-exist"
	broken := writeFix(t, "broken.pdb", "ATOM      1  CA  GLY A   1       X.XXX  Y.YYY  Z.ZZZ\n")

	st := New(contact.NewAnalyzer(4.5))
	if err := st.Load([]string{good, good2}); err != nil {
		t.Fatal(err)
	}
	before := st.Current()

	for _, bad := range []string{missing, broken} {
		err := st.Load([]string{good, bad})
		var re *ReloadError
		if !errors.As(err, &re) {
			t.Fatalf("bad input %s: want *ReloadError, got %T: %v", bad, err, err)
		}
		if st.Current() != before {
			t.Fatalf("bad input %s: previous snapshot must be retained", bad)
		}
		t.Logf("bad input %s => rollback, retained v%d comparison %s",
			bad, re.Version, st.Current().Comparison.Fingerprint)
	}
	if st.Current().Comparison != before.Comparison {
		t.Fatal("retained snapshot must keep the previous comparison object")
	}
}

func TestDegradedEntryExcludedFromComparisonButListed(t *testing.T) {
	// ok.pdb: healthy A-B interface. deg.pdb: ALA missing CB => degraded but
	// still parsed; its interface still counts (degraded != excluded).
	ok := writeFix(t, "ok.pdb", pdbChains("AB", 3))
	degContent := ""
	for _, a := range []struct {
		n, el string
		dx    float64
	}{{"N", "N", -1.4}, {"CA", "C", 0}, {"C", "C", 1.4}, {"O", "O", 1.4}} {
		dy := 0.0
		if a.n == "O" {
			dy = 1.1
		}
		degContent += fmt.Sprintf(
			"ATOM  %5d %-4s%1s%-3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n",
			1, " "+a.n, " ", "ALA", "A", 1, a.dx, dy, 0.0, 1.0, 0.0, a.el)
	}
	deg := writeFix(t, "deg.pdb", degContent+"END\n")

	st := New(contact.NewAnalyzer(4.5))
	if err := st.Load([]string{ok, deg}); err != nil {
		t.Fatalf("degraded entry must not fail the batch: %v", err)
	}
	snap := st.Current()
	if snap.State != StateDegraded {
		t.Fatalf("batch with degraded entry must be degraded, got %s", snap.State)
	}
	if snap.Entries[ok].State != StateOK || snap.Entries[deg].State != StateDegraded {
		t.Fatalf("per-entry states wrong: ok=%s deg=%s",
			snap.Entries[ok].State, snap.Entries[deg].State)
	}
	// The degraded file has no usable residues (ALA missing CB) => no
	// interface there; the comparison must reflect that per-source.
	if len(snap.Comparison.Common) != 0 {
		t.Fatalf("degraded file has no usable residues; want no common interface, got %+v",
			snap.Comparison.Common)
	}
	if len(snap.Comparison.Partial) != 1 || len(snap.Comparison.Partial[0].PresentIn) != 1 {
		t.Fatalf("A-B interface should be partial (only in ok.pdb), got %+v", snap.Comparison.Partial)
	}
	t.Logf("degraded entry kept but distinguishable: state=%s; A-B partial presentIn=%v",
		snap.Entries[deg].State, snap.Comparison.Partial[0].PresentIn)
}
