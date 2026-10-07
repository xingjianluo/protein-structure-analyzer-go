package store

import (
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"reflect"
	"sync"
	"testing"

	"github.com/xingjianluo/protein-structure-analyzer-go/compare"
	"github.com/xingjianluo/protein-structure-analyzer-go/contact"
)

// pdbWithHeteroChain builds a two-chain PDB where chain B is a HETATM
// (modified residue) record 3 A from chain A.
func pdbWithHeteroChain() string {
	var b []byte
	mk := func(rec, res, chain string, serial int, y float64) {
		for i, a := range []struct {
			n, el string
			dx    float64
		}{
			{"N", "N", -1.4}, {"CA", "C", 0}, {"C", "C", 1.4}, {"O", "O", 1.4},
		} {
			b = append(b, fmt.Sprintf(
				"%-6s%5d %-4s%1s%-3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n",
				rec, serial*4+i, " "+a.n, " ", res, chain, 1, a.dx, y, 0.0, 1.0, 0.0, a.el)...)
		}
	}
	mk("ATOM", "GLY", "A", 0, 0)
	mk("HETATM", "MSE", "B", 1, 3)
	return string(b) + "END\n"
}

func loadPair(t *testing.T, st *Store, contentA, contentB string) (string, string) {
	t.Helper()
	pa := writeFix(t, "a.pdb", contentA)
	pb := writeFix(t, "b.pdb", contentB)
	if err := st.Load([]string{pa, pb}); err != nil {
		t.Fatalf("load pair: %v", err)
	}
	return pa, pb
}

func TestOptionsChangeInvalidatesOldConclusions(t *testing.T) {
	st := New(contact.NewAnalyzer(4.5))
	pa, pb := loadPair(t, st, pdbWithHeteroChain(), pdbWithHeteroChain())

	snap1 := st.Current()
	if len(snap1.Comparison.Interfaces) != 0 {
		t.Fatalf("default options must exclude HETATM chain, got %+v", snap1.Comparison.Interfaces)
	}
	if snap1.Options.IncludeHetero {
		t.Fatal("snapshot must record IncludeHetero=false")
	}
	t.Logf("v%d: options=%+v interfaces=%d (hetero chain excluded by default)",
		snap1.Version, snap1.Options, len(snap1.Comparison.Interfaces))

	// Flip the switch and reload the very same files.
	st.SetOptions(compare.Options{IncludeHetero: true})
	if err := st.Load([]string{pa, pb}); err != nil {
		t.Fatalf("reload with new options: %v", err)
	}
	snap2 := st.Current()
	if snap2.Version != snap1.Version+1 {
		t.Fatalf("want version %d, got %d", snap1.Version+1, snap2.Version)
	}
	if !snap2.Options.IncludeHetero {
		t.Fatal("new snapshot must record IncludeHetero=true")
	}
	if len(snap2.Comparison.Interfaces) != 1 || !snap2.Comparison.Interfaces[0].Shared {
		t.Fatalf("with hetero included the A-B interface must appear: %+v", snap2.Comparison.Interfaces)
	}
	// The store must no longer serve the old conclusion.
	if st.Current() == snap1 {
		t.Fatal("store still serves the pre-switch snapshot")
	}
	// The old snapshot object itself stays coherent for its holders.
	if snap1.Options.IncludeHetero || len(snap1.Comparison.Interfaces) != 0 {
		t.Fatal("old snapshot mutated after option change")
	}
	t.Logf("v%d: options=%+v interfaces=%d; old conclusion replaced, never mixed",
		snap2.Version, snap2.Options, len(snap2.Comparison.Interfaces))
}

func TestTypedFailuresAndRollbackPreservesComparison(t *testing.T) {
	st := New(contact.NewAnalyzer(4.5))
	pa, pb := loadPair(t, st, pdbChains("AB", 3), pdbChains("AB", 3))
	good := st.Current()
	if good.Comparison == nil || good.Comparison.Compared != 2 {
		t.Fatalf("fixture: want comparison over 2 structures, got %+v", good.Comparison)
	}

	missing := filepath.Join(t.TempDir(), "does-not-exist.pdb")
	garbage := writeFix(t, "garbage.pdb", "HEADER broken\nno atoms here\nEND\n")

	err := st.Load([]string{pa, pb, missing, garbage})
	var re *ReloadError
	if !errors.As(err, &re) {
		t.Fatalf("expected *ReloadError, got %T: %v", err, err)
	}
	kinds := map[string]FailureKind{}
	for _, f := range re.Failures {
		kinds[f.Path] = f.Kind
		t.Logf("typed failure: %s kind=%s msg=%s", f.Path, f.Kind, f.Message)
	}
	if kinds[missing] != FailureIO {
		t.Fatalf("unopenable file must be FailureIO, got %s", kinds[missing])
	}
	if kinds[garbage] != FailureParse {
		t.Fatalf("unparsable content must be FailureParse, got %s", kinds[garbage])
	}
	if !re.Retained || re.Version != good.Version {
		t.Fatalf("rollback metadata wrong: %+v", re)
	}
	after := st.Current()
	if after != good {
		t.Fatal("failed reload must retain the previous snapshot")
	}
	if !reflect.DeepEqual(after.Comparison, good.Comparison) {
		t.Fatal("retained comparison result changed after failed reload")
	}
	t.Logf("reload rolled back; retained v%d comparison with %d interface(s)",
		after.Version, len(after.Comparison.Interfaces))
}

func TestDegradedFileExcludedFromComparison(t *testing.T) {
	// ALA missing its CB atom parses but is degraded.
	degraded := ""
	for i, a := range []struct {
		n, el string
		dx    float64
	}{
		{"N", "N", -1.4}, {"CA", "C", 0}, {"C", "C", 1.4}, {"O", "O", 1.4},
	} {
		degraded += fmt.Sprintf(
			"ATOM  %5d %-4s%1s%-3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n",
			i+1, " "+a.n, " ", "ALA", "A", 1, a.dx, 0.0, 0.0, 1.0, 0.0, a.el)
	}
	degraded += "END\n"

	st := New(contact.NewAnalyzer(4.5))
	pOK1 := writeFix(t, "ok1.pdb", pdbChains("AB", 3))
	pOK2 := writeFix(t, "ok2.pdb", pdbChains("AB", 3))
	pDeg := writeFix(t, "deg.pdb", degraded)
	if err := st.Load([]string{pOK1, pOK2, pDeg}); err != nil {
		t.Fatalf("degraded input must not fail the load: %v", err)
	}
	snap := st.Current()
	if snap.State != StateDegraded {
		t.Fatalf("snapshot state must be degraded, got %s", snap.State)
	}
	if snap.Entries[pDeg].State != StateDegraded {
		t.Fatalf("entry state must be degraded, got %s", snap.Entries[pDeg].State)
	}
	comp := snap.Comparison
	if comp.Compared != 2 {
		t.Fatalf("degraded file must be excluded from comparison, compared=%d", comp.Compared)
	}
	classes := map[string]compare.FileClass{}
	for _, f := range comp.Files {
		classes[f.ID] = f.Class
	}
	if classes[pDeg] != compare.ClassDegraded {
		t.Fatalf("compare must classify the degraded file by type, got %s", classes[pDeg])
	}
	for _, it := range comp.Interfaces {
		for _, id := range it.PresentIn {
			if id == pDeg {
				t.Fatalf("degraded file leaked into comparison: %+v", it)
			}
		}
	}
	t.Logf("state=%s compared=%d; degraded file classified %s and excluded",
		snap.State, comp.Compared, classes[pDeg])
}

func TestFingerprintTracksDataNotOrder(t *testing.T) {
	st := New(contact.NewAnalyzer(4.5))
	pa := writeFix(t, "fa.pdb", pdbChains("AB", 3))
	pb := writeFix(t, "fb.pdb", pdbChains("AB", 50))

	if err := st.Load([]string{pa, pb}); err != nil {
		t.Fatal(err)
	}
	fp1 := st.Current().Fingerprint
	// Same data, different submission order: identical fingerprint.
	st2 := New(contact.NewAnalyzer(4.5))
	if err := st2.Load([]string{pb, pa}); err != nil {
		t.Fatal(err)
	}
	if st2.Current().Fingerprint != fp1 {
		t.Fatal("fingerprint must not depend on submission order")
	}
	// Same paths, changed content: fingerprint must change.
	if err := os.WriteFile(pb, []byte(pdbChains("AB", 3)), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := st.Load([]string{pa, pb}); err != nil {
		t.Fatal(err)
	}
	fp2 := st.Current().Fingerprint
	if fp2 == fp1 {
		t.Fatal("fingerprint must change when file content changes")
	}
	t.Logf("fingerprint stable across order, changed on data change: %s... -> %s...", fp1[:12], fp2[:12])
}

func TestConcurrentReadersSeeConsistentComparison(t *testing.T) {
	p := writeFix(t, "race.pdb", pdbWithHeteroChain())
	st := New(contact.NewAnalyzer(4.5))
	if err := st.Load([]string{p}); err != nil {
		t.Fatal(err)
	}

	const readers = 8
	var wg sync.WaitGroup
	stop := make(chan struct{})
	for r := 0; r < readers; r++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for {
				select {
				case <-stop:
					return
				default:
					snap := st.Current()
					if snap == nil || snap.Comparison == nil {
						t.Error("reader saw nil snapshot or comparison")
						return
					}
					// A committed snapshot must be self-consistent: the
					// comparison options always match the snapshot options,
					// and the interface count matches the switch.
					if !reflect.DeepEqual(snap.Options, snap.Comparison.Options) {
						t.Errorf("half-updated snapshot: options=%+v comparison.options=%+v",
							snap.Options, snap.Comparison.Options)
						return
					}
					want := 0
					if snap.Options.IncludeHetero {
						want = 1
					}
					if len(snap.Comparison.Interfaces) != want {
						t.Errorf("inconsistent comparison: options=%+v interfaces=%d",
							snap.Options, len(snap.Comparison.Interfaces))
						return
					}
				}
			}
		}()
	}

	wg.Add(1)
	go func() {
		defer wg.Done()
		for i := 0; i < 100; i++ {
			st.SetOptions(compare.Options{IncludeHetero: i%2 == 1})
			if err := st.Load([]string{p}); err != nil {
				t.Errorf("reload %d: %v", i, err)
				return
			}
		}
	}()

	for i := 0; i < 200; i++ {
		st.Current()
	}
	close(stop)
	wg.Wait()
	t.Log("concurrent readers only ever saw complete, self-consistent snapshots")
}
