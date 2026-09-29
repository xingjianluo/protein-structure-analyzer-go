package store

import (
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"sync"
	"testing"

	"github.com/xingjianluo/protein-structure-analyzer-go/contact"
)

func pdbChains(chainIDs string, gap float64) string {
	// each chain: one GLY residue with N,CA,C,O; chain k placed at y=k*gap
	var b []byte
	for k := 0; k < len(chainIDs); k++ {
		id := string(chainIDs[k])
		y := float64(k) * gap
		for _, a := range []struct {
			n, el string
			dx    float64
		}{
			{"N", "N", -1.4}, {"CA", "C", 0}, {"C", "C", 1.4}, {"O", "O", 1.4},
		} {
			dy := 0.0
			if a.n == "O" {
				dy = 1.1
			}
			b = append(b, fmt.Sprintf(
				"ATOM  %5d %-4s%1s%-3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n",
				k*4+1, " "+a.n, " ", "GLY", id, 1, a.dx, y+dy, 0.0, 1.0, 0.0, a.el)...)
		}
	}
	return string(b) + "END\n"
}

func writeFix(t *testing.T, name, content string) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), name)
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatal(err)
	}
	t.Logf("fixture %s content:\n%s", name, content)
	return p
}

func TestLoadOKAndGet(t *testing.T) {
	p := writeFix(t, "ok.pdb", pdbChains("AB", 3))
	st := New(contact.NewAnalyzer(4.5))
	if err := st.Load([]string{p}); err != nil {
		t.Fatalf("load: %v", err)
	}
	snap := st.Current()
	if snap == nil || snap.State != StateOK || snap.Version != 1 {
		t.Fatalf("want ok v1 snapshot, got %+v", snap)
	}
	rep, state, ok := st.Get(p)
	if !ok || state != StateOK || rep == nil {
		t.Fatalf("get failed: ok=%v state=%v rep=%v", ok, state, rep)
	}
	if len(rep.Interfaces) != 1 {
		t.Fatalf("want 1 interface, got %d", len(rep.Interfaces))
	}
	t.Logf("loaded v%d state=%s interface=%s-%s contacts=%d",
		snap.Version, snap.State,
		rep.Interfaces[0].ChainA, rep.Interfaces[0].ChainB, rep.Interfaces[0].ContactCount)
}

func TestReloadFailureRollsBack(t *testing.T) {
	good := writeFix(t, "good.pdb", pdbChains("AB", 3))
	bad := writeFix(t, "bad.pdb", "HEADER broken\nno atoms here\nEND\n")

	st := New(contact.NewAnalyzer(4.5))
	if err := st.Load([]string{good}); err != nil {
		t.Fatalf("initial load: %v", err)
	}
	before := st.Current()

	err := st.Load([]string{good, bad})
	var re *ReloadError
	if !errors.As(err, &re) {
		t.Fatalf("expected *ReloadError, got %T: %v", err, err)
	}
	if !re.Retained || re.Version != before.Version {
		t.Fatalf("rollback metadata wrong: %+v", re)
	}
	after := st.Current()
	if after != before {
		t.Fatal("previous snapshot must be retained on failed reload")
	}
	if len(after.Reports) != 1 {
		t.Fatalf("retained snapshot should have 1 report, got %d", len(after.Reports))
	}
	t.Logf("reload failed distinctly (%T) and retained v%d: %v", err, after.Version, err)
}

func TestDegradedStateDistinguished(t *testing.T) {
	// ALA missing CB => degraded but still loaded.
	content := ""
	for _, a := range []struct {
		n, el string
		dx    float64
	}{
		{"N", "N", -1.4}, {"CA", "C", 0}, {"C", "C", 1.4}, {"O", "O", 1.4},
	} {
		dy := 0.0
		if a.n == "O" {
			dy = 1.1
		}
		content += fmt.Sprintf(
			"ATOM  %5d %-4s%1s%-3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n",
			1, " "+a.n, " ", "ALA", "A", 1, a.dx, dy, 0.0, 1.0, 0.0, a.el)
	}
	p := writeFix(t, "deg.pdb", content+"END\n")
	st := New(contact.NewAnalyzer(4.5))
	if err := st.Load([]string{p}); err != nil {
		t.Fatalf("degraded input should still load: %v", err)
	}
	if st.Current().State != StateDegraded {
		t.Fatalf("want degraded state, got %s", st.Current().State)
	}
	_, state, _ := st.Get(p)
	if state != StateDegraded {
		t.Fatalf("want degraded entry state, got %s", state)
	}
	t.Logf("degraded load state=%s, distinguished from error and ok", state)
}

func TestReloadReplacesOldConclusions(t *testing.T) {
	contactP := writeFix(t, "contact.pdb", pdbChains("AB", 3))
	separatedP := writeFix(t, "separated.pdb", pdbChains("AB", 50))

	st := New(contact.NewAnalyzer(4.5))
	if err := st.Load([]string{contactP}); err != nil {
		t.Fatal(err)
	}
	oldRep, _, _ := st.Get(contactP)
	if len(oldRep.Interfaces) != 1 {
		t.Fatalf("fixture: want 1 interface, got %d", len(oldRep.Interfaces))
	}

	// Reload with content where the same path key is replaced by far chains:
	// overwrite the file on disk, then reload.
	if err := os.WriteFile(contactP, []byte(pdbChains("AB", 50)), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := st.Load([]string{contactP}); err != nil {
		t.Fatalf("reload: %v", err)
	}
	newSnap := st.Current()
	if newSnap.Version != 2 {
		t.Fatalf("want version 2, got %d", newSnap.Version)
	}
	newRep, _, _ := st.Get(contactP)
	if len(newRep.Interfaces) != 0 {
		t.Fatalf("after reload want 0 interfaces, got %d (stale conclusion?)", len(newRep.Interfaces))
	}
	// Old report pointer must no longer be reachable via the store.
	cur, _, _ := st.Get(contactP)
	if cur == oldRep {
		t.Fatal("store still returns the old report after reload")
	}
	t.Logf("reload v%d replaced old contact conclusion (interfaces: %d -> 0)",
		newSnap.Version, len(newRep.Interfaces))

	// Loading a different file set removes the old key entirely.
	if err := st.Load([]string{separatedP}); err != nil {
		t.Fatal(err)
	}
	if _, _, ok := st.Get(contactP); ok {
		t.Fatal("old path key must disappear after loading a different file set")
	}
	t.Log("old keys removed after loading a new file set")
}

func TestConcurrentReadWhileReload(t *testing.T) {
	p := writeFix(t, "race.pdb", pdbChains("AB", 3))
	st := New(contact.NewAnalyzer(4.5))
	if err := st.Load([]string{p}); err != nil {
		t.Fatal(err)
	}

	const readers = 16
	const rounds = 200
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
					if snap == nil {
						t.Error("reader saw nil snapshot after initial load")
						return
					}
					rep, state, ok := st.Get(p)
					if !ok || state == StateError || rep == nil {
						t.Errorf("reader observed inconsistent get: ok=%v state=%v rep=%v", ok, state, rep)
						return
					}
					// Every committed snapshot must be self-consistent.
					if len(snap.Reports) != len(snap.Entries) {
						t.Errorf("half-updated snapshot: reports=%d entries=%d",
							len(snap.Reports), len(snap.Entries))
						return
					}
				}
			}
		}()
	}

	wg.Add(1)
	go func() {
		defer wg.Done()
		for i := 0; i < rounds; i++ {
			gap := 3.0
			if i%2 == 1 {
				gap = 50
			}
			if err := os.WriteFile(p, []byte(pdbChains("AB", gap)), 0o644); err != nil {
				t.Error(err)
				return
			}
			if err := st.Load([]string{p}); err != nil {
				t.Errorf("reload %d: %v", i, err)
				return
			}
		}
	}()

	// Failed reloads interleave too and must not disturb readers.
	wg.Add(1)
	go func() {
		defer wg.Done()
		bad := filepath.Join(t.TempDir(), "bad.pdb")
		_ = os.WriteFile(bad, []byte("nothing\n"), 0o644)
		for i := 0; i < 100; i++ {
			if err := st.Load([]string{p, bad}); err == nil {
				t.Error("expected failed reload error")
				return
			}
		}
	}()

	// Let readers run until reloads settle, then stop.
	for i := 0; i < rounds; i++ {
		st.Current()
	}
	close(stop)
	wg.Wait()
	t.Logf("concurrent test finished: %d successful reloads interleaved with %d readers",
		rounds, readers)
}
