package compare

import (
	"fmt"
	"reflect"
	"testing"
	"time"

	"github.com/xingjianluo/protein-structure-analyzer-go/structure"
)

// makeResidue builds one usable GLY residue at the given origin.
func makeResidue(seq int, o structure.Vec3) *structure.Residue {
	mk := func(name, el string, dx float64) structure.Atom {
		return structure.Atom{Name: name, Element: el, Coord: structure.Vec3{X: o.X + dx, Y: o.Y, Z: o.Z}, Occupancy: 1}
	}
	return &structure.Residue{
		Seq: seq, Name: "GLY", IsStandardAA: true,
		Atoms: []structure.Atom{mk("N", "N", -1.4), mk("CA", "C", 0), mk("C", "C", 1.4), mk("O", "O", 1.4)},
	}
}

// makeChain builds a chain of n GLY residues spaced 10 A apart along x,
// starting at origin. Intra-chain distances stay well above the cutoff.
func makeChain(id string, n int, o structure.Vec3) *structure.Chain {
	ch := &structure.Chain{ID: id}
	for i := 0; i < n; i++ {
		ch.Residues = append(ch.Residues, makeResidue(i+1, structure.Vec3{X: o.X + float64(i)*10, Y: o.Y, Z: o.Z}))
	}
	return ch
}

// makeStructure assembles an in-memory structure with an effective model.
func makeStructure(id string, chains ...*structure.Chain) *structure.Structure {
	m := &structure.Model{Serial: 1, Chains: chains}
	return &structure.Structure{Source: id, Models: []*structure.Model{m}, Effective: m}
}

// makeHeteroChain builds a chain whose single residue is a non-standard
// (HETATM) record, e.g. a modified residue or metal ion.
func makeHeteroChain(id string, o structure.Vec3) *structure.Chain {
	r := makeResidue(1, o)
	r.Name = "MSE"
	r.Hetero = true
	r.IsStandardAA = false
	return &structure.Chain{ID: id, Residues: []*structure.Residue{r}}
}

func inputIDs(rep *Report) []string {
	ids := make([]string, 0, len(rep.Files))
	for _, f := range rep.Files {
		ids = append(ids, f.ID)
	}
	return ids
}

func TestDeterministicRegardlessOfOrder(t *testing.T) {
	// Three related structures: AB interface everywhere, AC only in two.
	s1 := makeStructure("s1", makeChain("A", 3, structure.Vec3{}), makeChain("B", 3, structure.Vec3{Y: 3}), makeChain("C", 3, structure.Vec3{Y: 3}))
	s2 := makeStructure("s2", makeChain("A", 3, structure.Vec3{}), makeChain("B", 3, structure.Vec3{Y: 3}), makeChain("C", 3, structure.Vec3{Y: 50}))
	s3 := makeStructure("s3", makeChain("A", 3, structure.Vec3{}), makeChain("B", 3, structure.Vec3{Y: 3}), makeChain("C", 3, structure.Vec3{Y: 3}))
	base := []Input{{ID: "s1", Struct: s1}, {ID: "s2", Struct: s2}, {ID: "s3", Struct: s3}}

	c := NewComparer(4.5)
	ref, err := c.Compare(base, Options{})
	if err != nil {
		t.Fatalf("compare: %v", err)
	}
	t.Logf("inputs=%v compared=%d interfaces=%+v", inputIDs(ref), ref.Compared, ref.Interfaces)

	// Every permutation and every repetition must yield an identical report.
	perms := [][]int{{0, 1, 2}, {2, 1, 0}, {1, 2, 0}, {2, 0, 1}, {0, 2, 1}, {1, 0, 2}}
	for _, p := range perms {
		for round := 0; round < 3; round++ {
			var in []Input
			for _, i := range p {
				in = append(in, base[i])
			}
			got, err := c.Compare(in, Options{})
			if err != nil {
				t.Fatalf("compare: %v", err)
			}
			if !reflect.DeepEqual(ref, got) {
				t.Fatalf("permutation %v round %d diverged:\nref=%+v\ngot=%+v", p, round, ref, got)
			}
		}
	}
	t.Log("all 6 permutations x 3 rounds produced identical reports")

	// Sanity of the content itself: AB shared, AC partial.
	var ab, ac *InterfacePresence
	for i := range ref.Interfaces {
		it := &ref.Interfaces[i]
		switch {
		case it.ChainA == "A" && it.ChainB == "B":
			ab = it
		case it.ChainA == "A" && it.ChainB == "C":
			ac = it
		}
	}
	if ab == nil || !ab.Shared || len(ab.PresentIn) != 3 {
		t.Fatalf("AB must be shared by all three: %+v", ab)
	}
	if ac == nil || ac.Shared || !reflect.DeepEqual(ac.PresentIn, []string{"s1", "s3"}) ||
		!reflect.DeepEqual(ac.AbsentFrom, []string{"s2"}) {
		t.Fatalf("AC must be partial (s1,s3 only): %+v", ac)
	}
	t.Logf("rationale: AB present in %v (shared), AC present in %v absent from %v",
		ab.PresentIn, ac.PresentIn, ac.AbsentFrom)
}

func TestHeteroExcludedByDefaultAndIncludedBySwitch(t *testing.T) {
	// Chain B is a HETATM (modified residue) 3 A away from chain A.
	mk := func(id string) *structure.Structure {
		return makeStructure(id, makeChain("A", 2, structure.Vec3{}), makeHeteroChain("B", structure.Vec3{Y: 3}))
	}
	in := []Input{{ID: "wt", Struct: mk("wt")}, {ID: "mut", Struct: mk("mut")}}
	c := NewComparer(4.5)

	def, err := c.Compare(in, Options{})
	if err != nil {
		t.Fatal(err)
	}
	if len(def.Interfaces) != 0 {
		t.Fatalf("default must exclude hetero records, got %+v", def.Interfaces)
	}
	if def.Options.IncludeHetero {
		t.Fatal("report must echo IncludeHetero=false")
	}
	t.Logf("default: hetero chain B excluded -> %d interfaces (options=%+v)", len(def.Interfaces), def.Options)

	inc, err := c.Compare(in, Options{IncludeHetero: true})
	if err != nil {
		t.Fatal(err)
	}
	if len(inc.Interfaces) != 1 || !inc.Interfaces[0].Shared {
		t.Fatalf("with IncludeHetero the A-B interface must appear shared: %+v", inc.Interfaces)
	}
	if !inc.Options.IncludeHetero {
		t.Fatal("report must echo IncludeHetero=true")
	}
	t.Logf("switch on: A-B interface present in %v (options=%+v)", inc.Interfaces[0].PresentIn, inc.Options)
}

func TestChainSelection(t *testing.T) {
	mk := func(id string) *structure.Structure {
		return makeStructure(id,
			makeChain("A", 2, structure.Vec3{}),
			makeChain("B", 2, structure.Vec3{Y: 3}),
			makeChain("C", 2, structure.Vec3{Y: 6}))
	}
	in := []Input{{ID: "s1", Struct: mk("s1")}, {ID: "s2", Struct: mk("s2")}}
	c := NewComparer(4.5)

	all, err := c.Compare(in, Options{})
	if err != nil {
		t.Fatal(err)
	}
	// Chains at Y=0,3,6: A-B (3 A) and B-C (3 A) are in range, A-C is not.
	if len(all.Interfaces) != 2 {
		t.Fatalf("all chains: want 2 interfaces, got %+v", all.Interfaces)
	}

	sel, err := c.Compare(in, Options{Chains: []string{"B", "A"}}) // unordered on purpose
	if err != nil {
		t.Fatal(err)
	}
	if len(sel.Interfaces) != 1 || sel.Interfaces[0].ChainA != "A" || sel.Interfaces[0].ChainB != "B" {
		t.Fatalf("selection A,B must keep only A-B: %+v", sel.Interfaces)
	}
	if !reflect.DeepEqual(sel.Options.Chains, []string{"A", "B"}) {
		t.Fatalf("normalized selection must be sorted, got %v", sel.Options.Chains)
	}
	t.Logf("selection %v -> only interface %s-%s; normalized options=%+v",
		[]string{"B", "A"}, sel.Interfaces[0].ChainA, sel.Interfaces[0].ChainB, sel.Options)
}

func TestDegradedAndExcludedClassification(t *testing.T) {
	ok1 := makeStructure("ok1", makeChain("A", 2, structure.Vec3{}), makeChain("B", 2, structure.Vec3{Y: 3}))
	ok2 := makeStructure("ok2", makeChain("A", 2, structure.Vec3{}), makeChain("B", 2, structure.Vec3{Y: 3}))
	degraded := makeStructure("deg", makeChain("A", 2, structure.Vec3{}), makeChain("B", 2, structure.Vec3{Y: 3}))
	degraded.Warnings = []structure.Issue{{Kind: structure.IssueMissingAtom, Message: "ALA 1: missing CB"}}
	empty := &structure.Structure{Source: "empty"}

	c := NewComparer(4.5)
	rep, err := c.Compare([]Input{
		{ID: "ok1", Struct: ok1}, {ID: "deg", Struct: degraded},
		{ID: "empty", Struct: empty}, {ID: "ok2", Struct: ok2},
	}, Options{})
	if err != nil {
		t.Fatal(err)
	}
	classes := map[string]FileClass{}
	for _, f := range rep.Files {
		classes[f.ID] = f.Class
	}
	t.Logf("classification: ok1=%s ok2=%s deg=%s empty=%s (degraded reason: %v)",
		classes["ok1"], classes["ok2"], classes["deg"], classes["empty"], rep.Files[1].Reasons)

	if classes["ok1"] != ClassCompared || classes["ok2"] != ClassCompared {
		t.Fatalf("healthy inputs must be compared: %v", classes)
	}
	if classes["deg"] != ClassDegraded {
		t.Fatalf("degraded input must be ClassDegraded, got %s", classes["deg"])
	}
	if classes["empty"] != ClassExcluded {
		t.Fatalf("empty input must be ClassExcluded, got %s", classes["empty"])
	}
	if rep.Compared != 2 {
		t.Fatalf("only the 2 healthy inputs may be compared, got %d", rep.Compared)
	}
	// The degraded structure must not leak into any presence list.
	for _, it := range rep.Interfaces {
		for _, id := range it.PresentIn {
			if id == "deg" || id == "empty" {
				t.Fatalf("degraded/excluded input %q participated in %+v", id, it)
			}
		}
	}
	if len(rep.Interfaces) != 1 || !rep.Interfaces[0].Shared {
		t.Fatalf("A-B must be shared among the 2 compared structures: %+v", rep.Interfaces)
	}
	t.Logf("degraded/excluded inputs kept out of comparison; compared=%d", rep.Compared)
}

func TestContactChangeDetection(t *testing.T) {
	near := makeStructure("near", makeChain("A", 2, structure.Vec3{}), makeChain("B", 2, structure.Vec3{Y: 3}))
	far := makeStructure("far", makeChain("A", 2, structure.Vec3{}), makeChain("B", 2, structure.Vec3{Y: 4.4}))
	gone := makeStructure("gone", makeChain("A", 2, structure.Vec3{}), makeChain("B", 2, structure.Vec3{Y: 50}))
	c := NewComparer(4.5)

	rep, err := c.Compare([]Input{
		{ID: "near", Struct: near}, {ID: "far", Struct: far}, {ID: "gone", Struct: gone},
	}, Options{})
	if err != nil {
		t.Fatal(err)
	}
	if len(rep.Changes) != 1 {
		t.Fatalf("want 1 chain-pair change record, got %+v", rep.Changes)
	}
	ch := rep.Changes[0]
	if !ch.Changed {
		t.Fatalf("contact count/presence differs across structures, must be flagged: %+v", ch)
	}
	for _, st := range ch.Stats {
		t.Logf("  %s: present=%v contacts=%d minDist=%.2f", st.ID, st.Present, st.Contacts, st.MinDist)
		if st.ID == "gone" && st.Present {
			t.Fatal("gone structure must show Present=false")
		}
	}

	// Identical structures must report no change.
	same, err := c.Compare([]Input{
		{ID: "r1", Struct: makeStructure("r1", makeChain("A", 2, structure.Vec3{}), makeChain("B", 2, structure.Vec3{Y: 3}))},
		{ID: "r2", Struct: makeStructure("r2", makeChain("A", 2, structure.Vec3{}), makeChain("B", 2, structure.Vec3{Y: 3}))},
	}, Options{})
	if err != nil {
		t.Fatal(err)
	}
	if len(same.Changes) != 1 || same.Changes[0].Changed {
		t.Fatalf("identical inputs must show no change: %+v", same.Changes)
	}
	t.Logf("identical replicates: Changed=%v (contacts=%d)", same.Changes[0].Changed, same.Changes[0].Stats[0].Contacts)
}

func TestDuplicateIDRejected(t *testing.T) {
	s := makeStructure("x", makeChain("A", 1, structure.Vec3{}))
	c := NewComparer(4.5)
	_, err := c.Compare([]Input{{ID: "dup", Struct: s}, {ID: "dup", Struct: s}}, Options{})
	if _, ok := err.(*DuplicateIDError); !ok {
		t.Fatalf("want *DuplicateIDError, got %T: %v", err, err)
	}
	t.Logf("duplicate IDs rejected by type: %v", err)
}

func TestBatchScaleCompletes(t *testing.T) {
	// Dozens of structures, each with two chains of 1000 residues.
	const nStructs = 30
	const nRes = 1000
	inputs := make([]Input, 0, nStructs)
	for i := 0; i < nStructs; i++ {
		id := fmt.Sprintf("s%02d", i)
		// Interface at 3 A for even structures, 50 A (none) for odd ones.
		gap := 3.0
		if i%2 == 1 {
			gap = 50
		}
		inputs = append(inputs, Input{ID: id, Struct: makeStructure(id,
			makeChain("A", nRes, structure.Vec3{}),
			makeChain("B", nRes, structure.Vec3{Y: gap}),
		)})
	}
	c := NewComparer(4.5)
	start := time.Now()
	rep, err := c.Compare(inputs, Options{})
	el := time.Since(start)
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("batch of %d structures x 2 chains x %d residues compared in %v", nStructs, nRes, el)
	if rep.Compared != nStructs {
		t.Fatalf("all %d structures must be compared, got %d", nStructs, rep.Compared)
	}
	if len(rep.Interfaces) != 1 {
		t.Fatalf("want the single A-B pair, got %+v", rep.Interfaces)
	}
	if got := len(rep.Interfaces[0].PresentIn); got != nStructs/2 {
		t.Fatalf("A-B must be present in %d structures, got %d", nStructs/2, got)
	}
	if el > 30*time.Second {
		t.Fatalf("batch comparison too slow: %v", el)
	}
}
