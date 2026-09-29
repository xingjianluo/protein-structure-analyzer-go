package contact

import (
	"testing"

	"github.com/xingjianluo/protein-structure-analyzer-go/structure"
)

// buildStructure builds a structure directly from atom specs:
// each spec is {chain, seq, resName, atomName, element, x, y, z}.
func buildStructure(t *testing.T, specs [][]any) *structure.Structure {
	t.Helper()
	type rkey struct {
		chain string
		seq   int
		name  string
	}
	type ckey = string
	_ = ckey("")
	model := &structure.Model{Serial: 1}
	chainIdx := map[string]*structure.Chain{}
	resIdx := map[rkey]*structure.Residue{}
	addChain := func(id string) *structure.Chain {
		if c, ok := chainIdx[id]; ok {
			return c
		}
		c := &structure.Chain{ID: id}
		chainIdx[id] = c
		model.Chains = append(model.Chains, c)
		return c
	}
	for _, sp := range specs {
		chain := sp[0].(string)
		seq := sp[1].(int)
		resName := sp[2].(string)
		atomName := sp[3].(string)
		el := sp[4].(string)
		c := sp[5].([]float64)
		ch := addChain(chain)
		k := rkey{chain, seq, resName}
		r, ok := resIdx[k]
		if !ok {
			r = &structure.Residue{Seq: seq, Name: resName, IsStandardAA: structure.IsStandardAminoAcid(resName)}
			resIdx[k] = r
			ch.Residues = append(ch.Residues, r)
		}
		r.Atoms = append(r.Atoms, structure.Atom{
			Name: atomName, Element: el,
			Coord: structure.Vec3{X: c[0], Y: c[1], Z: c[2]}, Occupancy: 1,
		})
	}
	return &structure.Structure{Source: "mem", Effective: model}
}

func glyAtoms(chain string, seq int, ox, oy, oz float64) [][]any {
	return [][]any{
		{chain, seq, "GLY", "N", "N", []float64{ox - 1.4, oy, oz}},
		{chain, seq, "GLY", "CA", "C", []float64{ox, oy, oz}},
		{chain, seq, "GLY", "C", "C", []float64{ox + 1.4, oy, oz}},
		{chain, seq, "GLY", "O", "O", []float64{ox + 1.4, oy + 1.1, oz}},
	}
}

func TestAnalyzeInterfaceContact(t *testing.T) {
	var specs [][]any
	specs = append(specs, glyAtoms("A", 1, 0, 0, 0)...)
	specs = append(specs, glyAtoms("B", 1, 0, 3, 0)...) // 3 A from A's atoms => contact at 4.5
	s := buildStructure(t, specs)

	a := NewAnalyzer(4.5)
	rep, err := a.Analyze(s)
	if err != nil {
		t.Fatalf("analyze: %v", err)
	}
	if len(rep.Interfaces) != 1 {
		t.Fatalf("want 1 interface, got %d: %+v", len(rep.Interfaces), rep.Interfaces)
	}
	if rep.Interfaces[0].ChainA != "A" || rep.Interfaces[0].ChainB != "B" {
		t.Fatalf("unexpected chain pair %+v", rep.Interfaces[0])
	}
	if len(rep.Contacts) != 1 {
		t.Fatalf("want 1 residue contact, got %d", len(rep.Contacts))
	}
	t.Logf("contact A/GLY1-B/GLY1 minDist=%.3f, cutoff=%.2f => contact",
		rep.Contacts[0].MinDist, rep.Cutoff)
}

func TestAnalyzeNoInterfaceWhenFar(t *testing.T) {
	var specs [][]any
	specs = append(specs, glyAtoms("A", 1, 0, 0, 0)...)
	specs = append(specs, glyAtoms("B", 1, 0, 50, 0)...)
	s := buildStructure(t, specs)
	rep, err := NewAnalyzer(4.5).Analyze(s)
	if err != nil {
		t.Fatalf("analyze: %v", err)
	}
	if len(rep.Interfaces) != 0 || len(rep.Contacts) != 0 {
		t.Fatalf("want no interface, got %+v", rep.Interfaces)
	}
	t.Log("chains 50 A apart: no interface, as expected")
}

func TestAnalyzeDeterministicAndOrderIndependent(t *testing.T) {
	var specs [][]any
	specs = append(specs, glyAtoms("A", 1, 0, 0, 0)...)
	specs = append(specs, glyAtoms("B", 1, 0, 3, 0)...)
	specs = append(specs, glyAtoms("A", 2, 0, 0, 6)...)
	specs = append(specs, glyAtoms("B", 2, 0, 3, 6)...)
	s1 := buildStructure(t, specs)
	// Mirror: submit chains in reverse order (B then A) via rebuilding.
	var rev [][]any
	rev = append(rev, glyAtoms("B", 1, 0, 3, 0)...)
	rev = append(rev, glyAtoms("A", 1, 0, 0, 0)...)
	rev = append(rev, glyAtoms("B", 2, 0, 3, 6)...)
	rev = append(rev, glyAtoms("A", 2, 0, 0, 6)...)
	s2 := buildStructure(t, rev)

	a := NewAnalyzer(4.5)
	r1, err := a.Analyze(s1)
	if err != nil {
		t.Fatal(err)
	}
	r2, err := a.Analyze(s2)
	if err != nil {
		t.Fatal(err)
	}
	if len(r1.Contacts) != len(r2.Contacts) {
		t.Fatalf("contact count differs: %d vs %d", len(r1.Contacts), len(r2.Contacts))
	}
	for i := range r1.Contacts {
		if r1.Contacts[i] != r2.Contacts[i] {
			t.Fatalf("contact %d differs under chain reorder:\n%+v\n%+v", i, r1.Contacts[i], r2.Contacts[i])
		}
	}
	// Repeat calculation yields identical result.
	r3, _ := a.Analyze(s1)
	for i := range r1.Contacts {
		if r1.Contacts[i] != r3.Contacts[i] {
			t.Fatal("repeated analysis not stable")
		}
	}
	t.Logf("deterministic: %d contacts identical under reorder and repeat", len(r1.Contacts))
}

func TestAnalyzeSkipsResidueWithMissingAtoms(t *testing.T) {
	// A: GLY complete at origin; B: GLY with missing C/O only, N+CA near A.
	specs := glyAtoms("A", 1, 0, 0, 0)
	specs = append(specs,
		[]any{"B", 1, "GLY", "N", "N", []float64{0, 2, 0}},
		[]any{"B", 1, "GLY", "CA", "C", []float64{0, 3, 0}},
	)
	s := buildStructure(t, specs)
	b := s.Effective.Chains[1].Residues[0]
	b.Missing = []string{"C", "O"} // simulate parser-level missing atom finding
	if b.UsableForContact() {
		t.Fatal("fixture error: residue should be unusable")
	}
	rep, err := NewAnalyzer(4.5).Analyze(s)
	if err != nil {
		t.Fatal(err)
	}
	if len(rep.Contacts) != 0 {
		t.Fatalf("residue with missing atoms must be skipped, got %d contacts", len(rep.Contacts))
	}
	t.Log("residue with missing atoms excluded from contact calculation")
}

func TestAnalyzeEmptyStructure(t *testing.T) {
	if _, err := NewAnalyzer(4.5).Analyze(nil); err == nil {
		t.Fatal("expected error for nil structure")
	}
}
