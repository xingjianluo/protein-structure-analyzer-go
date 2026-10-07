package compare

import (
	"errors"
	"testing"
	"time"

	"github.com/xingjianluo/protein-structure-analyzer-go/contact"
	"github.com/xingjianluo/protein-structure-analyzer-go/structure"
)

// --- fixtures -------------------------------------------------------------

type atomSpec struct {
	chain   string
	seq     int
	res     string
	atom    string
	el      string
	x, y, z float64
	hetero  bool
}

// buildStructure assembles a Structure directly from atom specs.
func buildStructure(t *testing.T, source string, specs []atomSpec) *structure.Structure {
	t.Helper()
	type rkey struct {
		chain string
		seq   int
	}
	model := &structure.Model{Serial: 1}
	chainIdx := map[string]*structure.Chain{}
	resIdx := map[rkey]*structure.Residue{}
	for _, sp := range specs {
		ch, ok := chainIdx[sp.chain]
		if !ok {
			ch = &structure.Chain{ID: sp.chain}
			chainIdx[sp.chain] = ch
			model.Chains = append(model.Chains, ch)
		}
		k := rkey{sp.chain, sp.seq}
		r, ok := resIdx[k]
		if !ok {
			r = &structure.Residue{
				Seq: sp.seq, Name: sp.res, Hetero: sp.hetero,
				IsStandardAA: structure.IsStandardAminoAcid(sp.res) && !sp.hetero,
			}
			resIdx[k] = r
			ch.Residues = append(ch.Residues, r)
		}
		r.Atoms = append(r.Atoms, structure.Atom{
			Name: sp.atom, Element: sp.el,
			Coord: structure.Vec3{X: sp.x, Y: sp.y, Z: sp.z}, Occupancy: 1,
		})
	}
	return &structure.Structure{Source: source, Effective: model}
}

// gly places a complete GLY residue (N, CA, C, O) at offset (ox, oy, oz).
func gly(chain string, seq int, ox, oy, oz float64) []atomSpec {
	return []atomSpec{
		{chain: chain, seq: seq, res: "GLY", atom: "N", el: "N", x: ox - 1.4, y: oy, z: oz},
		{chain: chain, seq: seq, res: "GLY", atom: "CA", el: "C", x: ox, y: oy, z: oz},
		{chain: chain, seq: seq, res: "GLY", atom: "C", el: "C", x: ox + 1.4, y: oy, z: oz},
		{chain: chain, seq: seq, res: "GLY", atom: "O", el: "O", x: ox + 1.4, y: oy + 1.1, z: oz},
	}
}

func concat(specs ...[]atomSpec) []atomSpec {
	var out []atomSpec
	for _, s := range specs {
		out = append(out, s...)
	}
	return out
}

// --- tests ----------------------------------------------------------------

func TestCompareCommonAndPartialInterfaces(t *testing.T) {
	// s1: A-B and B-C interfaces. s2: A-B only. s3: A-B and B-C.
	mk := func(source string, withBC bool) *structure.Structure {
		specs := concat(
			gly("A", 1, 0, 0, 0),
			gly("B", 1, 0, 3, 0), // 3 A from A/1 => A-B contact
		)
		if withBC {
			specs = append(specs, gly("C", 1, 0, 6, 0)...) // 3 A from B/1 => B-C contact
		} else {
			specs = append(specs, gly("C", 1, 0, 60, 0)...) // far away: no B-C
		}
		return buildStructure(t, source, specs)
	}
	batch := []*structure.Structure{mk("s2", false), mk("s1", true), mk("s3", true)}

	res, err := New(nil).Compare(batch, contact.DefaultOptions())
	if err != nil {
		t.Fatalf("compare: %v", err)
	}
	t.Logf("inputs: s1(A-B,B-C) s2(A-B) s3(A-B,B-C); sources sorted=%v", res.Sources)

	if len(res.Common) != 1 || res.Common[0].ChainA != "A" || res.Common[0].ChainB != "B" {
		t.Fatalf("want common=[A-B], got %+v", res.Common)
	}
	if got := res.Common[0].PresentIn; len(got) != 3 {
		t.Fatalf("A-B must be present in all 3 structures, got %v", got)
	}
	if len(res.Partial) != 1 || res.Partial[0].ChainA != "B" || res.Partial[0].ChainB != "C" {
		t.Fatalf("want partial=[B-C], got %+v", res.Partial)
	}
	want := []string{"s1", "s3"}
	got := res.Partial[0].PresentIn
	if len(got) != 2 || got[0] != want[0] || got[1] != want[1] {
		t.Fatalf("B-C presence want %v, got %v", want, got)
	}
	t.Logf("rationale: A-B in all 3 => common; B-C only in %v => partial", got)
}

func TestCompareDeterministicUnderOrderAndRepeat(t *testing.T) {
	mk := func(source string, shift float64) *structure.Structure {
		return buildStructure(t, source, concat(
			gly("A", 1, 0, 0, 0),
			gly("B", 1, 0, 3+shift, 0),
			gly("A", 2, 0, 0, 8),
			gly("B", 2, 0, 3, 8),
		))
	}
	batch := []*structure.Structure{mk("a", 0), mk("b", 0.2), mk("c", 0)}
	reversed := []*structure.Structure{batch[2], batch[1], batch[0]}

	c := New(nil)
	r1, err := c.Compare(batch, contact.DefaultOptions())
	if err != nil {
		t.Fatal(err)
	}
	r2, err := c.Compare(reversed, contact.DefaultOptions())
	if err != nil {
		t.Fatal(err)
	}
	if r1.Fingerprint != r2.Fingerprint {
		t.Fatalf("fingerprint differs under reorder:\n%s\n%s", r1.Fingerprint, r2.Fingerprint)
	}
	if len(r1.Common) != len(r2.Common) || len(r1.Partial) != len(r2.Partial) || len(r1.Changes) != len(r2.Changes) {
		t.Fatal("result sizes differ under reorder")
	}
	for i := range r1.Changes {
		a, b := r1.Changes[i], r2.Changes[i]
		if a.ChainA != b.ChainA || a.SeqA != b.SeqA || a.MaxSpread != b.MaxSpread {
			t.Fatalf("change %d differs under reorder: %+v vs %+v", i, a, b)
		}
	}
	// Repeat with the same comparer: cached, identical pointer.
	r3, err := c.Compare(batch, contact.DefaultOptions())
	if err != nil {
		t.Fatal(err)
	}
	if r3 != r1 {
		t.Fatal("repeated identical compare should be served from cache")
	}
	t.Logf("deterministic: fingerprint=%s stable under reorder and repeat", r1.Fingerprint)
}

func TestCompareContactChangesReported(t *testing.T) {
	// Contact A1-B1 exists in all three but its distance in s2 is 1.0 A larger
	// (spread > tolerance). Contact A2-B2 exists only in s1.
	near := concat(gly("A", 1, 0, 0, 0), gly("B", 1, 0, 3, 0))
	s1 := buildStructure(t, "s1", concat(near, gly("A", 2, 0, 0, 8), gly("B", 2, 0, 3, 8)))
	s2 := buildStructure(t, "s2", concat(
		gly("A", 1, 0, 0, 0), gly("B", 1, 0, 4, 0), // distance grows by 1.0
		gly("A", 2, 0, 0, 8), gly("B", 2, 0, 60, 8), // A2-B2 absent
	))
	s3 := buildStructure(t, "s3", concat(near, gly("A", 2, 0, 0, 8), gly("B", 2, 0, 60, 8)))

	res, err := New(nil).Compare([]*structure.Structure{s1, s2, s3}, contact.DefaultOptions())
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("inputs: A1-B1 distances s1/s3=3.0 s2=4.0 (spread 1.0 > tol %.1f); A2-B2 only in s1",
		DistanceChangeTolerance)

	var distChange, presenceChange *ContactChange
	for i := range res.Changes {
		ch := &res.Changes[i]
		t.Logf("change: %s/%d-%s/%d presentIn=%v spread=%.2f",
			ch.ChainA, ch.SeqA, ch.ChainB, ch.SeqB, ch.PresentIn, ch.MaxSpread)
		if ch.SeqA == 1 {
			distChange = ch
		}
		if ch.SeqA == 2 {
			presenceChange = ch
		}
	}
	if distChange == nil {
		t.Fatal("A1-B1 distance change not reported")
	}
	if distChange.MaxSpread < DistanceChangeTolerance {
		t.Fatalf("spread %.2f below tolerance but reported as distance change", distChange.MaxSpread)
	}
	if len(distChange.PresentIn) != 3 {
		t.Fatalf("A1-B1 present everywhere, got %v", distChange.PresentIn)
	}
	if presenceChange == nil {
		t.Fatal("A2-B2 presence change not reported")
	}
	if len(presenceChange.PresentIn) != 1 || presenceChange.PresentIn[0] != "s1" {
		t.Fatalf("A2-B2 should be present only in s1, got %v", presenceChange.PresentIn)
	}
}

func TestCompareHeteroSwitch(t *testing.T) {
	// Chain B holds a ZN ion (HETATM) 3 A from A/GLY1 instead of a GLY.
	metal := []atomSpec{
		{chain: "B", seq: 1, res: "ZN", atom: "ZN", el: "ZN", x: 0, y: 3, z: 0, hetero: true},
	}
	specs := concat(gly("A", 1, 0, 0, 0), metal)
	s := buildStructure(t, "metal", specs)

	def := contact.DefaultOptions()
	res, err := New(nil).Compare([]*structure.Structure{s}, def)
	if err != nil {
		t.Fatal(err)
	}
	if len(res.Common) != 0 || len(res.Partial) != 0 {
		t.Fatalf("hetero must be excluded by default, got common=%+v partial=%+v", res.Common, res.Partial)
	}
	t.Log("default: ZN (HETATM) excluded => no interface")

	withHetero := contact.Options{Cutoff: contact.DefaultCutoff, IncludeHetero: true}
	res2, err := New(nil).Compare([]*structure.Structure{s}, withHetero)
	if err != nil {
		t.Fatal(err)
	}
	if len(res2.Common) != 1 {
		t.Fatalf("with IncludeHetero want 1 interface, got %+v", res2.Common)
	}
	t.Logf("IncludeHetero=true: ZN participates => interface %s-%s", res2.Common[0].ChainA, res2.Common[0].ChainB)
}

func TestCompareChainSelection(t *testing.T) {
	s := buildStructure(t, "tri", concat(
		gly("A", 1, 0, 0, 0),
		gly("B", 1, 0, 3, 0),
		gly("C", 1, 0, 6, 0),
	))
	all, err := New(nil).Compare([]*structure.Structure{s}, contact.DefaultOptions())
	if err != nil {
		t.Fatal(err)
	}
	if len(all.Common) != 2 {
		t.Fatalf("want A-B and B-C interfaces, got %+v", all.Common)
	}
	sel := contact.Options{Cutoff: contact.DefaultCutoff, Chains: []string{"A", "B"}}
	res, err := New(nil).Compare([]*structure.Structure{s}, sel)
	if err != nil {
		t.Fatal(err)
	}
	if len(res.Common) != 1 || res.Common[0].ChainB != "B" {
		t.Fatalf("chain selection A,B must drop B-C, got %+v", res.Common)
	}
	t.Logf("chains=[A B]: only interface %s-%s remains; options fingerprint=%s",
		res.Common[0].ChainA, res.Common[0].ChainB, res.Options.Fingerprint())
}

func TestCompareCacheInvalidation(t *testing.T) {
	s1 := buildStructure(t, "s1", concat(gly("A", 1, 0, 0, 0), gly("B", 1, 0, 3, 0)))
	s2 := buildStructure(t, "s2", concat(gly("A", 1, 0, 0, 0), gly("B", 1, 0, 3, 0)))
	c := New(nil)

	r1, err := c.Compare([]*structure.Structure{s1, s2}, contact.DefaultOptions())
	if err != nil {
		t.Fatal(err)
	}
	// Same batch again: cached.
	r1b, err := c.Compare([]*structure.Structure{s2, s1}, contact.DefaultOptions())
	if err != nil {
		t.Fatal(err)
	}
	if r1b != r1 {
		t.Fatal("identical batch must hit cache")
	}
	// Switch changed: old conclusion must not be reused.
	r2, err := c.Compare([]*structure.Structure{s1, s2}, contact.Options{Cutoff: 3.0})
	if err != nil {
		t.Fatal(err)
	}
	if r2 == r1 || r2.Fingerprint == r1.Fingerprint {
		t.Fatal("cutoff change must invalidate cached conclusion")
	}
	// Data changed: old conclusion must not be reused.
	s2b := buildStructure(t, "s2", concat(gly("A", 1, 0, 0, 0), gly("B", 1, 0, 4, 0)))
	r3, err := c.Compare([]*structure.Structure{s1, s2b}, contact.DefaultOptions())
	if err != nil {
		t.Fatal(err)
	}
	if r3.Fingerprint == r1.Fingerprint {
		t.Fatal("data change must invalidate cached conclusion")
	}
	t.Logf("fingerprints: base=%s cutoff3=%s moved=%s (all distinct)",
		r1.Fingerprint, r2.Fingerprint, r3.Fingerprint)
}

func TestCompareBatchErrorTypeAndCacheUntouched(t *testing.T) {
	good := buildStructure(t, "good", concat(gly("A", 1, 0, 0, 0), gly("B", 1, 0, 3, 0)))
	empty := &structure.Structure{Source: "empty"} // no effective model
	c := New(nil)

	prev, err := c.Compare([]*structure.Structure{good}, contact.DefaultOptions())
	if err != nil {
		t.Fatal(err)
	}
	_, err = c.Compare([]*structure.Structure{good, empty}, contact.DefaultOptions())
	var berr *BatchError
	if !errors.As(err, &berr) {
		t.Fatalf("want *BatchError, got %T: %v", err, err)
	}
	var ierr *contact.InputError
	if !errors.As(berr.Failures[0], &ierr) {
		t.Fatalf("failures must carry *contact.InputError, got %T", berr.Failures[0])
	}
	t.Logf("batch error typed: input %d source=%q err=%v", ierr.Index, ierr.Source, ierr.Err)

	// Cache must still serve the pre-failure conclusion for the old batch.
	again, err := c.Compare([]*structure.Structure{good}, contact.DefaultOptions())
	if err != nil {
		t.Fatal(err)
	}
	if again != prev {
		t.Fatal("failed batch must not clobber the cached good result")
	}
	t.Log("after failed batch, previous cached conclusion still served")
}

func TestCompareScalesToBatch(t *testing.T) {
	// 12 structures x 2 chains x 250 residues each; every B residue is 3 A
	// from the same-numbered A residue.
	const structsN, resN = 12, 250
	var batch []*structure.Structure
	for i := 0; i < structsN; i++ {
		var specs []atomSpec
		for r := 1; r <= resN; r++ {
			z := float64(r) * 6
			specs = append(specs, gly("A", r, 0, 0, z)...)
			specs = append(specs, gly("B", r, 0, 3, z)...)
		}
		batch = append(batch, buildStructure(t, "bulk-"+string(rune('a'+i)), specs))
	}
	start := time.Now()
	res, err := New(nil).Compare(batch, contact.DefaultOptions())
	elapsed := time.Since(start)
	if err != nil {
		t.Fatal(err)
	}
	if len(res.Common) != 1 || res.Common[0].ContactCounts[res.Sources[0]] != resN {
		t.Fatalf("want one common A-B interface with %d contacts, got %+v", resN, res.Common)
	}
	if len(res.Changes) != 0 {
		t.Fatalf("identical geometry must yield no changes, got %d", len(res.Changes))
	}
	t.Logf("scale: %d structures x 2 chains x %d residues compared in %v; contacts per structure=%d",
		structsN, resN, elapsed, resN)
}
