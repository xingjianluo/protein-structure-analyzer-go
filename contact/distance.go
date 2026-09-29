package contact

import (
	"math"

	"github.com/xingjianluo/protein-structure-analyzer-go/structure"
)

// Distance returns the Euclidean distance between two points in Angstrom.
func Distance(a, b structure.Vec3) float64 {
	dx := a.X - b.X
	dy := a.Y - b.Y
	dz := a.Z - b.Z
	return math.Sqrt(dx*dx + dy*dy + dz*dz)
}

func isHeavy(el string) bool {
	return el != "H" && el != "D" && el != ""
}

// ResidueMinDistance returns the minimum heavy-atom distance between two
// residues. Returns +Inf when either residue has no heavy atoms.
func ResidueMinDistance(r1, r2 *structure.Residue) float64 {
	best := math.Inf(1)
	for i := range r1.Atoms {
		a := &r1.Atoms[i]
		if !isHeavy(a.Element) {
			continue
		}
		for j := range r2.Atoms {
			b := &r2.Atoms[j]
			if !isHeavy(b.Element) {
				continue
			}
			if d := Distance(a.Coord, b.Coord); d < best {
				best = d
			}
		}
	}
	return best
}
