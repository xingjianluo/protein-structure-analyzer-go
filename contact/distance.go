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

// bbox is an axis-aligned bounding box over heavy-atom coordinates.
type bbox struct {
	min, max structure.Vec3
}

// residueBBox computes the bounding box of a residue's heavy atoms.
func residueBBox(r *structure.Residue) bbox {
	b := bbox{
		min: structure.Vec3{X: math.Inf(1), Y: math.Inf(1), Z: math.Inf(1)},
		max: structure.Vec3{X: math.Inf(-1), Y: math.Inf(-1), Z: math.Inf(-1)},
	}
	for i := range r.Atoms {
		a := &r.Atoms[i]
		if !isHeavy(a.Element) {
			continue
		}
		c := a.Coord
		if c.X < b.min.X {
			b.min.X = c.X
		}
		if c.Y < b.min.Y {
			b.min.Y = c.Y
		}
		if c.Z < b.min.Z {
			b.min.Z = c.Z
		}
		if c.X > b.max.X {
			b.max.X = c.X
		}
		if c.Y > b.max.Y {
			b.max.Y = c.Y
		}
		if c.Z > b.max.Z {
			b.max.Z = c.Z
		}
	}
	return b
}

func axisGap(min1, max1, min2, max2 float64) float64 {
	if max1 < min2 {
		return min2 - max1
	}
	if max2 < min1 {
		return min1 - max2
	}
	return 0
}

// bboxDist2 returns the squared minimum distance between two boxes. Any pair
// of points inside the boxes is at least this far apart.
func bboxDist2(b1, b2 bbox) float64 {
	dx := axisGap(b1.min.X, b1.max.X, b2.min.X, b2.max.X)
	dy := axisGap(b1.min.Y, b1.max.Y, b2.min.Y, b2.max.Y)
	dz := axisGap(b1.min.Z, b1.max.Z, b2.min.Z, b2.max.Z)
	return dx*dx + dy*dy + dz*dz
}
