package structure

import (
	"math"
	"sort"
)

// ContactOptions 控制接触判定参数。
type ContactOptions struct {
	Cutoff          float64 // 原子距离阈值（埃），默认 4.5
	IncludeHetAtoms bool    // 是否让非标准残基参与界面判定
}

// DefaultContactOptions 返回默认参数。
func DefaultContactOptions() ContactOptions {
	return ContactOptions{Cutoff: 4.5}
}

func (o ContactOptions) cutoff() float64 {
	if o.Cutoff <= 0 {
		return 4.5
	}
	return o.Cutoff
}

// ResidueContact 表示一对残基间接触。
type ResidueContact struct {
	ChainA, ChainB string
	SeqA, SeqB     int
	ICodeA, ICodeB string
	NameA, NameB   string
	MinDistance    float64
	AtomPairs      int
}

// Interface 表示两条链之间的界面。
type Interface struct {
	ChainA, ChainB string
	Contacts       []ResidueContact
	BuriedAreaHint float64 // 接触对数 * 平均贡献，仅作排序提示
}

// residueLess 定义残基的规范顺序，保证结果与输入顺序无关。
func residueLess(chainA string, seqA int, icoA, chainB string, seqB int, icoB string) bool {
	if chainA != chainB {
		return chainA < chainB
	}
	if seqA != seqB {
		return seqA < seqB
	}
	return icoA < icoB
}

// Contacts 计算模型内所有链间残基接触。
//
// 判定规则：两残基分属不同链，且各自至少一对原子距离 <= Cutoff。
// 非标准残基与缺原子残基默认不参与（可用 IncludeHetAtoms 放开非标准残基），
// 缺原子残基始终排除，避免不完整数据污染界面结论。
// 结果按 (链A, 链B, 序号A, 插入码A, 序号B, 插入码B) 排序，保证可复现。
func Contacts(m *Model, opts ContactOptions) []ResidueContact {
	if m == nil {
		return nil
	}
	cut := opts.cutoff()
	cut2 := cut * cut

	type resRef struct {
		chain string
		res   *Residue
	}
	var eligible []resRef
	for _, ch := range m.Chains {
		for _, r := range ch.Residues {
			if r.Incomplete {
				continue
			}
			if r.NonStandard && !opts.IncludeHetAtoms {
				continue
			}
			eligible = append(eligible, resRef{ch.ID, r})
		}
	}

	var out []ResidueContact
	for i := 0; i < len(eligible); i++ {
		for j := i + 1; j < len(eligible); j++ {
			a, b := eligible[i], eligible[j]
			if a.chain == b.chain {
				continue
			}
			min2 := math.Inf(1)
			pairs := 0
			for _, aa := range a.res.Atoms {
				for _, bb := range b.res.Atoms {
					dx, dy, dz := aa.X-bb.X, aa.Y-bb.Y, aa.Z-bb.Z
					d2 := dx*dx + dy*dy + dz*dz
					if d2 <= cut2 {
						pairs++
						if d2 < min2 {
							min2 = d2
						}
					}
				}
			}
			if pairs == 0 {
				continue
			}
			c := ResidueContact{
				ChainA: a.chain, SeqA: a.res.Seq, ICodeA: a.res.ICode, NameA: a.res.Name,
				ChainB: b.chain, SeqB: b.res.Seq, ICodeB: b.res.ICode, NameB: b.res.Name,
				MinDistance: math.Sqrt(min2),
				AtomPairs:   pairs,
			}
			// 规范方向：较小的残基标识在前
			if residueLess(c.ChainB, c.SeqB, c.ICodeB, c.ChainA, c.SeqA, c.ICodeA) {
				c.ChainA, c.SeqA, c.ICodeA, c.NameA, c.ChainB, c.SeqB, c.ICodeB, c.NameB =
					c.ChainB, c.SeqB, c.ICodeB, c.NameB, c.ChainA, c.SeqA, c.ICodeA, c.NameA
			}
			out = append(out, c)
		}
	}
	sort.Slice(out, func(i, j int) bool {
		a, b := out[i], out[j]
		if residueLess(a.ChainA, a.SeqA, a.ICodeA, b.ChainA, b.SeqA, b.ICodeA) {
			return true
		}
		if residueLess(b.ChainA, b.SeqA, b.ICodeA, a.ChainA, a.SeqA, a.ICodeA) {
			return false
		}
		return residueLess(a.ChainB, a.SeqB, a.ICodeB, b.ChainB, b.SeqB, b.ICodeB)
	})
	return out
}

// Interfaces 由接触结果聚合链间界面，按 (链A, 链B) 排序，保证可复现。
func Interfaces(m *Model, opts ContactOptions) []Interface {
	contacts := Contacts(m, opts)
	type pair struct{ a, b string }
	byPair := map[pair][]ResidueContact{}
	var order []pair
	for _, c := range contacts {
		p := pair{c.ChainA, c.ChainB}
		if _, ok := byPair[p]; !ok {
			order = append(order, p)
		}
		byPair[p] = append(byPair[p], c)
	}
	sort.Slice(order, func(i, j int) bool {
		if order[i].a != order[j].a {
			return order[i].a < order[j].a
		}
		return order[i].b < order[j].b
	})
	out := make([]Interface, 0, len(order))
	for _, p := range order {
		cs := byPair[p]
		hint := 0.0
		for _, c := range cs {
			hint += float64(c.AtomPairs)
		}
		out = append(out, Interface{
			ChainA: p.a, ChainB: p.b,
			Contacts:       cs,
			BuriedAreaHint: hint,
		})
	}
	return out
}
