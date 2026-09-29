package structure

// ContactOptions 控制接触判定参数。
type ContactOptions struct {
	Cutoff          float64 // 原子距离阈值（埃），默认 4.5
	IncludeHetAtoms bool    // 是否让非标准残基参与界面判定
}

// DefaultContactOptions 返回默认参数。
func DefaultContactOptions() ContactOptions {
	return ContactOptions{Cutoff: 4.5}
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

// Contacts 计算模型内所有链间残基接触，结果按规范顺序排序，保证可复现。
func Contacts(m *Model, opts ContactOptions) []ResidueContact {
	return nil
}

// Interfaces 由接触结果聚合链间界面，排序确定。
func Interfaces(m *Model, opts ContactOptions) []Interface {
	return nil
}
