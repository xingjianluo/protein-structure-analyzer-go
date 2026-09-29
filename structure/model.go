// Package structure 定义统一的链-残基-原子结构模型。
package structure

// Atom 表示一个原子记录。
type Atom struct {
	Serial    int
	Name      string  // 原子名，如 CA、CB
	AltLoc    string  // 可变构象标识，已按规则归一化
	ResName   string  // 残基名
	ChainID   string  // 链标识
	ResSeq    int     // 残基序号
	ICode     string  // 插入码
	X, Y, Z   float64 // 坐标（埃）
	Occupancy float64
	Element   string
	Het       bool // 是否来自 HETATM（非标准残基/配体）
}

// Residue 表示一个残基（由 链ID+序号+插入码 唯一标识）。
type Residue struct {
	ChainID     string
	Seq         int
	ICode       string
	Name        string
	Atoms       []*Atom
	NonStandard bool // 非标准残基（来自 HETATM）
	Incomplete  bool // 缺少关键骨架原子（CA/C 等）
}

// Chain 表示一条链。
type Chain struct {
	ID       string
	Residues []*Residue
}

// Model 表示一个结构模型（多模型文件取第一个模型，见解析规则）。
type Model struct {
	Chains []*Chain
}

// IssueKind 区分解析过程中各类可区分的处理结论。
type IssueKind string

const (
	IssueMalformedRecord IssueKind = "MALFORMED_RECORD" // 畸形记录，已跳过
	IssueTruncatedRecord IssueKind = "TRUNCATED_RECORD" // 截断记录，已跳过
	IssueMultiModel      IssueKind = "MULTI_MODEL"      // 多模型：仅采用第一个模型
	IssueAltLoc          IssueKind = "ALT_LOC"          // 可变构象：保留占用率最高者
	IssueInsertionCode   IssueKind = "INSERTION_CODE"   // 含插入码的残基
	IssueNonStandardRes  IssueKind = "NON_STANDARD_RES" // 非标准残基（HETATM）
	IssueMissingAtoms    IssueKind = "MISSING_ATOMS"    // 残基缺少关键原子
	IssueDuplicateAtom   IssueKind = "DUPLICATE_ATOM"   // 重复原子记录，保留首条
	IssueUnknownRecord   IssueKind = "UNKNOWN_RECORD"   // 未识别的记录类型，忽略
)

// Issue 描述一条解析期结论（不静默跳过任何异常输入）。
type Issue struct {
	Kind   IssueKind
	Line   int
	Detail string
}

// ParseReport 汇总一次解析的全部结论。
type ParseReport struct {
	Source string
	Issues []Issue
	Atoms  int
}

// HasKind 报告是否包含某类结论。
func (r *ParseReport) HasKind(k IssueKind) bool {
	for _, is := range r.Issues {
		if is.Kind == k {
			return true
		}
	}
	return false
}
