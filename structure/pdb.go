package structure

import "io"

// ParseError 表示解析层面的致命错误（文件整体不可用）。
type ParseError struct {
	Source string
	Reason string
}

func (e *ParseError) Error() string { return "parse " + e.Source + ": " + e.Reason }

// ParsePDB 从 reader 解析 PDB 格式结构文件，返回统一模型与解析报告。
// 畸形/截断记录记入报告并跳过；若完全无可用原子则返回 ParseError。
func ParsePDB(source string, r io.Reader) (*Model, *ParseReport, error) {
	return nil, nil, &ParseError{Source: source, Reason: "not implemented"}
}
