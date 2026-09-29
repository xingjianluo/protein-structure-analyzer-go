package structure

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"
)

// ParseError 表示解析层面的致命错误（文件整体不可用）。
type ParseError struct {
	Source string
	Reason string
}

func (e *ParseError) Error() string { return "parse " + e.Source + ": " + e.Reason }

// 关键骨架原子，用于判定残基是否缺原子。
const backboneCA = "CA"

// 标准氨基酸残基名集合。
var standardResidues = map[string]bool{
	"ALA": true, "ARG": true, "ASN": true, "ASP": true, "CYS": true,
	"GLN": true, "GLU": true, "GLY": true, "HIS": true, "ILE": true,
	"LEU": true, "LYS": true, "MET": true, "PHE": true, "PRO": true,
	"SER": true, "THR": true, "TRP": true, "TYR": true, "VAL": true,
}

// col 安全地截取定长列；越界返回空串。
func col(line string, start, end int) string {
	if len(line) < start {
		return ""
	}
	if len(line) < end {
		return strings.TrimSpace(line[start-1:])
	}
	return strings.TrimSpace(line[start-1 : end])
}

// ParsePDB 从 reader 解析 PDB 格式结构文件，返回统一模型与解析报告。
//
// 处理规则（每条规则都会在报告中留下可区分的结论）：
//   - 多模型（MODEL/ENDMDL）：仅采用第一个模型，其余模型记录 IssueMultiModel；
//   - 可变构象（altLoc）：同一原子位置保留占用率最高者（并列取 altLoc 字母序最小），
//     被丢弃的构象记录 IssueAltLoc；
//   - 插入码：并入残基唯一标识，记录 IssueInsertionCode；
//   - 非标准残基（HETATM 或非标准残基名）：保留于模型并标记 NonStandard，
//     记录 IssueNonStandardRes；
//   - 缺少 CA 骨架原子的残基：标记 Incomplete，记录 IssueMissingAtoms；
//   - 畸形/截断记录：记录 IssueMalformedRecord / IssueTruncatedRecord 并跳过；
//   - 重复原子记录：保留首条，记录 IssueDuplicateAtom。
//
// 若整个文件没有任何可用原子，返回 ParseError。
func ParsePDB(source string, r io.Reader) (*Model, *ParseReport, error) {
	rep := &ParseReport{Source: source}

	type resKey struct {
		chain string
		seq   int
		icode string
	}
	type atomKey struct {
		res  resKey
		name string
	}

	residues := map[resKey]*Residue{}
	var resOrder []resKey
	// altLoc 候选：key -> 已选原子（按规则可能整体替换）
	chosen := map[atomKey]*Atom{}
	seenAlt := map[resKey]bool{}
	seenICode := map[resKey]bool{}
	seenNonStd := map[resKey]bool{}

	inModel := false   // 是否处于 MODEL 块内
	modelSeen := false // 是否已遇到第一个模型
	skipModel := false // 当前处于被跳过的后续模型

	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 0, 1<<20), 1<<20)
	lineNo := 0
	for sc.Scan() {
		lineNo++
		line := strings.TrimRight(sc.Text(), "\r\n")
		rec := col(line, 1, 6)

		switch rec {
		case "MODEL":
			if modelSeen {
				skipModel = true
				rep.Issues = append(rep.Issues, Issue{
					Kind: IssueMultiModel, Line: lineNo,
					Detail: "additional MODEL discarded; only the first model is used",
				})
			} else {
				modelSeen = true
				inModel = true
			}
			continue
		case "ENDMDL":
			inModel = false
			skipModel = false
			continue
		case "ATOM", "HETATM":
			// 继续处理
		case "":
			if strings.TrimSpace(line) != "" {
				rep.Issues = append(rep.Issues, Issue{
					Kind: IssueUnknownRecord, Line: lineNo,
					Detail: "unrecognized record: " + truncate(line, 40),
				})
			}
			continue
		default:
			// HEADER/REMARK/TER/END 等已知但不参与结构计算的记录，静默忽略属正常行为
			continue
		}

		if skipModel {
			continue
		}
		if modelSeen && !inModel {
			// MODEL 块之外的坐标记录（畸形布局）
			rep.Issues = append(rep.Issues, Issue{
				Kind: IssueMalformedRecord, Line: lineNo,
				Detail: rec + " record outside MODEL block, skipped",
			})
			continue
		}

		// 截断判定：坐标列（54 列）不完整即视为截断
		if len(line) < 54 {
			rep.Issues = append(rep.Issues, Issue{
				Kind: IssueTruncatedRecord, Line: lineNo,
				Detail: rec + " record shorter than 54 columns, skipped",
			})
			continue
		}

		serial, err1 := strconv.Atoi(col(line, 7, 11))
		resSeq, err2 := strconv.Atoi(col(line, 23, 26))
		x, err3 := strconv.ParseFloat(col(line, 31, 38), 64)
		y, err4 := strconv.ParseFloat(col(line, 39, 46), 64)
		z, err5 := strconv.ParseFloat(col(line, 47, 54), 64)
		occ := 1.0
		if s := col(line, 55, 60); s != "" {
			if v, err := strconv.ParseFloat(s, 64); err == nil {
				occ = v
			}
		}
		if err1 != nil || err2 != nil || err3 != nil || err4 != nil || err5 != nil {
			rep.Issues = append(rep.Issues, Issue{
				Kind: IssueMalformedRecord, Line: lineNo,
				Detail: rec + " record has unparseable numeric fields, skipped",
			})
			continue
		}

		at := &Atom{
			Serial:    serial,
			Name:      col(line, 13, 16),
			AltLoc:    col(line, 17, 17),
			ResName:   col(line, 18, 20),
			ChainID:   col(line, 22, 22),
			ResSeq:    resSeq,
			ICode:     col(line, 27, 27),
			X:         x, Y: y, Z: z,
			Occupancy: occ,
			Element:   col(line, 77, 78),
			Het:       rec == "HETATM",
		}
		if at.Name == "" || at.ResName == "" {
			rep.Issues = append(rep.Issues, Issue{
				Kind: IssueMalformedRecord, Line: lineNo,
				Detail: rec + " record missing atom or residue name, skipped",
			})
			continue
		}

		rk := resKey{at.ChainID, at.ResSeq, at.ICode}
		res, ok := residues[rk]
		if !ok {
			res = &Residue{
				ChainID: at.ChainID, Seq: at.ResSeq, ICode: at.ICode,
				Name: at.ResName, NonStandard: at.Het || !standardResidues[at.ResName],
			}
			residues[rk] = res
			resOrder = append(resOrder, rk)
			if at.ICode != "" && !seenICode[rk] {
				seenICode[rk] = true
				rep.Issues = append(rep.Issues, Issue{
					Kind: IssueInsertionCode, Line: lineNo,
					Detail: fmt.Sprintf("residue %s:%d%s uses insertion code", at.ChainID, at.ResSeq, at.ICode),
				})
			}
			if res.NonStandard && !seenNonStd[rk] {
				seenNonStd[rk] = true
				rep.Issues = append(rep.Issues, Issue{
					Kind: IssueNonStandardRes, Line: lineNo,
					Detail: fmt.Sprintf("non-standard residue %s %s:%d kept but flagged", at.ResName, at.ChainID, at.ResSeq),
				})
			}
		}

		ak := atomKey{rk, at.Name}
		if prev, dup := chosen[ak]; dup {
			// 可变构象：保留占用率最高者，并列取 altLoc 字母序最小（空串视为最小）
			if betterAltLoc(at, prev) {
				rep.Issues = append(rep.Issues, Issue{
					Kind: IssueAltLoc, Line: lineNo,
					Detail: fmt.Sprintf("altLoc %q of %s %s:%d replaced by higher-occupancy conformer",
						prev.AltLoc, at.Name, at.ChainID, at.ResSeq),
				})
				replaceAtom(res, prev, at)
				chosen[ak] = at
			} else {
				kind := IssueAltLoc
				if prev.AltLoc == "" && at.AltLoc == "" {
					kind = IssueDuplicateAtom
				}
				rep.Issues = append(rep.Issues, Issue{
					Kind: kind, Line: lineNo,
					Detail: fmt.Sprintf("duplicate/altLoc %q of %s %s:%d discarded",
						at.AltLoc, at.Name, at.ChainID, at.ResSeq),
				})
			}
			continue
		}
		if at.AltLoc != "" && !seenAlt[rk] {
			seenAlt[rk] = true
		}
		chosen[ak] = at
		res.Atoms = append(res.Atoms, at)
	}
	if err := sc.Err(); err != nil {
		return nil, nil, &ParseError{Source: source, Reason: "read failed: " + err.Error()}
	}

	// 缺原子判定：无 CA 的氨基酸残基标记为 Incomplete
	for _, rk := range resOrder {
		res := residues[rk]
		if res.NonStandard {
			continue
		}
		hasCA := false
		for _, a := range res.Atoms {
			if a.Name == backboneCA {
				hasCA = true
				break
			}
		}
		if !hasCA {
			res.Incomplete = true
			rep.Issues = append(rep.Issues, Issue{
				Kind: IssueMissingAtoms,
				Detail: fmt.Sprintf("residue %s %s:%d%s lacks backbone CA atom",
					res.Name, res.ChainID, res.Seq, res.ICode),
			})
		}
	}

	// 规范排序：链按 ID，残基按 (序号, 插入码)，原子按序号 —— 保证可复现
	sort.Slice(resOrder, func(i, j int) bool {
		a, b := resOrder[i], resOrder[j]
		if a.chain != b.chain {
			return a.chain < b.chain
		}
		if a.seq != b.seq {
			return a.seq < b.seq
		}
		return a.icode < b.icode
	})
	chainIdx := map[string]*Chain{}
	var chainIDs []string
	m := &Model{}
	for _, rk := range resOrder {
		res := residues[rk]
		sort.Slice(res.Atoms, func(i, j int) bool { return res.Atoms[i].Serial < res.Atoms[j].Serial })
		ch, ok := chainIdx[rk.chain]
		if !ok {
			ch = &Chain{ID: rk.chain}
			chainIdx[rk.chain] = ch
			chainIDs = append(chainIDs, rk.chain)
			m.Chains = append(m.Chains, ch)
		}
		ch.Residues = append(ch.Residues, res)
		rep.Atoms += len(res.Atoms)
	}
	_ = chainIDs

	if rep.Atoms == 0 {
		return nil, rep, &ParseError{Source: source, Reason: "no usable atom records found"}
	}
	return m, rep, nil
}

// betterAltLoc 判定候选原子是否应替换已选原子（占用率更高，或并列时 altLoc 更小）。
func betterAltLoc(cand, cur *Atom) bool {
	if cand.Occupancy != cur.Occupancy {
		return cand.Occupancy > cur.Occupancy
	}
	return cand.AltLoc != "" && (cur.AltLoc == "" || cand.AltLoc < cur.AltLoc)
}

func replaceAtom(res *Residue, old, nu *Atom) {
	for i, a := range res.Atoms {
		if a == old {
			res.Atoms[i] = nu
			return
		}
	}
}

func truncate(s string, n int) string {
	if len(s) <= n {
		return s
	}
	return s[:n] + "..."
}
