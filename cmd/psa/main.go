// psa 是本地蛋白质结构界面判定命令行工具。
//
// 用法：
//
//	psa [-cutoff 4.5] [-include-het] file1.pdb [file2.pdb ...]
//
// 对每个输入文件打印解析结论（各类 Issue 计数）与链间界面判定结果。
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/xingjianluo/protein-structure-analyzer-go/structure"
)

func main() {
	cutoff := flag.Float64("cutoff", 4.5, "原子接触距离阈值（埃）")
	includeHet := flag.Bool("include-het", false, "让非标准残基参与界面判定")
	verbose := flag.Bool("v", false, "打印每条解析结论明细")
	flag.Parse()
	if flag.NArg() == 0 {
		fmt.Fprintln(os.Stderr, "usage: psa [-cutoff A] [-include-het] [-v] file.pdb ...")
		os.Exit(2)
	}

	opts := structure.ContactOptions{Cutoff: *cutoff, IncludeHetAtoms: *includeHet}
	exit := 0
	for _, path := range flag.Args() {
		if err := analyze(path, opts, *verbose); err != nil {
			fmt.Fprintf(os.Stderr, "ERROR %s: %v\n", path, err)
			exit = 1
		}
	}
	os.Exit(exit)
}

func analyze(path string, opts structure.ContactOptions, verbose bool) error {
	f, err := os.Open(path)
	if err != nil {
		return err
	}
	defer f.Close()

	model, rep, err := structure.ParsePDB(path, f)
	if err != nil {
		return err
	}

	fmt.Printf("== %s ==\n", path)
	fmt.Printf("atoms: %d, chains: %d\n", rep.Atoms, len(model.Chains))

	counts := map[structure.IssueKind]int{}
	for _, is := range rep.Issues {
		counts[is.Kind]++
	}
	if len(counts) == 0 {
		fmt.Println("parse issues: none")
	} else {
		fmt.Println("parse issues:")
		for _, k := range []structure.IssueKind{
			structure.IssueMultiModel, structure.IssueAltLoc, structure.IssueInsertionCode,
			structure.IssueNonStandardRes, structure.IssueMissingAtoms,
			structure.IssueDuplicateAtom, structure.IssueMalformedRecord,
			structure.IssueTruncatedRecord, structure.IssueUnknownRecord,
		} {
			if n := counts[k]; n > 0 {
				fmt.Printf("  %-18s %d\n", k, n)
			}
		}
	}
	if verbose {
		for _, is := range rep.Issues {
			fmt.Printf("    [%s] line %d: %s\n", is.Kind, is.Line, is.Detail)
		}
	}

	ifaces := structure.Interfaces(model, opts)
	if len(ifaces) == 0 {
		fmt.Println("interfaces: none")
		return nil
	}
	fmt.Println("interfaces:")
	for _, ifa := range ifaces {
		fmt.Printf("  %s <-> %s : %d residue contacts\n", ifa.ChainA, ifa.ChainB, len(ifa.Contacts))
		for _, c := range ifa.Contacts {
			fmt.Printf("    %s:%d%s(%s) - %s:%d%s(%s)  min=%.2fA pairs=%d\n",
				c.ChainA, c.SeqA, c.ICodeA, c.NameA,
				c.ChainB, c.SeqB, c.ICodeB, c.NameB,
				c.MinDistance, c.AtomPairs)
		}
	}
	return nil
}
