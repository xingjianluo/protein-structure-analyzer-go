package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/xingjianluo/protein-structure-analyzer-go/contact"
	"github.com/xingjianluo/protein-structure-analyzer-go/store"
	"github.com/xingjianluo/protein-structure-analyzer-go/structure"
)

func main() {
	cutoff := flag.Float64("cutoff", contact.DefaultCutoff, "heavy-atom contact cutoff in Angstrom")
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "usage: %s [--cutoff angstrom] <structure-file> [more files...]\n", os.Args[0])
		fmt.Fprintf(os.Stderr, "  Parses local PDB/mmCIF files and reports inter-chain residue interfaces.\n")
		flag.PrintDefaults()
	}
	flag.Parse()
	if flag.NArg() == 0 {
		flag.Usage()
		os.Exit(2)
	}

	analyzer := contact.NewAnalyzer(*cutoff)
	st := store.New(analyzer)
	paths := flag.Args()
	if err := st.Load(paths); err != nil {
		fmt.Fprintf(os.Stderr, "load failed: %v\n", err)
		os.Exit(1)
	}
	snap := st.Current()
	fmt.Printf("snapshot version %d, state=%s, cutoff=%.2f A\n", snap.Version, snap.State, *cutoff)
	for _, p := range snap.Paths {
		e := snap.Entries[p]
		s, err := structure.ParseFile(p)
		if err != nil {
			continue
		}
		fmt.Printf("\n== %s (state=%s, model=%d, heavy atoms=%d, usable residues=%d) ==\n",
			p, e.State, s.Effective.Serial, s.HeavyAtomCount(), s.ContactUsableResidueCount())
		for _, w := range e.Warnings {
			fmt.Printf("  [%s] %s\n", w.Kind, w.String())
		}
		rep := snap.Reports[p]
		if len(rep.Interfaces) == 0 {
			fmt.Println("  no inter-chain interfaces")
			continue
		}
		for _, it := range rep.Interfaces {
			fmt.Printf("  interface %s--%s: %d residue contact(s), min distance %.2f A\n",
				it.ChainA, it.ChainB, it.ContactCount, it.MinDist)
		}
		for _, c := range rep.Contacts {
			fmt.Printf("    contact %s/%s%d%c - %s/%s%d%c  d=%.2f A\n",
				c.ChainA, c.NameA, c.SeqA, printable(c.ICodeA),
				c.ChainB, c.NameB, c.SeqB, printable(c.ICodeB), c.MinDist)
		}
	}
}

func printable(b byte) byte {
	if b == 0 {
		return ' '
	}
	return b
}
