package main

import (
	"flag"
	"fmt"
	"os"
	"strings"

	"github.com/xingjianluo/protein-structure-analyzer-go/compare"
	"github.com/xingjianluo/protein-structure-analyzer-go/contact"
	"github.com/xingjianluo/protein-structure-analyzer-go/store"
	"github.com/xingjianluo/protein-structure-analyzer-go/structure"
)

func main() {
	cutoff := flag.Float64("cutoff", contact.DefaultCutoff, "heavy-atom contact cutoff in Angstrom")
	includeHetero := flag.Bool("include-hetero", false, "include non-standard records (HETATM: metals, modified residues) in interface analysis")
	chains := flag.String("chains", "", "comma-separated chain IDs to compare (default: all chains)")
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "usage: %s [flags] <structure-file> [more files...]\n", os.Args[0])
		fmt.Fprintf(os.Stderr, "  Parses local PDB/mmCIF files, reports inter-chain interfaces,\n")
		fmt.Fprintf(os.Stderr, "  and compares interfaces across the whole batch.\n")
		flag.PrintDefaults()
	}
	flag.Parse()
	if flag.NArg() == 0 {
		flag.Usage()
		os.Exit(2)
	}

	opts := compare.Options{IncludeHetero: *includeHetero}
	if *chains != "" {
		for _, c := range strings.Split(*chains, ",") {
			if c = strings.TrimSpace(c); c != "" {
				opts.Chains = append(opts.Chains, c)
			}
		}
	}

	analyzer := contact.NewAnalyzer(*cutoff)
	st := store.New(analyzer)
	st.SetOptions(opts)
	paths := flag.Args()
	if err := st.Load(paths); err != nil {
		fmt.Fprintf(os.Stderr, "load failed: %v\n", err)
		os.Exit(1)
	}
	snap := st.Current()
	fmt.Printf("snapshot version %d, state=%s, cutoff=%.2f A, fingerprint=%s\n",
		snap.Version, snap.State, *cutoff, snap.Fingerprint[:12])
	fmt.Printf("options: include-hetero=%v chains=%v\n",
		snap.Options.IncludeHetero, chainList(snap.Options.Chains))
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

	comp := snap.Comparison
	fmt.Printf("\n== batch comparison (%d of %d files compared) ==\n", comp.Compared, len(comp.Files))
	for _, f := range comp.Files {
		if f.Class != compare.ClassCompared {
			fmt.Printf("  file %s: %s (%d issue(s), excluded from comparison)\n", f.ID, f.Class, len(f.Reasons))
		}
	}
	if len(comp.Interfaces) == 0 {
		fmt.Println("  no interfaces to compare")
		return
	}
	for _, it := range comp.Interfaces {
		kind := "shared"
		if !it.Shared {
			kind = "partial"
		}
		fmt.Printf("  interface %s--%s [%s]: present in %v", it.ChainA, it.ChainB, kind, it.PresentIn)
		if len(it.AbsentFrom) > 0 {
			fmt.Printf(", absent from %v", it.AbsentFrom)
		}
		fmt.Println()
	}
	for _, ch := range comp.Changes {
		if !ch.Changed {
			continue
		}
		fmt.Printf("  changed contacts %s--%s:\n", ch.ChainA, ch.ChainB)
		for _, st := range ch.Stats {
			if st.Present {
				fmt.Printf("    %s: %d contact(s), min %.2f A\n", st.ID, st.Contacts, st.MinDist)
			} else {
				fmt.Printf("    %s: no interface\n", st.ID)
			}
		}
	}
}

func chainList(cs []string) string {
	if len(cs) == 0 {
		return "(all)"
	}
	return strings.Join(cs, ",")
}

func printable(b byte) byte {
	if b == 0 {
		return ' '
	}
	return b
}
