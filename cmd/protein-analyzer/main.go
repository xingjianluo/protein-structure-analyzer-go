package main

import (
	"flag"
	"fmt"
	"os"
	"strings"

	"github.com/xingjianluo/protein-structure-analyzer-go/contact"
	"github.com/xingjianluo/protein-structure-analyzer-go/store"
	"github.com/xingjianluo/protein-structure-analyzer-go/structure"
)

func main() {
	cutoff := flag.Float64("cutoff", contact.DefaultCutoff, "heavy-atom contact cutoff in Angstrom")
	includeHetero := flag.Bool("include-hetero", false, "include HETATM / non-standard residues (metal ions, modified residues) in interface determination")
	chains := flag.String("chains", "", "comma-separated chain IDs to include (default: all chains)")
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "usage: %s [flags] <structure-file> [more files...]\n", os.Args[0])
		fmt.Fprintf(os.Stderr, "  Parses local PDB/mmCIF files, reports inter-chain interfaces per file\n")
		fmt.Fprintf(os.Stderr, "  and a cross-structure comparison of the whole batch.\n")
		flag.PrintDefaults()
	}
	flag.Parse()
	if flag.NArg() == 0 {
		flag.Usage()
		os.Exit(2)
	}

	opts := contact.Options{Cutoff: *cutoff, IncludeHetero: *includeHetero}
	if *chains != "" {
		for _, c := range strings.Split(*chains, ",") {
			if c = strings.TrimSpace(c); c != "" {
				opts.Chains = append(opts.Chains, c)
			}
		}
	}

	analyzer := contact.NewAnalyzer(*cutoff)
	st := store.New(analyzer)
	paths := flag.Args()
	if err := st.LoadWithOptions(paths, opts); err != nil {
		fmt.Fprintf(os.Stderr, "load failed: %v\n", err)
		os.Exit(1)
	}
	snap := st.Current()
	fmt.Printf("snapshot version %d, state=%s, cutoff=%.2f A, include-hetero=%v, chains=%s\n",
		snap.Version, snap.State, *cutoff, opts.IncludeHetero, chainLabel(opts.Chains))
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

	cmp := snap.Comparison
	fmt.Printf("\n== batch comparison (fingerprint %s, %d structures) ==\n", cmp.Fingerprint, len(cmp.Sources))
	fmt.Printf("  common interfaces (present in all %d):\n", len(cmp.Sources))
	for _, it := range cmp.Common {
		fmt.Printf("    %s--%s: contacts per structure %v\n", it.ChainA, it.ChainB, it.ContactCounts)
	}
	fmt.Println("  partial interfaces (present in some):")
	for _, it := range cmp.Partial {
		fmt.Printf("    %s--%s: present in %v\n", it.ChainA, it.ChainB, it.PresentIn)
	}
	fmt.Println("  changed contacts:")
	for _, ch := range cmp.Changes {
		fmt.Printf("    %s/%s%d - %s/%s%d: present in %v, distance spread %.2f A\n",
			ch.ChainA, ch.NameA, ch.SeqA, ch.ChainB, ch.NameB, ch.SeqB, ch.PresentIn, ch.MaxSpread)
	}
}

func chainLabel(chains []string) string {
	if len(chains) == 0 {
		return "all"
	}
	return strings.Join(chains, ",")
}

func printable(b byte) byte {
	if b == 0 {
		return ' '
	}
	return b
}
