package structure

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"
)

// rawAtom is an unprocessed ATOM/HETATM record.
type rawAtom struct {
	line    int
	hetero  bool
	name    string
	altLoc  byte
	resName string
	chainID string
	resSeq  int
	iCode   byte
	x, y, z float64
	occ     float64
	element string
}

type resKey struct {
	seq   int
	iCode byte
}

type resBuilder struct {
	key   resKey
	name  string
	alt   []byte // alt loc labels encountered, in order
	atoms []rawAtom
}

type chainBuilder struct {
	id       string
	order    []resKey
	residues map[resKey]*resBuilder
}

type modelBuilder struct {
	serial int
	chains []*chainBuilder // insertion order
	index  map[string]*chainBuilder
}

func newModelBuilder(serial int) *modelBuilder {
	return &modelBuilder{serial: serial, index: map[string]*chainBuilder{}}
}

func (m *modelBuilder) chain(id string) *chainBuilder {
	if c, ok := m.index[id]; ok {
		return c
	}
	c := &chainBuilder{id: id, residues: map[resKey]*resBuilder{}}
	m.index[id] = c
	m.chains = append(m.chains, c)
	return c
}

func (c *chainBuilder) residue(k resKey, name string) *resBuilder {
	if r, ok := c.residues[k]; ok {
		return r
	}
	r := &resBuilder{key: k, name: name}
	c.residues[k] = r
	c.order = append(c.order, k)
	return r
}

// fixed-width helpers: PDB columns are 1-based; [lo,hi] inclusive.
func col(line string, lo, hi int) string {
	if len(line) < lo {
		return ""
	}
	if hi > len(line) {
		hi = len(line)
	}
	return line[lo-1 : hi]
}

func parseFloatField(s string) (float64, bool) {
	s = strings.TrimSpace(s)
	if s == "" {
		return 0, false
	}
	v, err := strconv.ParseFloat(s, 64)
	return v, err == nil
}

func parseIntField(s string) (int, bool) {
	s = strings.TrimSpace(s)
	if s == "" {
		return 0, false
	}
	v, err := strconv.Atoi(s)
	return v, err == nil
}

// guessElement derives the element from an atom name when columns 77-78 are
// absent. PDB atom names place the element right-justified in columns 13-14.
func guessElement(name string) string {
	name = strings.TrimSpace(name)
	for i := 0; i < len(name); i++ {
		ch := name[i]
		if (ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z') {
			el := string(ch)
			// Two-letter elements start with the digit-stripped first chars,
			// e.g. "CL" in "CL12" / "1CL?"; PDB right-justifies, so checking
			// the leading alpha run is sufficient for common cases.
			return strings.ToUpper(el)
		}
	}
	return ""
}

// parseAtomLine parses one fixed-width ATOM/HETATM line.
func parseAtomLine(line string, lineNo int, hetero bool) (rawAtom, Issue, bool) {
	ra := rawAtom{line: lineNo, hetero: hetero}
	// Coordinates end at column 54; shorter lines are truncated.
	if len(strings.TrimRight(line, "\r")) < 54 {
		return ra, Issue{Kind: IssueTruncatedRecord, Line: lineNo,
			Message: fmt.Sprintf("record has %d columns, want at least 54", len(strings.TrimRight(line, "\r")))}, false
	}
	ra.name = strings.TrimSpace(col(line, 13, 16))
	if a := strings.TrimSpace(col(line, 17, 17)); a != "" {
		ra.altLoc = a[0]
	}
	ra.resName = strings.TrimSpace(col(line, 18, 20))
	ra.chainID = col(line, 22, 22)
	seq, ok := parseIntField(col(line, 23, 26))
	if !ok {
		return ra, Issue{Kind: IssueMalformedRecord, Line: lineNo,
			Message: "unparsable residue sequence number"}, false
	}
	ra.resSeq = seq
	if ic := col(line, 27, 27); ic != " " && ic != "" {
		ra.iCode = ic[0]
	}
	var good bool
	if ra.x, good = parseFloatField(col(line, 31, 38)); !good {
		return ra, Issue{Kind: IssueMalformedRecord, Line: lineNo,
			Message: "unparsable X coordinate"}, false
	}
	if ra.y, good = parseFloatField(col(line, 39, 46)); !good {
		return ra, Issue{Kind: IssueMalformedRecord, Line: lineNo,
			Message: "unparsable Y coordinate"}, false
	}
	if ra.z, good = parseFloatField(col(line, 47, 54)); !good {
		return ra, Issue{Kind: IssueMalformedRecord, Line: lineNo,
			Message: "unparsable Z coordinate"}, false
	}
	if len(line) >= 60 {
		if occ, ok := parseFloatField(col(line, 55, 60)); ok {
			ra.occ = occ
		}
	}
	if len(line) >= 78 {
		ra.element = strings.TrimSpace(col(line, 77, 78))
	}
	if ra.element == "" {
		ra.element = guessElement(ra.name)
	}
	if ra.name == "" {
		return ra, Issue{Kind: IssueMalformedRecord, Line: lineNo,
			Message: "empty atom name"}, false
	}
	return ra, Issue{}, true
}

func (b *modelBuilder) ingest(ra rawAtom) {
	c := b.chain(ra.chainID)
	k := resKey{seq: ra.resSeq, iCode: ra.iCode}
	r := c.residue(k, ra.resName)
	r.atoms = append(r.atoms, ra)
	if ra.altLoc != ' ' && ra.altLoc != 0 {
		found := false
		for _, a := range r.alt {
			if a == ra.altLoc {
				found = true
				break
			}
		}
		if !found {
			r.alt = append(r.alt, ra.altLoc)
		}
	}
}

// chooseAltLoc deterministically selects one alternate conformation: highest
// total occupancy, ties broken by smallest ASCII label.
func chooseAltLoc(r *resBuilder) byte {
	if len(r.alt) == 0 {
		return ' '
	}
	labels := append([]byte{}, r.alt...)
	sort.Slice(labels, func(i, j int) bool { return labels[i] < labels[j] })
	sums := map[byte]float64{}
	for _, a := range r.atoms {
		if a.altLoc != ' ' && a.altLoc != 0 {
			sums[a.altLoc] += a.occ
		}
	}
	best := labels[0]
	bestSum := sums[best]
	for _, l := range labels[1:] {
		if sums[l] > bestSum {
			best, bestSum = l, sums[l]
		}
	}
	return best
}

func whereLabel(chainID string, seq int, iCode byte, name string) string {
	id := chainID
	if id == "" {
		id = " "
	}
	ic := ""
	if iCode != 0 {
		ic = string(iCode)
	}
	return fmt.Sprintf(" [%s/%s%d%s]", id, name, seq, ic)
}

// finalize converts the first model builder into a Model, applying alternate
// conformation selection and missing-atom/non-standard classification.
func finalize(b *modelBuilder, warns *[]Issue) *Model {
	m := &Model{Serial: b.serial}
	for _, cb := range b.chains {
		ch := &Chain{ID: cb.id}
		for _, k := range cb.order {
			rb := cb.residues[k]
			res := &Residue{Seq: k.seq, ICode: k.iCode, Name: rb.name}
			chosen := chooseAltLoc(rb)
			if chosen != ' ' {
				res.AltID = chosen
				*warns = append(*warns, Issue{Kind: IssueAlternateConformation,
					Where:   whereLabel(cb.id, k.seq, k.iCode, rb.name),
					Message: fmt.Sprintf("selected alternate location '%c'", chosen)})
			}
			seen := map[string]bool{}
			for _, ra := range rb.atoms {
				keep := ra.altLoc == ' ' || ra.altLoc == 0 || ra.altLoc == chosen
				if !keep || seen[ra.name] {
					continue
				}
				seen[ra.name] = true
				res.Atoms = append(res.Atoms, Atom{
					Name: ra.name, Element: ra.element, AltLoc: ra.altLoc,
					Coord: Vec3{ra.x, ra.y, ra.z}, Occupancy: ra.occ,
				})
			}
			firstHetero := false
			for _, ra := range rb.atoms {
				firstHetero = ra.hetero
				break
			}
			res.Hetero = firstHetero
			res.IsStandardAA = !firstHetero && IsStandardAminoAcid(rb.name)
			if res.IsStandardAA {
				present := map[string]bool{}
				for _, a := range res.Atoms {
					present[a.Name] = true
				}
				for _, want := range ExpectedHeavyAtoms(rb.name) {
					if !present[want] {
						res.Missing = append(res.Missing, want)
					}
				}
				if len(res.Missing) > 0 {
					*warns = append(*warns, Issue{Kind: IssueMissingAtom, Line: rb.atoms[0].line,
						Where:   whereLabel(cb.id, k.seq, k.iCode, rb.name),
						Message: "missing expected heavy atoms: " + strings.Join(res.Missing, ",")})
				}
			} else {
				*warns = append(*warns, Issue{Kind: IssueNonStandardResidue, Line: rb.atoms[0].line,
					Where:   whereLabel(cb.id, k.seq, k.iCode, rb.name),
					Message: "non-standard or hetero residue kept but flagged"})
			}
			if k.iCode != 0 {
				*warns = append(*warns, Issue{Kind: IssueInsertionCode, Line: rb.atoms[0].line,
					Where:   whereLabel(cb.id, k.seq, k.iCode, rb.name),
					Message: "insertion code honored as distinct residue"})
			}
			sort.Strings(res.Missing)
			ch.Residues = append(ch.Residues, res)
		}
		m.Chains = append(m.Chains, ch)
	}
	return m
}

// ParsePDB parses PDB-format data from r into the unified model.
func ParsePDB(r io.Reader, source string) (*Structure, error) {
	s := &Structure{Source: source}
	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 0, 64*1024), 1024*1024)
	var first *modelBuilder
	current := newModelBuilder(1)
	lineNo := 0
	atomCount := 0
	seenModel := false
	skipModel := false
	multiWarned := false
	finish := func() {
		if first == nil && len(current.chains) > 0 {
			first = current
		}
	}
	for scanner.Scan() {
		lineNo++
		raw := strings.TrimRight(scanner.Text(), "\r")
		if len(raw) < 6 {
			continue
		}
		rec := strings.TrimSpace(raw[:6])
		switch rec {
		case "MODEL":
			serial, ok := parseIntField(col(raw, 11, 14))
			if !ok {
				serial = lineNo
			}
			if seenModel {
				if !multiWarned {
					s.Warnings = append(s.Warnings, Issue{Kind: IssueMultiModel, Line: lineNo,
						Message: "multiple MODEL blocks present; only the first is used"})
					multiWarned = true
				}
				skipModel = true
			} else {
				current = newModelBuilder(serial)
			}
			seenModel = true
		case "ENDMDL":
			if !skipModel {
				finish()
			}
		case "ATOM", "HETATM":
			if skipModel {
				continue
			}
			atomCount++
			ra, iss, ok := parseAtomLine(raw, lineNo, rec == "HETATM")
			if !ok {
				s.Warnings = append(s.Warnings, iss)
				continue
			}
			current.ingest(ra)
		}
	}
	if err := scanner.Err(); err != nil {
		return nil, &ParseError{Source: source, Reason: "read error: " + err.Error(), Kind: ParseErrIO}
	}
	if !seenModel {
		finish()
	} else if first == nil && len(current.chains) > 0 {
		// First MODEL block without a trailing ENDMDL.
		first = current
	}
	if atomCount == 0 {
		return nil, &ParseError{Source: source, Reason: "no ATOM/HETATM records found", Kind: ParseErrContent}
	}
	if first == nil || len(first.chains) == 0 {
		return nil, &ParseError{Source: source, Reason: "no usable atomic records in first model", Kind: ParseErrContent}
	}
	s.Models = append(s.Models, finalize(first, &s.Warnings))
	s.Effective = s.Models[0]
	sort.SliceStable(s.Warnings, func(i, j int) bool {
		if s.Warnings[i].Line != s.Warnings[j].Line {
			return s.Warnings[i].Line < s.Warnings[j].Line
		}
		return s.Warnings[i].Kind < s.Warnings[j].Kind
	})
	return s, nil
}

// ParsePDBFile parses a local PDB file.
func ParsePDBFile(path string) (*Structure, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, &ParseError{Source: path, Reason: "cannot open file: " + err.Error(), Kind: ParseErrIO}
	}
	defer f.Close()
	data, err := io.ReadAll(f)
	if err != nil {
		return nil, &ParseError{Source: path, Reason: "read error: " + err.Error(), Kind: ParseErrIO}
	}
	return ParsePDB(bytes.NewReader(data), path)
}
