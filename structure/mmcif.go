package structure

import (
	"bytes"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"
)

// tokenizeMMCIF splits mmCIF content into whitespace tokens, honoring single
// and double quoted strings, semicolon-delimited text blocks and '#' comments.
// The tokens '.' and '?' are kept verbatim so callers can detect missing data.
func tokenizeMMCIF(data []byte) []string {
	var toks []string
	i := 0
	n := len(data)
	for i < n {
		c := data[i]
		switch {
		case c == '\n' || c == '\r' || c == ' ' || c == '\t':
			i++
		case c == '#':
			for i < n && data[i] != '\n' {
				i++
			}
		case c == ';':
			// Text block: starts at column 1 and ends at a line starting ';'.
			i++
			start := i
			for i < n {
				if data[i] == '\n' && i+1 < n && data[i+1] == ';' {
					break
				}
				i++
			}
			toks = append(toks, string(data[start:i]))
			for i < n && data[i] != '\n' {
				i++
			}
			if i < n {
				i++
			}
		case c == '\'' || c == '"':
			q := c
			i++
			start := i
			var b strings.Builder
			for i < n && data[i] != q {
				if data[i] == '\n' {
					b.WriteByte(' ')
				} else {
					b.WriteByte(data[i])
				}
				i++
			}
			_ = start
			if i < n {
				i++
			}
			toks = append(toks, b.String())
		default:
			start := i
			for i < n {
				ch := data[i]
				if ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' {
					break
				}
				i++
			}
			toks = append(toks, string(data[start:i]))
		}
	}
	return toks
}

func mmcifCol(headers map[string]int, names ...string) int {
	for _, nm := range names {
		if idx, ok := headers[nm]; ok {
			return idx
		}
	}
	return -1
}

// ParseMMCIF parses mmCIF data using only the _atom_site category.
func ParseMMCIF(r io.Reader, source string) (*Structure, error) {
	data, err := io.ReadAll(r)
	if err != nil {
		return nil, &ParseError{Source: source, Reason: "read error: " + err.Error()}
	}
	toks := tokenizeMMCIF(data)
	s := &Structure{Source: source}

	models := map[int]*modelBuilder{}
	var modelSerials []int
	atomCount := 0
	dropped := 0

	get := func(row []string, idx int) string {
		if idx < 0 || idx >= len(row) {
			return ""
		}
		v := row[idx]
		if v == "." || v == "?" {
			return ""
		}
		return strings.TrimSpace(v)
	}

	for p := 0; p < len(toks); {
		if toks[p] != "loop_" {
			p++
			continue
		}
		p++
		var headers []string
		for p < len(toks) && strings.HasPrefix(toks[p], "_") {
			headers = append(headers, strings.ToLower(toks[p]))
			p++
		}
		isAtomSite := false
		for _, h := range headers {
			if strings.HasPrefix(h, "_atom_site.") {
				isAtomSite = true
				break
			}
		}
		if !isAtomSite {
			// Skip this loop's data rows.
			for p < len(toks) && toks[p] != "loop_" && !strings.HasPrefix(toks[p], "_") {
				p++
			}
			continue
		}
		col := map[string]int{}
		for i, h := range headers {
			col[strings.TrimPrefix(h, "_atom_site.")] = i
		}
		cGroup := mmcifCol(col, "group_pdb")
		cAtom := mmcifCol(col, "auth_atom_id", "label_atom_id")
		cAlt := mmcifCol(col, "label_alt_id")
		cComp := mmcifCol(col, "label_comp_id", "auth_comp_id")
		cChain := mmcifCol(col, "auth_asym_id", "label_asym_id")
		cSeq := mmcifCol(col, "auth_seq_id", "label_seq_id")
		cIns := mmcifCol(col, "pdbx_pdb_ins_code")
		cX := mmcifCol(col, "cartn_x")
		cY := mmcifCol(col, "cartn_y")
		cZ := mmcifCol(col, "cartn_z")
		cOcc := mmcifCol(col, "occupancy")
		cEl := mmcifCol(col, "type_symbol")
		cModel := mmcifCol(col, "pdbx_pdb_model_num")

		width := len(headers)
		for p+width <= len(toks) {
			// Stop if the next tokens start a new construct.
			if toks[p] == "loop_" || strings.HasPrefix(toks[p], "_") || toks[p] == "data_" || strings.HasPrefix(toks[p], "data_") {
				break
			}
			row := toks[p : p+width]
			p += width
			atomCount++
			bad := func(msg string) {
				dropped++
				s.Warnings = append(s.Warnings, Issue{Kind: IssueMalformedRecord,
					Line:    0,
					Message: "atom_site row dropped: " + msg})
			}
			x, errx := strconv.ParseFloat(get(row, cX), 64)
			y, erry := strconv.ParseFloat(get(row, cY), 64)
			z, errz := strconv.ParseFloat(get(row, cZ), 64)
			if errx != nil || erry != nil || errz != nil {
				bad("unparsable coordinates")
				continue
			}
			seqStr := get(row, cSeq)
			if seqStr == "" {
				bad("missing residue sequence id")
				continue
			}
			seq, err := strconv.Atoi(seqStr)
			if err != nil {
				bad("unparsable residue sequence id")
				continue
			}
			name := get(row, cAtom)
			if name == "" {
				bad("missing atom id")
				continue
			}
			serial := 1
			if ms := get(row, cModel); ms != "" {
				if v, err := strconv.Atoi(ms); err == nil {
					serial = v
				}
			}
			mb, ok := models[serial]
			if !ok {
				mb = newModelBuilder(serial)
				models[serial] = mb
				modelSerials = append(modelSerials, serial)
			}
			occ := 0.0
			if os2 := get(row, cOcc); os2 != "" {
				occ, _ = strconv.ParseFloat(os2, 64)
			}
			var alt byte
			if a := get(row, cAlt); a != "" {
				alt = a[0]
			}
			var ic byte
			if v := get(row, cIns); v != "" {
				ic = v[0]
			}
			ra := rawAtom{
				hetero:  get(row, cGroup) == "HETATM",
				name:    name,
				altLoc:  alt,
				resName: get(row, cComp),
				chainID: get(row, cChain),
				resSeq:  seq,
				iCode:   ic,
				x:       x, y: y, z: z,
				occ:     occ,
				element: get(row, cEl),
			}
			mb.ingest(ra)
		}
	}
	if atomCount == 0 {
		return nil, &ParseError{Source: source, Reason: "no _atom_site rows found"}
	}
	sort.Ints(modelSerials)
	if len(modelSerials) > 1 {
		s.Warnings = append(s.Warnings, Issue{Kind: IssueMultiModel,
			Message: fmt.Sprintf("%d models present; only the first (serial %d) is used",
				len(modelSerials), modelSerials[0])})
	}
	first := models[modelSerials[0]]
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

// ParseMMCIFFile parses a local mmCIF file.
func ParseMMCIFFile(path string) (*Structure, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, &ParseError{Source: path, Reason: "cannot open file: " + err.Error()}
	}
	defer f.Close()
	data, err := io.ReadAll(f)
	if err != nil {
		return nil, &ParseError{Source: path, Reason: "read error: " + err.Error()}
	}
	return ParseMMCIF(bytes.NewReader(data), path)
}

// ParseFile auto-detects PDB vs mmCIF by extension and content.
func ParseFile(path string) (*Structure, error) {
	lower := strings.ToLower(path)
	if strings.HasSuffix(lower, ".cif") || strings.HasSuffix(lower, ".mmcif") {
		return ParseMMCIFFile(path)
	}
	return ParsePDBFile(path)
}
