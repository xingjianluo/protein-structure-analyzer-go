package structure

import (
	"bytes"
	"strings"
)

// ParseBytes parses structure content already held in memory, choosing the
// PDB or mmCIF parser the same way ParseFile does (by file extension).
func ParseBytes(data []byte, source string) (*Structure, error) {
	lower := strings.ToLower(source)
	if strings.HasSuffix(lower, ".cif") || strings.HasSuffix(lower, ".mmcif") {
		return ParseMMCIF(bytes.NewReader(data), source)
	}
	return ParsePDB(bytes.NewReader(data), source)
}
