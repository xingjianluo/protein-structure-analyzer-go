package store

import "strings"

// FailureKind distinguishes why one input failed, so callers can react by
// type instead of matching message text.
type FailureKind int

const (
	FailureUnknown FailureKind = iota
	// FailureIO: the file could not be opened or read.
	FailureIO
	// FailureParse: the content could not be parsed into a structure
	// (including files whose coordinate columns are entirely unreadable).
	FailureParse
	// FailureAnalysis: the file parsed but interface analysis failed.
	FailureAnalysis
)

func (k FailureKind) String() string {
	switch k {
	case FailureIO:
		return "io"
	case FailureParse:
		return "parse"
	case FailureAnalysis:
		return "analysis"
	default:
		return "unknown"
	}
}

// EntryFailure is one typed per-file failure inside a ReloadError.
type EntryFailure struct {
	Path    string
	Kind    FailureKind
	Message string
}

// ReloadError reports a failed reload. It is distinguishable from normal
// operation by type and carries both the reason and the retained version.
type ReloadError struct {
	AttemptedPaths []string
	Version        int
	Retained       bool
	Reason         string
	// Failures lists every failed input with a typed kind.
	Failures []EntryFailure
}

func (e *ReloadError) Error() string {
	base := e.Reason
	if e.Retained {
		base += " (retained snapshot version " + itoa(e.Version) + ")"
	} else {
		base += " (no previous snapshot available)"
	}
	if len(e.Failures) > 0 {
		parts := make([]string, 0, len(e.Failures))
		for _, f := range e.Failures {
			parts = append(parts, f.Path+"["+f.Kind.String()+"]: "+f.Message)
		}
		base += ": " + strings.Join(parts, "; ")
	}
	return base
}

func itoa(n int) string {
	if n == 0 {
		return "0"
	}
	neg := n < 0
	if neg {
		n = -n
	}
	var b [20]byte
	i := len(b)
	for n > 0 {
		i--
		b[i] = byte('0' + n%10)
		n /= 10
	}
	if neg {
		i--
		b[i] = '-'
	}
	return string(b[i:])
}
