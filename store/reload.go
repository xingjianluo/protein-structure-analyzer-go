package store

import "strings"

// ReloadError reports a failed reload. It is distinguishable from normal
// operation by type and carries both the reason and the retained version.
type ReloadError struct {
	AttemptedPaths []string
	Version        int
	Retained       bool
	Reason         string
	Failures       []string
}

func (e *ReloadError) Error() string {
	base := e.Reason
	if e.Retained {
		base += " (retained snapshot version " + itoa(e.Version) + ")"
	} else {
		base += " (no previous snapshot available)"
	}
	if len(e.Failures) > 0 {
		base += ": " + strings.Join(e.Failures, "; ")
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
