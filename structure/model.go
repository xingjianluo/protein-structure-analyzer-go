package structure

import "fmt"

// Vec3 is a Cartesian coordinate in Angstrom.
type Vec3 struct {
	X, Y, Z float64
}

// Atom is a single atom record within a residue.
type Atom struct {
	Name      string
	Element   string
	AltLoc    byte
	Coord     Vec3
	Occupancy float64
}

// Residue is a residue (standard or hetero) inside a chain.
type Residue struct {
	Seq    int
	ICode  byte
	Name   string
	Hetero bool
	// AltID is the selected alternate location identifier (' ' when none).
	AltID byte
	Atoms []Atom
	// Missing lists expected heavy-atom names that were absent.
	Missing []string
	// IsStandardAA is true for the 20 standard amino acids on ATOM records.
	IsStandardAA bool
}

// UsableForContact reports whether the residue is allowed to participate in
// contact calculation. Residues with missing expected atoms are excluded so
// that incomplete data is never treated as a normal structure.
func (r *Residue) UsableForContact() bool {
	return len(r.Missing) == 0 && len(r.Atoms) > 0
}

// Chain is a collection of residues sharing a chain ID.
type Chain struct {
	ID       string
	Residues []*Residue
}

// Model is one MODEL/ENDMDL block (model 1 when no MODEL records exist).
type Model struct {
	Serial int
	Chains []*Chain
}

// Structure is the unified chain-residue-atom model.
type Structure struct {
	Source   string
	Models   []*Model
	Warnings []Issue
	// Effective is the deterministically selected model used for calculation.
	Effective *Model
}

// IssueKind enumerates the distinguishable problems found while parsing.
type IssueKind int

const (
	IssueUnknown IssueKind = iota
	// IssueMultiModel: several MODEL blocks were present; only the first is used.
	IssueMultiModel
	// IssueAlternateConformation: altLoc records were collapsed to one choice.
	IssueAlternateConformation
	// IssueInsertionCode: an insertion code was present and honored.
	IssueInsertionCode
	// IssueNonStandardResidue: a hetero / non-standard residue was kept but flagged.
	IssueNonStandardResidue
	// IssueMissingAtom: an expected heavy atom of a standard residue is absent.
	IssueMissingAtom
	// IssueTruncatedRecord: a structural record was shorter than its fixed width.
	IssueTruncatedRecord
	// IssueMalformedRecord: a structural record had unparsable fields and was dropped.
	IssueMalformedRecord
	// IssueUnsupportedRecord: a recognized structural record type is not handled.
	IssueUnsupportedRecord
)

// String returns a stable identifier for the issue kind.
func (k IssueKind) String() string {
	switch k {
	case IssueMultiModel:
		return "multi-model"
	case IssueAlternateConformation:
		return "alternate-conformation"
	case IssueInsertionCode:
		return "insertion-code"
	case IssueNonStandardResidue:
		return "non-standard-residue"
	case IssueMissingAtom:
		return "missing-atom"
	case IssueTruncatedRecord:
		return "truncated-record"
	case IssueMalformedRecord:
		return "malformed-record"
	case IssueUnsupportedRecord:
		return "unsupported-record"
	default:
		return "unknown"
	}
}

// Issue describes one distinguishable parse outcome.
type Issue struct {
	Kind    IssueKind
	Line    int
	Where   string
	Message string
}

func (i Issue) String() string {
	return fmt.Sprintf("%s (line %d%s): %s", i.Kind, i.Line, i.Where, i.Message)
}

// ParseErrorKind distinguishes why parsing failed, so callers can react by
// type instead of matching message text.
type ParseErrorKind int

const (
	ParseErrUnknown ParseErrorKind = iota
	// ParseErrIO: the file could not be opened or read.
	ParseErrIO
	// ParseErrContent: the content could not be parsed into a structure
	// (e.g. no coordinate records at all).
	ParseErrContent
)

// ParseError is returned when input cannot be turned into a usable structure.
type ParseError struct {
	Source string
	Reason string
	Kind   ParseErrorKind
}

func (e *ParseError) Error() string {
	if e.Source == "" {
		return "structure parse error: " + e.Reason
	}
	return fmt.Sprintf("structure parse error in %s: %s", e.Source, e.Reason)
}

// HasIssue reports whether the structure collected an issue of the given kind.
func (s *Structure) HasIssue(k IssueKind) bool {
	for _, w := range s.Warnings {
		if w.Kind == k {
			return true
		}
	}
	return false
}

// ChainByID looks up a chain in the effective model.
func (s *Structure) ChainByID(id string) *Chain {
	if s.Effective == nil {
		return nil
	}
	for _, c := range s.Effective.Chains {
		if c.ID == id {
			return c
		}
	}
	return nil
}
