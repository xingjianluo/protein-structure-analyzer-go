package contact

import "fmt"

// ErrEmptyStructure is returned when there is no effective model to analyze.
var ErrEmptyStructure = &AnalysisError{Reason: "empty or unparsed structure"}

// AnalysisError marks a contact-analysis failure.
type AnalysisError struct{ Reason string }

func (e *AnalysisError) Error() string { return "contact analysis error: " + e.Reason }

// InputError identifies which input in a batch failed.
type InputError struct {
	Index  int
	Source string
	Err    error
}

func (e *InputError) Error() string {
	return fmt.Sprintf("input %d (%s): %v", e.Index, e.Source, e.Err)
}

func (e *InputError) Unwrap() error { return e.Err }
