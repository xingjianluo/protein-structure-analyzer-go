package store

import (
	"crypto/sha256"
	"encoding/hex"
	"errors"
	"os"
	"sort"
	"sync"

	"github.com/xingjianluo/protein-structure-analyzer-go/compare"
	"github.com/xingjianluo/protein-structure-analyzer-go/contact"
	"github.com/xingjianluo/protein-structure-analyzer-go/structure"
)

// LoadState distinguishes a healthy snapshot from degraded and failed loads.
type LoadState int

const (
	StateUnknown LoadState = iota
	StateOK
	StateDegraded
	StateError
)

func (s LoadState) String() string {
	switch s {
	case StateOK:
		return "ok"
	case StateDegraded:
		return "degraded"
	case StateError:
		return "error"
	default:
		return "unknown"
	}
}

// EntryState is the per-file state within a snapshot.
type EntryState struct {
	Source   string
	State    LoadState
	Warnings []structure.Issue
	Error    string
}

// Snapshot is an immutable load+analysis+comparison result.
type Snapshot struct {
	Version  int
	State    LoadState
	Paths    []string
	Reports  map[string]*contact.InterfaceReport
	Entries  map[string]*EntryState
	Warnings map[string][]structure.Issue
	// Options echoes the normalized comparison options used for this
	// snapshot, so callers can verify which switches a result belongs to.
	Options compare.Options
	// Comparison is the cross-structure interface comparison of this batch.
	Comparison *compare.Report
	// Fingerprint identifies the exact input set (paths + file content).
	// Any change to the data produces a different fingerprint.
	Fingerprint string
}

// Store keeps the latest successfully committed snapshot behind an RWMutex.
type Store struct {
	// loadMu serializes whole load attempts so concurrent reloads cannot race.
	loadMu   sync.Mutex
	mu       sync.RWMutex
	current  *Snapshot
	analyzer *contact.Analyzer
	comparer *compare.Comparer
	opts     compare.Options
	version  int
}

// New creates a Store backed by the given analyzer.
func New(a *contact.Analyzer) *Store {
	if a == nil {
		a = contact.NewAnalyzer(contact.DefaultCutoff)
	}
	return &Store{analyzer: a, comparer: compare.NewComparer(a.Cutoff())}
}

// SetOptions replaces the comparison options used by subsequent Load calls.
// It never mutates the committed snapshot: the new switches take effect only
// after a successful reload, so a committed result always matches its own
// Snapshot.Options.
func (st *Store) SetOptions(o compare.Options) {
	st.mu.Lock()
	st.opts = o
	st.mu.Unlock()
}

// Options returns the options that will be used by the next Load.
func (st *Store) Options() compare.Options {
	st.mu.RLock()
	defer st.mu.RUnlock()
	return st.opts
}

func classifyEntry(s *structure.Structure, parseErr error) (LoadState, []structure.Issue, string) {
	if parseErr != nil {
		return StateError, nil, parseErr.Error()
	}
	if s.IsDegraded() {
		return StateDegraded, s.Warnings, ""
	}
	return StateOK, s.Warnings, ""
}

// failureKind maps a parse/analysis error to a typed failure kind.
func failureKind(err error) FailureKind {
	var pe *structure.ParseError
	if errors.As(err, &pe) {
		switch pe.Kind {
		case structure.ParseErrIO:
			return FailureIO
		case structure.ParseErrContent:
			return FailureParse
		}
	}
	return FailureAnalysis
}

// Load parses, analyzes and compares the given local paths and atomically
// replaces the current snapshot. On any failure the previous snapshot is
// retained and a *ReloadError is returned; concurrent readers never see a
// partial update.
func (st *Store) Load(paths []string) error {
	// Serialize load attempts; readers are not blocked during parsing.
	st.loadMu.Lock()
	defer st.loadMu.Unlock()

	st.mu.RLock()
	opts := st.opts
	st.mu.RUnlock()

	// Build the entire new snapshot off-lock, then swap once.
	sortedPaths := append([]string{}, paths...)
	sort.Strings(sortedPaths)

	h := sha256.New()
	reports := map[string]*contact.InterfaceReport{}
	entries := map[string]*EntryState{}
	warningsBySource := map[string][]structure.Issue{}
	structs := map[string]*structure.Structure{}
	worst := StateOK
	var failures []EntryFailure

	for _, path := range sortedPaths {
		data, rerr := os.ReadFile(path)
		if rerr != nil {
			failures = append(failures, EntryFailure{
				Path: path, Kind: FailureIO, Message: rerr.Error(),
			})
			entries[path] = &EntryState{Source: path, State: StateError, Error: rerr.Error()}
			continue
		}
		// Fingerprint depends only on the (sorted) set of path+content.
		h.Write([]byte(path))
		h.Write([]byte{0})
		h.Write(data)
		h.Write([]byte{0})

		s, perr := structure.ParseBytes(data, path)
		state, warns, errMsg := classifyEntry(s, perr)
		if perr != nil {
			entries[path] = &EntryState{Source: path, State: StateError, Error: errMsg}
			failures = append(failures, EntryFailure{
				Path: path, Kind: failureKind(perr), Message: errMsg,
			})
			continue
		}
		rep, aerr := st.analyzer.Analyze(s)
		if aerr != nil {
			entries[path] = &EntryState{Source: path, State: StateError, Error: aerr.Error()}
			failures = append(failures, EntryFailure{
				Path: path, Kind: FailureAnalysis, Message: aerr.Error(),
			})
			continue
		}
		if state == StateDegraded && worst != StateError {
			worst = StateDegraded
		}
		entries[path] = &EntryState{Source: path, State: state, Warnings: warns}
		reports[path] = rep
		warningsBySource[path] = warns
		structs[path] = s
	}

	if len(failures) > 0 {
		// Rollback: do not touch current.
		st.mu.RLock()
		prevVersion := 0
		if st.current != nil {
			prevVersion = st.current.Version
		}
		st.mu.RUnlock()
		return &ReloadError{
			AttemptedPaths: sortedPaths,
			Version:        prevVersion,
			Retained:       st.current != nil,
			Reason:         "one or more inputs failed; previous snapshot retained",
			Failures:       failures,
		}
	}

	// Cross-structure comparison over every parsed input. Degraded inputs
	// are classified and excluded inside compare, never silently compared.
	inputs := make([]compare.Input, 0, len(structs))
	for _, p := range sortedPaths {
		if s, ok := structs[p]; ok {
			inputs = append(inputs, compare.Input{ID: p, Struct: s})
		}
	}
	comparison, cerr := st.comparer.Compare(inputs, opts)
	if cerr != nil {
		st.mu.RLock()
		prevVersion := 0
		if st.current != nil {
			prevVersion = st.current.Version
		}
		st.mu.RUnlock()
		return &ReloadError{
			AttemptedPaths: sortedPaths,
			Version:        prevVersion,
			Retained:       st.current != nil,
			Reason:         "comparison failed; previous snapshot retained",
			Failures:       []EntryFailure{{Path: "", Kind: FailureAnalysis, Message: cerr.Error()}},
		}
	}

	st.mu.Lock()
	st.version++
	snap := &Snapshot{
		Version:     st.version,
		State:       worst,
		Paths:       sortedPaths,
		Reports:     reports,
		Entries:     entries,
		Warnings:    warningsBySource,
		Options:     comparison.Options,
		Comparison:  comparison,
		Fingerprint: hex.EncodeToString(h.Sum(nil)),
	}
	st.current = snap
	st.mu.Unlock()
	return nil
}

// Current returns the committed snapshot (never a partially built one).
func (st *Store) Current() *Snapshot {
	st.mu.RLock()
	defer st.mu.RUnlock()
	return st.current
}

// Get returns the report and state for one source path.
func (st *Store) Get(source string) (*contact.InterfaceReport, LoadState, bool) {
	st.mu.RLock()
	defer st.mu.RUnlock()
	if st.current == nil {
		return nil, StateUnknown, false
	}
	if e, ok := st.current.Entries[source]; ok {
		if e.State == StateError {
			return nil, StateError, true
		}
	}
	rep, ok := st.current.Reports[source]
	if !ok {
		return nil, st.current.State, false
	}
	return rep, st.current.Entries[source].State, true
}
