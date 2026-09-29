package store

import (
	"sort"
	"sync"

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

// Snapshot is an immutable load+analysis result.
type Snapshot struct {
	Version  int
	State    LoadState
	Paths    []string
	Reports  map[string]*contact.InterfaceReport
	Entries  map[string]*EntryState
	Warnings map[string][]structure.Issue
}

// Store keeps the latest successfully committed snapshot behind an RWMutex.
type Store struct {
	// loadMu serializes whole load attempts so concurrent reloads cannot race.
	loadMu   sync.Mutex
	mu       sync.RWMutex
	current  *Snapshot
	analyzer *contact.Analyzer
	version  int
}

// New creates a Store backed by the given analyzer.
func New(a *contact.Analyzer) *Store {
	if a == nil {
		a = contact.NewAnalyzer(contact.DefaultCutoff)
	}
	return &Store{analyzer: a}
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

// Load parses and analyzes the given local paths and atomically replaces the
// current snapshot. On any failure the previous snapshot is retained and a
// *ReloadError is returned; concurrent readers never see a partial update.
func (st *Store) Load(paths []string) error {
	// Serialize load attempts; readers are not blocked during parsing.
	st.loadMu.Lock()
	defer st.loadMu.Unlock()
	// Build the entire new snapshot off-lock, then swap once.
	sortedPaths := append([]string{}, paths...)
	sort.Strings(sortedPaths)

	reports := map[string]*contact.InterfaceReport{}
	entries := map[string]*EntryState{}
	warningsBySource := map[string][]structure.Issue{}
	worst := StateOK
	var failures []string

	for _, path := range sortedPaths {
		s, perr := structure.ParseFile(path)
		state, warns, errMsg := classifyEntry(s, perr)
		if perr != nil {
			entries[path] = &EntryState{Source: path, State: StateError, Error: errMsg}
			failures = append(failures, path+": "+errMsg)
			worst = StateError
			continue
		}
		rep, aerr := st.analyzer.Analyze(s)
		if aerr != nil {
			entries[path] = &EntryState{Source: path, State: StateError, Error: aerr.Error()}
			failures = append(failures, path+": "+aerr.Error())
			worst = StateError
			continue
		}
		if state == StateDegraded && worst != StateError {
			worst = StateDegraded
		}
		entries[path] = &EntryState{Source: path, State: state, Warnings: warns}
		reports[path] = rep
		warningsBySource[path] = warns
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

	st.mu.Lock()
	st.version++
	snap := &Snapshot{
		Version:  st.version,
		State:    worst,
		Paths:    sortedPaths,
		Reports:  reports,
		Entries:  entries,
		Warnings: warningsBySource,
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
