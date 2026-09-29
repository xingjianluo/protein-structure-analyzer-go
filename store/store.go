// Package store 提供支持并发查询与在线重载的结构数据存储。
package store

import (
	"errors"
	"fmt"
	"os"
	"sort"
	"sync"

	"github.com/xingjianluo/protein-structure-analyzer-go/structure"
)

// State 表示存储当前状态。
type State string

const (
	StateEmpty    State = "EMPTY"    // 尚未成功加载任何数据
	StateReady    State = "READY"    // 数据可用
	StateDegraded State = "DEGRADED" // 最近一次重载失败，保留上一份可用数据
)

// Entry 是一个已加载的结构及其解析报告。
type Entry struct {
	Name   string
	Model  *structure.Model
	Report *structure.ParseReport
}

// Status 描述存储当前状态快照。
type Status struct {
	State      State
	Generation uint64
	Entries    []string
	LastError  string
}

// ErrNotFound 表示查询的结构不存在。
var ErrNotFound = errors.New("structure not loaded")

// Store 线程安全地持有当前数据快照。
//
// 重载流程：先在锁外完整解析全部输入并构建新快照，成功后在写锁内原子替换；
// 任一输入失败则整体回滚，保留上一份可用快照，状态转为 DEGRADED 并记录失败原因。
// 每次成功重载递增 Generation，查询结果携带 Generation，
// 调用方可据此识别结论所基于的数据版本，避免误用旧结论。
type Store struct {
	mu         sync.RWMutex
	generation uint64
	entries    map[string]*Entry
	state      State
	lastErr    string
}

// New 创建空存储。
func New() *Store {
	return &Store{state: StateEmpty, entries: map[string]*Entry{}}
}

// snapshot 返回当前快照的只读视图（map 不复制，调用方不得修改）。
func (s *Store) snapshot() (map[string]*Entry, uint64) {
	s.mu.RLock()
	defer s.mu.RUnlock()
	return s.entries, s.generation
}

// Load 解析给定文件并原子替换当前快照；任一文件失败则整体回滚。
// 同名文件以后出现者为准。传入空列表等价于清空。
func (s *Store) Load(paths []string) error {
	entries := make(map[string]*Entry, len(paths))
	for _, p := range paths {
		f, err := os.Open(p)
		if err != nil {
			s.markFailure(fmt.Errorf("load %s: %w", p, err))
			return err
		}
		model, rep, err := structure.ParsePDB(p, f)
		f.Close()
		if err != nil {
			s.markFailure(err)
			return err
		}
		entries[p] = &Entry{Name: p, Model: model, Report: rep}
	}

	s.mu.Lock()
	s.entries = entries
	s.generation++
	s.state = StateReady
	s.lastErr = ""
	s.mu.Unlock()
	return nil
}

// markFailure 记录一次失败的重载：保留旧快照，仅更新状态与失败原因。
func (s *Store) markFailure(err error) {
	s.mu.Lock()
	defer s.mu.Unlock()
	if s.state == StateReady {
		s.state = StateDegraded
	}
	s.lastErr = err.Error()
}

// Interfaces 返回指定结构的链间界面及数据版本号；数据不存在时返回 ErrNotFound。
func (s *Store) Interfaces(name string, opts structure.ContactOptions) ([]structure.Interface, uint64, error) {
	entries, gen := s.snapshot()
	e, ok := entries[name]
	if !ok {
		return nil, gen, ErrNotFound
	}
	// Model 加载后不可变，可在锁外安全计算
	return structure.Interfaces(e.Model, opts), gen, nil
}

// Contacts 返回指定结构的残基接触及数据版本号。
func (s *Store) Contacts(name string, opts structure.ContactOptions) ([]structure.ResidueContact, uint64, error) {
	entries, gen := s.snapshot()
	e, ok := entries[name]
	if !ok {
		return nil, gen, ErrNotFound
	}
	return structure.Contacts(e.Model, opts), gen, nil
}

// Report 返回指定结构的解析报告。
func (s *Store) Report(name string) (*structure.ParseReport, uint64, error) {
	entries, gen := s.snapshot()
	e, ok := entries[name]
	if !ok {
		return nil, gen, ErrNotFound
	}
	return e.Report, gen, nil
}

// Status 返回当前状态。
func (s *Store) Status() Status {
	s.mu.RLock()
	defer s.mu.RUnlock()
	names := make([]string, 0, len(s.entries))
	for n := range s.entries {
		names = append(names, n)
	}
	sort.Strings(names)
	return Status{
		State:      s.state,
		Generation: s.generation,
		Entries:    names,
		LastError:  s.lastErr,
	}
}
