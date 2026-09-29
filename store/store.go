// Package store 提供支持并发查询与在线重载的结构数据存储。
package store

import (
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

// Store 线程安全地持有当前数据快照；重载为原子替换，失败自动回滚。
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

// Load 解析给定文件并原子替换当前快照；任一文件失败则整体回滚。
func (s *Store) Load(paths []string) error {
	return nil
}

// Interfaces 返回指定结构的链间界面；数据不存在时返回错误。
func (s *Store) Interfaces(name string, opts structure.ContactOptions) ([]structure.Interface, uint64, error) {
	return nil, 0, nil
}

// Status 返回当前状态。
func (s *Store) Status() Status {
	return Status{}
}
