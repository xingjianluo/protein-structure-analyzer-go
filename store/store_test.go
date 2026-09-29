package store

import (
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"sync"
	"testing"

	"github.com/xingjianluo/protein-structure-analyzer-go/structure"
)

// atomLine 构造一条 PDB ATOM 记录（与 structure 包测试相同的列格式）。
func atomLine(serial int, name, resName, chain string, seq int, x, y, z float64) string {
	return fmt.Sprintf("ATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f",
		serial, name, resName, chain, seq, x, y, z, 1.0)
}

// writePDB 写入临时结构文件，返回路径。gap 控制两条链的间距。
func writePDB(t *testing.T, dir, name string, gap float64) string {
	t.Helper()
	content := atomLine(1, "CA", "GLY", "A", 1, 0, 0, 0) + "\n" +
		atomLine(2, "CA", "SER", "B", 1, gap, 0, 0) + "\n"
	p := filepath.Join(dir, name)
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatal(err)
	}
	return p
}

func TestLoadAndQuery(t *testing.T) {
	dir := t.TempDir()
	p := writePDB(t, dir, "a.pdb", 3.0)
	s := New()
	if st := s.Status(); st.State != StateEmpty {
		t.Fatalf("new store must be EMPTY, got %s", st.State)
	}
	if err := s.Load([]string{p}); err != nil {
		t.Fatalf("load: %v", err)
	}
	ifaces, gen, err := s.Interfaces(p, structure.DefaultContactOptions())
	if err != nil {
		t.Fatalf("interfaces: %v", err)
	}
	t.Logf("gen=%d interfaces=%+v", gen, ifaces)
	if gen != 1 || len(ifaces) != 1 || ifaces[0].ChainA != "A" || ifaces[0].ChainB != "B" {
		t.Fatalf("unexpected result gen=%d ifaces=%+v", gen, ifaces)
	}
	if st := s.Status(); st.State != StateReady || st.LastError != "" {
		t.Fatalf("status: %+v", st)
	}
}

func TestReloadReplacesStaleResults(t *testing.T) {
	dir := t.TempDir()
	near := writePDB(t, dir, "near.pdb", 3.0)
	far := writePDB(t, dir, "far.pdb", 100.0)
	s := New()
	if err := s.Load([]string{near}); err != nil {
		t.Fatal(err)
	}
	ifaces1, gen1, _ := s.Interfaces(near, structure.DefaultContactOptions())
	if len(ifaces1) != 1 {
		t.Fatalf("expected interface before reload, got %+v", ifaces1)
	}
	// 用无界面的数据重载：旧结论不得继续返回
	if err := s.Load([]string{far}); err != nil {
		t.Fatal(err)
	}
	if _, _, err := s.Interfaces(near, structure.DefaultContactOptions()); !errors.Is(err, ErrNotFound) {
		t.Fatalf("stale entry must be gone, err=%v", err)
	}
	ifaces2, gen2, err := s.Interfaces(far, structure.DefaultContactOptions())
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("before: gen=%d ifaces=%d; after: gen=%d ifaces=%d", gen1, len(ifaces1), gen2, len(ifaces2))
	if gen2 <= gen1 {
		t.Fatalf("generation must increase: %d -> %d", gen1, gen2)
	}
	if len(ifaces2) != 0 {
		t.Fatalf("expected no interfaces after reload, got %+v", ifaces2)
	}
}

func TestFailedReloadRollsBack(t *testing.T) {
	dir := t.TempDir()
	good := writePDB(t, dir, "good.pdb", 3.0)
	bad := filepath.Join(dir, "bad.pdb")
	if err := os.WriteFile(bad, []byte("HEADER    NO ATOMS HERE\nEND\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	missing := filepath.Join(dir, "missing.pdb")

	s := New()
	if err := s.Load([]string{good}); err != nil {
		t.Fatal(err)
	}
	before, genBefore, _ := s.Interfaces(good, structure.DefaultContactOptions())

	// 失败原因一：文件不存在（I/O 错误）
	err := s.Load([]string{missing})
	if err == nil {
		t.Fatal("expected error for missing file")
	}
	var pe *structure.ParseError
	if errors.As(err, &pe) {
		t.Fatalf("missing file must be an I/O error, not ParseError: %v", err)
	}
	st := s.Status()
	t.Logf("after missing-file reload: state=%s lastErr=%q", st.State, st.LastError)
	if st.State != StateDegraded || st.LastError == "" {
		t.Fatalf("expected DEGRADED with reason, got %+v", st)
	}

	// 失败原因二：文件存在但无可用原子（ParseError，可与 I/O 错误区分）
	err = s.Load([]string{bad})
	if !errors.As(err, &pe) {
		t.Fatalf("expected *structure.ParseError, got %T: %v", err, err)
	}
	t.Logf("parse failure reason: %s", pe.Reason)

	// 回滚验证：旧数据与旧结论仍然可用，generation 不变
	after, genAfter, err := s.Interfaces(good, structure.DefaultContactOptions())
	if err != nil {
		t.Fatalf("old snapshot must survive failed reloads: %v", err)
	}
	if genAfter != genBefore || len(after) != len(before) {
		t.Fatalf("rollback violated: gen %d->%d, ifaces %d->%d", genBefore, genAfter, len(before), len(after))
	}

	// 成功重载后状态恢复 READY 且清除失败原因
	if err := s.Load([]string{good}); err != nil {
		t.Fatal(err)
	}
	if st := s.Status(); st.State != StateReady || st.LastError != "" {
		t.Fatalf("expected READY after successful reload, got %+v", st)
	}
}

func TestQueryNotFound(t *testing.T) {
	s := New()
	_, _, err := s.Interfaces("nope.pdb", structure.DefaultContactOptions())
	if !errors.Is(err, ErrNotFound) {
		t.Fatalf("expected ErrNotFound, got %v", err)
	}
}

func TestConcurrentQueryAndReload(t *testing.T) {
	dir := t.TempDir()
	near := writePDB(t, dir, "near.pdb", 3.0)
	far := writePDB(t, dir, "far.pdb", 100.0)
	s := New()
	if err := s.Load([]string{near, far}); err != nil {
		t.Fatal(err)
	}

	const workers = 8
	const rounds = 200
	var wg sync.WaitGroup
	errs := make(chan error, workers*rounds)

	// 读方：并发查询两个结构
	for w := 0; w < workers; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for i := 0; i < rounds; i++ {
				ifaces, _, err := s.Interfaces(near, structure.DefaultContactOptions())
				if err != nil && !errors.Is(err, ErrNotFound) {
					errs <- fmt.Errorf("query: %w", err)
					return
				}
				// 结果要么不存在（重载间隙不存在——本实现原子替换故不应发生），
				// 要么恰好 0 或 1 个界面，绝不可能是半更新状态
				if err == nil && len(ifaces) != 1 {
					errs <- fmt.Errorf("near.pdb must have exactly 1 interface, got %d", len(ifaces))
					return
				}
				_ = s.Status()
			}
		}()
	}
	// 写方：交替重载 单文件/双文件 快照
	for i := 0; i < 20; i++ {
		var err error
		if i%2 == 0 {
			err = s.Load([]string{near, far})
		} else {
			err = s.Load([]string{near})
		}
		if err != nil {
			t.Fatalf("reload %d: %v", i, err)
		}
	}
	wg.Wait()
	close(errs)
	for err := range errs {
		t.Fatal(err)
	}
	st := s.Status()
	t.Logf("final state=%s generation=%d entries=%v", st.State, st.Generation, st.Entries)
	if st.State != StateReady {
		t.Fatalf("expected READY, got %+v", st)
	}
}
