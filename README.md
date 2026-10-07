# 蛋白质结构解析与复合物界面判定（本地、可复现）

纯 Go 实现，零外部依赖，不访问任何网络服务或在线数据库。输入本地结构文件，
解析为统一的「链–残基–原子」模型，并基于重原子距离判定残基接触与链间界面。

## 能力概览

- 本地解析 **PDB**（固定列 `ATOM`/`HETATM`）与 **mmCIF**（`_atom_site` 类别）。
- 统一模型：`Structure → Model → Chain → Residue → Atom`。
- 对多模型、可变构象、插入码、非标准残基、缺原子、畸形/截断记录给出
  **可区分的 `Issue` 结论**，绝不静默跳过，也不让坏数据参与计算。
- 残基接触与链间界面判定，结果只依赖输入内容：重复计算、乱序提交结论一致；
  底层数据重新加载后旧结论不再可达。
- **批量跨结构对比**：一组结构一起提交，直接给出共有界面、部分存在的界面、
  以及同一对链之间发生明显变化的残基接触；结论按规范序输出，与提交顺序无关。
- 非标准记录（金属离子、修饰残基等 `HETATM`）**默认不参与**界面判定，
  可用开关显式纳入；也可只挑选若干条链参与对比。所选开关与范围记录在结果中，
  更换开关或数据后旧结论自动失效（按选项+输入指纹缓存）。
- 并发安全的存储层：读写并发无半更新状态、无数据竞争（`go test -race` 验证）；
  重载失败整体回滚，保留上一份可用快照，并以独立错误类型返回失败原因。

## 包结构

| 包 | 文件 | 职责 |
| --- | --- | --- |
| `structure` | `model.go` | 统一模型、`IssueKind`、`ParseError` |
| | `standard.go` | 20 种标准氨基酸及其期望重原子表 |
| | `pdb.go` | PDB 固定列解析 |
| | `mmcif.go` | mmCIF `_atom_site` 解析与格式自动识别（`ParseFile`） |
| | `analyze.go` | 状态分级、重原子统计 |
| `contact` | `distance.go` | 欧氏距离、残基间最小重原子距离 |
| | `contact.go` | 接触/界面判定、批量分析、确定性排序 |
| | `options.go` | `Options`（cutoff、杂原子开关、链选择）及其指纹 |
| | `errors.go` | 可断言的错误类型 |
| `compare` | `compare.go` | 跨结构批量对比、结果指纹、带失效语义的缓存 |
| `store` | `store.go` | 快照、并发读写、原子加载、失败回滚 |
| | `reload.go` | `ReloadError`（区分失败原因与保留版本） |
| `cmd/protein-analyzer` | `main.go` | 本地命令行入口 |

## 解析规则

### 统一模型

- 一条输入 = 一个 `Structure`；PDB 无 `MODEL` 记录时视为 1 号模型。
- 链按出现顺序保留；残基键为 `(链, 残基序号, 插入码)`，插入码不同视为不同残基。

### 各类异常的明确结论（`structure.IssueKind`）

| 场景 | 结论 | 对计算的影响 |
| --- | --- | --- |
| 多模型 `IssueMultiModel` | 只使用**第一个模型**（mmCIF 取最小 model 序号），并记录该问题 | 其余模型不参与计算 |
| 可变构象 `IssueAlternateConformation` | 每个残基独立选择：选中构象原子占据率之和最大；并列时取 ASCII 最小标签 | 仅选中构象参与计算，结果确定 |
| 插入码 `IssueInsertionCode` | 作为独立残基保留并标记 | 正常参与计算 |
| 非标准/杂原子残基 `IssueNonStandardResidue` | `HETATM` 或非 20 种标准氨基酸：保留但显式标记 | **默认不参与**界面判定；设置 `IncludeHetero` 后按重原子距离正常参与 |
| 缺原子 `IssueMissingAtom` | 标准氨基酸缺失任何期望重原子（骨架 N/CA/C/O + 标准侧链）时，缺失名进入 `Residue.Missing` | 该残基**不参与**接触判定 |
| 截断记录 `IssueTruncatedRecord` | 记录短于坐标列（第 54 列）：丢弃该行并记录 | 不计入模型 |
| 畸形记录 `IssueMalformedRecord` | 坐标/残基序号不可解析、原子名为空等：丢弃该行并记录 | 不计入模型 |
| 无任何原子 | 返回 `*structure.ParseError`（可 `errors.As` 断言） | 加载失败 |

解析后状态分三级（`store.LoadState`）：

- `ok`：无数据质量问题；
- `degraded`：解析成功但存在缺原子、丢弃畸形/截断记录或构象折叠等；
- `error`：文件无法解析或分析失败——与正常状态明确区分。

### 接触与界面判定

- 两个**不同链**上的残基，任一重原子（排除 H/D 氢/氘）对之间最小距离
  `≤ cutoff`（默认 **4.5 Å**，可通过 `NewAnalyzer(cutoff)` 或 `--cutoff` 调整）
  即记为一次残基接触。
- `UsableForContact()` 为假的残基（缺期望原子或无原子）被排除。
- **非标准残基默认被排除**：`HETATM` 记录或非 20 种标准氨基酸（金属离子、
  修饰残基、配体等）不参与界面判定；`contact.Options.IncludeHetero=true`
  （CLI `--include-hetero`）可显式纳入。
- **链选择**：`contact.Options.Chains`（CLI `--chains A,B`）非空时，只有列出的
  链参与判定；为空表示全部链。
- 界面 = 两链之间全部接触的聚合：接触数 + 最小距离。
- 输出顺序是**规范序**：链对按链 ID 排序（每对内部 `ChainA ≤ ChainB`），
  接触按 `(ChainA, ChainB, SeqA, ICodeA, SeqB, ICodeB, NameA)` 排序。
  因此同一份输入无论以何种顺序提交，报告逐字节一致。
- 浮点距离直接比较 `d <= cutoff`，无任何随机或外部输入，重复计算结果相同。

### 批量跨结构对比（`compare` 包）

- `Comparer.Compare(structures, opts)` 对一批结构在同一组 `Options` 下逐一分析，
  再聚合为一份 `compare.Result`：
  - `Common`：在**每一份**结构中都存在的链间界面；
  - `Partial`：只在**部分**结构中存在的界面，`PresentIn` 按来源名排序；
  - `Changes`：残基级接触的变化——在某些结构中缺失，或距离极差
    （max−min）超过 `DistanceChangeTolerance`（**0.5 Å**）。
- 结果排序固定：界面按 `(ChainA, ChainB)`，变化按
  `(ChainA, ChainB, SeqA, ICodeA, SeqB, ICodeB)`；来源列表统一按名字排序。
  同一批文件**无论以什么顺序提交、重复算多少次，结论逐字段一致**。
- **指纹与缓存**：`Result.Fingerprint` 由选项指纹 + 全部输入的分析内容哈希
  （FNV-1a）构成。`Comparer` 以该指纹缓存上一次结论：完全相同的批次直接命中；
  **只要换了开关、链选择或任何一份数据，指纹必然不同，绝不复用旧结论**。
- 批量中任一输入无法分析时返回 `*compare.BatchError`（可用 `errors.As` 断言），
  其 `Failures` 逐项为 `*contact.InputError`（含序号、来源、底层错误），
  且缓存中的上一份结论不受影响。
- 性能：聚合为 O(总接触数)；12 份 × 2 链 × 250 残基的批量在单测中约 0.04 s
  完成（见 `TestCompareScalesToBatch`）。

### 并发与在线重新加载

- `Store.Load(paths)` / `Store.LoadWithOptions(paths, opts)` 先在锁外构建
  **完整新快照**，再一次性替换指针；读路径只持有 `RLock`，任何时刻读到的
  快照都是自洽完整的（要么完整的旧结果、要么完整的新结果）。
- 并发 `Load` 由独立的加载锁串行化；重载期间查询永不阻塞于解析过程。
- **任何一个输入失败则整批不提交**：当前快照原样保留，返回 `*store.ReloadError`，
  其中包含 `Reason`、逐文件 `Failures`、被保留快照的 `Version` 与 `Retained` 标志，
  可与成功/降级状态断言区分。失败类型可按来源区分：打不开的文件与坐标列
  无法解析的文件都产生 `*structure.ParseError`，逐文件状态见 `Snapshot.Entries`
  （`ok` / `degraded` / `error`）。
- 每个快照携带 `Options` 与 `Comparison`（含指纹）：换了开关或数据后提交的
  新快照指纹必然不同，旧结论随旧快照一起不再可达。
- 快照不可变；加载新文件集合后，旧路径键与旧报告指针不再可达。

## 支持的输入范围

- PDB：`ATOM`/`HETATM`（含 `MODEL`/`ENDMDL`、altLoc、插入码、占据率、元素列）。
- mmCIF：`loop_` + `_atom_site.*`（支持单/双引号、`;` 文本块、`#` 注释、`.`/`?` 缺失值）；
  列名优先取 `auth_*`，缺省回退 `label_*`。
- 仅使用本地文件，按扩展名（`.cif`/`.mmcif`）或回退 PDB 解析。
- 不解析二级结构、化学键、生物大分子装配（BIOMT）等；不联网。

## 本地验证

```bash
# 构建
go build ./...

# 全部自动化测试（含竞态检测）
go test -race ./...

# 静态检查
go vet ./...
go fmt 已在提交前通过（gofmt -l . 无输出）
```

测试覆盖：多模型、可变构象确定性选择、插入码、非标准残基、缺原子排除、
截断与畸形记录区分、空文件失败、重复计算与乱序一致性、重载替换旧结论、
重载失败回滚、16 读协程 × 200 次成功重载 + 100 次失败重载的竞态场景；
批量对比的共有/部分界面归类、接触变化（缺失或距离极差 > 0.5 Å）识别、
排序与乱序/重复一致性、杂原子开关与链选择、开关或数据变化后缓存失效、
坏文件（打不开 / 坐标列不可解析）触发整体回滚且保留旧对比结论、
降级条目与正常/错误条目的逐文件区分、以及 12 份 × 500 残基的批量规模验证。
单测通过 `t.Logf` 打印关键夹具输入与判定依据（所选构象、最小距离、保留版本、
对比指纹等），`go test -v` 可见。

## 命令行用法

```bash
go run ./cmd/protein-analyzer [--cutoff 4.5] [--include-hetero] [--chains A,B] file1.pdb file2.cif ...
```

输出示例：

```
snapshot version 1, state=ok, cutoff=4.50 A, include-hetero=false, chains=all

== /path/two.pdb (state=ok, model=1, heavy atoms=8, usable residues=2) ==
  interface A--B: 1 residue contact(s), min distance 1.90 A
    contact A/GLY1  - B/GLY1   d=1.90 A

== batch comparison (fingerprint cmp-..., 2 structures) ==
  common interfaces (present in all 2):
    A--B: contacts per structure map[...]
  partial interfaces (present in some):
  changed contacts:
```

任一文件解析失败时进程以非零码退出，打印 `ReloadError`（含逐文件原因），
不产生任何部分结果。
