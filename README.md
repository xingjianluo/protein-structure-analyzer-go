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
- **批量跨结构对比**：一组结构一次提交，直接给出哪些链间界面是共有、
  哪些只出现在部分结构中，以及同一对链的接触在各结构间是否发生明显变化；
  非标准记录（金属离子、修饰残基等 HETATM）默认不参与，可用开关显式纳入，
  也可只挑选关心的若干条链。
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
| `contact` | `distance.go` | 欧氏距离、残基间最小重原子距离、包围盒预过滤 |
| | `contact.go` | 接触/界面判定、批量分析、确定性排序 |
| | `errors.go` | 可断言的错误类型 |
| `compare` | `compare.go` | 跨结构批量对比：共有/部分界面、接触变化、输入分级 |
| | `filter.go` | 对比选项归一化、链选择与 HETATM 过滤 |
| `store` | `store.go` | 快照、并发读写、原子加载、失败回滚、对比选项与输入指纹 |
| | `reload.go` | `ReloadError` 与类型化的逐文件失败（`FailureKind`） |
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
| 非标准/杂原子残基 `IssueNonStandardResidue` | `HETATM` 或非 20 种标准氨基酸：保留但显式标记 | 单结构分析按重原子距离正常参与；**批量对比默认排除**，需 `IncludeHetero` 开关显式纳入 |
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
- 界面 = 两链之间全部接触的聚合：接触数 + 最小距离。
- 输出顺序是**规范序**：链对按链 ID 排序（每对内部 `ChainA ≤ ChainB`），
  接触按 `(ChainA, ChainB, SeqA, ICodeA, SeqB, ICodeB, NameA)` 排序。
  因此同一份输入无论以何种顺序提交，报告逐字节一致。
- 浮点距离直接比较 `d <= cutoff`，无任何随机或外部输入，重复计算结果相同。

### 并发与在线重新加载

- `Store.Load(paths)` 先在锁外构建**完整新快照**，再一次性替换指针；
  读路径只持有 `RLock`，任何时刻读到的快照都是自洽完整的。
- 并发 `Load` 由独立的加载锁串行化；重载期间查询永不阻塞于解析过程。
- **任何一个输入失败则整批不提交**：当前快照原样保留，返回 `*store.ReloadError`，
  其中包含 `Reason`、逐文件 `Failures`、被保留快照的 `Version` 与 `Retained` 标志，
  可与成功/降级状态断言区分。
- 快照不可变；加载新文件集合后，旧路径键与旧报告指针不再可达。

### 批量跨结构对比（`compare` 包）

一组相关结构（同一复合物的不同条件、野生型 + 突变体等）一次提交，
`Store.Load` 在生成快照时自动完成跨结构对比，结果挂在 `Snapshot.Comparison`：

- **输入分级（按类型区分，不是比文字）**：每个文件在 `Report.Files` 中带有
  `compare.FileClass`——`compared`（健康，参与对比）、`degraded`（解析成功但有
  缺原子/畸形记录等质量问题，**被排除在对比之外**，绝不拿残缺数据当正常结构）、
  `excluded`（完全不可用）。完全打不开的文件不会进入对比报告，而是让整次加载
  以 `*store.ReloadError` 回滚，其 `Failures` 用 `FailureKind` 区分
  `io`（打不开/读不出）与 `parse`（内容无法解析，包括坐标列整段缺失）。
- **共有 / 部分界面**：`Report.Interfaces` 列出全部出现过的链对，
  `Shared=true` 表示在所有参与对比的结构中都存在；`PresentIn`/`AbsentFrom`
  逐个列出在哪些结构里出现/缺席。
- **接触变化**：`Report.Changes` 对每对链给出逐结构的接触数与最小距离；
  当存在性不同、接触数不同、或最小距离相差超过 `compare.SignificantDistDelta`
  （0.5 Å）时 `Changed=true`。
- **规范顺序**：输入先按 ID（文件路径）排序再计算，所有输出列表按
  `(ChainA, ChainB)` 排序。同一批文件无论以什么顺序提交、重复算多少次，
  报告逐字段一致（`TestDeterministicRegardlessOfOrder` 用全排列验证）。
- **规模**：接触判定带包围盒预过滤，几十份 × 每份上千残基的批量对比
  在亚秒级完成（`TestBatchScaleCompletes`：30 份 × 2 链 × 1000 残基）。

### 对比开关与选择范围

`compare.Options`（通过 `store.SetOptions` 设置，下一次 `Load` 生效）：

- `IncludeHetero`（默认 `false`）：金属离子、修饰残基等非标准记录
  （`HETATM`）**默认不参与**界面判定；置 `true` 显式纳入。
- `Chains`（默认空 = 全部链）：只让列出的链参与对比；归一化后排序去重。

每个快照都记录自己实际使用的 `Snapshot.Options`（与 `Comparison.Options` 一致）
和输入指纹 `Snapshot.Fingerprint`（排序后的路径 + 文件内容的 SHA-256）。
换了开关或换了数据必须重新 `Load`：新结果只会出现在新快照里，
旧快照整体被替换、绝不混用；指纹只随数据变化，与提交顺序无关。

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
批量对比的排列/重复一致性、HETATM 默认排除与开关纳入、链选择归一化、
降级/排除输入的类型化分类、开关切换后旧结论失效、类型化失败（io/parse）
与回滚保留对比结果、指纹随数据变化、并发读取快照自洽、30 份 × 千残基规模冒烟。
单测通过 `t.Logf` 打印关键夹具输入与判定依据（所选构象、最小距离、保留版本等）。

## 命令行用法

```bash
go run ./cmd/protein-analyzer [--cutoff 4.5] [--include-hetero] [--chains A,B] file1.pdb file2.cif ...
```

输出示例（两份结构批量对比）：

```
snapshot version 1, state=ok, cutoff=4.50 A, fingerprint=a80ae1da32a8
options: include-hetero=false chains=(all)

== /path/wt.pdb (state=ok, model=1, heavy atoms=8, usable residues=2) ==
  interface A--B: 1 residue contact(s), min distance 3.00 A
    contact A/GLY1  - B/GLY1   d=3.00 A

== batch comparison (2 of 2 files compared) ==
  interface A--B [partial]: present in [/path/wt.pdb], absent from [/path/mut.pdb]
  changed contacts A--B:
    /path/mut.pdb: no interface
    /path/wt.pdb: 1 contact(s), min 3.00 A
```

任一文件解析失败时进程以非零码退出，打印 `ReloadError`（含逐文件原因），
不产生任何部分结果。
