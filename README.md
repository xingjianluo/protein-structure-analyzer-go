# protein-structure-analyzer-go

本地可运行的蛋白质结构解析与复合物界面判定工具。以本地 PDB 结构文件为输入，
不依赖任何外部服务或在线数据库，判定结论确定、可复现。

## 功能概览

- **结构解析**（`structure` 包）：把 PDB 文件解析为统一的 链-残基-原子 模型，
  对多模型、可变构象、插入码、非标准残基、缺原子、畸形/截断记录给出
  明确且可区分的处理结论（`ParseReport.Issues`），不静默跳过。
- **界面判定**（`structure` 包）：基于原子距离的残基接触与链间界面判定，
  结果规范排序，重复计算或改变输入顺序结论一致。
- **并发存储**（`store` 包）：支持并发查询与在线重载；重载为原子替换，
  失败自动回滚并保留上一份可用结论，失败原因与正常状态可区分。
- **命令行工具**（`cmd/psa`）：本地验证入口。

## 输入范围

- 仅支持 **PDB 格式**（定长列）的本地文件；识别 `ATOM` / `HETATM` / `MODEL` /
  `ENDMDL` 记录，其余记录（HEADER、REMARK、TER 等）不参与计算。
- 坐标单位：埃（Å）。文件不含任何可用原子记录时解析失败（`*structure.ParseError`）。

## 解析规则与可区分的处理结论

每条结论以 `Issue{Kind, Line, Detail}` 记录于 `ParseReport`，`Kind` 可断言区分：

| IssueKind | 触发场景 | 处理决策 |
|---|---|---|
| `MULTI_MODEL` | 文件含多个 `MODEL` 块 | 仅采用第一个模型，其余整体丢弃 |
| `ALT_LOC` | 同一原子位置有多个可变构象 | 保留占用率最高者，并列取 altLoc 字母序最小者 |
| `INSERTION_CODE` | 残基带插入码 | 插入码并入残基唯一标识（链+序号+插入码） |
| `NON_STANDARD_RES` | `HETATM` 或非标准残基名 | 保留于模型并标记 `NonStandard`，默认不参与界面判定 |
| `MISSING_ATOMS` | 标准残基缺少骨架 CA 原子 | 标记 `Incomplete`，不参与界面判定 |
| `DUPLICATE_ATOM` | 完全重复的原子记录 | 保留首条 |
| `MALFORMED_RECORD` | 数值字段不可解析、记录位于 MODEL 块外等 | 跳过该记录 |
| `TRUNCATED_RECORD` | 记录短于坐标列（54 列） | 跳过该记录 |
| `UNKNOWN_RECORD` | 无法识别的非空行 | 忽略并记录 |

致命错误（文件无任何可用原子、读取失败）返回 `*structure.ParseError`，
与单条记录级 Issue 区分。

## 界面判定规则

- 两残基分属不同链，且至少一对原子距离 ≤ 阈值（默认 4.5 Å，可用
  `ContactOptions.Cutoff` 调整）即构成一个**残基接触**。
- 同一对链之间的全部接触聚合为一个**链间界面**。
- `Incomplete`（缺原子）残基始终排除；`NonStandard` 残基默认排除，
  可经 `ContactOptions.IncludeHetAtoms` 放开。
- 结果按 （链， 残基序号， 插入码） 规范排序，接触对方向规范化：
  同一输入重复计算、或同一批结构按不同顺序提交，结论完全一致。

## 并发与重载语义（store 包）

- `Store.Load` 先在锁外完整解析全部输入，成功后在写锁内**原子替换**快照；
  任一输入失败则整体回滚，旧快照继续可用。
- 每次成功重载递增 `Generation`；查询接口返回结果时携带 Generation，
  调用方可识别结论所基于的数据版本，避免误用旧结论。
- 状态机：`EMPTY`（从未成功加载）→ `READY`（可用）→ `DEGRADED`
  （最近一次重载失败，保留旧数据，`Status.LastError` 记录失败原因）；
  再次成功重载后恢复 `READY` 并清空失败原因。
- 失败原因可区分：文件不存在等为 I/O 错误，文件存在但无可用原子为
  `*structure.ParseError`（可用 `errors.As` 断言）。

## 本地验证方法

```bash
# 运行全部测试（含竞态检测）
go test ./... -race

# 查看关键输入与判定依据日志
go test ./structure/ ./store/ -v

# 命令行端到端验证
go run ./cmd/psa -v testdata/demo.pdb
```

测试覆盖：多模型、可变构象、插入码、非标准残基、缺原子、畸形/截断记录、
空文件错误、原子顺序无关性、重复计算一致性、重载替换旧结论、
失败回滚（区分 I/O 错误与解析错误）、并发查询与并发重载（`-race` 验证）。

## 项目结构

```
structure/   结构模型、PDB 解析、接触与界面判定
store/       并发安全的数据存储与在线重载
cmd/psa/     命令行工具
testdata/    演示用结构文件
```
