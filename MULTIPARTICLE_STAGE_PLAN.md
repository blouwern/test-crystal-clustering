# 多粒子聚类测试与算法改进 — 阶段计划

**版本**：2026-10-29  
**状态**：Stage 1 完成；Stage 2 待开始；Stage 6 暂缓  
**目标**：用多个单粒子样本拼接出多粒子样本，系统测试当前聚类算法的多粒子分辨能力，并为后续算法改进建立可比较的基线。

> 本工作不会在一次对话内完成。Stage 2 起，每完成一个阶段结束一次对话，供人工审查；每次对话只做一个阶段。

---

## 1. 背景与约束

当前 `reconstruction.py` 的两阈值波前聚类已在三个**单粒子**样本（`ep_iso` / `gamma_iso` / `mu_rest`）上充分扫描：最优工作点 `Nth=0, Eth_seed=7 MeV, Eth_wave=0.02–0.1 MeV`，exact-match ≈ 0.868，能量捕获 ≈ 0.976。单粒子指标已接近能量接受度地板（≈0.875），无法再区分算法优劣；单粒子样本也**看不到多粒子最主要的失效模式——两个相邻簇射被合并**。

因此需要用多个单粒子事例拼接出多粒子事例，测试：

1. 算法能否给出正确的团簇数 N；
2. 每个真实簇射能否被单独找回（效率/纯度）；
3. 两个簇射在多近的距离内会被合并；
4. 能量是否守恒、是否被错误分配。

**约束**：

* 样本**不需要时间信息**；只使用晶体能量图。
* 现有 processed 文件只有 `Edeps`（568 维 `vector<float>`），单粒子真值由样本构造保证；processed 文件没有 `NTrack`。
* 原始 MC（`data/raw/*.root`）的 `G4Run0/ECALSimHit` 含 `ModID` 与三维击中位置 `x`，可用于计算晶体中心与簇射质心（Stage 3 需要）。
* 拼接假设探测器响应线性、两个粒子互不影响；该假设对"能量图聚类"测试足够，但**不能用于物理产额/效率**。

---

## 2. 仓库结构

几何位置生成、粒子拼接、评估、参数扫描是**同级**的核心模块；评估按关注点拆分为多个源文件：

```text
test-crystal-clustering/
├── MULTIPARTICLE_STAGE_PLAN.md          # 本文件
├── data/
│   ├── processed/                        # 单粒子样本（LFS）
│   ├── raw/                              # 原始 MC（含 x 位置）
│   └── utilities/
│       ├── ecal_neighbor_info.root       # 拓扑（568 晶体）
│       └── ecal_module_positions.npz     # Stage 3 生成：晶体中心
└── pytrial/
    ├── reconstruction.py                 # 现有聚类算法
    ├── scoring.py                        # 现有单粒子评估
    ├── geometry.py                       # 【新】晶体位置 / 质心 / 距离
    ├── stitching.py                      # 【新】多粒子拼接核心
    ├── evaluation/                       # 【新】多粒子评估（按关注点拆分）
    │   ├── __init__.py
    │   ├── counts.py                     # 计数类指标（Stage 1 已实现）
    │   ├── matching.py                   # 真值匹配 / 效率 / 纯度 / 合并 / 分裂（Stage 2/3）
    │   └── energy.py                     # 能量类指标（Stage 2）
    ├── sweeps/                           # 现有参数扫描（与上面同级）
    ├── smoke_tests/
    │   ├── check_evt_info/               # 现有
    │   ├── pre_ana_duplicate_mod_hit/    # 现有
    │   └── simple_evt_stitch/            # 【新】Stage 1 冒烟测试
    │       ├── README.md                 # 测试报告
    │       ├── stitch_smoke.py           # 拼接算法
    │       ├── run_smoke.py              # 重建 + 评估
    │       ├── meta.json                 # 生成配置
    │       ├── data/                     # 测试数据（拼接后的 ROOT）
    │       └── results/                  # 指标、团簇输出、图
    └── results/                          # 现有扫描结果
```

**分层原则**：

* `geometry.py`、`stitching.py`、`evaluation/`、`sweeps/` 平级，互不嵌套；
* `evaluation/` 按“计数 / 匹配 / 能量”拆文件，避免单文件膨胀；
* 冒烟测试只属于 `smoke_tests/`，不污染核心模块；
* 每个阶段的数据、配置、结果放在对应阶段自己的目录下。

---

## 3. 阶段总览

| 阶段 | 方案 | 内容 | 数据 | 主要指标 | 状态 |
|:--|:--|:--|:--|:--|:--|
| **Stage 1** | **A** | 最简计数拼接冒烟测试：N=2，2k 事例 | `ep_iso`，pool 20k（正能量） | 计数类 + I/O/能量守恒 | **本次完成** |
| **Stage 2** | **B** | 逐初级能量矩阵；N=2/3/4；同种 + 混合；各 20k | `ep_iso`/`gamma_iso`/`mu_rest` | 计数 + 效率/纯度/合并/分裂 + 能量 | 待开始 |
| **Stage 3** | **C** | 分离度分箱；N=2；6 个距离箱，每箱 ≥2k | 同上 + 晶体位置 | 合并率 vs 距离 | 待开始 |
| **Stage 4** | **D** | 晶体集不相交对照；N=2 | 同上 | 无重叠上限 | 待开始 |
| **Stage 5** | **E** | Poisson(2–3) + 混合种类，真实感基准 | 同上 | 综合性能 | 待开始 |
| **Stage 6** | **F** | 固定总能量/能量归一，系统检查 | — | 阈值移动检查 | **暂缓，不执行** |

**审查流程**：Stage 2 起，每完成一个阶段结束一次对话；对话中给出该阶段总结、关键结果、复现方式，以及下一阶段的具体内容。Stage 6 暂不排期。

---

## 4. Stage 1：方案 A 冒烟测试（已完成）

### 4.1 目的

只验证“多粒子拼接 → 当前重建 → 计数评估”这条链路能跑通，并检查：

* 拼接后的 ROOT 结构正确（`Edeps` 长度 568、`NTrack=2`、`SourceIdx`）；
* 逐事例能量守恒（拼接后能量 = 两个源事例能量之和）；
* 当前重建能处理多粒子输入；
* 计数指标（exact-2、欠/过分裂、MAE、RMSE 等）能正确计算；
* 结果可复现（固定 RNG、固定 pool、记录源事例索引）。

**不验证**：逐初级真值、效率/纯度、合并/分裂、合并率 vs 距离。这些留给 Stage 2/3。

### 4.2 设置

| 项 | 值 |
|:--|:--|
| 样本 | `ep_iso`（电磁簇射，最容易） |
| 源 pool | 前 20 000 事例中总能量 >0 的事例 |
| 拼接重数 | N = 2（每个事例内不放回抽 2 个源事例） |
| 事例数 | 2 000 |
| 随机种子 | 20261018 |
| 重建参数 | `Nth=0`, `Eth_seed=7.0 MeV`, `Eth_wave=0.05 MeV`（单粒子推荐点） |
| 输出 | `pytrial/smoke_tests/simple_evt_stitch/data/stitched_ep_iso_n2_2k.root` |

### 4.3 算法

```text
读取 ep_iso 前 20k 事例的 Edeps
过滤：只保留 sum(Edeps) > 0 的事例
for 每个拼接事例 (2k):
    idx = rng.choice(pool, size=2, replace=False)
    Edeps_stitched = Edeps[idx[0]] + Edeps[idx[1]]
    NTrack = 2
    SourceIdx = idx
写 ROOT：EdepOfEachEvt/Edeps + NTrack + SourceIdx
```

### 4.4 测试报告

位置：`pytrial/smoke_tests/simple_evt_stitch/README.md`。必须包含：

* 测试内容与原理；
* 数据与参数；
* 复现命令；
* 测试结果（计数指标表、分布图、检查项通过情况）；
* 局限（无逐初级真值，不能测合并/分裂）。

### 4.5 通过标准

* 2 000 个事例全部写入、全部读出；
* 每个事例 `NTrack == 2`、`len(Edeps) == 568`；
* 逐事例能量守恒最大绝对误差 < 1e-3 MeV；
* 重建正常结束，2 000 个事例的团簇数全部可读；
* 计数指标有限、在 [0,1] 内；
* 结果可用固定种子复现。

---

## 5. Stage 2：方案 B — 逐初级能量矩阵

### 5.1 目的

在 Stage 1 的基础上加入**逐初级真值**，使评估从“只看总数”升级为“能测每个簇射是否被单独找回”。

### 5.2 设置

| 项 | 值 |
|:--|:--|
| 重数 | N = 2, 3, 4 |
| 种类 | 同种（ep、gamma、mu 各一组）+ 混合（每初级随机种类） |
| 事例数 | 每组 20 000 |
| 源 pool | 每样本前 200 000 事例（正能量；同种/混合按需抽样） |
| 空间选择 | 随机（不控制重叠），保留自然分离度分布 |
| 输出 | `stitched_N{n}_{species}.root`，每文件含 `Edeps`、`NTrack`、逐初级稀疏真值 |

### 5.3 真值格式（建议）

| 分支 | 类型 | 含义 |
|:--|:--|:--|
| `Edeps` | `vector<float>[568]` | 相加后的总能量（重建用） |
| `NTrack` | `int` | 初级数 N |
| `Species` | `vector<int>[N]` | 0=ep, 1=gamma, 2=mu |
| `SourceIdx` | `vector<int>[N]` | 源事例索引 |
| `TruthCrystal` | `vector<int>` | 有贡献的晶体（稀疏） |
| `TruthPrimary` | `vector<int>` | 该贡献属于第几个初级 |
| `TruthEnergy` | `vector<float>` | 该初级在该晶体的能量 |
| `TruthCentroid` | `vector<float>[3N]` | 逐初级能量加权质心 |
| `Separation` | `float` | 两两质心距离的最小值（mm） |
| `OverlapEnergyFraction` | `float` | 被 ≥2 个初级共享的晶体能量占比 |

### 5.4 评估

* 计数：exact-N、欠/过分裂、MAE、RMSE、bias、mean/std/max；
* 匹配：
  * 构造重叠矩阵 `M[i,j]`（重建团簇 i 中来自真值初级 j 的能量）；
  * 匈牙利算法最大化总重叠；
  * 效率 `ε_j`、纯度 `p_i`；
  * 合并率、分裂率、匹配率；
* 能量：每个真值簇的能量回收、总捕获、未聚类能量。

### 5.5 交付物

* `pytrial/stitching.py` 的 `stitch_with_truth` 实现；
* `pytrial/evaluation/matching.py`、`energy.py` 实现；
* 拼接样本、评估脚本、结果表/图；
* 阶段报告（同 Stage 1 格式）。

### 5.6 通过标准

* 每个文件真值自洽：`Σ Edeps == Σ TruthEnergy`；
* 匹配矩阵每行/列能量守恒；
* 计数与匹配指标有限、在 [0,1]；
* 给出 N=2/3/4 的基线表与图；
* 复现命令写入报告。

---

## 6. Stage 3：方案 C — 分离度分箱（N=2）

### 6.1 目的

把“多粒子分辨”变成一条可比较的曲线：**合并率 vs 两簇射质心距离**。

### 6.2 设置

* 用 `geometry.py` 从原始 MC 的 `G4Run0/ECALSimHit`（`ModID`, `x`）生成 `568×3` 晶体中心表；
* 对每个拼接 N=2 事例，算两个初级的能量加权质心，得到距离（mm）；
* 距离分箱（以晶距 ≈46 mm 为单位）：

| 箱 | 范围（晶距） | 物理含义 |
|:--|:--|:--|
| 0 | 0–0.5 | 几乎完全重叠 |
| 1 | 0.5–1 | 同一晶体邻域 |
| 2 | 1–2 | 相邻晶体 |
| 3 | 2–3 | 隔一块晶体 |
| 4 | 3–5 | 明显分开 |
| 5 | >5 | 基本独立 |

* 每箱 ≥2 000 事例；同种 + 混合各一组；
* 用 Stage 2 的真值格式与匹配指标。

### 6.3 交付物

* `geometry.py` + `data/utilities/ecal_module_positions.npz`；
* `stitching.py` 的 separation-bin 模式；
* 每箱样本、`merge_rate_vs_separation.png`、拟合的“合并半径”；
* 阶段报告。

### 6.4 通过标准

* 位置表覆盖 568 晶体，晶距与邻接表一致（≈46 mm）；
* 每箱统计量 ≥2 000；
* 合并率曲线单调（距离越大合并率越低），并给出 10%/50% 合并对应的距离；
* 复现命令写入报告。

---

## 7. Stage 4：方案 D — 晶体集不相交对照（N=2）

### 7.1 目的

给出“完全不重叠”时的算法上限，作为 Stage 2/3 的对照；不做独立物理结论。

### 7.2 设置

* 要求两个源事例的点亮晶体集合互不重叠（或共享能量 < 10⁻³ MeV）；
* N=2，20 000 事例；同种 + 混合；
* 真值格式与 Stage 2 相同（此方案下自然无共享）。

### 7.3 交付物

* `stitching.py` 的 disjoint 模式；
* 对照样本、计数 + 匹配指标表、与 Stage 2 随机样本的对比图；
* 阶段报告。

### 7.4 通过标准

* 所有事例的晶体集交集为空（或共享能量低于阈值）；
* 给出与随机重叠样本的差异，量化“重叠造成的性能损失”；
* 复现命令写入报告。

---

## 8. Stage 5：方案 E — Poisson 重数 + 混合种类（真实感基准）

### 8.1 目的

在受控的 N=2/3/4 之外，给出一个更接近真实物理的基准：多重数服从 Poisson(λ)，种类混合。

### 8.2 设置

* N ~ Poisson(λ)，λ = 2、3（两个样本）；
* 每个初级随机从 `ep_iso`/`gamma_iso`/`mu_rest` 抽；
* 每个样本 20 000 事例；
* 真值格式与 Stage 2 相同（N 可变，`NTrack` 逐事例）；
* 评估：计数（exact-N 按事例 N 计算）、匹配、能量；按 N 分组的性能表。

### 8.3 交付物

* `stitching.py` 的 Poisson/mixed 模式；
* 样本、按 N 分组的指标表、综合图；
* 阶段报告。

### 8.4 通过标准

* `NTrack` 分布与目标 Poisson 一致（KS 检验不拒绝）；
* 给出按 N 分组的计数/匹配/能量指标；
* 复现命令写入报告。

---

## 9. Stage 6：方案 F — 固定总能量/能量归一（暂缓）

**暂不执行。** 目的：区分“重数增加”与“总能量增加”对阈值的影响；把每个单粒子事例能量乘以 1/N，使拼接后总能量与单粒子相当。若后续需要，再单独排期。

---

## 10. 评估指标定义（跨阶段统一）

### 10.1 计数类

| 指标 | 定义 |
|:--|:--|
| `exact_N` | `n_clusters == N` 的事例比例 |
| `under_rate` | `n_clusters < N` 的比例 |
| `over_rate` | `n_clusters > N` 的比例 |
| `MAE` | `mean(|n_clusters − N|)` |
| `RMSE` | `sqrt(mean((n_clusters − N)²))` |
| `bias` | `mean(n_clusters) − N` |
| `mean_clusters` / `std` / `max` | 分布统计 |
| `histogram` | 团簇数 0..max 的计数 |

### 10.2 匹配类（Stage 2+）

对每个事例：

```text
M[i,j] = 重建团簇 i 中来自真值初级 j 的能量
```

用匈牙利算法最大化总重叠，然后：

| 指标 | 定义 |
|:--|:--|
| 效率 `ε_j` | 初级 j 被匹配到的能量 / 初级 j 总能量 |
| 纯度 `p_i` | 团簇 i 中属于其匹配初级的能量 / 团簇 i 总能量 |
| 合并率 | 两个初级的最佳匹配团簇相同的事例/初级对比例 |
| 分裂率 | 一个初级的能量被分到 ≥2 个团簇（最佳团簇 <90%）的比例 |
| 匹配率 | `ε>50%` 且 `p>50%` 的初级比例 |

### 10.3 能量类（Stage 2+）

总能量捕获、每个真值簇的能量回收、收束率、未聚类能量。

### 10.4 统计

* 每组 ≥20 000 事例；分箱 ≥2 000；
* 计数指标报告二项标准误；
* 不同拼接事例若复用源事例，误差可能相关；主结果使用分块 pool（每块独立）。

---

## 11. 复现与数据

### 11.1 数据可用性

* processed ROOT 由 Git LFS 跟踪；若工作区只有指针文件：

```bash
git lfs fetch --include="data/processed/**" origin
git lfs checkout
```

* 原始 MC（`data/raw/*.root`）用于 Stage 3 的晶体位置；只读 `G4Run0/ECALSimHit` 的 `ModID` 与 `x`。

### 11.2 各阶段命令

Stage 1：

```bash
cd pytrial/smoke_tests/simple_evt_stitch
python3 stitch_smoke.py
python3 run_smoke.py
```

Stage 2 起，命令写入各阶段自己的 README / 阶段报告。

### 11.3 数据与结果存放

* 每阶段样本：`pytrial/smoke_tests/<stage>/data/` 或 `pytrial/results/multiparticle/<stage>/data/`；
* 指标/图：同阶段的 `results/`；
* 生成配置：`meta.json`（seed、pool、N、种类、分箱、代码版本）。

> 约定：Stage 2 起的正式样本放在 `pytrial/results/multiparticle/` 下，按阶段分子目录；Stage 1 冒烟测试按用户要求放在 `pytrial/smoke_tests/simple_evt_stitch/`。

---

## 12. 风险与假设

1. **线性叠加假设**：只对能量图聚类有效，不用于物理产额。
2. **共享晶体**：必须逐初级记账；否则真值错误（Stage 1 不涉及）。
3. **阈值偏移**：多粒子总能量 ≈ N×，最优 `Eth_seed`/`Eth_wave` 可能移动；Stage 2 起在拼接样本上重扫。
4. **合并 vs 漏掉**：软簇射可能因低于种子阈值而未成簇，不能全部算作合并；用效率/纯度区分。
5. **统计相关**：拼接事例可能复用源事例；主结果用分块 pool。
6. **无时间信息**：同一晶体上的重叠簇射在能量上不可区分；这是本测试的边界，不是缺陷。
7. **边缘效应**：质心靠近探测器边缘的事例接受度不同；必要时加 fiducial cut 或单独报告。

---

## 13. 当前状态

* **Stage 1（A）**：✅ 已完成（2026-10-29）。2 000 个 N=2 拼接事例；**exact-2 = 0.8575**，under = 0.1260，over = 0.0165，MAE = 0.1430；能量守恒最大误差 1.53×10⁻⁵ MeV；所有检查通过。报告：`pytrial/smoke_tests/simple_evt_stitch/README.md`。
* **Stage 2（B）**：待开始。下一次对话的内容：实现逐初级真值拼接 + 匹配/能量评估，生成 N=2/3/4 同种/混合样本，给出第一版多粒子基线。
* **Stage 3–5**：待排期。
* **Stage 6（F）**：暂缓。

---

*本计划随阶段推进更新；每个阶段结束时在 §13 记录状态与关键结果。*
