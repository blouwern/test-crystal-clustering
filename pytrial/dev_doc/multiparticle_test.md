# 多粒子测试：运行流程与数据格式

**文档版本**：2026-11-12  
**对应代码**：`pytrial/reconstruction.py`、`pytrial/stitching.py`、`pytrial/geometry.py`、`pytrial/evaluation/`、`pytrial/results/multiparticle/stage2/`、`pytrial/results/multiparticle/stage3/`  
**对应报告**：[Stage 2 报告](../results/multiparticle/stage2/README.md)、[Stage 3 报告](../results/multiparticle/stage3/README.md)、[多粒子阶段计划](../../MULTIPARTICLE_STAGE_PLAN.md)

本文说明多粒子测试从**单粒子样本**到**重建**、再到**评估结果**的完整运行流程，并逐个给出数据文件、ROOT 分支、内存结构和函数接口的格式。

## 目录

0. 一句话概括
1. 全景流程
2. 重建算法的输入契约
3. 数据层 1：单粒子样本
4. 数据层 2：拼接样本
5. Stage 2：样本生成流程
6. Stage 2：评估流程
7. Stage 3：样本生成流程
8. Stage 3：评估流程
9. 重建算法内部流程
10. 评估接口速查
11. 复现与运行
12. 不变量、验收检查与常见坑
13. 完整示例：一个 2γ 事例从输入到指标
14. 关键文件清单
15. 参考结果与量级
16. 小结

---

## 0. 一句话概括

多粒子测试的本质是：

```text
单粒子 Edeps[568]  ──按初级相加──▶  多粒子 Edeps[568] + 逐初级稀疏真值 Truth*
                                              │
                                              ▼
                              reconstruction.reconstruction()
                              只读 Edeps + 邻接表，输出“每个团簇包含哪些 ModID”
                                              │
                                              ▼
                              评估层用 Truth* 反算每个团簇中的
                              信号能量 / 本底能量，再做计数、匹配、合并/分裂、能量与本底指标
```

因此，**重建算法只依赖两样输入**：`Edeps`（568 维晶体能量向量）和 `ecal_neighbor_info.root`（晶体邻接表）。`NTrack`、`Species`、`Truth*`、`Bkg*` 等分支都是评估层使用的真值/记账信息，`reconstruction.py` 完全不读。

---

## 1. 全景流程

### 1.1 数据流总图

```text
data/raw/{ep_iso,gamma_iso}.root
  G4Run0/ECALSimHit: EvtID, ModID, Edep, t, x[3]
        │
        │  data/read_edep_of_each_evt.C
        │  （按 (EvtID, ModID) 累加 Edep，逐事件填成长度 568 的向量）
        ▼
data/processed/edep_of_each_evt_{ep_iso,gamma_iso}.root
  EdepOfEachEvt/Edeps : vector<float>[568]
        │
        ├──────────────────────────────────────────────┐
        │ Stage 2                                      │ Stage 3
        │ generate.py                                  │ generate.py
        │  · 随机抽 N 个源事例逐晶体相加                 │  · 先算每个 gamma 事例的能量加权质心
        │  · 每个拼接事例抽一次 γ 本底                   │  · 按质心距离分 6 个 separation bin
        │  · 写 _sig / _bkg 两个版本                    │  · 每箱 2500 对，写 _sig / _bkg
        ▼                                              ▼
pytrial/results/multiparticle/stage2/data/      pytrial/results/multiparticle/stage3/data/
  stitched_2gamma_{sig,bkg}.root                  stage3_2gamma_{sig,bkg}.root
  stitched_ep2gamma_{sig,bkg}.root                （额外有 Separation / SepBin）
        │                                              │
        └──────────────────────┬───────────────────────┘
                               │
                               ▼
                 <work_dir>/processed/<sample>.root   （符号链接）
                 <work_dir>/utilities/ecal_neighbor_info.root
                               │
                               ▼
                 reconstruction.reconstruction(...)
                 读 Edeps + neighbors，波前传播聚类
                               │
                               ▼
                 results/clusters_<tag>.root
                 ECALClusters：EvtID/I + ModIDList(vector<int>)
                               │
                               ▼
                 evaluation.io.load_clusters()
                 evaluation.io.load_multiparticle_sample()
                               │
          ┌────────────────────┼────────────────────┐
          ▼                    ▼                    ▼
   matching.py           background.py          energy.py
   （重叠矩阵/指派/      （假团簇/纯度/         （信号捕获/
     效率/合并/分裂）      本底捕获）             收束/逐初级响应）
          │                    │                    │
          └────────────────────┼────────────────────┘
                               ▼
                    counts.py / separation.py
                               │
                               ▼
        stage2_metrics.csv + points/*.json + best_points.json + ...
        stage3_metrics.csv + fits.json + merge_rate_*.png
```

一次完整运行按时间顺序是：

1. 上游：`data/read_edep_of_each_evt.C` 从原始 MC 生成 `data/processed/edep_of_each_evt_{ep_iso,gamma_iso}.root`；`ecal_neighbor_info.root` 由外部几何流程提供；
2. Stage 2 生成：`stage2/generate.py` 读这两个单粒子样本，构造 `2gamma` / `ep2gamma`，抽本底，写 `stitched_*_{sig,bkg}.root`；
3. Stage 2 评估：`stage2/evaluate.py` 建 `work/` 符号链接，跑 112 个重建点，读回 `ECALClusters` + `Truth*`，写 CSV/JSON/图；
4. Stage 3 生成：`stage3/generate.py` 用晶体位置表算质心，按分离度箱构造 `2gamma` 对，写 `stage3_2gamma_{sig,bkg}.root`；
5. Stage 3 评估：`stage3/evaluate.py` 跑 28 个重建点，逐箱统计合并率，拟合 50%/10% 半径。

### 1.2 三个阶段的分工

| 项目 | Stage 1（冒烟测试） | Stage 2（主样本） | Stage 3（分离度曲线） |
|:--|:--|:--|:--|
| 位置 | `pytrial/smoke_tests/simple_evt_stitch/` | `pytrial/results/multiparticle/stage2/` | `pytrial/results/multiparticle/stage3/` |
| 目标 | 验证“拼接 → 重建 → 计数”链路 | 测多粒子总体性能、本底影响 | 测合并率随两簇射距离的变化 |
| 信号 | `ep_iso`，N=2 | `2gamma`（N=2）、`ep2gamma`（N=3） | `2gamma`（N=2） |
| 真值 | 只有 `NTrack` + `SourceIdx`，无逐初级信息 | `Truth*` 逐初级 + `Bkg*` | 同 Stage 2，另加 `Separation` / `SepBin` |
| 事例选择 | 源 pool 只保留 `ΣEdeps > 0` | **不做任何沉积能量筛选** | 源事例要求 `ΣEdeps > 0`（质心需要） |
| 每个拼接事例的本底 | 无 | 抽一次 | 抽一次 |
| 评估 | `evaluation/counts.py` | 计数 + 匹配 + 本底 + 能量 | 逐分离度箱的合并/分裂/计数 |
| 重建点数 | 1 | 112 | 28 |

> Stage 1 是历史冒烟测试，代码仍保留在 `stitching.py` 的 `stitch_count_only()` / `write_stitched()` 中；Stage 2 起改用 `write_multiparticle_sample()` 的逐初级真值格式。本文重点讲 Stage 2 / Stage 3。

### 1.3 `pytrial/results/multiparticle/` 目录

```text
pytrial/results/multiparticle/
├── stage2/
│   ├── generate.py                 # 生成 2γ / e+2γ 的 sig/bkg 样本
│   ├── evaluate.py                 # 112 点重建 + 评估
│   ├── run_stage2.sh               # generate.py → evaluate.py
│   ├── meta.json                   # 生成参数、能量守恒检查
│   ├── data/
│   │   ├── stitched_2gamma_sig.root
│   │   ├── stitched_2gamma_bkg.root
│   │   ├── stitched_ep2gamma_sig.root
│   │   └── stitched_ep2gamma_bkg.root
│   ├── work/                       # 运行时生成：processed/ + utilities/ 符号链接
│   └── results/
│       ├── stage2_metrics.csv      # 112 行 × 35 列
│       ├── points/*.json           # 每个参数点的 row + hist
│       ├── best_points.json
│       ├── background_impact.json
│       ├── clusters_*.root         # 运行时生成；当前 results/ 未保留
│       └── *.png
└── stage3/
    ├── generate.py                 # 分离度分箱样本
    ├── evaluate.py                 # 28 点重建 + 逐箱评估
    ├── run_stage3.sh
    ├── meta.json                   # pitch、箱边界、能量/本底统计
    ├── data/
    │   ├── stage3_2gamma_sig.root
    │   └── stage3_2gamma_bkg.root
    ├── work/                       # 运行时生成
    └── results/
        ├── stage3_metrics.csv      # 168 行 × 18 列
        ├── fits.json               # 28 条合并半径
        ├── clusters_*.root         # 运行时生成；当前 results/ 未保留
        └── *.png
```

> `clusters_*.root` 是 `evaluate.py` 运行时写出的中间产物；仓库里的 `results/` 只保留了指标、JSON 和图。重新运行 `run_stage2.sh` / `run_stage3.sh` 会重新生成它们。

---


## 2. 重建算法的输入契约（最重要的接口）

### 2.1 目录布局

`reconstruction.reconstruction()` 只接受一个数据目录 `data_dir`，并在其中按固定相对路径找文件：

```text
<data_dir>/
├── processed/
│   └── <filename_edep_of_each_evt>     # ROOT 文件，树 EdepOfEachEvt
│       └── Edeps : vector<float>       # 长度 nMod（当前 568），单位 MeV
└── utilities/
    └── <neighbor_file_name>            # 默认 ecal_neighbor_info.root
        └── 树 ECALCrystalNeighbors
            ├── neighbors        : vector<int>   # 必须
            └── second_neighbors : vector<int>   # 可选，当前仅加载不使用
```

对应代码：

```python
input_path = Path(data_dir) / "processed" / filename_edep_of_each_evt
neighbor_file = ROOT.TFile.Open(
    str(Path(data_dir) / "utilities" / neighbor_file_name), "READ"
)
```

`data_dir` 的默认值是环境变量 `ECAL_CLUSTERING_DATA_DIR`，未设置时为 `"."`：

```python
if data_dir is None:
    data_dir = os.environ.get("ECAL_CLUSTERING_DATA_DIR", ".")
```

### 2.2 `Edeps` 的格式

| 项 | 值 |
|:--|:--|
| 树名 | `EdepOfEachEvt` |
| 分支名 | `Edeps` |
| ROOT 类型 | `vector<float>`（`std::vector<float>`） |
| 长度 | `nMod`，当前 568；第 `i` 个元素对应 `ModID = i` 的晶体 |
| 单位 | MeV |
| 索引 | 0-based，直接等于 `G4Run0/ECALSimHit.ModID` |
| 缺失晶体 | 用 0.0 填充 |
| 同一晶体多次击中 | 在写 `Edeps` 前已经按 `(EvtID, ModID)` 累加 |
| 零沉积事例 | 保留，`Edeps` 全 0；重建时不会产生任何团簇 |

`reconstruction` 内部会把输入向量补齐/截断到 `nMod`：

```python
nMod = neighbor_tree.GetEntries()       # 不是硬编码 568
edeps_full = [0.0] * nMod
length = min(len(edeps), nMod)
for i in range(length):
    edeps_full[i] = edeps[i]
```

因此：

- 输入向量比 `nMod` 短：尾部补 0；
- 输入向量比 `nMod` 长：尾部截断；
- `nMod` 由邻接表的条目数决定，而不是由 `Edeps` 的长度决定。

### 2.3 邻接表的格式

| 项 | 值 |
|:--|:--|
| 树名 | `ECALCrystalNeighbors` |
| 条目数 | 568（即 `nMod`） |
| 条目 `i` 的含义 | 晶体 `i` 的邻居列表 |
| `neighbors` | `vector<int>`，直接相邻的 `ModID` 列表 |
| `second_neighbors` | `vector<int>`，二阶邻居；当前文件存在，但重建中的 gap-jumping 回退被注释掉，因此不使用 |
| 单位/坐标系 | 不涉及坐标；只用拓扑 |

`reconstruction` 读取方式：

```python
neighbor_tree = neighbor_file.Get("ECALCrystalNeighbors")
nMod = neighbor_tree.GetEntries()
neighbor_vec = vector("int")()
neighbor_tree.SetBranchAddress("neighbors", neighbor_vec)

# second_neighbors 是可选的：没有该分支时不崩溃，只打印一条 INFO
has_second = any(b.GetName() == "second_neighbors"
                 for b in neighbor_tree.GetListOfBranches())
```

> 当前 `data/utilities/ecal_neighbor_info.root` 同时有 `neighbors` 和 `second_neighbors`，568 条。`second_neighbors` 的生成脚本在 `data/utilities/add_secondary_neighbors_deprecated.py`，已标记 deprecated；重建中的回退逻辑也被注释，实际不参与聚类。

### 2.4 重建输出契约

| 项 | 值 |
|:--|:--|
| 文件名 | 由 `output_path` 指定；`None` 时在当前工作目录生成 `ecal_clusters_Eth_s-..._Eth_w-..._Nth-....root` |
| 树名 | `ECALClusters` |
| 分支 `EvtID` | `int`，输入事件序号（0-based） |
| 分支 `ModIDList` | `vector<int>`，该团簇包含的晶体 `ModID` |
| 条目粒度 | **一个团簇一条 entry**，不是“一个事件一条” |
| 零团簇事件 | **不写任何 entry**；读取时必须显式提供总事件数 |
| 团簇内能量 | **不存储**；评估时用 `Truth*` 对该团簇的 `ModIDList` 求和得到 |
| 团簇顺序 | 按种子发现顺序；`ModIDList[0]` 是种子，后续是波前扩展顺序，不是按能量排序 |

输出 schema：

```text
ECALClusters
├── EvtID/I            event index
└── ModIDList          vector<int>  crystal IDs in this cluster
```

这个契约是 Stage 2 / Stage 3 评估能够工作的关键：**重建只负责“哪些晶体属于同一个团簇”**，不负责“这个团簇里有多少信号能量、多少本底能量”。后者由评估层用逐初级真值计算。

---

## 3. 数据层 1：单粒子样本

### 3.1 原始 MC 文件

原始文件位于：

```text
data/raw/ep_iso.root
data/raw/gamma_iso.root
data/raw/mu_rest.root
data/raw/mu2egg_rest.root
```

多粒子测试只用到 `ep_iso` 和 `gamma_iso`。相关树/分支：

| 树 | 分支 | 类型 | 用途 |
|:--|:--|:--|:--|
| `G4Run0/ECALSimHit` | `EvtID` | `Int_t` | 事件号 |
| `G4Run0/ECALSimHit` | `ModID` | `Short_t` | 晶体编号（0-based） |
| `G4Run0/ECALSimHit` | `Edep` | `Float_t` | 该次 hit 的沉积能量（MeV） |
| `G4Run0/ECALSimHit` | `t` | `Double_t` | 击中时间（ns） |
| `G4Run0/ECALSimHit` | `x[3]` | `Float_t[3]` | 击中位置（mm），Stage 3 用来算晶体中心/质心 |
| `G4Run0/SimPrimaryVertex` | `PDGID` | `Int_t` | 初级粒子种类：`-11` = e+，`22` = γ |
| `G4Run0/SimPrimaryVertex` | `Ek0` | `Float_t` | 初级初始动能；`ep_iso` / `gamma_iso` 均为 53.085 MeV |

### 3.2 `raw → processed`：C++ 宏

当前正式使用的转换脚本是 [data/read_edep_of_each_evt.C](../../data/read_edep_of_each_evt.C)。它做四件事：

1. 用 `ROOT::RDataFrame` 读 `G4Run0/ECALSimHit`；
2. `nMod = max(ModID) + 1`，`nEvt = max(EvtID) + 1`；
3. 用 `map<pair<evtID, modID>, float>` 把同一事件、同一晶体的所有 hit 累加；
4. 对每个事件填一个长度 `nMod` 的 `std::vector<float>`，写入 `EdepOfEachEvt/Edeps`。

关键代码：

```cpp
auto nMod{ static_cast<short>(*df.Max("ModID")+1) };
auto nEvt{ static_cast<int>(*df.Max("EvtID")+1) };

std::map<std::pair<int, short>, float> evt_mod_edep_map;
df.Foreach([&](int evtID, short modID, float edep){
    evt_mod_edep_map[{evtID, modID}] += edep;
}, {"EvtID","ModID","Edep"});
```

输出：

```text
EdepOfEachEvt
└── Edeps : vector<float>[568]      # 每个事件一条，单位 MeV
```

特点：

- **没有** `NTrack`、`Species`、`Truth*`：单粒子样本本身没有多粒子真值，真值由拼接阶段构造；
- 每个事件都有一条长度 568 的向量，零沉积事件是 568 个 0；
- 同一晶体的多个 hit 已合并；
- 不做时间筛选。

> `pytrial/_deprecated_data_processing/read_evt_edeps.py` 是旧版 Python 实现，额外做了 `filter_timeout_hit`（同一晶体重复 hit 的时间挑选）。该目录已标记 deprecated，当前 Stage 2/3 的输入不是它生成的。

### 3.3 `data/processed/` 的现状

当前目录：

```text
data/processed/
├── edep_of_each_evt_ep_iso.root        # 1,000,000 事件
├── edep_of_each_evt_gamma_iso.root     # 1,000,000 事件
├── edep_of_each_evt_mu_rest.root       # 1,000,000 事件
├── edep_of_each_evt_ep_iso_bkg.root
├── edep_of_each_evt_gamma_iso_bkg.root
└── edep_of_each_evt_mu_rest_bkg.root
```

其中：

- 无 `_bkg` 的文件：只有 `Edeps`，是 Stage 2/3 的**信号来源**；
- `_bkg` 文件：由 [bkg_est_add/add_background.py](../bkg_est_add/add_background.py) 生成，`Edeps` 已包含本底，且多一个 `NBkg` 分支。**Stage 2/3 不使用这些单粒子 `_bkg` 文件**，原因见 §5.6。

单粒子 `_bkg` 文件的 schema 是：

```text
EdepOfEachEvt
├── Edeps      : vector<float>[568]   # 信号 + 本底
├── NBkg       : Int_t                # 注入的本底 hit 数
└── BkgEdeps   : vector<float>[568]   # 可选（--store-bkg 时才有）；纯本底分量
```

文件内还有一个名为 `bkg_model` 的 `TNamed`，内容是 JSON 字符串，记录本底模型参数、随机种子和统计摘要。Stage 2/3 不用它，只从 `bkg_estimation.json` 读 `n_bkg`，再用 `DepositionSpectrum` + `draw_injections` 在拼接事例上重新抽样。

Stage 2/3 用到的 `ep_iso` / `gamma_iso` processed 文件，`Edeps` 长度都是 568。实际检查：

```text
data/processed/edep_of_each_evt_ep_iso.root    entries=1,000,000  Edeps len=568
data/processed/edep_of_each_evt_gamma_iso.root entries=1,000,000  Edeps len=568
```

> 注意：`data/processed/edep_of_each_evt_gamma_iso.root` 内部有同名 `EdepOfEachEvt` 的多个 cycle；`TFile::Get("EdepOfEachEvt")` 返回最高 cycle。消费者不需要关心 cycle 号。

### 3.4 `E_sig`：单粒子沉积信号能量

对每个单粒子事例：

```text
E_sig = Σ_i Edeps[i]        # MeV
```

它不等于初级初始能量 53.085 MeV；它包含探测器泄漏/沉积涨落：

- `gamma_iso` 事例 0：`ΣEdeps = 0.000 MeV`（零沉积）；
- `gamma_iso` 事例 1：`ΣEdeps = 53.085 MeV`（完全沉积，13 个晶体非零）；
- `gamma_iso` 事例 2：`ΣEdeps = 53.085 MeV`（9 个晶体非零）。

Stage 2 的“all-event”结果保留这些零沉积事例；Stage 3 因为要算质心，只保留 `E_sig > 0` 的源事例。

### 3.5 晶体位置表（Stage 3 专用）

`reconstruction` 不使用位置；Stage 3 需要它来算质心。文件：

```text
data/utilities/ecal_module_positions.npz
├── centre : float64[568, 3]    # 每个晶体的平均 MC hit 位置（mm）
└── counts : int64[568]         # 每个晶体用于平均的 hit 数；0 表示从未被击中
```

由 [geometry.py](../geometry.py) 的 `read_module_positions()` 生成：

```python
df = ROOT.RDataFrame("G4Run0/ECALSimHit", raw_path).Range(0, n_hits)
data = df.AsNumpy(["ModID", "x"])
mod = np.asarray(data["ModID"], dtype=np.int64)
pts = np.stack([np.asarray(v, dtype=np.float64) for v in data["x"]])

acc = np.zeros((n_mod, 3), dtype=np.float64)
cnt = np.zeros(n_mod, dtype=np.int64)
np.add.at(acc, mod, pts)
np.add.at(cnt, mod, 1)
centre = np.divide(acc, cnt[:, None], out=np.zeros_like(acc),
                   where=cnt[:, None] > 0)
```

当前 `n_mod=568`，`n_hits=500_000`，568/568 个模块都有 hit；`estimate_pitch()` 给出最近邻距离中位数 **46.2409 mm**。

---

## 4. 数据层 2：拼接样本

Stage 2 / Stage 3 写出的文件使用同一个 schema，由 [stitching.write_multiparticle_sample()](../stitching.py) 写入。树名统一为 `EdepOfEachEvt`。

### 4.1 分支总表

| 分支 | ROOT 类型 | 长度/形状 | 写入者 | 读取者 | 含义 |
|:--|:--|:--|:--|:--|:--|
| `Edeps` | `vector<float>` | 568 | `write_multiparticle_sample` | `reconstruction` | 信号 + 本底的总能量图；重建唯一使用的分支 |
| `NTrack` | `Int_t` | 1 | 同上 | `io.load_multiparticle_sample` | 信号粒子数 N（2 或 3） |
| `Species` | `vector<int>` | N | 同上 | `io` | 每个初级的种类：0=e+，1=γ |
| `SourceIdx` | `vector<int>` | N | 同上 | `io`（记账） | 每个初级对应的源事例索引 |
| `PrimaryEnergy` | `vector<float>` | N | 同上 | `io` | 每个初级的初始能量（固定 53.085 MeV） |
| `SourceEnergy` | `vector<float>` | N | 同上 | `io` | 每个初级在探测器中的沉积能量 `E_sig`（可为 0） |
| `SignalEnergy` | `Float_t` | 1 | 同上 | `io` | 所有初级沉积能量之和 |
| `BkgEnergyTotal` | `Float_t` | 1 | 同上 | `io` | 该事例注入的本底总能量 |
| `NBkg` | `Int_t` | 1 | 同上 | `io` | 本底 hit 数（Poisson 抽样） |
| `TruthCrystal` | `vector<int>` | 稀疏，长度不定 | 同上 | `io` | 逐初级真值的晶体编号 |
| `TruthPrimary` | `vector<int>` | 与 `TruthCrystal` 等长 | 同上 | `io` | 该条真值属于哪个初级（0..N-1） |
| `TruthEnergy` | `vector<float>` | 与 `TruthCrystal` 等长 | 同上 | `io` | 该初级在该晶体的信号能量 |
| `BkgCrystal` | `vector<int>` | 稀疏，长度 = NBkg | 同上 | `io` | 本底 hit 的晶体编号 |
| `BkgEnergy` | `vector<float>` | 与 `BkgCrystal` 等长 | 同上 | `io` | 本底 hit 的能量 |
| `Separation` | `Float_t` | 1（仅 Stage 3） | 同上 | `io` | 两个初级质心距离（mm） |
| `SepBin` | `Int_t` | 1（仅 Stage 3） | 同上 | `io` | 分离度箱编号 0..5 |

### 4.2 稀疏真值的语义

`TruthCrystal` / `TruthPrimary` / `TruthEnergy` 是**三个平行数组**，每条 entry 表示“某个初级在某个晶体上的信号贡献”：

```text
TruthCrystal[k] = ModID
TruthPrimary[k] = primary index (0..N-1)
TruthEnergy[k]  = energy (MeV)
```

例如：

```text
TruthCrystal = [64, 78, 172, 280, 523, ..., 449, 55, ...]
TruthPrimary = [ 0,  0,   0,   0,   0, ...,   1,  1, ...]
TruthEnergy  = [0.93, 0.40, 1.09, 0.35, 13.42, ..., ...]
```

关键点：

- **共享晶体**会以多条 entry 出现：晶体 `c` 同时被初级 0 和 1 击中时，会分别有 `(c, 0, E0)` 和 `(c, 1, E1)` 两条；这是合并/分裂诊断的基础。
- `TruthEnergy` 对同一初级求和应等于 `SourceEnergy[primary]`；对所有 entry 求和应等于 `SignalEnergy`。
- 零沉积的初级没有任何 entry（`SourceEnergy=0`）。
- `BkgCrystal` / `BkgEnergy` 是独立的本底真值，允许同一晶体重复出现；`io` 会把重复项累加。

对任意晶体 `c`，三种表示满足：

```text
Edeps[c] ≈ Σ_j T[j, c] + B[c]          # 重建看到的唯一能量图（float32 舍入内）
T[j, c]  = 初级 j 在晶体 c 的信号能量   # 由 Truth* 稀疏三元组展开
B[c]     = 本底在晶体 c 的能量          # 由 Bkg* 稀疏二元组展开
```

重建只看到 `Edeps`，因此**同一晶体上的信号和本底在重建阶段不可区分**；评估层用 `Truth*` / `Bkg*` 事后把每个团簇的能量分解成信号和本底。这是整个多粒子测试能够同时报告“合并/分裂”和“本底污染”的前提。

### 4.3 `SourceIdx` 的池语义

`SourceIdx` 不是全局唯一 ID，而是“在对应源样本前 100,000 个事件中的下标”：

- `2gamma`：两个 `SourceIdx` 都指向 `gamma_iso` 的 pool；
- `ep2gamma`：`SourceIdx[0]` 指向 `ep_iso` 的 pool，`SourceIdx[1]` / `[2]` 指向 `gamma_iso` 的 pool；
- 不同场景、不同文件之间的相同数字没有可比性；
- 跨拼接事例允许复用同一个源事例（有放回抽样）；同一拼接事例内部，两个 gamma 源保证不同（重采样直到不相等）。

### 4.4 文件级不变量

对每个拼接文件应满足：

```text
len(Edeps) = 568
len(Species) = len(SourceIdx) = len(PrimaryEnergy) = len(SourceEnergy) = NTrack
Σ TruthEnergy            ≈ SignalEnergy                    （float32/float64 容差 ~1e-5 MeV）
Σ SourceEnergy           ≈ SignalEnergy
Σ_{primary=p} TruthEnergy ≈ SourceEnergy[p]
Σ Edeps (sig)             ≈ SignalEnergy
Σ Edeps (bkg)             ≈ SignalEnergy + BkgEnergyTotal
len(BkgCrystal) = len(BkgEnergy) = NBkg
```

Stage 2 的 `meta.json` 记录了实测最大偏差：

| 场景 | `ΣTruthEnergy − SignalEnergy` 最大绝对偏差 | `ΣEdeps_sig − SignalEnergy` 最大绝对偏差 | `ΣEdeps_bkg − (SignalEnergy+BkgEnergyTotal)` 最大绝对偏差 |
|:--|--:|--:|--:|
| `2gamma` | 7.57×10⁻⁶ MeV | 1.53×10⁻⁵ MeV | 3.05×10⁻⁵ MeV |
| `ep2gamma` | 7.63×10⁻⁶ MeV | 3.05×10⁻⁵ MeV | 3.05×10⁻⁵ MeV |

这些偏差来自 float32 累加，远小于任何物理分辨率，不影响指标。

### 4.5 一个具体事例（Stage 2，`2gamma`，事件 0）

从 `pytrial/results/multiparticle/stage2/data/stitched_2gamma_bkg.root` 读出：

| 分支 | 值 |
|:--|:--|
| `NTrack` | 2 |
| `Species` | `[1, 1]` |
| `SourceIdx` | `[89149, 72614]` |
| `PrimaryEnergy` | `[53.085, 53.085]` MeV |
| `SourceEnergy` | `[45.7241, 53.0849]` MeV |
| `SignalEnergy` | 98.8090 MeV |
| `NBkg` | 2 |
| `BkgEnergyTotal` | 9.4273 MeV |
| `Σ Edeps` | 108.2363 MeV |
| `Edeps` 非零晶体数 | 14 |
| `Truth` entry 数 | 12 |
| `Σ TruthEnergy` | 98.8090 MeV |
| 前 5 条 truth | `(64,0,0.9336)`, `(78,0,0.4018)`, `(172,0,1.0866)`, `(280,0,0.3452)`, `(523,0,13.4171)` |

`TruthCrystal` 中初级 0 的晶体和初级 1 的晶体没有重叠；实际用 seed=13、wave=0.05、Nth=0 重建后得到两个团簇：

```text
cluster 0: [524, 64, 78, 172, 523, 280]   c_sig = 45.7241 MeV
cluster 1: [449, 55, 82, 163, 298, 564]   c_sig = 53.0849 MeV
```

`c_bkg` 两个团簇都是 0：该事例的两个本底 hit 没有落在任何团簇内。评估结果见 §13。

---

## 5. Stage 2：样本生成流程

### 5.1 入口与常量

```text
pytrial/results/multiparticle/stage2/
├── generate.py
├── evaluate.py
├── run_stage2.sh          # python3 generate.py && python3 evaluate.py
├── data/                  # 生成的 ROOT
├── meta.json              # 生成参数与自洽检查
└── results/               # 指标、点、图
```

`generate.py` 的关键常量：

```python
HERE = Path(__file__).resolve().parent
PYTRIAL = HERE.parents[2]          # .../pytrial
REPO = HERE.parents[3]             # 仓库根目录
POOL_EVENTS = 100_000
N_EVENTS = 20_000
SEED = 20261106
PRIMARY_ENERGY_MEV = 53.085
N_MOD = 568
```

### 5.2 运行步骤

```text
run_stage2.sh
  └─ python3 generate.py
       ├─ 读 bkg_est_add/bkg_estimation.json 得到 n_bkg
       ├─ 构造 DepositionSpectrum（N(5,2²) 截断于 0）
       ├─ stitching.load_edeps() 读 ep_iso / gamma_iso 前 100k
       ├─ 场景 2gamma：抽两个不同的 gamma 源，逐晶体相加
       ├─ 场景 ep2gamma：抽一个 e+ 源 + 两个不同的 gamma 源
       ├─ 每个场景抽一次本底，写 _sig / _bkg 两个文件
       └─ 写 meta.json
  └─ python3 evaluate.py
       ├─ 建立 work/processed + work/utilities 符号链接
       ├─ 读 _bkg 样本的 Truth*
       ├─ 对 2×2×7×2×2 = 112 个点调用 reconstruction
       ├─ 逐点评估并写 points/*.json
       └─ 汇总 CSV / best_points / background_impact / 图
```

### 5.3 读源样本

```python
ep = stitching.load_edeps(
    REPO / "data/processed/edep_of_each_evt_ep_iso.root", n_events=POOL_EVENTS)
ga = stitching.load_edeps(
    REPO / "data/processed/edep_of_each_evt_gamma_iso.root", n_events=POOL_EVENTS)
e_ep = stitching.event_energies(ep)
e_ga = stitching.event_energies(ga)
```

`stitching.load_edeps()` 的实现：

```python
with ROOT.TFile.Open(str(path), "READ") as f:
    tree = f.Get(tree_name)
    total = tree.GetEntries()
    n = total if n_events is None else min(int(n_events), total)
    edeps = vector("float")()
    tree.SetBranchAddress("Edeps", edeps)
    out = np.zeros((n, n_mod), dtype=np.float32)
    for i in range(n):
        tree.GetEntry(i)
        m = min(edeps.size(), n_mod)
        out[i, :m] = np.asarray(edeps, dtype=np.float32)[:m]
return out
```

返回：

```text
ep: float32[100000, 568]
ga: float32[100000, 568]
e_ep: float64[100000]      # ΣEdeps per event
e_ga: float64[100000]
```

内存量：`100000 × 568 × 4 B ≈ 227 MB` 每个矩阵；两个矩阵约 454 MB。Stage 2 生成阶段是内存峰值最高的地方。

### 5.4 场景 `2gamma`（N=2）

```python
rng = np.random.default_rng([SEED, 1])
pool = np.arange(len(e_ga), dtype=np.int64)

i1 = rng.choice(pool, size=N_EVENTS, replace=True)
i2 = rng.choice(pool, size=N_EVENTS, replace=True)
same = i1 == i2
while same.any():
    i2[same] = rng.choice(pool, size=int(same.sum()), replace=True)
    same = i1 == i2

sig_2g = ga[i1] + ga[i2]                     # float32[20000,568]
species_2g = np.ones((N_EVENTS, 2), dtype=np.int32)          # [1,1]
src_2g = np.stack([i1, i2], axis=1).astype(np.int32)
src_e_2g = np.stack([e_ga[i1], e_ga[i2]], axis=1).astype(np.float32)
sig_e_2g = src_e_2g.sum(axis=1)              # float32[20000]
truth_2g = build_truth([ga[i1], ga[i2]])
```

要点：

- `pool` 是**全部 100k** `gamma_iso` 事例，不做 `E_sig` 筛选；
- 跨事例有放回；同一事例内 `i1 != i2`（`while` 重采样）；
- `sig_2g` 是两个源事例的逐晶体和；线性叠加假设；
- `species` 全 1（gamma）；
- `src_e_2g` 的每一列是该源的 `E_sig`，允许为 0；
- `build_truth` 生成逐初级稀疏真值。

### 5.5 场景 `ep2gamma`（N=3）

```python
rng = np.random.default_rng([SEED, 2])
ie = rng.choice(ep_pool, size=N_EVENTS, replace=True)
j1 = rng.choice(ga_pool, size=N_EVENTS, replace=True)
j2 = rng.choice(ga_pool, size=N_EVENTS, replace=True)
# 同样重采样直到 j1 != j2

sig_ep2g = ep[ie] + ga[j1] + ga[j2]
species_ep2g = np.tile(np.array([0, 1, 1], dtype=np.int32), (N_EVENTS, 1))
src_ep2g = np.stack([ie, j1, j2], axis=1).astype(np.int32)
src_e_ep2g = np.stack([e_ep[ie], e_ga[j1], e_ga[j2]], axis=1).astype(np.float32)
sig_e_ep2g = src_e_ep2g.sum(axis=1)
truth_ep2g = build_truth([ep[ie], ga[j1], ga[j2]])
```

要点：

- `SourceIdx[0]` 是 `ep_iso` 的下标，`SourceIdx[1]` / `[2]` 是 `gamma_iso` 的下标；
- 不要求 `ie` 与 `j1` / `j2` 数值不同（它们来自不同文件，语义不同）；
- `species = [0,1,1]`：0 = 正电子，1 = gamma；
- 不做任何总能量窗口（上一版 100–125 MeV 窗口已废弃）。

### 5.6 为什么不用单粒子的 `_bkg` 文件？

`data/processed/*_bkg.root` 是“每个单粒子事例注入一次本底”。如果直接把 N 个 `_bkg` 事例相加，一个拼接事例会得到 **N 次独立本底抽样**，相当于把本底率放大 N 倍；而物理上多个粒子来自**同一个事件窗口**，本底只应抽一次。

因此 Stage 2 / Stage 3 的做法是：

1. 只读**无本底**的单粒子样本；
2. 在拼接完成后，对每个拼接事例调用一次 `draw_injections()`；
3. 把本底加到拼接后的 `Edeps` 上。

### 5.7 本底抽样

`draw_injections()` 来自 [bkg_est_add/add_background.py](../bkg_est_add/add_background.py)：

```python
def draw_injections(n_events, n_crystals, n_bkg, spectrum, seed, stream):
    rng = np.random.default_rng([int(seed), int(stream)])
    counts = rng.poisson(n_bkg, size=int(n_events)).astype(np.int64)
    n_hits = int(counts.sum())
    mods = rng.integers(0, int(n_crystals), size=n_hits, dtype=np.int64).astype(np.int32)
    energies = spectrum.rvs(n_hits, rng).astype(np.float32)
    return counts, mods, energies
```

返回：

```text
counts: int64[n_events]       # 每个事例的本底 hit 数 k ~ Poisson(n_bkg)
mods:   int32[Σcounts]        # 按事件顺序展开的晶体编号；有放回，允许重复
energies: float32[Σcounts]    # 截断正态 N(5,2²) 的能量
```

当前参数：

| 量 | 值 |
|:--|--:|
| `n_bkg` | 2.424660396855316 hit/event |
| 单 hit 能量均值 | 5.001089 MeV |
| 单 hit 能量 sigma | 1.414214 MeV |
| 有本底事例比例 | ≈ 0.91 |
| Stage 2 `2gamma` 的平均本底能量 | 12.1074 MeV |
| Stage 2 `ep2gamma` 的平均本底能量 | 12.1404 MeV |

`stitching.background_vectors()` 把扁平数组拆成：

```python
vectors: float32[n_events, 568]   # 稠密的本底能量图
mods_list: list[np.ndarray]       # 每个事例的晶体编号
energies_list: list[np.ndarray]   # 每个事例的 hit 能量
```

同一事例内重复的 `ModID` 用 `np.add.at(vectors[i], m, e)` 累加。

### 5.8 写出 `_sig` / `_bkg`

`generate.py` 的 `write_pair()` 对每个场景调用两次 `write_multiparticle_sample()`：

```python
# signal-only
stitching.write_multiparticle_sample(
    sig_path, edeps_sig, ntrack, species, source_idx, primary_energy,
    source_energy, signal_energy, zero_bkg, zero_nbkg,
    tc, tp, te, empty, empty)

# signal + background
stitching.write_multiparticle_sample(
    bkg_path, edeps_bkg, ntrack, species, source_idx, primary_energy,
    source_energy, signal_energy, bkg_total, bkg["counts"],
    tc, tp, te, bkg_mods, bkg_energies)
```

两次调用的区别只有：

| 参数 | `_sig` | `_bkg` |
|:--|:--|:--|
| `edeps` | `sig_2g` / `sig_ep2g` | `sig + bkg_vecs` |
| `bkg_energy_total` | 0 | `bkg_vecs.sum(axis=1)` |
| `nbkg` | 0 | `counts` |
| `bkg_crystal` / `bkg_energy` | 空列表 | `mods_list` / `energies_list` |

`Truth*`、`SourceIdx`、`Species`、`PrimaryEnergy`、`SourceEnergy`、`SignalEnergy` 在两个文件中**完全相同**，因此可以做配对比较（同一批信号，只差本底）。

### 5.9 写出函数 `write_multiparticle_sample()`

签名：

```python
def write_multiparticle_sample(
    path,
    edeps, ntrack, species, source_idx, primary_energy,
    source_energy, signal_energy, bkg_energy_total, nbkg,
    truth_crystal, truth_primary, truth_energy,
    bkg_crystal, bkg_energy,
    separation=None, sep_bin=None,
    tree_name="EdepOfEachEvt",
)
```

内部行为：

1. 把 numpy 数组转成 float32/int32（如果还不是）；
2. `separation` 和 `sep_bin` 必须同时给或同时不给，否则抛 `ValueError`；
3. 创建 `TTree("EdepOfEachEvt")`；
4. 按固定顺序 `Branch()` 所有分支；
5. 逐事件 `clear()` 每个 `std::vector`，`push_back` 当前事件的值，`Fill()`；
6. `Write()` + `Close()`，返回 `Path`。

注意：这个函数**不做长度/自洽校验**；长度和能量守恒由调用方（`generate.py`）保证，并在 `meta.json` 中记录检查结果。

### 5.10 Stage 2 生成结果

| 文件 | 事例数 | 内容 |
|:--|--:|:--|
| `data/stitched_2gamma_sig.root` | 20,000 | 2γ 信号，无本底 |
| `data/stitched_2gamma_bkg.root` | 20,000 | 同信号 + 一次 γ 本底 |
| `data/stitched_ep2gamma_sig.root` | 20,000 | e+2γ 信号，无本底 |
| `data/stitched_ep2gamma_bkg.root` | 20,000 | 同信号 + 一次 γ 本底 |

`meta.json` 记录：

```json
{
  "stage": "2",
  "seed": 20261106,
  "pool_events": 100000,
  "n_events": 20000,
  "primary_energy_mev": 53.085,
  "background": {
    "n_bkg_per_event": 2.424660396855316,
    "spectrum_mean_mev": 5.001089363774569,
    "spectrum_sigma_mev": 1.4142135623730951,
    "draw": "one Poisson draw per stitched event"
  },
  "scenarios": {
    "2gamma": {
      "signal_energy_mean_mev": 90.8483,
      "signal_energy_min_mev": 0.0,
      "signal_energy_max_mev": 106.17,
      "bkg_energy_mean_mev": 12.1074,
      "nbkg_mean": 2.4223,
      "frac_events_with_bkg": 0.91265
    },
    "ep2gamma": {
      "signal_energy_mean_mev": 137.2886,
      "signal_energy_min_mev": 0.0,
      "signal_energy_max_mev": 160.18,
      "bkg_energy_mean_mev": 12.1404,
      "nbkg_mean": 2.4235,
      "frac_events_with_bkg": 0.9075
    }
  }
}
```

---

## 6. Stage 2：评估流程

### 6.1 入口与常量

`evaluate.py`：

```python
SCENARIOS = {
    "2gamma":   {"ntrack": 2, "bkg": DATA / "stitched_2gamma_bkg.root"},
    "ep2gamma": {"ntrack": 3, "bkg": DATA / "stitched_ep2gamma_bkg.root"},
}
VERSIONS = ("sig", "bkg")
SEEDS = (5, 6, 7, 8, 10, 13, 14)
WAVES = (0.02, 0.05)
NTHS = (0, 1)
N_EVENTS = 20_000
FAKE_THRESHOLD = 1.0
MIN_PRIMARY_ENERGY = 20.0
```

网格大小：

```text
2 scenarios × 2 versions × 7 seeds × 2 waves × 2 Nth = 112 次重建
```

> `evaluate.py` 顶部 docstring 写的是 “seeds {7,8,10,13,14} ... 80 points”，这是旧版注释；以代码常量 `SEEDS` 和 `total_points` 为准（112 点）。

### 6.2 `prepare_workdir()`：符号链接桥

`reconstruction` 只认 `<data_dir>/processed/...` 和 `<data_dir>/utilities/...`，而 Stage 2 的样本在 `stage2/data/`。`evaluate.py` 不复制大文件，而是建符号链接：

```python
def prepare_workdir():
    proc = WORK / "processed"
    util = WORK / "utilities"
    proc.mkdir(parents=True, exist_ok=True)
    util.mkdir(parents=True, exist_ok=True)
    for scenario in SCENARIOS:
        for version in VERSIONS:
            src = DATA / f"stitched_{scenario}_{version}.root"
            dst = proc / src.name
            if dst.is_symlink() or dst.exists():
                dst.unlink()
            dst.symlink_to(src.resolve())
    topo = REPO / "data" / "utilities" / "ecal_neighbor_info.root"
    dst = util / topo.name
    if dst.is_symlink() or dst.exists():
        dst.unlink()
    dst.symlink_to(topo.resolve())
```

建好后的目录：

```text
stage2/work/
├── processed/
│   ├── stitched_2gamma_sig.root   -> <stage2>/data/stitched_2gamma_sig.root
│   ├── stitched_2gamma_bkg.root   -> <stage2>/data/stitched_2gamma_bkg.root
│   ├── stitched_ep2gamma_sig.root -> <stage2>/data/stitched_ep2gamma_sig.root
│   └── stitched_ep2gamma_bkg.root -> <stage2>/data/stitched_ep2gamma_bkg.root
└── utilities/
    └── ecal_neighbor_info.root    -> <repo>/data/utilities/ecal_neighbor_info.root
```

这样调用重建时：

```python
reconstruction.reconstruction(
    filename_edep_of_each_evt=sample.name,   # 如 "stitched_2gamma_bkg.root"
    energy_threshold_seed=float(seed),
    energy_threshold_wave=float(wave),
    n_neighbor_threshold=int(nth),
    data_dir=str(WORK),
    output_path=str(cluster_path),
)
```

### 6.3 读真值：`io.load_multiparticle_sample()`

```python
truth = {
    sc: io.load_multiparticle_sample(cfg["bkg"])
    for sc, cfg in SCENARIOS.items()
}
```

注意：**`sig` 和 `bkg` 共用同一份 truth**，`truth` 从 `_bkg` 文件读；`_sig` 文件的 `Truth*` 与它完全相同。

返回的字典：

```python
{
  "path": str,
  "n_events": int,                              # 20000
  "ntrack": int32[n_events],                    # 每事例的 N
  "species": list[int32[N]],                    # 每事例一个数组
  "source_idx": list[int32[N]],
  "primary_energy": list[float32[N]],
  "source_energy": list[float32[N]],
  "signal_energy": float32[n_events],
  "bkg_energy_total": float32[n_events],
  "nbkg": int32[n_events],
  "signal_arrays": list[float32[N, 568]],       # 逐初级信号能量矩阵
  "bkg_arrays": list[float32[568]],             # 稠密本底能量图
  "separation": float32[n_events],              # 无该分支时全 NaN
  "sep_bin": int32[n_events],                   # 无该分支时全 -1
}
```

`signal_arrays[i]` 的构造：

```python
sig = np.zeros((nt, n_mod), dtype=np.float32)
for c, p, e in zip(tc, tp, te):
    sig[int(p), int(c)] += float(e)
```

`bkg_arrays[i]` 的构造：

```python
bkg = np.zeros(n_mod, dtype=np.float32)
for c, e in zip(bc, be):
    bkg[int(c)] += float(e)
```

因此评估层拿到的是**稠密的逐初级矩阵**：`truth["signal_arrays"][i][j, c]` = 初级 `j` 在晶体 `c` 的信号能量。共享晶体在这里自然表现为同一列在两个初级行上都有值。

> `load_multiparticle_sample()` **不读 `Edeps`**。它只读 `NTrack`、`Species`、`SourceIdx`、`PrimaryEnergy`、`SourceEnergy`、`SignalEnergy`、`BkgEnergyTotal`、`NBkg`、`Truth*`、`Bkg*`（以及可选的 `Separation` / `SepBin`）。`Edeps` 只给重建用。

### 6.4 逐点循环

```python
for scenario, cfg in SCENARIOS.items():
    n_truth = cfg["ntrack"]
    for version in VERSIONS:
        sample = DATA / f"stitched_{scenario}_{version}.root"
        for seed in SEEDS:
            for wave in WAVES:
                for nth in NTHS:
                    tag = f"{scenario}_{version}_seed{seed:g}_wave{wave:g}_nth{nth}"
                    cluster_path = RESULTS / f"clusters_{tag}.root"
                    reconstruction.reconstruction(
                        filename_edep_of_each_evt=sample.name,
                        energy_threshold_seed=float(seed),
                        energy_threshold_wave=float(wave),
                        n_neighbor_threshold=int(nth),
                        data_dir=str(WORK),
                        output_path=str(cluster_path),
                    )
                    clusters = io.load_clusters(cluster_path, N_EVENTS)
                    row, detail = evaluate_point(truth[scenario], clusters, version, n_truth)
                    ...
```

`io.load_clusters()`：

```python
clusters = [[] for _ in range(n_events)]
with ROOT.TFile.Open(str(cluster_path), "READ") as f:
    tree = f.Get("ECALClusters")
    for entry in tree:
        evt = int(entry.EvtID)
        if 0 <= evt < len(clusters):
            clusters[evt].append([int(m) for m in entry.ModIDList])
return clusters
```

返回：

```python
clusters[event] = [[mod_id, ...], [mod_id, ...], ...]
```

因为 `clusters` 预分配了 `n_events`，零团簇事件自然保留为空列表；这与 `ECALClusters` 中“零团簇事件没有 entry”互补。

### 6.5 单事例评估 `evaluate_point()`

对每个事件 `i`：

```python
sig = truth["signal_arrays"][i]        # float32[N,568]
bkg = truth["bkg_arrays"][i] if version == "bkg" else zero_bkg
e_sig = sig.sum()
e_bkg = bkg.sum() if version == "bkg" else 0.0
cluster_mods = clusters[i]             # list[list[int]]

if cluster_mods:
    c_sig = np.array([sig[:, mods].sum() for mods in cluster_mods])
    c_bkg = np.array([bkg[mods].sum() for mods in cluster_mods])
else:
    c_sig = np.zeros(0)
    c_bkg = np.zeros(0)

n_total = len(cluster_mods)
n_signal = int(np.sum(c_sig >= FAKE_THRESHOLD))
```

然后计算：

```python
overlap = matching.build_overlap_matrix(cluster_mods, sig)
assignment = matching.assign(overlap)
eff, match_purity = matching.efficiency_purity(
    overlap, assignment, truth["source_energy"][i])
ms = matching.merge_split(overlap, truth["source_energy"][i])
bg = background.event_background(
    c_sig, c_bkg, e_sig, e_bkg, FAKE_THRESHOLD)
en = energy.event_energy(
    c_sig, c_bkg, e_sig, overlap, assignment,
    truth["source_energy"][i], MIN_PRIMARY_ENERGY)
```

### 6.6 重叠矩阵与匹配

对每个事例，`overlap` 是 `n_clusters × n_primary` 矩阵：

```text
M[k, j] = Σ_{c ∈ cluster_k} T[j, c]
```

其中 `T[j, c]` 是初级 `j` 在晶体 `c` 的信号能量（`truth["signal_arrays"][i][j, c]`）。

`matching.assign(M)` 返回 `assignment[k] = j` 或 `-1`：

- `n_clusters * n_primary <= 20` 时，穷举“簇子集 × 初级排列”，取总重叠最大；等价于匈牙利；
- 更大的矩阵回退到 `scipy.optimize.linear_sum_assignment`。

`efficiency_purity()`：

```text
ε_j = M[k_j, j] / E_j          # k_j 是指派给初级 j 的团簇；E_j = SourceEnergy[j]
p_k = M[k, a(k)] / Σ_j M[k, j]
```

`merge_split()` 不用匈牙利，而是看**每个初级的最佳团簇**：

```text
best_k(j) = argmax_k M[k, j]
best_eff(j) = M[best_k(j), j] / E_j
merge: 两个或更多初级的 best_k 相同
split: best_eff(j) < 0.9
```

这也是 Stage 3 合并率/分裂率的定义基础。

### 6.7 本底与能量指标

`background.event_background()` 对单个事例：

```text
n_fake        = #{k : c_sig[k] < 1 MeV}
signal_purity = Σ c_sig / (Σ c_sig + Σ c_bkg)
bkg_capture   = Σ c_bkg / e_bkg                   # e_bkg=0 时为 NaN
lead          = argmax_k (c_sig[k] + c_bkg[k])
lead_bkg_frac = c_bkg[lead] / (c_sig[lead] + c_bkg[lead])
```

`energy.event_energy()`：

```text
signal_capture      = Σ c_sig / e_sig
lead_containment    = c_sig[lead] / e_sig
energy_bias(j)      = (M[a(j), j] - E_j) / E_j     # 仅 E_j >= 20 MeV 的已匹配初级
```

### 6.8 汇总指标定义

`evaluate_point()` 对 20,000 个事例累加后写出以下字段（`stage2_metrics.csv` 的列）：

| 字段 | 定义 / 分母 |
|:--|:--|
| `n_events` | 20,000 |
| `ntrack` | 2 或 3 |
| `exact_total` / `under_total` / `over_total` | 总团簇数 `n_total` 与 N 比较的比例 |
| `mae_total` | `mean(abs(n_total − N))` |
| `mean_clusters_total` | `mean(n_total)` |
| `exact_signal` / `under_signal` / `over_signal` | 信号簇数 `n_signal`（`c_sig >= 1 MeV`）与 N 比较 |
| `mae_signal` | `mean(abs(n_signal − N))` |
| `mean_clusters_signal` | `mean(n_signal)` |
| `fake_clusters_per_event` | `Σ n_fake / n_events` |
| `frac_events_with_fake` | 至少有一个假团簇的事例比例 |
| `signal_purity_energy` | `Σ c_sig / (Σ c_sig + Σ c_bkg)`（所有事例所有团簇的能量总和） |
| `signal_purity_mean` | 只在有团簇的事例上平均 `bg["signal_purity"]` |
| `bkg_capture` | `Σ c_bkg / Σ e_bkg`；`_sig` 版本为 NaN |
| `signal_capture` | `Σ c_sig / Σ e_sig` |
| `lead_signal_containment` | 只在有团簇且 `e_sig > 0` 的事例上平均 `c_sig[lead]/e_sig` |
| `lead_bkg_fraction` | 只在有团簇的事例上平均 `c_bkg[lead]/c_tot[lead]` |
| `efficiency_mean` | 对**所有初级**平均 `ε_j`；未匹配/零沉积初级贡献 0 |
| `match_purity_mean` | 对**信号团簇**平均 `p_k`；假团簇不计入 |
| `matched_primary_fraction` | 满足 `ε > 0.5` 且 `p > 0.5` 的簇-初级指派数 / 总初级数 |
| `merge_rate_primaries` | `Σ merged_primaries / 总初级数` |
| `merge_rate_events` | 至少有一个合并的事例 / `n_events` |
| `split_rate_primaries` | `Σ split_primaries / 总初级数` |
| `energy_bias_mean` | 已匹配且 `E_j >= 20 MeV` 的初级的 `(M[k,j]−E_j)/E_j` 均值 |
| `energy_rms` | 同一组 bias 的 `sqrt(mean(bias²))` |
| `n_energy` | 进入能量统计的初级数 |

> **分母陷阱**：`efficiency_mean`、`merge_rate_primaries`、`split_rate_primaries` 的分母都是“所有初级”，包括 `E_sig=0` 的零沉积初级。零沉积初级既不能被匹配，也不会被判为 split/merge（`best_energy=0` 时 `valid=False`），所以这些指标会被零沉积事例稀释。报告中的“条件指标”（所有源 `E_sig>0`）就是为了排除这种稀释。

### 6.9 输出文件

```text
stage2/results/
├── stage2_metrics.csv          # 112 行 × 35 列，每个参数点一行
├── points/<tag>.json           # 每个点一个：{"row": {...}, "hist_total": {...}, "hist_signal": {...}}
├── best_points.json            # 4 个最佳点（2 场景 × sig/bkg）
├── background_impact.json      # bkg 最佳点与同参数 sig 点的配对差
├── exact_signal_vs_seed_<scenario>.png
├── fake_purity_vs_seed_<scenario>.png
├── merge_split_best.png
└── count_hist_best.png
```

`best_points.json` 的选择规则：

```python
def pick_best(rows):
    return max(rows, key=lambda r: (
        r["exact_signal"],
        r["signal_purity_energy"],
        -r["fake_clusters_per_event"],
    ))
```

`background_impact.json` 对每个场景取 `version="bkg"` 的最佳点，再找同 `seed/wave/nth` 的 `version="sig"` 行，计算 `exact`、`merge`、`split`、`efficiency` 的配对差。

当前最佳点（报告口径）：

| 场景 / 版本 | seed | wave | Nth | exact_signal | fake/event | merge(events) | split(primaries) |
|:--|--:|--:|--:|--:|--:|--:|--:|
| `2gamma` bkg | 13 | 0.05 | 0 | 0.7098 | 0.00015 | 0.0544 | 0.0326 |
| `2gamma` sig | 13 | 0.05 | 0 | 0.7100 | 0.0000 | 0.0542 | 0.0329 |
| `ep2gamma` bkg | 10 | 0.05 | 0 | 0.53915 | 0.0032 | 0.1572 | 0.0354 |
| `ep2gamma` sig | 10 | 0.05 | 0 | 0.53965 | 0.0000 | 0.1565 | 0.0356 |

### 6.10 运行时间

报告给出的单进程运行时间约 **30 分钟**（112 点 × 20,000 事例）。主要成本是 `reconstruction` 的 Python 逐事件循环，以及每个参数点 `evaluate_point()` 的逐事例循环；`io.load_multiparticle_sample()` 只在评估开始时展开一次真值，但它会在内存中保留 20,000 个 `(N,568)` 矩阵和 20,000 个 `(568,)` 本底数组。

---

## 7. Stage 3：样本生成流程

### 7.1 目标

Stage 3 把“多粒子分辨”变成一条曲线：

```text
合并率 vs 两个簇射的能量加权质心距离（单位：crystal pitch）
```

只做 `2gamma`，N=2。每个分离度箱 2,500 事例，共 6 箱，15,000 事例。

### 7.2 晶体位置与 pitch

首先生成位置表（只需一次）：

```bash
python3 pytrial/geometry.py \
  --raw data/raw/ep_iso.root \
  --out data/utilities/ecal_module_positions.npz \
  --n-hits 500000
```

`geometry.py` 的关键函数：

| 函数 | 作用 | 返回 |
|:--|:--|:--|
| `read_module_positions(raw_path, n_hits, n_mod)` | 从 `G4Run0/ECALSimHit` 按 `ModID` 平均 `x[3]` | `centre[568,3]`, `counts[568]` |
| `save_positions(path, centre, counts)` | 写压缩 `.npz` | `Path` |
| `load_positions(path)` | 读 `.npz` | `(centre, counts)` |
| `energy_weighted_centroid(edeps, centre)` | `Σ E_i x_i / Σ E_i` | `(3,)` 或 NaN |
| `pairwise_distances(centroids)` | 两两欧氏距离 | `(n,n)`，对角线 NaN |
| `min_pairwise_distance(centroids)` | 最小两两距离 | `float` / NaN |
| `estimate_pitch(centre, counts)` | 被击中模块的最近邻距离中位数 | `float`（mm） |

当前：

```text
568/568 模块有 hit
pitch = 46.24094384427734 mm
```

### 7.3 质心

对每个 `E_sig > 0` 的 `gamma_iso` 事例：

```python
centroids = (ga[pool] @ centre) / e_ga[pool, None]
```

- `ga[pool]`：`(m, 568)`；
- `centre`：`(568, 3)`；
- 结果 `(m, 3)`，单位 mm；
- `pool = np.flatnonzero(e_ga > 0)`，当前 `87,747 / 100,000`（12.3% 零沉积被排除，因为质心无定义）。

### 7.4 分箱与配对采样

箱边界（单位 pitch）：

```python
EDGES_PITCH = np.array([0.0, 0.5, 1.0, 2.0, 3.0, 5.0, np.inf])
N_BINS = 6
EDGES_MM = EDGES_PITCH * pitch
N_PER_BIN = 2500
```

| bin | 范围（pitch） | 范围（mm，当前 pitch=46.2409） | 采样方式 |
|--:|:--|:--|:--|
| 0 | 0–0.5 | 0–23.12 | `sample_knn_pairs`，k=200 |
| 1 | 0.5–1 | 23.12–46.24 | `sample_knn_pairs`，k=200 |
| 2 | 1–2 | 46.24–92.48 | `sample_knn_pairs`，k=500 |
| 3 | 2–3 | 92.48–138.72 | `sample_random_pairs` |
| 4 | 3–5 | 138.72–231.20 | `sample_random_pairs` |
| 5 | >5 | >231.20 | `sample_random_pairs` |

`sample_knn_pairs()`：

1. 随机选一个种子事例 `s`；
2. `cKDTree.query(centroids[s], k=k)`；
3. 取距离在 `[lo, hi)` 的候选；
4. 随机选一个候选 `p`，用 `(min(s,p), max(s,p))` 去重；
5. 直到该箱收集到 2,500 个唯一对。

`sample_random_pairs()`：

1. 随机抽 `i, j`（有放回）；
2. 计算距离，保留 `i != j` 且 `d ∈ [lo, hi)` 的对；
3. 用 `(min(i,j), max(i,j))` 去重；
4. 直到 2,500 个唯一对。

采样结果按 bin 顺序拼接：前 2,500 个是 bin 0，接着 2,500 个是 bin 1，依此类推。因此**文件中的事件顺序是“按箱块排列”**，不是全局随机打乱。

> `stage3/generate.py` 的 `EDGES_PITCH` 与 `evaluation/separation.py` 的 `DEFAULT_EDGES_PITCH` 必须保持一致。评估时 `truth["sep_bin"]` 直接来自文件；`DEFAULT_EDGES_PITCH` 主要用于输出行的 `sep_lo_pitch` / `sep_hi_pitch` 标签。如果修改分箱，必须重新生成 Stage 3 样本。

### 7.5 组装与写出

```python
src_a = np.array([pool[a] for a, _ in all_pairs], dtype=np.int32)
src_b = np.array([pool[c] for _, c in all_pairs], dtype=np.int32)
e_a = e_ga[src_a].astype(np.float32)
e_b = e_ga[src_b].astype(np.float32)
signal_energy = (e_a + e_b).astype(np.float32)
species = np.ones((total, 2), dtype=np.int32)       # 两个 gamma
src = np.stack([src_a, src_b], axis=1).astype(np.int32)
src_e = np.stack([e_a, e_b], axis=1).astype(np.float32)
primary_e = np.full((total, 2), PRIMARY_ENERGY_MEV, dtype=np.float32)

edeps_sig = (ga[src_a] + ga[src_b]).astype(np.float32)
truth = build_truth([ga[src_a], ga[src_b]])
bg_draw = draw_injections(total, N_MOD, n_bkg, spectrum, SEED, 3)
bg_vecs, bg_mods, bg_energies = stitching.background_vectors(*bg_draw)
edeps_bkg = (edeps_sig + bg_vecs).astype(np.float32)
```

写出：

```python
stitching.write_multiparticle_sample(
    sig_path, edeps_sig, ntrack, species, src, primary_e, src_e,
    signal_energy, zero_bkg, zero_nbkg, *truth, empty, empty,
    separation=sep_mm, sep_bin=sep_bin)

stitching.write_multiparticle_sample(
    bkg_path, edeps_bkg, ntrack, species, src, primary_e, src_e,
    signal_energy, bkg_total, bg_draw[0], *truth, bg_mods, bg_energies,
    separation=sep_mm, sep_bin=sep_bin)
```

与 Stage 2 的差别：

- 只做 `2gamma`；
- 源 pool 只保留 `E_sig>0`（质心需要）；
- 不是完全随机抽样，而是按分离度箱定向抽样；
- 额外写 `Separation`（mm）和 `SepBin`（0..5）。

### 7.6 Stage 3 生成结果

| 文件 | 事例数 | 内容 |
|:--|--:|:--|
| `data/stage3_2gamma_sig.root` | 15,000 | 6 箱 × 2,500，signal-only |
| `data/stage3_2gamma_bkg.root` | 15,000 | 同信号 + 一次 γ 本底 |

`meta.json` 记录：

```json
{
  "stage": "3",
  "pitch_mm": 46.24094384427734,
  "edges_pitch": [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, Infinity],
  "edges_mm": [0.0, 23.120..., 46.240..., 92.481..., 138.722..., 231.204..., Infinity],
  "bin_counts": [2500, 2500, 2500, 2500, 2500, 2500],
  "bin_median_mm": [14.516, 26.772, 46.746, 119.376, 189.836, 467.257],
  "signal_energy_mean_mev": 103.1257,
  "signal_energy_min_mev": 2.044,
  "signal_energy_max_mev": 106.170,
  "source_energy_mean_by_primary_mev": [51.5777, 51.5481],
  "background": {
    "n_bkg_per_event": 2.424660396855316,
    "spectrum_mean_mev": 5.001089363774569,
    "frac_events_with_bkg": 0.9104
  },
  "max_abs_truth_minus_signal_mev": 7.57e-06,
  "max_abs_edeps_bkg_minus_signal_plus_bkg_mev": 3.05e-05
}
```

> `meta.json` 由 Python `json.dumps` 写出，`np.inf` 会写成非标准的 `Infinity`；Python 的 `json.loads` 可以正常读回。

> **条件样本**：每个 bin 等统计量，不是物理分离度分布。因此不能把 6 个 bin 的平均当作物理平均；物理平均要看 Stage 2 的随机 2γ 样本。

---

## 8. Stage 3：评估流程

### 8.1 入口与常量

```python
VERSIONS = ("sig", "bkg")
SEEDS = (5, 6, 7, 8, 10, 13, 14)
NTHS = (0, 1)
WAVE = 0.05
N_EVENTS = 15_000
FAKE_THRESHOLD = 1.0
BINS = separation.DEFAULT_EDGES_PITCH
N_BINS = 6
```

网格：

```text
2 versions × 7 seeds × 2 Nth = 28 次重建
```

每个点输出 6 行（每个分离度箱一行），所以：

```text
stage3_metrics.csv = 28 × 6 = 168 行
fits.json          = 28 条拟合记录
```

### 8.2 工作目录

与 Stage 2 相同，`prepare_workdir()` 把：

```text
stage3/work/processed/stage3_2gamma_{sig,bkg}.root  -> <stage3>/data/...
stage3/work/utilities/ecal_neighbor_info.root       -> <repo>/data/utilities/...
```

然后调用 `reconstruction`。

### 8.3 读真值

```python
truth = io.load_multiparticle_sample(DATA / "stage3_2gamma_bkg.root")
pitch = float(json.loads((HERE / "meta.json").read_text())["pitch_mm"])
```

因为文件有 `Separation` / `SepBin`，返回的 `truth` 中：

```python
truth["separation"]  # float32[15000]，单位 mm
truth["sep_bin"]     # int32[15000]，0..5
```

### 8.4 逐箱评估 `evaluate_point()`

```python
stats = {
    b: {"n": 0, "merge": 0, "split": 0, "exact2": 0,
        "eff_sum": 0.0, "eff_n": 0, "mp_sum": 0.0, "mp_n": 0,
        "fake": 0, "sig_in": 0.0, "bkg_in": 0.0, "bkg_total": 0.0}
    for b in range(N_BINS)
}

for i in range(n_events):
    b = int(truth["sep_bin"][i])
    if b < 0 or b >= N_BINS:
        continue
    st = stats[b]
    ...
```

单事例的计算与 Stage 2 相同：

```python
c_sig = [sig[:, mods].sum() for mods in cluster_mods]
c_bkg = [bkg[mods].sum() for mods in cluster_mods]
n_signal = int(np.sum(c_sig >= FAKE_THRESHOLD))

overlap = matching.build_overlap_matrix(cluster_mods, sig)
assignment = matching.assign(overlap)
eff, match_purity = matching.efficiency_purity(overlap, assignment, source_energy)
ms = matching.merge_split(overlap, source_energy)
bg = background.event_background(c_sig, c_bkg, e_sig, e_bkg, FAKE_THRESHOLD)
```

累加：

```python
st["merge"] += int(ms["merge_event"])
st["split"] += int(ms["split_primaries"])
st["exact2"] += int(n_signal == 2)
st["eff_sum"] += float(eff.sum())
st["eff_n"] += int(eff.size)          # 每个事例 2 个初级
st["fake"] += int(bg["n_fake"])
st["sig_in"] += float(c_sig.sum())
st["bkg_in"] += float(c_bkg.sum())
st["bkg_total"] += e_bkg
```

### 8.5 每箱指标定义

| 字段 | 定义 / 分母 |
|:--|:--|
| `bin` | 0..5 |
| `sep_lo_pitch` / `sep_hi_pitch` | 箱边界（pitch） |
| `median_sep_mm` / `median_sep_pitch` | 该箱内 `Separation` 的中位数 |
| `n_events` | 该箱事例数（2500） |
| `merge_rate` | 至少一个合并的事例 / `n_events` |
| `split_rate` | `Σ split_primaries / st["eff_n"]`；分母是初级数（2 × n_events） |
| `exact2` | `n_signal == 2` 的事例 / `n_events` |
| `efficiency_mean` | `Σ eff / st["eff_n"]`；分母是初级数 |
| `match_purity_mean` | 只在信号团簇上平均 `p_k` |
| `fake_clusters_per_event` | `Σ n_fake / n_events` |
| `signal_purity` | `Σ c_sig / (Σ c_sig + Σ c_bkg)` |
| `bkg_capture` | `Σ c_bkg / Σ e_bkg` |

> Stage 3 的源事例都满足 `E_sig>0`，所以这里的 `efficiency_mean` / `merge_rate` / `split_rate` **没有** Stage 2 的零沉积稀释；这也是 Stage 3 条件曲线可以和 Stage 2 “条件指标”对比的原因。

### 8.6 合并半径 `fits.json`

对每个 `version × seed × Nth`：

```python
x = [r["median_sep_pitch"] for r in sel]   # 6 个箱的中位分离度
y = [r["merge_rate"] for r in sel]        # 6 个箱的合并率

fits.append({
    "version": version,
    "seed": seed,
    "nth": nth,
    "r50_pitch": separation.crossing_radius(x, y, 0.5),
    "r10_pitch": separation.crossing_radius(x, y, 0.1),
})
```

`crossing_radius(x, y, level)` 的逻辑：

1. 从 `i=1` 开始找第一个满足 `y[i-1] >= level > y[i]` 的区间；
2. 在该区间内线性插值：
   ```text
   t = (level - y[i-1]) / (y[i] - y[i-1])
   x_cross = x[i-1] + t * (x[i] - x[i-1])
   ```
3. 如果曲线一开始就低于 `level`，返回 `x[0]`；
4. 如果从不向下穿过 `level`，返回 NaN。

因为只有 6 个箱，插值是粗粒度的；半径的不确定度主要来自箱边界，而不是每箱 2,500 的统计量。

当前结果（seed=13, Nth=0）：

| 版本 | 50% 半径（pitch） | 50% 半径（mm） | 10% 半径（pitch） | 10% 半径（mm） |
|:--|--:|--:|--:|--:|
| signal-only | 3.410 | 157.7 | 7.430 | 343.6 |
| +background | 3.418 | 158.1 | 7.477 | 345.8 |

### 8.7 输出文件

```text
stage3/results/
├── stage3_metrics.csv          # 168 行 × 18 列
├── fits.json                   # 28 条：r50/r10 pitch + mm
├── merge_rate_vs_separation.png
└── merge_rate_best.png
```

`merge_rate_vs_separation.png` 画 `merge_rate` vs `median_sep_pitch`，覆盖 seed 5/7/10/14 与 Nth 0/1；`merge_rate_best.png` 使用 `pick_stage2_best()` 从 Stage 2 CSV 中选出的最佳 2γ bkg 点（当前 seed=13, Nth=0）。

---

## 9. 重建算法内部流程

这一节把 `reconstruction.reconstruction()` 的内部数据流展开，便于把“输入格式”和“输出格式”对上。

### 9.1 函数签名

```python
def reconstruction(
    filename_edep_of_each_evt,
    energy_threshold_seed,
    energy_threshold_wave,
    n_neighbor_threshold,
    n_fallback_chance=1,
    neighbor_file_name="ecal_neighbor_info.root",
    data_dir=None,
    output_path=None,
):
```

| 参数 | 类型 | 含义 |
|:--|:--|:--|
| `filename_edep_of_each_evt` | `str` | `data_dir/processed/` 下的文件名 |
| `energy_threshold_seed` | `float` | 种子阈值（MeV）；低于此值的晶体不能当种子 |
| `energy_threshold_wave` | `float` | 波前阈值（MeV）；低于此值的晶体不能加入团簇 |
| `n_neighbor_threshold` | `int` | 种子所需的“点亮邻居”数 |
| `n_fallback_chance` | `int` | 旧版 gap-jumping 参数；当前未使用 |
| `neighbor_file_name` | `str` | `data_dir/utilities/` 下的邻接表文件名 |
| `data_dir` | `str 或 None` | 数据根目录；`None` 时读环境变量 |
| `output_path` | `str/Path 或 None` | 输出 ROOT；`None` 时自动命名到当前工作目录 |

返回：`Path(output_path)`。

### 9.2 算法步骤

```text
Step 1  读邻接表：neighbors[i] = 晶体 i 的邻居列表；nMod = 条目数
Step 2  准备输出树：ECALClusters / EvtID / ModIDList
Step 3  逐事件 cluster_event(edeps, evt_id)
Step 4  遍历输入树的所有 entry
Step 5  写输出树并关闭
```

`cluster_event()` 的内部：

```text
edeps_full = 长度 nMod 的 float 列表（不足补 0，超出截断）
states     = 长度 nMod 的 int 列表，0 = 未聚类

sorted_idx = 按 edeps_full 从大到小排序的晶体编号

is_lit_on(c)       = edeps_full[c] >= energy_threshold_wave
is_seedable(c)     = states[c] == 0
                     and edeps_full[c] >= energy_threshold_seed
                     and #{nb in neighbors[c] : is_lit_on(nb)} >= n_neighbor_threshold
is_propagatable(c) = states[c] == 0 and is_lit_on(c)

seed_ptr = 0
while True:
    从 seed_ptr 向后找第一个满足 is_seedable(sorted_idx[seed_ptr]) 的晶体
    如果找不到：break
    seed_id = sorted_idx[seed_ptr]
    states[seed_id] = 1
    cluster_list = [seed_id]
    wave = neighbors[seed_id]
    next_wave = []
    crystal_state = 1

    while wave:
        crystal_state += 1
        for cur_id in wave:
            if is_propagatable(cur_id):
                states[cur_id] = crystal_state
                cluster_list.append(cur_id)
                for nb in neighbors[cur_id]:
                    if nb not in next_wave:
                        next_wave.append(nb)
        wave, next_wave = next_wave, []

    把 cluster_list 写入 ModIDList，EvtID=event_id，Fill()
```

### 9.3 关键语义

1. **`Nth` 只约束种子，不约束波前扩展。**
   - 一个晶体一旦被种子选中，后续加入团簇只要求 `states == 0` 且 `E >= Eth_wave`；
   - 它不需要自己有 `Nth` 个点亮邻居。
   - 这与旧版 `pipeline.md` 的描述不同（旧文档说“crystal joins if ... has at least n_neighbor_threshold lit neighbors”），以当前代码为准。

2. **种子按能量降序扫描。**
   - `sorted_idx` 按 `edeps_full` 降序；
   - `seed_ptr` 单调前进；已经被聚类的晶体 `states != 0`，会被跳过；
   - 因此高能晶体优先成为种子。

3. **波前是 BFS。**
   - 从种子的邻居开始，逐层向外；
   - `next_wave` 用 `not in` 去重；
   - `states` 用 `crystal_state`（2,3,...）标记层级，但后续只检查 `== 0`。

4. **每个团簇一条输出 entry。**
   - `EvtID` 相同、`ModIDList` 不同；
   - 一个事件产生 K 个团簇就有 K 条 entry；
   - 零团簇事件没有 entry。

5. **`second_neighbors` / gap-jumping 未启用。**
   - `sec_neighbors` 会被读入，但对应回退代码被注释；
   - `n_fallback_chance` 参数保留但无效果。

6. **阈值关系。**
   - 代码只在 `energy_threshold_seed <= energy_threshold_wave` 时打印警告，不阻止运行；
   - 物理上应该 `seed >= wave`；Stage 2/3 的网格都满足。

### 9.4 CLI

`reconstruction.py` 的 `__main__` 入口：

```bash
python3 reconstruction.py <filename> \
    [Eth_seed=10] [Eth_wave=5] [Nth=3] [output_name] [expected_lambda] \
    [--event-by-event] [--data-path path] [--data-dir path]
```

它会先调用 `reconstruction()`，再调用 `scoring.run()`。注意：

- **Stage 2 / Stage 3 不用这个 CLI**；它们直接 import `reconstruction.reconstruction()`；
- `scoring.run()` 是单粒子/计数评估，不能区分“合并”和“漏掉”，也不是多粒子评估；
- 旧 `pipeline.md` 里的 CLI 顺序（`[Eth] [Nth]`）和函数签名已经过时，以当前代码为准。

---

## 10. 评估接口速查

### 10.1 生产 / 拼接

| 函数 | 位置 | 输入 | 输出 |
|:--|:--|:--|:--|
| `stitching.load_edeps(path, n_events, tree_name, n_mod)` | [stitching.py](../stitching.py) | 单粒子 ROOT | `float32[n,568]` |
| `stitching.event_energies(edeps)` | 同上 | `[n,568]` | `float64[n]` |
| `stitching.positive_energy_mask(edeps, eps=0)` | 同上 | `[n,568]` | `bool[n]` |
| `stitching.stitch_count_only(edeps, n_primary, n_events, rng)` | 同上 | 单粒子矩阵 | `(stitched, ntrack, source_idx)`（Stage 1） |
| `stitching.write_stitched(path, edeps, ntrack, source_idx)` | 同上 | 计数拼接结果 | `Edeps/NTrack/SourceIdx` ROOT（Stage 1） |
| `stitching.background_vectors(counts, mods, energies)` | 同上 | 扁平本底抽样 | `(vectors, mods_list, energies_list)` |
| `stitching.write_multiparticle_sample(...)` | 同上 | 拼接数组 + 真值 | 多粒子 ROOT（Stage 2/3） |

### 10.2 拓扑 / 几何

| 函数 | 位置 | 输入 | 输出 |
|:--|:--|:--|:--|
| `geometry.read_module_positions(raw_path, n_hits, n_mod)` | [geometry.py](../geometry.py) | 原始 MC | `(centre[568,3], counts[568])` |
| `geometry.save_positions(path, centre, counts)` | 同上 | 位置表 | `.npz` |
| `geometry.load_positions(path)` | 同上 | `.npz` | `(centre, counts)` |
| `geometry.energy_weighted_centroid(edeps, centre)` | 同上 | 一个事例的 `Edeps` | `(3,)` / NaN |
| `geometry.pairwise_distances(centroids)` | 同上 | `(n,3)` | `(n,n)` |
| `geometry.min_pairwise_distance(centroids)` | 同上 | `(n,3)` | `float` / NaN |
| `geometry.estimate_pitch(centre, counts)` | 同上 | 位置表 | `float`（mm） |

### 10.3 评估

| 函数 | 位置 | 输入 | 输出 |
|:--|:--|:--|:--|
| `io.load_multiparticle_sample(path, n_mod=568)` | [evaluation/io.py](../evaluation/io.py) | 多粒子 ROOT | truth 字典（§6.3） |
| `io.load_clusters(cluster_path, n_events)` | 同上 | `ECALClusters` ROOT | `clusters[event] = [[mods], ...]` |
| `counts.read_cluster_counts(cluster_path, n_events)` | [evaluation/counts.py](../evaluation/counts.py) | 团簇 ROOT | `int64[n_events]` |
| `counts.signal_cluster_counts(cluster_sig_list, n_events, threshold=1.0)` | 同上 | 每事例团簇信号能量 | `int64[n_events]` |
| `counts.count_metrics(counts, n_truth)` | 同上 | 计数数组 + N | exact/under/over/MAE/RMSE/直方图 |
| `matching.build_overlap_matrix(cluster_mods, signal_arrays)` | [evaluation/matching.py](../evaluation/matching.py) | 团簇 + 逐初级矩阵 | `M[n_clusters, n_primary]` |
| `matching.assign(matrix)` | 同上 | `M` | `assignment[n_clusters]`（-1 = 未指派） |
| `matching.efficiency_purity(M, a, primary_totals)` | 同上 | `M`, 指派, `E_j` | `(eff[n_primary], purity[n_clusters])` |
| `matching.merge_split(M, primary_totals, split_threshold=0.9)` | 同上 | `M`, `E_j` | best_cluster / merge / split 字典 |
| `background.event_background(c_sig, c_bkg, e_sig, e_bkg, fake_threshold=1.0)` | [evaluation/background.py](../evaluation/background.py) | 每团簇信号/本底能量 | 假团簇/纯度/捕获字典 |
| `energy.event_energy(c_sig, c_bkg, e_sig, M, a, primary_totals, min_primary_energy=20.0)` | [evaluation/energy.py](../evaluation/energy.py) | 能量与匹配 | capture/containment/bias 字典 |
| `separation.assign_bins(sep_pitch, edges)` | [evaluation/separation.py](../evaluation/separation.py) | 分离度 | bin index |
| `separation.crossing_radius(x, y, level=0.5)` | 同上 | 中位分离度 + 合并率 | 半径（线性插值） |

### 10.4 本底模型

| 函数/类 | 位置 | 用途 |
|:--|:--|:--|
| `bkg_model.estimate()` | [bkg_est_add/bkg_model.py](../bkg_est_add/bkg_model.py) | 计算 `n_bkg`、几何、能谱 |
| `bkg_model.DepositionSpectrum(mean_mev=5.0, var_mev2=2.0)` | 同上 | 截断正态本底能谱 |
| `add_background.draw_injections(n_events, n_crystals, n_bkg, spectrum, seed, stream)` | [bkg_est_add/add_background.py](../bkg_est_add/add_background.py) | 抽 `counts/mods/energies` |

---

## 11. 复现与运行

### 11.1 环境

需要：

```text
Python 3
ROOT（PyROOT，当前 6.40.04）
numpy
scipy
matplotlib
```

数据文件走 Git LFS；如果只有指针文件：

```bash
git lfs fetch --include="data/processed/**" origin
git lfs checkout
```

### 11.2 Stage 2

```bash
cd pytrial/results/multiparticle/stage2
bash run_stage2.sh
```

等价于：

```bash
python3 generate.py
python3 evaluate.py
```

输出：

```text
data/stitched_2gamma_sig.root
data/stitched_2gamma_bkg.root
data/stitched_ep2gamma_sig.root
data/stitched_ep2gamma_bkg.root
meta.json
results/stage2_metrics.csv
results/points/*.json
results/best_points.json
results/background_impact.json
results/*.png
```

### 11.3 Stage 3

```bash
cd pytrial/results/multiparticle/stage3

# 1) 首次运行：生成晶体位置表（只需一次）
python3 ../../../geometry.py \
  --raw ../../../../data/raw/ep_iso.root \
  --out ../../../../data/utilities/ecal_module_positions.npz \
  --n-hits 500000

# 2) 生成样本 + 评估
bash run_stage3.sh
```

输出：

```text
data/stage3_2gamma_sig.root
data/stage3_2gamma_bkg.root
meta.json
results/stage3_metrics.csv
results/fits.json
results/merge_rate_vs_separation.png
results/merge_rate_best.png
```

### 11.4 单点手动运行

如果想只跑一个参数点：

```bash
# 1) 建 work 目录（或直接复用 evaluate.py 的 prepare_workdir）
mkdir -p /tmp/ecalrun/processed /tmp/ecalrun/utilities
ln -sf "$PWD/pytrial/results/multiparticle/stage2/data/stitched_2gamma_bkg.root" \
       /tmp/ecalrun/processed/
ln -sf "$PWD/data/utilities/ecal_neighbor_info.root" \
       /tmp/ecalrun/utilities/

# 2) 调用重建
cd pytrial
python3 - <<'PY'
import reconstruction
reconstruction.reconstruction(
    filename_edep_of_each_evt="stitched_2gamma_bkg.root",
    energy_threshold_seed=13.0,
    energy_threshold_wave=0.05,
    n_neighbor_threshold=0,
    data_dir="/tmp/ecalrun",
    output_path="/tmp/ecalrun/clusters.root",
)
PY
```

---

## 12. 不变量、验收检查与常见坑

### 12.1 生成阶段应检查

1. `len(Edeps) == 568`；
2. `len(Species) == len(SourceIdx) == len(PrimaryEnergy) == len(SourceEnergy) == NTrack`；
3. `Σ TruthEnergy == SignalEnergy`（float32 容差）；
4. `Σ Edeps_sig == SignalEnergy`；
5. `Σ Edeps_bkg == SignalEnergy + BkgEnergyTotal`；
6. `len(BkgCrystal) == len(BkgEnergy) == NBkg`；
7. `sig` 和 `bkg` 的 `SourceIdx` / `Truth*` 完全一致（配对比较的前提）。

Stage 2 的 `meta.json` 和 Stage 3 的 `meta.json` 已记录 3–5 的实测最大偏差。

### 12.2 评估阶段应检查

1. `clusters` 的长度等于 `N_EVENTS`（不是 `ECALClusters.GetEntries()`）；
2. `EvtID` 在 `[0, n_events)` 内；
3. `overlap` 的行和 = 该团簇的总信号能量（≤ 团簇总能量，后者还包含本底）；
4. `overlap` 的列和 ≤ `SourceEnergy[j]`（只有该初级的信号晶体全部被聚类时才取等号）；
5. 所有比例在 `[0,1]`，NaN 只出现在本底/零能量相关的预期位置；
6. `hist_total` / `hist_signal` 的计数之和等于 `n_events`。

### 12.3 常见坑

| 坑 | 说明 |
|:--|:--|
| 用 `ECALClusters.GetEntries()` 当事件数 | 零团簇事件没有 entry；必须用输入事件数 `N_EVENTS`。`io.load_clusters` / `counts.read_cluster_counts` 都要求显式 `n_events`。 |
| 以为 `ECALClusters` 存了能量 | 只存 `ModIDList`；团簇信号/本底能量必须用 `Truth*` 反算。 |
| 把 `data/processed/*_bkg.root` 当 Stage 2/3 输入 | Stage 2/3 只用无 `_bkg` 的单粒子文件，然后对拼接事例抽一次本底；直接用单粒子 `_bkg` 会把本底重复 N 次。 |
| 把 `SourceIdx` 当全局 ID | 它只是在对应源样本 pool 中的下标；`ep2gamma` 的第 0 个和第 1/2 个来自不同文件。 |
| 用 `pipeline.md` 的旧 API | 旧文档的 `reconstruction(filename, Eth, Nth)` 和 CLI 顺序已过时；当前需要 `Eth_seed`、`Eth_wave`、`Nth`。 |
| 以为 `Nth` 约束波前扩展 | 实际只约束种子；波前扩展只要求 `E >= Eth_wave` 且未聚类。 |
| 以为 `second_neighbors` 生效 | 回退代码被注释，`n_fallback_chance` 无效；二阶邻居当前不参与聚类。 |
| 对 `_sig` 版本算 `bkg_capture` | `e_bkg=0`，结果为 NaN；这是预期行为。 |
| 用 `exact_total` 代替 `exact_signal` | 有本底假团簇时两者不同；信号物理量通常看 `exact_signal`。 |
| 把 Stage 3 六箱平均当物理平均 | Stage 3 是等统计量条件样本；物理平均看 Stage 2 的随机 2γ 样本。 |
| 忽略零沉积稀释 | Stage 2 的 `efficiency_mean` / `merge_rate_primaries` / `split_rate_primaries` 分母包含零沉积初级；条件指标才是“可重建子样本”的性能。 |
| 以为 `merge` 和 `efficiency` 用同一套指派 | `merge/split` 用每个初级的 `argmax` 最佳团簇；`efficiency/purity` 用匈牙利指派。两者在复杂事例中可能不同。 |
| 以为 Stage 2 docstring 的 80 点是准的 | 代码是 7 个 seed、112 点；顶部 docstring 是旧版。 |

### 12.4 已知文档不一致

1. [pipeline.md](../pipeline.md) 描述的是旧版单粒子 API：
   - 旧签名 `reconstruction(filename, energy_threshold, n_neighbor_threshold, ...)`；
   - 旧 CLI `[Eth] [Nth]`；
   - 旧描述“加入团簇需要 `Nth` 个点亮邻居”。
   当前代码以 `reconstruction.py` 为准：`Eth_seed` + `Eth_wave` + `Nth`，`Nth` 只用于种子。

2. [stage2/evaluate.py](../results/multiparticle/stage2/evaluate.py) 顶部 docstring 写的是 “seeds {7,8,10,13,14} ... 80 points”；实际常量是 7 个 seed、112 点。

3. `docs/deprecated_files.md` 中提到的默认邻接表文件名（`ecal_neighbor_info_new_dedup.root`）与当前代码默认值（`ecal_neighbor_info.root`）不同；当前 Stage 2/3 的 `work/utilities/` 链接的是 `ecal_neighbor_info.root`。

---

## 13. 完整示例：一个 2γ 事例从输入到指标

以 Stage 2 的 `2gamma`、事件 0、参数 `seed=13, wave=0.05, Nth=0` 为例。

### 13.1 输入真值

```text
NTrack          = 2
Species         = [1, 1]
SourceIdx       = [89149, 72614]
PrimaryEnergy   = [53.085, 53.085] MeV
SourceEnergy    = [45.7241, 53.0849] MeV
SignalEnergy    = 98.8090 MeV
NBkg            = 2
BkgEnergyTotal  = 9.4273 MeV
Σ Edeps         = 108.2363 MeV
```

`Truth` 的 12 条稀疏 entry 中，前 5 条：

```text
(TruthCrystal=64,  TruthPrimary=0, TruthEnergy=0.9336)
(TruthCrystal=78,  TruthPrimary=0, TruthEnergy=0.4018)
(TruthCrystal=172, TruthPrimary=0, TruthEnergy=1.0866)
(TruthCrystal=280, TruthPrimary=0, TruthEnergy=0.3452)
(TruthCrystal=523, TruthPrimary=0, TruthEnergy=13.4171)
```

`TruthPrimary=1` 的条目对应第二个 gamma。`io.load_multiparticle_sample()` 把它们展开成：

```python
signal_arrays[0].shape == (2, 568)
signal_arrays[0][0, 64]  == 0.9336
signal_arrays[0][0, 523] == 13.4171
signal_arrays[0][1, ...] == 第二个 gamma 的贡献
```

### 13.2 重建输出

`reconstruction` 产生：

```text
EvtID = 0, ModIDList = [524, 64, 78, 172, 523, 280]
EvtID = 0, ModIDList = [449, 55, 82, 163, 298, 564]
```

两个团簇分别对应两个 gamma；`524` 和 `449` 是各自的种子，列表中第一个元素是种子。

### 13.3 评估层反算

```python
c_sig = [45.7241, 53.0849]   # 两个团簇的信号能量
c_bkg = [0.0, 0.0]           # 本底 hit 不在这两个团簇内
```

重叠矩阵：

```text
        primary 0   primary 1
cluster 0   45.7241       0.0
cluster 1   0.0          53.0849
```

指派：

```text
assignment = [0, 1]          # cluster 0 -> primary 0, cluster 1 -> primary 1
efficiency = [1.0, 1.0]
purity     = [1.0, 1.0]
merge      = False
split      = [False, False]
```

本底：

```text
n_fake          = 0
signal_purity   = 1.0
bkg_capture     = 0.0
lead_bkg_frac   = 0.0
```

能量：

```text
signal_capture        = 1.0000
lead_signal_containment = 0.5372   # 领头团簇是 cluster 1，53.0849 / 98.8090
energy_bias           = [0.0, 0.0] # 两个初级都匹配且 E_j > 20 MeV
```

### 13.4 这个例子说明了什么

- `Edeps` 只给重建提供“哪些晶体有能量”，重建输出只给“哪些晶体属于同一个团簇”；
- 每个团簇里到底有多少信号/本底，是评估层用 `Truth*` 对 `ModIDList` 求和得到的；
- 由于两个 gamma 的晶体集不重叠，重叠矩阵是对角的，匹配和合并判断都很干净；
- 如果两个 gamma 共享晶体，`TruthCrystal` 会出现同一晶体、不同 `TruthPrimary` 的两条 entry，重叠矩阵的同一列/行就会有多个非零，合并/分裂指标才有意义。

---

## 14. 关键文件清单

### 14.1 核心代码

| 文件 | 作用 |
|:--|:--|
| [reconstruction.py](../reconstruction.py) | 波前传播聚类；输入 `Edeps` + 邻接表，输出 `ECALClusters` |
| [stitching.py](../stitching.py) | 单粒子读取、拼接、本底稠密化、多粒子样本写出 |
| [geometry.py](../geometry.py) | 晶体位置、质心、距离、pitch |
| [evaluation/io.py](../evaluation/io.py) | 多粒子样本与团簇文件的读取 |
| [evaluation/matching.py](../evaluation/matching.py) | 重叠矩阵、指派、效率/纯度、合并/分裂 |
| [evaluation/background.py](../evaluation/background.py) | 假团簇、信号纯度、本底捕获 |
| [evaluation/energy.py](../evaluation/energy.py) | 信号捕获、领头簇收束、逐初级能量响应 |
| [evaluation/counts.py](../evaluation/counts.py) | 计数指标 |
| [evaluation/separation.py](../evaluation/separation.py) | 分离度分箱、合并半径插值 |
| [bkg_est_add/bkg_model.py](../bkg_est_add/bkg_model.py) | 本底模型与 `n_bkg` 计算 |
| [bkg_est_add/add_background.py](../bkg_est_add/add_background.py) | 本底抽样 `draw_injections()` |

### 14.2 Stage 2 驱动

| 文件 | 作用 |
|:--|:--|
| [stage2/generate.py](../results/multiparticle/stage2/generate.py) | 生成 2γ / e+2γ 的 sig/bkg 样本 |
| [stage2/evaluate.py](../results/multiparticle/stage2/evaluate.py) | 112 点重建 + 评估 |
| [stage2/run_stage2.sh](../results/multiparticle/stage2/run_stage2.sh) | 顺序执行 generate + evaluate |
| [stage2/meta.json](../results/multiparticle/stage2/meta.json) | 生成参数与自洽检查 |
| [stage2/README.md](../results/multiparticle/stage2/README.md) | Stage 2 报告与物理解读 |

### 14.3 Stage 3 驱动

| 文件 | 作用 |
|:--|:--|
| [stage3/generate.py](../results/multiparticle/stage3/generate.py) | 按分离度箱生成 2γ 样本 |
| [stage3/evaluate.py](../results/multiparticle/stage3/evaluate.py) | 28 点重建 + 逐箱评估 + 半径拟合 |
| [stage3/run_stage3.sh](../results/multiparticle/stage3/run_stage3.sh) | 顺序执行 generate + evaluate |
| [stage3/meta.json](../results/multiparticle/stage3/meta.json) | pitch、箱边界、能量/本底统计 |
| [stage3/README.md](../results/multiparticle/stage3/README.md) | Stage 3 报告与物理解读 |

### 14.4 数据与结果

| 路径 | 内容 |
|:--|:--|
| `data/raw/{ep_iso,gamma_iso}.root` | 原始 MC；`G4Run0/ECALSimHit` |
| `data/processed/edep_of_each_evt_{ep_iso,gamma_iso}.root` | 单粒子 `Edeps` |
| `data/utilities/ecal_neighbor_info.root` | 568 晶体邻接表 |
| `data/utilities/ecal_module_positions.npz` | 568 晶体中心 + hit 数 |
| `stage2/data/stitched_*_{sig,bkg}.root` | Stage 2 拼接样本 |
| `stage3/data/stage3_2gamma_{sig,bkg}.root` | Stage 3 分离度样本 |
| `stage2/results/stage2_metrics.csv` | 112 行指标 |
| `stage2/results/points/*.json` | 每点明细 + 直方图 |
| `stage2/results/best_points.json` | 最佳点 |
| `stage2/results/background_impact.json` | 本底配对影响 |
| `stage3/results/stage3_metrics.csv` | 168 行（28 点 × 6 箱） |
| `stage3/results/fits.json` | 28 条 50%/10% 半径 |

---

## 15. 参考结果与量级

### 15.1 样本量级

| 样本 | 事例数 | 信号 | N | 平均 `E_sig` |
|:--|--:|:--|--:|--:|
| `stitched_2gamma_*` | 20,000 | 2γ | 2 | 90.85 MeV |
| `stitched_ep2gamma_*` | 20,000 | e+2γ | 3 | 137.29 MeV |
| `stage3_2gamma_*` | 15,000 | 2γ | 2 | 103.13 MeV（条件样本） |

### 15.2 Stage 2 最佳点（bkg）

| 场景 | seed | wave | Nth | exact_signal | fake/event | merge(events) | split(primaries) |
|:--|--:|--:|--:|--:|--:|--:|--:|
| `2gamma` | 13 | 0.05 | 0 | 0.7098 | 0.00015 | 5.44% | 3.26% |
| `ep2gamma` | 10 | 0.05 | 0 | 0.5392 | 0.0032 | 15.72% | 3.54% |

### 15.3 Stage 3 主曲线（seed=13, Nth=0）

| 箱（pitch） | 中位分离 | 合并率 | exact-2 | 效率 | 匹配纯度 |
|:--|--:|--:|--:|--:|--:|
| 0–0.5 | 0.314 | 0.9992 | 0.0008 | 0.4962 | 0.5154 |
| 0.5–1 | 0.579 | 0.9984 | 0.0012 | 0.4957 | 0.5175 |
| 1–2 | 1.011 | 0.9920 | 0.0056 | 0.4955 | 0.5318 |
| 2–3 | 2.582 | 0.8328 | 0.1660 | 0.5755 | 0.6525 |
| 3–5 | 4.105 | 0.2268 | 0.7672 | 0.8671 | 0.9346 |
| >5 | 10.105 | 0.0012 | 0.9892 | 0.9758 | 0.9987 |

50% 合并半径 **3.410 pitch ≈ 157.7 mm**（signal-only）；+background 为 **3.418 pitch ≈ 158.1 mm**。10% 半径 **7.430 pitch ≈ 343.6 mm**（+background 7.477 pitch ≈ 345.8 mm）。

---

## 16. 小结

- **重建的输入契约很简单**：`processed/<file>` 的 `Edeps` + `utilities/ecal_neighbor_info.root`；重建只输出团簇成员，不输出能量。
- **多粒子测试的复杂度在真值记账**：`TruthCrystal/TruthPrimary/TruthEnergy` 用稀疏三元组保存“哪个初级在哪块晶体沉积了多少”，使得共享晶体、合并、分裂都能在评估层被精确分解。
- **Stage 2 回答“总体性能如何”**：2γ / e+2γ、有/无本底、112 个参数点，指标覆盖计数、匹配、本底、能量。
- **Stage 3 回答“多近会合并”**：按质心距离分箱，给出合并率曲线和 50%/10% 合并半径；阈值在 5–14 MeV、Nth 0/1 内几乎不改变曲线，限制来自几何重叠。
- **看结果时注意口径**：all-event 包含零沉积；条件指标只统计所有源 `E_sig>0`；`exact_signal` 用 `e_sig≥1 MeV` 区分信号簇和假团簇；`merge/split` 用最佳团簇，`efficiency/purity` 用匈牙利指派。

---

*文档完。代码与数据路径均以当前仓库为准；若后续 Stage 4/5 增加新的样本模式，应在本文件 §4 的分支表、§9 的重建契约和 §12 的检查项中同步扩展。*
