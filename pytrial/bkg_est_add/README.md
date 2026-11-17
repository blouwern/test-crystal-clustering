# 本底估计与注入（bkg_est_add）

本目录分两套独立的本底模型，方法一致（不直接模拟实验室本底，而是用「源项 + 解析物理 + 泊松统计」的保守数值方法），分别估计 **γ** 与 **中子** 本底，并注入 `data/processed/` 的预处理数据。

```
pytrial/bkg_est_add/
├── README.md          # 本文件：两套本底的对照与使用说明
├── gamma/             # γ 本底（30 μSv/h 剂量率 → n_bkg = 2.4247）
└── neutron/           # 中子本底（200 n/cm²/s 通量 → n = 0.1621，最保守上界）
```

## 命名约定

| 文件 | 含义 |
|:---|:---|
| `edep_of_each_evt_<dataset>.root` | 原始信号 |
| `edep_of_each_evt_<dataset>_bkg_g.root` | 注入 **γ** 本底 |
| `edep_of_each_evt_<dataset>_bkg_n.root` | 注入**中子**本底 |

两种注入文件都从原始信号文件独立生成，不级联；若要同时含两种本底，把两个文件的 `Edeps` 相加即可（随机流独立，互不影响）。

```bash
# γ：生成 _bkg_g
cd pytrial/bkg_est_add/gamma && python3 add_background.py

# 中子：生成 _bkg_n
cd pytrial/bkg_est_add/neutron && python3 add_neutron_background.py
```

## 两套本底的对照

| 项 | γ（`gamma/`）| 中子（`neutron/`）|
|:---|:---|:---|
| 源项 | 环境剂量率 30 μSv/h | 大厅中子通量 <200 n cm⁻² s⁻¹ |
| 归一化 | `P_total = D·M` | `r = Φ·S/4·ε` |
| 沉积模型 | 截断正态 N(5 MeV, 2 MeV²) | 俘获 γ 级联，`Q_γ·f_dep` |
| 空间分布 | 均匀随机 | 均匀随机（表层俘获，单晶注入）|
| 时间结构 | 事例窗内瞬时 | 连续束、稳态 |
| 每事例 hit 数 | **2.4247** | **0.1621** |
| 每事例能量 | **12.13 MeV** | **0.7121 MeV** |
| 假团簇 | 主要威胁（0.003/事例 @ 10 MeV 阈值）| **0**（单次上限 6.86 MeV < 7 MeV）|
| 分支 | `NBkg` | `NNeutron` |
| 后缀 | `_bkg_g` | `_bkg_n` |
| 随机种子 | 20261105 | 20261116 |
| 文档 | [`gamma/README.md`](gamma/README.md) | [`neutron/README.md`](neutron/README.md) |

**一句话**：γ 本底是"数得清但会造假簇"的本底；中子本底在纯 CsI 里只剩辐射俘获通道，是一个 **0.71 MeV/事例（≤2.2%）的能量污染，且不产生假团簇**。

## 口径提醒

- 两套模型都只处理**本底能量沉积**，不模拟时间信息；`Edeps` 是沉积能（现有 MC 的 `nOptPho` 全为 0，没有可见能）。
- γ 模型把剂量当量率直接当作 CsI 吸收剂量率（O(1) 系统误差）；中子模型则完全不用剂量系数，只依赖通量上限与俘获物理。
- 中子模型是**单边保守上界**：假设全部中子都是热中子、俘获概率为 1。真实值只会更小。
- 参数扫描代码（`pytrial/sweeps/config.py`）的 `BKG_SUFFIX` 已同步为 `_bkg_g`。
