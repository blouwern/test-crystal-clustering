# 已废弃文件清单（`_deprecated`）

> 整理时间：2026-09 · 原因：邻接表的 ModID 对齐问题被查清后，一批基于错误拓扑的产物与陈旧代码不再可信。
> 约定：这些文件**保留在磁盘上**供追溯，路径带 `_deprecated` 标记，并已通过 `.gitignore` 停止跟踪。

## 忽略规则（`.gitignore`）

```
**/*_deprecated*
**/*_deprecated/**
_deprecated/
```

## 1. 邻接表及其生成脚本

| 废弃文件 | 原因 |
|:--|:--|
| `data/utilities/ecal_neighbor_info_deprecated.root` | 原始邻接表。在 `reconstruction.py` 自己的索引假设下，边长中位数 **434.8 mm（9.4 倍晶距）**，与真实最近邻重合度仅 7.4%，0/568 完全吻合 |
| `data/utilities/ecal_neighbor_info_added_deprecated.root` | 同上，多一个 `second_neighbors` 分支；曾是 `reconstruction.py` 的默认值 |
| `data/utilities/ecal_neighbor_info_geo_deprecated.root` | 我临时用 MC 真值位置重建的**裁判用**邻接表（k=6 最近邻 + 二阶邻居）。物理自洽，但只是过渡产物，不应作为正式数据 |
| `data/utilities/add_secondary_neighbors_deprecated.py` | 为错误的 `ecal_neighbor_info.root` 生成二阶邻居 |
| `data/utilities/build_neighbor_info_from_geometry_deprecated.py` | 生成上面那张临时裁判表的脚本 |

**替代品**：`data/utilities/ecal_neighbor_info_new.root`（几何构建方对齐模拟参数后的新表，568 条目）
与其去重版 `ecal_neighbor_info_new_dedup.root`（`reconstruction.py` 当前默认值）。

## 2. 陈旧的算法副本

| 废弃目录 | 原因 |
|:--|:--|
| `algorithm_deprecated/` | 整个目录。`cpp/Reconstruction.C` 有 2 处语法错误**无法编译**（L16 多余的 `l`、L43 多余的 `g`），且选种循环 `while (state != 0) ++seed_rank;` 无越界保护；`cpp/CryInfo.H` 只服务于它；`python/reconstruction.py` 与 `pytrial/reconstruction.py` 已分叉 18 行（缺少 gap-jumping），`python/scoring.py` 与 `pytrial/scoring.py` 重复 |

**当前实现**：`pytrial/reconstruction.py` + `pytrial/scoring.py`。

## 3. 基于错误拓扑得到的结果

| 废弃路径 | 原因 |
|:--|:--|
| `pytrial/results/_deprecated/Eth10_Nth1{,_v2}/`、`Eth10_Nth2{,_v2}/` | 早期 1M 事例汇总，使用错误的邻接表 |
| `pytrial/results/_deprecated/SWEEP_ANALYSIS.md` | 基于错误拓扑的 108 点扫描分析 |
| `pytrial/results/_deprecated/sweep_26-09-14/` | 上述扫描的全部图表与 CSV |
| `pytrial/results/_deprecated/sweep_26-09-16_geo/` | 两阈值扫描（Nth=3），跑在临时几何表上 |
| `pytrial/results/_deprecated/nth_scan_geo/` | Nth 放宽扫描，跑在临时几何表上 |
| `pytrial/results/_deprecated/nth_scan_added/` | Nth 放宽扫描，跑在错误的邻接表上 |

**注意**：这些数字**不能**用于论文或开题报告。新表修好后需要重跑。

## 4. 保留（未废弃）

| 路径 | 说明 |
|:--|:--|
| `data/utilities/verify_neighbor_graph.py` | **文件检查程序**，T0–T6 检验套件 |
| `docs/neighbor_graph_verification.md` | 邻接表检验的完整流程与逻辑 |
| `docs/deprecated_files.md` | 本文件 |
| `pytrial/results/sweep/26-09-16/README.md`、`NTH_SCAN.md` | 扫描说明与结论（方法仍然有效） |
| `pytrial/results/sweep/26-09-16/topology_comparison.csv` | 新旧拓扑同参数对照 —— 这是定位该 bug 的**关键证据**，保留 |
| `pytrial/results/sweep/26-09-16/neighbor_graph_verification.json` | 检验结果原始数据 |
| `data/read_edep_of_each_evt.C` | 已逐元素验证能量写入正确（max\|diff\| = 0，全域相对偏差 2.3e-14） |
| `pytrial/_defective_algorithm/`、`pytrial/_deprecated_data_processing/` | 历史上已标记，维持现状 |

## 5. 本次同时做的兼容性修改（`pytrial/reconstruction.py`）

1. **默认邻接表**改为 `ecal_neighbor_info_new_dedup.root`（原来指向已废弃的 `_added`）。
2. **`second_neighbors` 改为可选**：新表只有 `neighbors` 一个分支，而原代码无条件调用
   `SetBranchAddress("second_neighbors", ...)`，会抛
   `ReferenceError: attempt to access a null-pointer`。现在先检查分支是否存在，缺失时打印
   `[INFO] ... gap-jumping fallback disabled` 并继续运行。

## 6. 从 git 恢复

`_deprecated` 内容仍在提交历史中，需要时：

```bash
git checkout <commit> -- <path>          # 取回单个文件
git log --all --oneline -- <path>        # 查历史
```
