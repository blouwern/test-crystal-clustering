# ECAL 聚类算法改进方向建议

> 目的：参考使用闪烁体/晶体量能器的粒子物理实验（ATLAS、CMS、BESIII、Belle II、LHCb 等）的物理事例重建算法，结合本仓库 `pytrial/results/SWEEP_ANALYSIS.md` 中参数扫描暴露出的具体短板，给出可落地的算法改进方向。
>
> 每条技术都给出 **① 文章标题与链接 ② 与该技术相关的原文片段（逐字引用） ③ 与本项目的对应关系**。
>
> 引用说明：所有引文均来自调研过程中实际抓取的页面（arXiv HTML / ar5iv / 官方文档）。抓取工具无法读取 PDF，因此个别论文的某些小节未能取到，凡未能亲自核实的都显式标注 **UNVERIFIED**。文末附参考文献总表。

---

## 0. 摘要：优先级排序

扫描结论（详见 `pytrial/results/SWEEP_ANALYSIS.md`）指向三个结构性缺陷：**过度分裂、阈值能量损失、簇间无能量共享**。据此给出优先级：

| 优先级 | 方向 | 解决什么问题 | 是否需要 ML | 预计收益 |
|:--:|:--|:--|:--:|:--|
| **P0** | A. 局部极大值分裂 + 簇间能量共享 | `Eth=5` 时 54% 事例 ≥2 簇、15% 能量在副簇 | 否 | 同时改善多重数与能量 |
| **P0** | B. 噪声显著度阈值取代全局 `Eth` | 全局阈值对 e/γ/μ 无适应性；μ 在高阈值崩溃 | 否 | 提高对粒子种类的鲁棒性 |
| **P1** | C. 时间信息入判据 | 无法区分同簇内的无关沉积 | 否 | 抑制噪声/堆积 |
| **P1** | D. 径迹-簇匹配与 split-off 判别 | 无法区分"真中性"与"带电簇碎片" | 否 | 粒子级指派 |
| **P1** | E. 评价指标升级（V-score + splitting rate） | 只看簇数，测不出分裂/合并 | 否 | **让后续所有改进可测量** |
| **P2** | F. 分层/几何/图结构重构 | 丢弃了位置、方向、簇形状 | 否 | 能量修正与形状判别的前置 |
| **P3** | G. 机器学习（GNN / 边分类器） | 需要更高上限时的路线 | 是 | 上限高，成本高 |

**建议顺序：E → A → B → C/D → F → G。** 先换指标（E），否则无法证明 A/B 有效；再做分裂与共享（A），这是当前最大的一块损失。

---

## 1. 现状与成熟实验的差距对照

| 能力 | 本仓库现状 | ATLAS | CMS ECAL | BESIII | Belle II | LHCb |
|:--|:--|:--|:--|:--|:--|:--|
| 种子定义 | 能量降序 + 全局阈值 | 显著度 ζ=E/σ > 4 | 局部极大值 + E_T 阈值 | 局部极大值 | 局部极大值 | 局部极大值（E_T>50 MeV） |
| 生长判据 | ≥Nth 个已点亮邻居 | 4σ/2σ/0σ 三级 + 拓扑连接 | 邻接 + 扫描方向单调性 | 连通区 | 邻居 >0.5 MeV，节点 >10 MeV 才继续 | 有向边指向种子 |
| 分裂 | **无** | 局部极大值分裂 | 一个连通区多个种子 | — | 每个局部极大值一个簇 | 弱连通分量 |
| 能量共享 | **无**（硬指派） | 边界单元加权共享 | 高斯分数共享 | CG 加权 | 对数加权 + 迭代收敛 | 重叠单元按能量分数 |
| 噪声模型 | **无**（靠 `Eth`） | 逐单元 σ_noise | 单晶 spike 抑制 E4/E1、E6/E2 | EINC 标定 | 事件级噪声估计 | 阈值 |
| 时间 | **未用** | 单元时间判据 | ±2 ns | <700 ns | 125 ns 窗口 | — |
| 径迹匹配 | **未用** | PF 单元级减除 | PF link | 径迹-簇匹配 | — | — |
| 评价指标 | 簇数 vs 粒子数 | — | — | — | — | — |

---

## 2. 方向 A（P0）：局部极大值分裂 + 簇间能量共享

当前算法把一个连通区强行指派成 1 个簇（每个晶体只能属于一个簇），这是"过度分裂 + 能量损失"的直接原因。成熟实验的共识是：**一个局部极大值对应一个粒子，晶体能量可以按距离加权分给多个簇。**

### A1. Belle II ECL：从局部极大值生长连通区 + 对数加权迭代共享（与本项目最接近的基线）

- **来源**：F. Wemmer et al. (Belle II), *Photon Reconstruction in the Belle II Calorimeter Using Graph Neural Networks*, arXiv:2306.04179, §5.1 "Baseline" — <https://arxiv.org/abs/2306.04179> / <https://ar5iv.labs.arxiv.org/html/2306.04179>
- **原文片段**：

  > "The clustering is performed in three steps. In the first step, all crystals are grouped into a connected set of crystals, so-called connected regions starting with LMs, as defined previously. In an iterative procedure all direct neighbors with energies above 0.5 MeV are added to this LM, and the process is continued if any neighbor itself has energy above 10 MeV. Overlapping connected regions are merged into one."

  > "If there is more than one LM in a connected region, the energy in each crystal of the connected region is assigned a distance-dependent weight and can be shared between different clusters. The distance is calculated from the cluster centroid to each crystal center, where the cluster centroid is updated iteratively using logarithmic energy weights. This process is repeated until all cluster centroids in a connected region are stable within 1 mm."

- **与本项目的对应**：这是把我们的"生长"直接升级为"分裂 + 共享"的模板。我们的 `Eth` 对应它的 0.5 MeV 加入阈值，`Nth` 对应它"节点能量 >10 MeV 才继续传播"的规则——**把"Nth 个已点亮邻居"换成"当前晶体是否足够亮（显著度）以继续波前"，就同时获得了分裂能力和一个更合理的生长判据**。共享权重只需晶体坐标（622 个模块的几何是现成的），完全不需要 ML。

### A2. CMS 粒子流聚类：一个连通区多个种子 + 高斯分数共享

- **来源**：CMS Collaboration, *Electron and photon reconstruction and identification with the CMS experiment at the CERN LHC*, JINST 16 (2021) P05014, arXiv:2012.06888, §0.4.1–0.4.2 — <https://arxiv.org/abs/2012.06888> / <https://ar5iv.labs.arxiv.org/html/2012.06888>
- **原文片段**：

  > "Energy deposits in several ECAL channels are clustered under the assumption that each local maximum above a certain energy threshold (1 GeV) corresponds to a single particle incident on the detector. An ECAL energy deposit may be shared between overlapping clusters, and a Gaussian shower profile is used to determine the fraction of the energy deposit to be assigned to each of the clusters."

  > "The energy reconstruction algorithm starts with the formation of clusters by grouping together crystals with energies exceeding a predefined threshold (typically ∼80 MeV in EB and ∼300 MeV in EE), which is generally 2 or 3 times bigger than the electronic noise expected for these crystals. A seed cluster is then defined as the one containing most of the energy deposited in any specific region, with a minimum transverse energy (E_T^seed) above 1 GeV."

- **补充（算法实现细节）**：CMS 官方文档 *SWGuideParticleFlowClustering* — <https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideParticleFlowClustering>
  > "It is possible to have several seeds in the same topological cluster." … "In case there are several seeds in the topological cluster, each of them will give a PFCluster. The energy of PFRecHits will give a contribution to each cluster. For each PFRecHit a fraction of its energy attributed to a given PFCluster is computed. … The shower shape is assumed to be gaussian. … The algorithm is iterative, starting from seed energy and position. At each step updated PFCluster is used to re-compute PFRecHitFractions. Iterations end when PFClusters positions are stable."
- **与本项目的对应**：**这是修"不分裂 / 不共享"这两个缺陷最直接的模板**。改动很小：允许一个连通区里有多个种子；每个模块按到各种子的高斯距离给出分数权重，迭代到簇位置稳定；簇数自然等于局部极大值数——而局部极大值数正是我们评价指标想衡量的量。所需数据只有模块位置与能量。

### A3. LHCb Graph Clustering：图结构化的分裂、重叠单元与能量分数

- **来源**：N. Valls Canudas et al. (LHCb), *Graph Clustering: a graph-based clustering algorithm for the electromagnetic calorimeter in LHCb*, LHCb-DP-2022-003, arXiv:2212.11061, §3 — <https://arxiv.org/abs/2212.11061> / <https://ar5iv.labs.arxiv.org/html/2212.11061>
- **原文片段**：

  > "By design, the target nodes of all edges in the graph are the seeds of the reconstructed clusters, where a seed is defined as a local maximum energy digit in the calorimeter grid over a threshold of 50 MeV in transverse energy. With this, the cluster seeds can be easily identified as nodes with only incoming edges."

  > "Furthermore, a node can be linked to more than one seed if it is susceptible to have energy deposits from more than one particle. These particular cases are called overlap cells."

  > "It searches for overlap vertices, identified by having two or more output edges, and accumulates the energy of all the connected nodes on all the clusters involved in the overlap, including the energy of the overlapping node equally fractioned for every involved cluster. Then, a weight is computed for every overlapping edge as the fraction between the energy of the target cluster and the sum of all the clusters involved in the overlap."

  > "It outperforms the previously used method by 65.4% in terms of computational time on average, with an equivalent efficiency and resolution."

- **与本项目的对应**：把"波前传播"形式化为"以种子为汇的有向图 + 弱连通分量"，就免费得到了分裂、重叠单元分数化和后续 ML 的挂载点。**LHCb 报告同等物理性能下快 65.4%**，说明这种重构不带来物理风险。我们的输出本来就是 `ModIDList`，改成"每模块带权重"只需扩展输出 schema。

### A4. 合并 π⁰ 情形：用两个最亮单元的比值决定是否扩大簇

- **来源**：同 A3，arXiv:2212.11061, §3.2.1
- **原文片段**：
  > "That is why the Graph Clustering algorithm adapts the shape of potential merged π⁰ candidates, expanding the cluster up to the neighbours of the second most energetic digit in the cluster. … we have studied the relation between the two most energetic digits as a ratio labeled R1 to define a threshold on which to make the cluster expansion. … the threshold for R1 is set to 25. The chosen value ensures that the residual energy left outside the cluster is less than 9% for the studied π⁰ samples …"
- **与本项目的对应**：给出了"局部极大值无法处理（两个极大值挨得太近）"时的具体处理规则，更重要的是给出了**参数整定方法论**：用带标签的 MC 扫描判别变量，把阈值定在"残留能量 < 某个目标"上，而不是凭经验挑。我们目前 `Eth`/`Nth` 的选取正缺少这种目标函数。

---

## 3. 方向 B（P0）：用噪声显著度取代全局能量阈值

扫描显示单个全局 `Eth` 无法同时服务 e⁺/γ/μ⁻。成熟实验的做法是把阈值表达为**以该通道自身噪声为单位的显著度**。

### B1. ATLAS 拓扑聚类：ζ = E/σ_noise 与 4σ/2σ/0σ 三级阈值

- **来源**：ATLAS Collaboration, *Improving topological cluster reconstruction using calorimeter cell timing in ATLAS*, Eur. Phys. J. C 84 (2024) 455, arXiv:2310.16497, §4 — <https://arxiv.org/abs/2310.16497> / <https://arxiv.org/html/2310.16497v2>
- **原文片段**：
  > "Topo-clusters are formed by a growing-volume algorithm, configured by three threshold parameters {S, N, P}:" — 主种子阈值 \|ζ_cell^EM\| > S、生长控制阈值 \|ζ_cell^EM\| > N、主单元过滤 \|ζ_cell^EM\| > P — "set to {S=4, N=2, P=0} in Runs 1, 2 and 3. The algorithm is seeded by cells … whose signal significance exceeds a threshold S … Seed cells are then sorted from highest to lowest energy significance and topo-clusters are grown by adding all neighbouring cells that satisfy Eq. (5)."

- **补充来源**：ATLAS Collaboration, *Jet reconstruction and performance using particle flow with the ATLAS Detector*, Eur. Phys. J. C77 (2017) 466, arXiv:1703.10485, §5 — <https://ar5iv.labs.arxiv.org/html/1703.10485>
  > "Topo-clusters of calorimeter cells are seeded by cells whose absolute energy measurements |E| exceed the expected noise by four times its standard deviation. The expected noise includes both electronic noise and the average contribution from pile-up … The topo-clusters are then expanded both laterally and longitudinally in two steps, first by iteratively adding all adjacent cells with absolute energies two standard deviations above noise, and finally adding all cells neighbouring the previous set."

- **与本项目的对应**：把单一的 `Eth` 换成**逐模块的 σ_noise + 三级显著度**：≥4σ 才能当种子，≥2σ 才能继续生长，≥0σ 做最后一圈收集。622 个模块的噪声表可以从标定/噪声跑数或 MC 得到，成本极低。**这直接解决"一个阈值服务不了三种粒子"的问题**，也是唯一能让"低能但成簇"的沉积（如 μ 的 MIP）保留下来而不引入噪声的机制。

### B2. ATLAS 拓扑噪声抑制的原理

- **来源**：ATLAS Collaboration, *Topological cell clustering in the ATLAS calorimeters and its performance in LHC Run 1*, Eur. Phys. J. C77 (2017) 490, arXiv:1603.02934, Abstract — <https://arxiv.org/abs/1603.02934>
- **原文片段**：
  > "The cluster formation follows cell signal-significance patterns generated by electromagnetic and hadronic showers. In this, the clustering algorithm implicitly performs a topological noise suppression by removing cells with insignificant signals which are not in close proximity to cells with significant signals."
  > "Calorimeter cells with insignificant signals found to not be connected to neighbouring cells with significant signals are considered noise and discarded from further jet, particle and missing transverse momentum reconstruction."
- **与本项目的对应**：**噪声抑制应该是拓扑的（孤立 + 不显著 = 噪声），而不是纯能量的。** 我们现在的 `Eth` 把"低能但紧邻主簇"的晶体也一并砍掉——这正是第 4 节量化的 12–41% 能量损失来源之一。

### B3. BESIII：噪声水平必须与簇形状一起标定（EINC）

- **来源**：V. Prasad, C. Liu, X. Ji, W. Li, H. Liu, X. Lou (BESIII), *Study of BESIII electromagnetic calorimeter performance with radiative lepton pair events*, arXiv:1606.00248, §4.1 — <https://arxiv.org/abs/1606.00248> / <https://ar5iv.labs.arxiv.org/html/1606.00248>
- **原文片段**：
  > "A number of variables have been developed to study the shower shapes of the different particles in the EMC, such as second moment and lateral moment (LAT). These two variables are used to quantify the transverse shower shape of the cluster and separate the EM showers from the hadronic showers."
  > "The EM showers tend to deposit a large fraction of their energy in one or two crystals, whereas the hadronic showers tend to be more spread out."
  > "The discrepancy between data and MC has been overcome while simulating a new MC sample with the EINC value of 0.27 MeV."
- **与本项目的对应**：两点可直接借用：**(a) 用横向矩（lateral moment / second moment）作为簇级判别变量**，这是区分电磁簇与 μ 类沉积的廉价手段；**(b) 噪声常数本身必须用簇形状数据反向标定**（BESIII 从 0.20 调到 0.27 MeV），否则显著度阈值会被系统性带偏。我们做显著度阈值之前，必须先有这一步标定。

### B4. CMS ECAL：单晶噪声/spike 的邻域一致性抑制

- **来源**：CMS Collaboration, *Particle-flow reconstruction and global event description with the CMS detector*, JINST 12 (2017) P10003, arXiv:1706.04965, §0.2.3 — <https://arxiv.org/abs/1706.04965> / <https://ar5iv.labs.arxiv.org/html/1706.04965>
- **原文片段**：
  > "Since these spikes mostly affect a single crystal and more rarely two neighbouring crystals, they are rejected by requiring the energy deposits to be compatible with arising from a particle shower: the ratios E4/E1 and E6/E2 should exceed 5% and 10% respectively, where E1 (E2) is the energy collected in the considered crystal (crystal pair) and E4 (E6) is the energy collected in the four (six) adjacent crystals."
  > "The timing of the energy deposits in excess of 1 GeV is also required to be compatible with the beam crossing time to better than ±2 ns."
- **与本项目的对应**：**一个下午就能加上的改动**：在聚类之前先剔除"能量与其 4/6 邻居不相容"的孤立亮晶（spike/噪声），这会直接减少 `Eth=5` 时那 54% 的过度分裂（噪声种子）。注意它同时给了时间窗口 ±2 ns 的用法。

### B5. CLUE：用局部能量密度替代"能量 + 邻居计数"

- **来源**：M. Rovere, Z. Chen, A. Di Pilato, F. Pantaleo, C. Seez, *CLUE: A Fast Parallel Clustering Algorithm for High Granularity Calorimeters in High Energy Physics*, Front. Big Data 3 (2020) 591315, arXiv:2001.09761, §2.2 — <https://arxiv.org/abs/2001.09761> / <https://pmc.ncbi.nlm.nih.gov/articles/PMC8080903/>
- **原文片段**：
  > "CLUE requires the following four parameters: dc is the cut-off distance in the calculation of local density; ρc is the minimum density to promote a point as a seed or the maximum density to demote a point as an outlier; δc and δo are the minimum separation requirements for seeds and outliers, respectively. The choice of these four parameters can be based on physics: for example, dc can be chosen based on the shower size and the lateral granularity of detectors; ρc can be chosen to exclude noise; δc and δo can be chosen based on the shower sizes and separations."
  > "Outliers and their descendant followers are guaranteed not to receive any cluster indices from seeds, which grants a noise rejection as shown in Figure 3."
- **补充（参数与探测器粒度的关系）**：E. Brondolin, M. Rovere, F. Pantaleo, *The k4Clue package…*, arXiv:2311.03089, §4.2 — <https://ar5iv.labs.arxiv.org/html/2311.03089>
  > "the critical distance, dc, … is chosen to take into account the granularity of the detector's geometry … Therefore, this is set to approximately two times the size of a single cell, dc = 40.0 mm."
  > "We then applied a pre-filter on the calorimeter cell energy at 2σ_noise, similar to the one applied in the clusterization process of the topological algorithm."
- **与本项目的对应**：把"Nth 个已点亮邻居"这种**计数式**判据换成**局部能量密度式**判据（在半径 dc 内求和），参数由簇的横向尺度与晶体粒度定出。对我们的 622 个模块，密度和"最近更高密度点"都是 O(N) 可算。这也解释了为什么 `Nth` 在我们的扫描里几乎无效：**邻居计数本身携带的物理信息太少。**

---

## 4. 方向 C（P1）：把时间信息引入聚类判据

### C1. ATLAS：聚类内部的单元时间判据

- **来源**：arXiv:2310.16497, Abstract 与 §2.1 — <https://arxiv.org/html/2310.16497v2>
- **原文片段**：
  > "The preferred version is found to reduce the out-of-time pile-up jet multiplicity by ∼50% for jet pT ∼20 GeV and by ∼80% for jet pT ≳50 GeV, while not disrupting the reconstruction of hadronic signals of interest, and improving the jet energy resolution by up to 5% for 20 < pT < 30 GeV."
  > "The cell time is only measured if the detected energy is above a certain configurable threshold. Typically, threshold values equal to three times the cell noise (3σ_noise) are used. If the reconstructed energy is below threshold, then the cell time is not computed and t = 0 is stored."
  > "In Run 2 and after applying calibration constants obtained from data, the constant term p0 was found to reach ∼200 ps, while the noise term p1 was O(1 GeV ns)."
- **与本项目的对应**：**只在 >3σ_noise 时才测量时间**，然后把时间与事例起始时间不相容的模块标记/剔除；对高能沉积再关闭该判据（避免误杀）。我们的晶体量能器本来就数字化脉冲波形，原始数据里也有 `t` 分支——这是"未用时间信息"这一短板的低成本解法。

### C2. HGCAL：簇级时间（单点时间不够用）

- **来源**：A. Lobanov, *Precision timing calorimetry with the CMS HGCAL*, JINST 15 (2020) C07003, arXiv:2005.13324, §2 — <https://arxiv.org/abs/2005.13324> / <https://ar5iv.labs.arxiv.org/html/2005.13324>
- **原文片段**：
  > "With the studies presented in the HGCAL TDR, it was shown that rejecting hits in the tails provides robust pileup rejection. Further, taking only cells within 2 cm of the shower axis, an efficiency of 100 % for photons is achievable with a resolution of 20 ps for pT>2 GeV."
  > "In order to efficiently reject particles originating from pileup, precision timing information of the order of 30 ps will be of great benefit. … the HGCAL will provide timing measurements for individual hits with signals above 12 fC … such that clusters resulting from particles with pT > 5 GeV should have a timing resolution better than 30 ps."
- **与本项目的对应**：**用整个簇（或只取靠近簇轴/种子的核心模块）做加权时间平均，而不是逐模块判断**；并对时间分布的尾部做剔除。这为将来处理多粒子/堆积事例提供了现成配方。

### C3. CMS ECAL：晶体量能器的时间分辨率量级

- **来源**：CMS Collaboration, *Time Reconstruction and Performance of the CMS Electromagnetic Calorimeter*, JINST 5 (2010) T03011, arXiv:0911.4044, Abstract 与 §0.6 — <https://arxiv.org/abs/0911.4044> / <https://ar5iv.labs.arxiv.org/html/0911.4044>
- **原文片段**：
  > "The resulting time resolution measured by lead tungstate crystals is better than 100 ps for energy deposits larger than 10 GeV."
  > "Such backgrounds are cosmic rays, beam halo muons, electronic noise, and out-of-time proton-proton interactions."
  > "As an example, 1 GeV energy deposits in the ECAL barrel have a time resolution of 1.5 ns."
- **与本项目的对应**：给出了晶体量能器时间分辨率的现实量级（低能 ~1.5 ns，>10 GeV <100 ps），可用于判断本探测器的时间判据是否可行，以及阈值该设在什么能量以上。

---

## 5. 方向 D（P1）：径迹-簇匹配与 split-off 判别

### D1. Pandora SDK：先过分裂、再按拓扑规则合并；用动量-能量一致性触发重聚类

- **来源**：J. S. Marshall, M. A. Thomson, *The Pandora Software Development Kit for Pattern Recognition*, Eur. Phys. J. C75 (2015) 439, arXiv:1506.05348, §1/§3/§11.1 — <https://arxiv.org/abs/1506.05348> / <https://ar5iv.labs.arxiv.org/html/1506.05348>
- **原文片段**：
  > "This design promotes an approach using many decoupled algorithms, each addressing specific topologies."
  > "The algorithms that address the problem must be able to build clusters of space-points and should be able to manipulate clusters by splitting them up or merging them together."
  > "The Clustering algorithm is configured so that it tends to split CaloHits from individual particles into multiple Clusters, rather than risk merging energy deposits from multiple particles into single Clusters. The Clusters are instead carefully merged together by a series of algorithms implementing well-defined topological rules."
  > "The compatibility of associated Tracks and Clusters is assessed, via comparison of Track momentum with associated Cluster energies. Significant discrepancies indicate pattern recognition problems and the reclustering approach described in Section 9 is used to improve the clustering."
  > "Clusters without associated tracks are examined to assess whether they genuinely represent electrically neutral particles, or whether they are more likely to be fragments of any nearby track-associated Clusters, representing charged particles."
- **与本项目的对应**：四个可直接借用的思想：**(i) 宁可先过度分裂，再用明确定义的拓扑规则合并**（我们现在的策略恰好相反：生长保守、靠阈值防合并）；**(ii) 用"径迹动量 vs 簇能量"的不一致触发重聚类**——我们有 e⁺/μ⁻ 样本，立刻可做；**(iii) "真中性 vs 带电簇碎片"的判别**就是 split-off 抑制的一般形式；**(iv) 注意 Pandora 明确禁止一个 hit 被两个簇共用**，而 Belle II/LHCb 允许分数共享——这是一个必须显式做出的设计选择。

### D2. BESIII：中性簇的定义与径迹匹配窗口

- **来源**：BESIII Collaboration, arXiv:2502.20821, §3 — <https://arxiv.org/abs/2502.20821> / <https://arxiv.org/html/2502.20821v2>
- **原文片段**：
  > "Neutral showers are reconstructed in the EMC. Showers not associated with any charged track are identified as photon candidates if they satisfy two additional criteria: (1) an energy deposition in the EMC of E_dep > 25 MeV in the barrel region corresponding to the polar angle |cos θ| < 0.8, while E_dep > 50 MeV in the end-cap region corresponding to 0.86 < |cos θ| < 0.92. (2) the EMC time difference from the event start time is required to be less than 700 ns to suppress electronic noise and showers unrelated to the event."
- **补充（匹配窗口）**：arXiv:1606.00248, §4.2
  > "The regions of θ_γ,γ_pred < 0.5 radian and E_γ/E_γ_pred ∈ [0.15,1.4] are considered to be the detected region of the photon in the EMC."
- **与本项目的对应**：给出了一套**完整、廉价**的粒子指派方案：与径迹匹配的簇算带电粒子；其余满足"区域相关能量阈值 + 时间窗口"的算光子。注意它是**区域相关阈值**（barrel 25 MeV / endcap 50 MeV）——这正是我们单一全局 `Eth` 的对照。

---

## 6. 方向 E（P1）：评价指标升级（**建议最先做**）

### E1. 能量加权 V-score（homogeneity + completeness）

- **来源**：M. Al Halabi et al., *Machine Learning Power Week 2023: Clustering in Hadronic Calorimeters*, arXiv:2508.09938, §2.4 — <https://arxiv.org/abs/2508.09938> / <https://ar5iv.labs.arxiv.org/html/2508.09938>
- **原文片段**：
  > "Two possible metrics are the homogeneity, which measures the extent to which each cluster contains hits from a single particle, and the completeness, which measures the extent to which all hits from a given particle are captured in a single cluster."
  > "We choose to weight hits linearly in the energy deposited in the hit, such that incorrectly clustering low-energy hits insubstantially affects the weighted Vscore, compared with clustering high-energy hits from different particles or splitting the hits in a single high-energy particle shower."
  > "This metric was chosen as it requires that a single algorithm must have good homogeneity and completeness in order to perform well - that is an homogeneity or completeness of 0 gives a V-score of 0, with a maximum possible validity score of 1."
- **与本项目的对应**：我们的"簇数 vs 粒子数"完全看不出**模块去了哪个簇**。能量加权的 V-score 把"过度分裂"（completeness 低）与"错误合并"（homogeneity 低）分开度量，且只需真值级的"粒子→模块"关联（单粒子 MC 里所有沉积都属同一粒子，可直接算 completeness）。**这是让 A/B/C/D 每一项改进变得可测量的前提。**

### E2. splitting rate：显式度量"一个光子被重建为两个"

- **来源**：Y. Maidannyk, F. Couderc, J. Malclès, M. Ö. Sahin, *Reconstruction of overlapping electromagnetic showers in calorimeters using Transformers*, Eur. Phys. J. C 86, 873 (2026), arXiv:2603.18172, Abstract — <https://arxiv.org/abs/2603.18172>
- **原文片段**：
  > "Performance is evaluated using efficiency, energy and position resolutions and splitting rate - reconstruction of two objects for a single photon. … For boosted π⁰ → γγ, the attention-based model retains di-photon mass reconstruction capability, where the standard algorithm becomes inefficient."
- **与本项目的对应**：我们恰好缺这个指标。扫描中 `Eth=5` 时 54% 事例 ≥2 簇，但"是过度分裂还是真有两个粒子"无法区分——splitting rate 正是为此设计的。

---

## 7. 方向 F（P2）：几何 / 分层 / 图结构

### F1. Arbor：把簇当成"树"，完全不使用能量

- **来源**：M. Ruan, *Arbor, a new approach of the Particle Flow Algorithm*, arXiv:1403.4784 — <https://arxiv.org/abs/1403.4784> / <https://ar5iv.labs.arxiv.org/html/1403.4784>
- **原文片段**：
  > "Arbor algorithm, as its name suggests, is inspired by the fact that the shower spatial development follows the topology of a tree. … Arbor reconstruct the long travelling charged particles generated at hadronic shower as the branches of the tree, while the tree structure is ensured with the constrain that no loop structure is allowed."
  > "In its current configuration, Arbor use only geometry information of the hits, i.e, no energy information is used."
  > "For all the connectors end at this hit, Arbor kept at most one connector that has the minimal angle to the reference direction. Therefore, no loop structure will be kept after the cleaning, and a tree structure based on the connectors emerges."
  > "Since the seed and the tree has one one correspondence, Arbor is a powerful algorithm to separate nearby showers, which is highly appreciated at Particle Flow Algorithm."
- **与本项目的对应**：一条**完全绕开 `Eth`/`Nth`** 的路线：只在相邻模块间建有向连接，每个模块只保留"与局部参考方向夹角最小"的那条连接，剩下的每棵树就是一个簇。每个簇天然带方向（种子→叶），正好补上我们缺失的簇形状/方向信息，而且对 622 个模块是 O(1)/节点。

### F2. HGCAL：逐层 2D 聚类 + 跨层关联 + PCA 簇变量

- **来源**：T. Cuisset on behalf of CMS, *Machine Learning for Event Reconstruction in the CMS Phase-2 High Granularity Calorimeter Endcap*, arXiv:2510.01851, §2 — <https://arxiv.org/abs/2510.01851>
- **原文片段**：
  > "To reduce computational complexity, rechits are first clustered within each calorimeter layer using the CLUE algorithm [7], forming two-dimensional clusters corresponding to transverse shower slices. These are then combined by the CLUE3D algorithm into tracksters, three-dimensional objects that represent full particle showers."
  > "Each trackster stores variables such as total energy, position, timing, and the results of a Principal Component Analysis (PCA) applied to its constituent 2D clusters, providing estimates of the shower's direction, length, and radius."
- **与本项目的对应**：如果我们的 ECAL 有纵向分层，逐层 2D 聚类再跨层关联比单次 3D 波前更抗纵向噪声，并且顺带产出能量修正所需的簇矩（轴、长度、半径）。

### F3. CMS ECAL 的 hybrid / island 聚类：固定条 + 动态扫描

- **来源**：CMS Collaboration, *Performance of photon reconstruction and identification with the CMS detector in pp collisions at sqrt(s) = 8 TeV*, JINST 10 (2015) P08010, arXiv:1502.02702, §0.4.2 — <https://arxiv.org/abs/1502.02702> / <https://ar5iv.labs.arxiv.org/html/1502.02702>
- **原文片段**：
  > "Clusters are built starting from a “seed crystal”: one containing a signal corresponding to a transverse energy greater than those of all its immediate neighbours and above a predefined threshold."
  > "In the barrel, where the crystals are arranged in an (η,ϕ) grid, the clusters have a fixed width of five crystals centred on the seed crystal, in the η direction. In the ϕ direction, adjacent strips of five crystals are added if their summed energy is above another predefined threshold."
  > "The RNINE variable is defined as the energy sum of the 3×3 crystals centred on the most energetic crystal in the supercluster divided by the energy of the supercluster."
- **补充（island 扫描的单向单调判据）**：CMS 官方文档 *SWGuideEcalRecoClustering* — <https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideEcalRecoClustering>
  > "For a crystal to be added to a cluster the following must be true: The crystal contains a rechit with positive energy. The crystal has not been assigned to another cluster already. The previous crystal added (in the same direction) has higher energy."
  > "Approximately 94% of the incident energy of a single electron or photon is contained in 3x3 crystals, and 97% in 5x5 crystals."
- **与本项目的对应**：**"沿扫描方向的前一个晶体能量必须更高"是一条极便宜的、确定性的生长终止/分界规则**，不需要调全局阈值就能在簇边界自动停下，也能在同能双峰间断开。另外 94%/97% 的横向包容数字可作为我们评估"收束率"上限的参照。

---

## 8. 方向 G（P3）：机器学习路线（分阶段）

### G1. 最小 ML 升级：把"邻居是否加入"变成一个可学习的边分类器

- **来源**：arXiv:2508.09938, §4.4 — <https://ar5iv.labs.arxiv.org/html/2508.09938>
- **原文片段**：
  > "The algorithm works in two stages: a) seed-to-seed classification, then b) seed-to-nonseed segmentation." … "any edge classifier, whether a simple multi-layer perceptron (MLP), GNN or transformer, can be used to prune spurious connections between a set of high-energy seeds and the rest of the energy deposits."
  > "The core component of this model, the edge-classifier MLP, achieved a best efficiency of 0.956, purity of 0.956, and an area under the receiver operator characteristic curve (AUC) of 0.998."
- **与本项目的对应**：**保留我们的种子与连通分量生长，只把"该邻居是否加入"这个布尔判断换成学到的边分数。** 这是从纯规则最小代价走向 ML 的路径，且不改变输出格式。

### G2. 端到端：GNN + object condensation

- **来源**：S. R. Qasim et al., *End-to-end multi-particle reconstruction in high occupancy imaging calorimeters with graph neural networks*, Eur. Phys. J. C82 (2022), arXiv:2204.01681 — <https://arxiv.org/abs/2204.01681> / <https://ar5iv.labs.arxiv.org/html/2204.01681>
- **原文片段**：
  > "The algorithm exploits a distance-weighted graph neural network, trained with object condensation, a graph segmentation technique. Through a single-shot approach, the reconstruction task is paired with energy regression."
  > "For instance, when two particles (p1 and p2) are maximally overlapping, we merge them into a single particle in the ground truth since disentangling such two clusters would be an impractical and imprecise task."
- **与本项目的对应**：给出"分裂/合并没有阈值"的终极形态（每 hit 学一个吸引/排斥势 + seediness），并且**最后一句对我们的评价体系是重要警告**：完全重叠的两个粒子在真值层面本就不可分，评价时必须定义"不可分辨即合并"的判据，否则会惩罚算法去做不可能的事。

### G3. Transformer 聚类与 overlapping shower

- **来源**：arXiv:2603.18172, Abstract — <https://arxiv.org/abs/2603.18172>（同 E2）
- **原文片段**：
  > "The study includes a two-step strategy in which candidate seed windows are identified and then jointly processed via distance-weighted message passing or attention mechanism and a single-step graph transformer, ClusTEX, which performs candidate selection and reconstruction in one inference stage."
- **与本项目的对应**：保留了"经典种子窗口"作为前端，说明**混合式（规则前端 + 学习后端）是现实路线**，适合我们在框架期先做规则改进、再逐步引入学习模块。

---

## 9. 落地路线图（对应到本仓库）

| 阶段 | 动作 | 触及文件 | 验收方式 |
|:--|:--|:--|:--|
| **0** | 换指标：能量加权 V-score + splitting rate；接受度修正（`E_tot > 阈值` 才算有效事例） | `pytrial/scoring.py`、`pytrial/sweeps/run_task.py` | 重跑 108 点扫描，指标可区分分裂/合并 |
| **1** | 加 spike/噪声抑制（CMS E4/E1、E6/E2）+ 显著度阈值原型（ATLAS 4σ/2σ/0σ，先假定统一 σ） | `pytrial/reconstruction.py` | μ 样本在高阈值下不再崩溃；阈值损失下降 |
| **2** | 局部极大值分裂 + 距离加权能量共享（Belle II / CMS PF 方案），输出改为 `ModIDList + weight` | `pytrial/reconstruction.py`、输出 schema、`pytrial/pipeline.md` | splitting rate 下降；`Eth=5` 的副簇能量损失显著降低 |
| **3** | 逐模块 σ_noise 标定（BESIII EINC 方法）+ 时间判据（ATLAS 3σ 触发、±2 ns 窗口、簇级加权时间） | 标定脚本 + `pytrial/reconstruction.py` | 显著度阈值可跨粒子种类通用 |
| **4** | 径迹-簇匹配 + 动量/能量一致性重聚类 + split-off 判别（Pandora / BESIII） | 新增 matching 模块 | 得到粒子级指派与 split-off 拒绝率 |
| **5**（可选） | 边分类器 / GNN（先 MLP 边分类，后 object condensation） | 新增 ML 模块 | 相对阶段 4 的 V-score 增益 |

---

## 10. 参考文献总表

| # | 文献 | 链接 | 关联方向 |
|:--|:--|:--|:--|
| 1 | ATLAS, *Topological cell clustering in the ATLAS calorimeters and its performance in LHC Run 1*, EPJC 77 (2017) 490 | [arXiv:1603.02934](https://arxiv.org/abs/1603.02934) | B |
| 2 | ATLAS, *Improving topological cluster reconstruction using calorimeter cell timing in ATLAS*, EPJC 84 (2024) 455 | [arXiv:2310.16497](https://arxiv.org/abs/2310.16497) | B, C |
| 3 | ATLAS, *Jet reconstruction and performance using particle flow with the ATLAS Detector*, EPJC 77 (2017) 466 | [arXiv:1703.10485](https://arxiv.org/abs/1703.10485) | B, D |
| 4 | ATLAS, *Electron and photon performance measurements … 2015–2017*, JINST 14 (2019) P12006 | [arXiv:1908.00005](https://arxiv.org/abs/1908.00005) | F |
| 5 | V. Prasad et al. (BESIII EMC), *Study of BESIII electromagnetic calorimeter performance with radiative lepton pair events* | [arXiv:1606.00248](https://arxiv.org/abs/1606.00248) | B, D |
| 6 | BESIII Collaboration, *Improved measurement of absolute branching fraction of the inclusive decay Λc⁺ → K_S⁰ X* | [arXiv:2502.20821](https://arxiv.org/abs/2502.20821) | D |
| 7 | Belle II, *Photon Reconstruction in the Belle II Calorimeter Using Graph Neural Networks* | [arXiv:2306.04179](https://arxiv.org/abs/2306.04179) | A |
| 8 | LHCb, *Graph Clustering: a graph-based clustering algorithm for the electromagnetic calorimeter in LHCb* | [arXiv:2212.11061](https://arxiv.org/abs/2212.11061) | A |
| 9 | *CLUE: A Fast Parallel Clustering Algorithm for High Granularity Calorimeters*, Front. Big Data 3 (2020) 591315 | [arXiv:2001.09761](https://arxiv.org/abs/2001.09761) | B |
| 10 | *The k4Clue package* | [arXiv:2311.03089](https://arxiv.org/abs/2311.03089) | B |
| 11 | CMS, *Particle-flow reconstruction and global event description with the CMS detector*, JINST 12 (2017) P10003 | [arXiv:1706.04965](https://arxiv.org/abs/1706.04965) | B, D |
| 12 | CMS, *Electron and photon reconstruction and identification … at the CERN LHC*, JINST 16 (2021) P05014 | [arXiv:2012.06888](https://arxiv.org/abs/2012.06888) | A |
| 13 | CMS, *Performance of photon reconstruction and identification … 8 TeV*, JINST 10 (2015) P08010 | [arXiv:1502.02702](https://arxiv.org/abs/1502.02702) | F |
| 14 | CMS, *Time Reconstruction and Performance of the CMS Electromagnetic Calorimeter*, JINST 5 (2010) T03011 | [arXiv:0911.4044](https://arxiv.org/abs/0911.4044) | C |
| 15 | CMS HGCAL, *Precision timing calorimetry with the CMS HGCAL*, JINST 15 (2020) C07003 | [arXiv:2005.13324](https://arxiv.org/abs/2005.13324) | C |
| 16 | CMS HGCAL, *Machine Learning for Event Reconstruction in the CMS Phase-2 HGC Endcap* | [arXiv:2510.01851](https://arxiv.org/abs/2510.01851) | F |
| 17 | *End-to-end multi-particle reconstruction … with graph neural networks*, EPJC 82 (2022) | [arXiv:2204.01681](https://arxiv.org/abs/2204.01681) | G |
| 18 | *Deep learning techniques for energy clustering in the CMS ECAL* | [arXiv:2204.10277](https://arxiv.org/abs/2204.10277) | F, G |
| 19 | *Machine Learning Power Week 2023: Clustering in Hadronic Calorimeters* | [arXiv:2508.09938](https://arxiv.org/abs/2508.09938) | E, G |
| 20 | *Reconstruction of overlapping electromagnetic showers in calorimeters using Transformers* (ClusTEX), EPJC 86, 873 (2026) | [arXiv:2603.18172](https://arxiv.org/abs/2603.18172) | E, G |
| 21 | J. S. Marshall, M. A. Thomson, *The Pandora Software Development Kit for Pattern Recognition*, EPJC 75 (2015) 439 | [arXiv:1506.05348](https://arxiv.org/abs/1506.05348) | D |
| 22 | M. Ruan, *Arbor, a new approach of the Particle Flow Algorithm* | [arXiv:1403.4784](https://arxiv.org/abs/1403.4784) | F |
| 23 | CMS 官方文档 *SWGuideEcalRecoClustering*（hybrid / island / RNINE） | [twiki](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideEcalRecoClustering) | F |
| 24 | CMS 官方文档 *SWGuideParticleFlowClustering*（多种子 + 高斯能量共享） | [twiki](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideParticleFlowClustering) | A |

### 调研局限（诚实声明）

- 抓取工具无法解析 PDF，因此部分会议报告/学位论文/期刊 PDF 无法逐字引用；引文改用 arXiv HTML / ar5iv 镜像。凡无法核实的片段已在正文中标注 **UNVERIFIED**。
- **ATLAS 拓扑聚类的分数共享解析公式**（arXiv:1603.02934 §3.1.3）未能取到（HTML 在 §3 前截断、Springer 需登录、CDS 有反爬）：分裂机制的存在与结构是已核实的，公式本身未核实。
- **CMS 粒子流 link 算法的数值判据**（arXiv:1706.04965 §0.4.1）未能取到，正文中改用 CMS 官方文档（其自身承认是未完成的 stub）。
- 未取得可引用一手来源的技术（如"sliding window"在晶体量能器语境下的原文）宁可不写，也不用二手转述冒充引文。
- **文献编号已逐个核验**：本文引用的全部 arXiv 编号都逐一访问过 `arxiv.org/abs/<ID>`（或 `export.arxiv.org` 镜像）核对标题、作者与期刊出处。核验中发现并已更正的 4 处书目错误：CMS 粒子流论文的期刊出处（JHEP → **JINST 12 (2017) P10003**）、Arbor 的作者（**仅 Manqi Ruan**，H. Videau 不在作者列表中）、BESIII EMC 性能论文的署名（**V. Prasad 等个人作者**，非 BESIII 合作组署名）、以及 Λc⁺ 论文标题漏掉的 "Improved"。
