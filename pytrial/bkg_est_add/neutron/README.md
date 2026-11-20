# 中子本底估计与本底注入（bkg_est_add/neutron）

本目录沿用 `../gamma/` 的同一套方法（不直接模拟大厅中子场，而是用「通量 + 俘获深度分布 + γ 衰减 + 泊松统计」的数值方法），在**本轮给定的最保守假设**下估计每个 100 ns 事例中的中子本底，并把能量沉积注入 `data/processed/` 的预处理数据。

**最保守假设**：所有入射中子都是**热中子**，且**俘获概率取 1**。由于每个中子最多只能俘获一次，这给出的是中子本底的**硬上界**，而不是最佳估计。

---

## 结果摘要

| 量 | 符号 | 数值 |
|:---|:---:|:---|
| 大厅中子通量（上限） | $\Phi$ | **200 n cm⁻² s⁻¹** |
| 外壳面积 / 有效面积 | $S_{\rm out}$, $S/4$ | 36 643.5 / 9 160.9 cm² |
| 几何覆盖率 | $\varepsilon$ | 0.8847（568/642）|
| 分子数密度 | $N_{\rm mol}$ | 1.0454×10²² cm⁻³ |
| 入射中子率 | $r_{\rm in}$ | 1.8322×10⁶ s⁻¹ |
| **俘获事例率** | $r_{\rm cap}$ | **1.6210×10⁶ s⁻¹** |
| 热俘获截面（每分子） | $\sigma_{\rm cap}$ | 35.2 b（Cs 29.0 + I 6.2）|
| 宏观俘获截面 | $\Sigma_{\rm cap}$ | 0.3680 cm⁻¹ |
| 俘获平均自由程 | $\lambda_c$ | **2.718 cm** |
| γ 线性衰减系数（6.88 MeV） | $\mu$ | 0.1691 cm⁻¹ |
| γ 平均自由程 | $\lambda_\gamma$ | **5.913 cm** |
| 每次俘获的 γ 总能量 | $Q_\gamma$ | 6.8798 MeV |
| **γ 沉积份额** | $f_{\rm dep}$ | **0.6494** |
| **单次俘获沉积能** | $E_\gamma$ | **4.4674 MeV** |
| **每事例俘获数** | $n$ | **0.16210** |
| **每事例沉积能** | $E_{\rm evt}$ | **0.7242 MeV** |
| 每事例逃逸能量 | — | 0.3910 MeV |
| 泊松概率 | P(0), P(1), P(2) | 0.8504, 0.1378, 0.0112 |

即：**约 15% 的事例含至少一次中子俘获**，每次俘获沉积 4.47 MeV，每事例平均沉积 0.724 MeV——占 32 MeV 信号的 2.26%、46 MeV 信号的 1.57%。

---

## 生成的数据

```
data/processed/edep_of_each_evt_ep_iso_bkg_n.root
data/processed/edep_of_each_evt_gamma_iso_bkg_n.root
data/processed/edep_of_each_evt_mu_rest_bkg_n.root
```

- 保留 `EdepOfEachEvt/Edeps` 原有 schema，可直接作为 `reconstruction.py` 的输入；
- 新增 `NNeutron` 分支，记录该事例注入的俘获数；
- 文件内另有 `TNamed("bkg_model_neutron", JSON)`，保存全部参数、随机种子与模型结果；
- 可选 `--store-bkg` 会额外写 `BkgEdeps`（本底分量），便于逐晶体分解。

**命名约定**（本轮统一）：

| 文件 | 含义 |
|:---|:---|
| `edep_of_each_evt_<dataset>.root` | 原始信号 |
| `edep_of_each_evt_<dataset>_bkg_g.root` | 注入 **γ** 本底（`../gamma/`，原 `_bkg` 已重命名）|
| `edep_of_each_evt_<dataset>_bkg_n.root` | 注入**中子**本底（本目录）|

两者都从原始信号文件独立生成，不级联；若要同时含两种本底，把两个文件的 `Edeps` 相加即可（随机流独立，互不影响）。

---

## 目录内容

| 文件 | 作用 |
|:---|:---|
| `neutron_model.py` | 第 1–6 步的模型：几何、俘获深度分布、γ 衰减、$f_{\rm dep}$ 积分、能量。可被其他脚本 import。 |
| `estimate_neutron.py` | 第 1–6 步的命令行入口：打印计算过程、写出 `neutron_estimation.json` 和深度分布图。 |
| `neutron_estimation.json` | **计算结果**（所有中间量、假设、解释）。 |
| `neutron_depth.png` | 俘获深度分布 $p(d)$ 与沉积份额 $g(d)$。 |
| `add_neutron_background.py` | **第 7 步**：泊松抽样俘获数 → 随机晶体 → 按俘获深度抽能量 → 累加到 `Edeps`。 |
| `verify_neutron_background.py` | 逐事例回读校验注入结果（可复现性、泊松统计、晶体编号均匀性、能量上限）。 |
| `neutron_verification.json` | 校验结果（三个数据集全部 PASS）。 |
| `PLAN.md` | 方案与推导（LaTeX 公式版），本 README 的完整推导版本。 |
| `REPORT.md` | **工作报告**：从中子能量估算、反应通道、最保守假设到事例率与能量沉积的完整分析与数据核验方法。 |

---

## 方法

### 第 1 步：几何

球壳内半径 260 mm、外半径 540 mm。各向同性通量下，中子进入凸体的速率是 $\Phi \cdot S/4$；对球体它恰好等于单向束的 $\Phi \cdot \pi R^2$，所以**角分布不影响一阶估计**。

$$S_{\rm out}=4\pi R_{\rm out}^2=36\,643.5\ \mathrm{cm^2},\qquad \frac{S_{\rm out}}{4}=9160.9\ \mathrm{cm^2}$$

理想 GP(8,0) 有 642 个晶胞，探测器上下游各一个缺口、实际 568 块：

$$\varepsilon=\frac{568}{642}=0.8847,\qquad N_{\rm mol}=\frac{\rho N_A}{M}=1.0454\times10^{22}\ \mathrm{cm^{-3}}$$

### 第 2 步：俘获事例率（保守取 $P_{\rm cap}=1$）

$$r_{\rm in}=\Phi\frac{S_{\rm out}}{4}=1.8322\times10^{6}\ \mathrm{s^{-1}}$$

$$r_{\rm cap}=r_{\rm in}\,\varepsilon\times1=1.6210\times10^{6}\ \mathrm{s^{-1}}$$

$$n=r_{\rm cap}\,\Delta t=1.6210\times10^{6}\times100\times10^{-9}=0.16210$$

这是全部通道的总上界：每个入射中子最多俘获一次。

### 第 3 步：俘获平均自由程与深度分布

天然 Cs 是 100% ¹³³Cs、天然 I 是 100% ¹²⁷I，分子截面直接相加（ENDF/B-VIII.0，0.0253 eV）：

$$\sigma_{\rm cap}=\sigma_\gamma(^{133}{\rm Cs})+\sigma_\gamma(^{127}{\rm I})=29.0+6.2=35.2\ \mathrm{b}$$

$$\Sigma_{\rm cap}=N_{\rm mol}\sigma_{\rm cap}=0.3680\ \mathrm{cm^{-1}},\qquad \lambda_c=\frac{1}{\Sigma_{\rm cap}}=2.718\ \mathrm{cm}$$

俘获深度服从指数分布：

$$p(d)=\frac{1}{\lambda_c}e^{-d/\lambda_c},\qquad 0\le d\le L=28\ \mathrm{cm}$$

**63.2% 的俘获发生在表面 2.72 cm 内，95.0% 在 8.15 cm 内。**

### 第 4 步：俘获 γ 的平均自由程

Cs 的俘获 γ 总能量 6.89 MeV、I 是 6.82 MeV，按丰度加权：

$$Q_\gamma=0.5115\times6.89+0.4885\times6.82=6.8798\ \mathrm{MeV}$$

$$\left.\frac{\mu}{\rho}\right|_{6.88\ \mathrm{MeV}}=0.0375\ \mathrm{cm^2/g}\;\Rightarrow\;\mu=4.51\times0.0375=0.1691\ \mathrm{cm^{-1}}$$

$$\lambda_\gamma=\frac{1}{\mu}=5.913\ \mathrm{cm}$$

$\lambda_\gamma$ 与晶体横向尺寸（3.6–7.5 cm）同量级，所以一个级联通常跨 1–2 块晶体。

### 第 5 步：沉积份额 $f_{\rm dep}$

俘获深度服从 $p(d)$；γ 各向同性，向内一半穿 $L-d$、向外一半穿 $d$ 后离开阵列：

$$g(d)=\frac{1}{2}\left[1-e^{-b(L-d)}\right]+\frac{1}{2}\left[1-e^{-bd}\right],\qquad b=\mu=0.1691\ \mathrm{cm^{-1}}$$

$$f_{\rm dep}=\int_0^L g(d)\,p(d)\,\mathrm{d}d=\frac{1}{2}\int_0^L\left[1-e^{-b(L-d)}\right]a e^{-ad}\,\mathrm{d}d+\frac{1}{2}\int_0^L\left[1-e^{-bd}\right]a e^{-ad}\,\mathrm{d}d$$

其中 $a=1/\lambda_c=0.3680\ \mathrm{cm^{-1}}$。两个初等积分：

$$\int_0^L e^{-b(L-d)}\,a e^{-ad}\,\mathrm{d}d=\frac{a}{a-b}e^{-bL}\left[1-e^{-(a-b)L}\right]$$

$$\int_0^L e^{-bd}\,a e^{-ad}\,\mathrm{d}d=\frac{a}{a+b}\left[1-e^{-(a+b)L}\right]$$

闭式：

$$f_{\rm dep}=\frac{1}{2}\left[1-\frac{a}{a-b}e^{-bL}\left(1-e^{-(a-b)L}\right)\right]+\frac{1}{2}\left[1-\frac{a}{a+b}\left(1-e^{-(a+b)L}\right)\right]$$

分项数值：

| 项 | 计算 | 结果 |
|:---|:---|:---|
| 向内半球 | $\frac{a}{a-b}=1.851$，$e^{-bL}=0.0088$，$1-e^{-(a-b)L}=0.9962$ | 0.4919 |
| 向外半球 | $\frac{a}{a+b}=0.6851$，$1-e^{-(a+b)L}=1-4\times10^{-7}$ | 0.1575 |
| **合计** | | **0.6494** |

数值积分（20 万点梯形）给出 0.6494，与闭式一致到 5×10⁻⁵。

> **读法**：向内那一半穿 ~25 cm（4 个自由程），几乎全沉积（0.98）；向外那一半只穿 ~2.7 cm（0.46 个自由程），只沉积 0.32。两半各占 50%，所以 $f_{\rm dep}=0.65$——**35% 的俘获 γ 从外表面逃掉了**。

### 第 6 步：能量

$$E_\gamma=Q_\gamma f_{\rm dep}=6.8798\times0.6494=4.4674\ \mathrm{MeV}\quad\text{（每次俘获）}$$

$$E_{\rm evt}=n E_\gamma=0.16210\times4.4674=0.7242\ \mathrm{MeV}\quad\text{（每事例）}$$

对照：每事例放出的总 γ 能量 $nQ_\gamma=1.1152$ MeV，其中 0.7242 MeV 沉积、0.3910 MeV 逃逸。

---

## 为什么"全部热中子 + 俘获概率 1"是最保守的

1. **把 $n$ 推到几何上界**：每个入射中子最多俘获一次，所以任何能谱下俘获率都不超过入射率。取 $f_{\rm th}=1$、$P_{\rm cap}=1$，等于假定每个进入阵列的中子都贡献一次俘获——这是 $n$ 的绝对上界。
2. **把单次沉积能量推到最大**：三个通道里只有俘获释放**无淬灭的电磁能量**（6.88 MeV）；弹性散射是高度淬灭的核反冲，非弹性只是 0.1–0.6 MeV 的分立 γ。
3. **沉积份额不额外放大**：$f_{\rm dep}=0.6494$ 是按真实深度分布积分的。若要取绝对上界，可令 $f_{\rm dep}=1$，则每事例 $nQ_\gamma=1.1152$ MeV；本模型采用 0.724 MeV，比绝对上界低 35%。

### 其它通道的量化比较

结论先给：**按每个入射中子算，俘获通道最多把 4.47 MeV 变成可见能；非弹性最多约 0.15 MeV，弹性淬灭后只有 keV 量级。** 下面逐条给依据，数值取"对其它通道最有利"的一端（1 MeV 快中子、2–3 次碰撞）。

#### 1. 弹性散射 $(n,n)$

两体运动学（质心系各向同性散射）：

$$E_R=\frac{2A}{(1+A)^2}E_n\left(1-\cos\theta_{\rm cm}\right),\qquad E_R^{\max}=\frac{4A}{(1+A)^2}E_n$$

对 ¹²⁷I（$A=127$）与 ¹³³Cs（$A=133$），$4A/(1+A)^2$ 分别是 0.031 与 0.030：**1 MeV 中子的最大反冲 31 keV，各向同性平均 15.5 keV。**

晶体内的碰撞次数：1 MeV 处 $\sigma_s\approx10$ b/分子（ENDF 量级），$\lambda_s=9.6$ cm，28 cm 内平均 1–3 次。每次平均转移 $0.0155E_n$，所以**每个 1 MeV 中子总沉积约 30 keV**。

淬灭：CsI 对 Cs/I 重反冲的闪烁效率远低于电子。CsI(Tl) 的核反冲淬灭因子在 10–100 keV 区间测得为 0.03–0.1（Saclay、UKDMC、KIMS、TEXONO 的测量，见文献 [5][6]）。纯 CsI 的光产额约为 CsI(Tl) 的 1/15（GDML：3500 vs 54 000 ph/MeV），其核反冲淬灭因子缺少系统测量；这里直接沿用 CsI(Tl) 的量级 0.03–0.1（取到 0.1 的上端），结论对该取值不敏感——即使 QF 取 0.1，也只有 3 keV/中子。

$$\text{可见能}=(30\ \mathrm{keV})\times(0.03\text{–}0.1)=0.9\text{–}3\ \mathrm{keV/中子}$$

#### 2. 非弹性散射 $(n,n'\gamma)$

阈能：¹³³Cs 第一激发态 81 keV、¹²⁷I 203 keV（ENSDF），所以 0.2 MeV 以上的中子就能开道。

截面：1 MeV 处 $\sigma_{\rm inel}\approx0.5$–2 b/分子（量级；精确值查 ENDF 的 MT=4），约为 $\sigma_s$（10 b/分子）的 5–20%。

每个中子：1–3 次碰撞 × 5–20% ≈ **0.05–0.6 次非弹性**。

γ 能量：低激发态退激给出 0.08–0.6 MeV 的分立线（Cs 81/161 keV、I 203/362 keV 起）；$\lambda_\gamma(0.2\text{–}0.6\ \mathrm{MeV})\approx2.3$–4.4 cm，远小于晶体尺寸，沉积份额约 0.8。

$$E_{\rm inel}=(0.05\text{–}0.6)\times(0.08\text{–}0.6\ \mathrm{MeV})\times0.8=3\text{–}300\ \mathrm{keV/中子}$$

#### 3. 其它通道

$(n,p)$、$(n,\alpha)$ 的阈能在 MeV 量级、截面 mb 量级（ENDF），比弹性低 2–3 个量级；产物是重带电粒子，淬灭与弹性同量级，可忽略。

#### 汇总

| 通道 | 单次可见能 | 每中子可见能 | 每事例（$n=0.162$）| 相对俘获 |
|:---|---:|---:|---:|---:|
| **辐射俘获** | 4.47 MeV | **4.47 MeV** | **724 keV** | 1 |
| 非弹性 $(n,n'\gamma)$ | ≤0.6 MeV | 3–300 keV | 0.5–50 keV | 1/15–1/1500 |
| 弹性 $(n,n)$ | ≤31 keV（淬灭前）| 0.9–3 keV | 0.15–0.5 keV | 1/1500–1/5000 |

**1–3 个数量级的差距来自三件事**：

1. **能量尺度**：俘获把整个中子分离能 6.88 MeV 一次放出；弹性每次只转移 $0.0155E_n$（1 MeV 时 15 keV），非弹性只有 0.1–0.6 MeV 的分立线。
2. **淬灭**：弹性通道的 30 keV 沉积再乘 0.03–0.1，落到 keV 量级。
3. **次数**：一个中子在 28 cm 里只散射 1–3 次，两个小通道都来不及累积。

**与阈值的关系**：只有俘获通道的单次能量（平均 4.47 MeV、上限 6.88 MeV）能逼近 $E_{\rm th,seed}=7$ MeV；非弹性单次 ≤0.6 MeV、弹性单次 ≤31 keV，离阈值差 1–2 个量级。

> **为什么这个比较不依赖谱**：4.47 MeV 与 ≤0.15 MeV 是**逐个中子**的比较——同一个中子如果被俘获，贡献 4.47 MeV；如果只发生弹性和非弹性，最多贡献 0.15 MeV。所以把热中子换成快中子只会让总可见能下降，全热中子（本模型的假设）确实是总可见能的**最大值**。若谱全是快的，其它通道加起来也只有 ~25 keV/事例，仍远低于 724 keV 的上界。

---

## 注入方法（第 7 步）

对每个事例 $i$：

1. 抽俘获数 $k\sim\mathrm{Poisson}(n)$（平均 0.162，85% 的事例 $k=0$）；
2. 在 568 块晶体中**均匀、有放回**地抽 $k$ 个编号（通量各向同性、阵列是球壳，每块晶体被从外击中的概率相同）；
3. 对每次俘获：从俘获深度分布抽 $d\sim\mathrm{Exp}(\lambda_c)$，按 $E=Q_\gamma\,g(d)$ 得到沉积能（均值 4.4674 MeV，范围 3.41–6.24 MeV）；
4. `Edeps[mod] += E`。

**单晶注入**：一次俘获的全部沉积能放在**一块**晶体里。这比物理情况（级联通常跨 1–2 块）更保守——它最大化了单晶能量，即便如此单晶也只在 3.4–6.1 MeV，仍低于 $E_{\rm th,seed}=7$ MeV。

随机数按 `numpy.random.default_rng([seed, crc32(dataset)])` 生成，种子默认 `20261116`（与 γ 的 `20261105` 独立，因此两个本底的随机流互不干扰），结果完全可复现。

```bash
cd pytrial/bkg_est_add/neutron
python3 add_neutron_background.py                          # 默认：三个数据集、全部 1e6 事例
python3 add_neutron_background.py --datasets mu_rest --max-events 10000   # 快速测试
python3 verify_neutron_background.py                       # 逐事例回读校验
```

---

## 注入结果

| 数据集 | 事例数 | 注入俘获数 | 注入总能量 (MeV) | 平均每次俘获 (MeV) |
|:---|---:|---:|---:|---:|
| `ep_iso` | 1,000,000 | 161,747 | 723 307 | 4.4718 |
| `gamma_iso` | 1,000,000 | 162,692 | 726 947 | 4.4682 |
| `mu_rest` | 1,000,000 | 161,513 | 721 735 | 4.4686 |

（平均每次俘获与模型 4.4674 MeV 一致；总能量取自各文件的 `bkg_model_neutron` 元数据。）

---

## 校验（`neutron_verification.json`）

`verify_neutron_background.py` 用输出文件里保存的种子**重新生成**全部随机数，再逐事例比较原始文件与注入文件：

- `NNeutron == 0` 的事例必须与原始文件逐位相同；
- `NNeutron > 0` 的事例的差值必须与重新生成的沉积能完全一致；
- `NNeutron` 分布与 `Poisson(n)` 的一致性（χ² 检验）；
- 晶体编号使用频率的均匀性（χ² 检验）；
- 平均每次俘获能量与 $Q_\gamma f_{\rm dep}$ 的一致性；
- **最大单次沉积 ≤ $Q_\gamma$**（硬上限）。

| 数据集 | 注入俘获 | 含俘获事例 | Poisson p | 均匀性 p | 平均能量 (MeV) | 最大单次 (MeV) | 结论 |
|:---|---:|---:|---:|---:|---:|---:|:---:|
| `ep_iso` | 161,747 | 149,389 (14.94%) | 0.719 | 0.444 | 4.4718 ± 0.0018 | 6.235 | PASS |
| `gamma_iso` | 162,692 | 150,016 (15.00%) | 0.213 | 0.322 | 4.4682 ± 0.0018 | 6.235 | PASS |
| `mu_rest` | 161,513 | 149,269 (14.93%) | 0.029 | 0.911 | 4.4686 ± 0.0018 | 6.235 | PASS |

三个数据集的逐事例最大偏差均为 **0.000 MeV**。

---

## 对聚类的影响（冒烟测试）

在 `mu_rest` 的前 2000 个事例上，用 `reconstruction.py`（`Eth_seed=7`、`Eth_wave=0.05`、`Nth=0`）分别对原始文件与 `_bkg_n` 文件各跑一遍：

| 文件 | 团簇总数 | 每事例团簇数分布 |
|:---|---:|:---|
| 原始 `mu_rest` | 1756 | {1: 1719, 2: 17, 3: 1} |
| `mu_rest_bkg_n` | 1756 | {1: 1719, 2: 17, 3: 1} |

**逐事例团簇数完全一致**（0 个事例改变）。这与"单次俘获沉积能 ≤6.88 MeV < $E_{\rm th,seed}=7$ MeV"的结论一致：中子只给事例加能量，不改变聚类结果。

## 假设与口径

1. **谱**：全部为热中子。这是上界；真实屏蔽大厅的谱里热中子只占一部分，结果只会更小。
2. **俘获概率**：1。热中子的真实值约 0.9（背散射逃逸），所以这一项高估约 10%。
3. **强子簇射**：不考虑（按本轮决定）。它是唯一能产生假团簇的通道，若将来测到硬谱需要回补。
4. **弹性/非弹性**：不建模。淬灭后的可见能只有 keV–几十 keV 量级，比俘获通道低 1–3 个量级。
5. **几何**：使用理想 GP(8,0) 的 642 个晶胞与实际 568 块晶体；$S/4$ 对球壳严格成立。
6. **事例窗**：100 ns（任务给定），连续束、稳态。
7. **单晶注入**：保守，最大化了单晶能量。
8. **能量口径**：`Edeps` 是**沉积能**（现有 MC 的 `nOptPho` 全为 0，没有可见能信息）。俘获 γ 是电磁过程，沉积能≈可见能，因此这一口径误差可忽略；核反冲才需要淬灭，而它们不在本通道里。

---

## 复现方式

```bash
cd pytrial/bkg_est_add/neutron

# 第 1–6 步：打印计算过程、写 neutron_estimation.json 和深度分布图
python3 estimate_neutron.py --plot neutron_depth.png

# 第 7 步：注入（约 30 s / 3 个文件）
python3 add_neutron_background.py

# 校验（约 70 s）
python3 verify_neutron_background.py --json neutron_verification.json
```

> 依赖：CERN ROOT 的 PyROOT、numpy、scipy、matplotlib。
> `data/processed/*.root` 受 `.gitattributes` 的 Git LFS 规则管理；`*_bkg_n.root` 为派生数据，可随时用上面的命令重新生成。

---

## 与 γ 本底的关系

| 项 | γ（`../gamma/`）| 中子（本目录）|
|:---|:---|:---|
| 源 | 30 μSv/h 剂量率 | 200 n cm⁻² s⁻¹ 通量 |
| 归一化 | $P_{\rm total}=D\cdot M$ | $r=\Phi S/4\,\varepsilon$ |
| 沉积谱 | N(5 MeV, 2 MeV²) 截断 | 俘获 γ 级联，$Q_\gamma f_{\rm dep}$ |
| 空间 | 均匀随机 | 均匀随机 + 表层俘获（单晶注入）|
| 时间 | 窗内瞬时 | 稳态（连续束）|
| 每事例 hit 数 | 2.4247 | **0.1621** |
| 每事例能量 | 12.13 MeV | **0.7242 MeV** |
| 假簇 | 主要威胁（0.003/事例 @ 10 MeV 阈值）| **0（单次上限 6.88 MeV < 阈值）** |
| 注入后缀 | `_bkg_g` | `_bkg_n` |
| 分支 | `NBkg` | `NNeutron` |
| 随机种子 | 20261105 | 20261116 |

**一句话**：γ 本底是"数得清但会造假簇"的本底；中子本底在纯 CsI 里只剩俘获通道，**是一个 0.724 MeV/事例（≤2.26%）的能量污染，且不产生假团簇**。

---

## 参考文献

1. 俘获截面、$Q_\gamma$ 与中子截面：ENDF/B-VIII.0（Cs、I），[NNDC](https://www.nndc.bnl.gov/endf/) 或 [IAEA Nuclear Data Services](https://nds.iaea.org/)。
2. 俘获 γ 线强（若日后细化 $f_{\rm dep}$）：IAEA PGAA 数据库，[nds.iaea.org/pgaa](https://nds.iaea.org/pgaa/)。
3. 质量衰减系数：NIST XCOM，[physics.nist.gov/xcom](https://physics.nist.gov/PhysRefData/Xcom/html/xcom1.html)。
4. ¹³³Cs、¹²⁷I 的能级（第一激发态 81 / 203 keV）：ENSDF，[nndc.bnl.gov/ensdf](https://www.nndc.bnl.gov/ensdf/)。
5. CsI(Tl) 核反冲淬灭因子（Saclay / UKDMC / KIMS / TEXONO 测量，0.03–0.1）：综述见 *New J. Phys.* **12** (2010) 075003，[doi:10.1088/1367-2630/12/7/075003](https://doi.org/10.1088/1367-2630/12/7/075003)。
6. *CsI(Tl) for WIMP dark matter searches*, arXiv:[hep-ex/0005031](https://arxiv.org/abs/hep-ex/0005031)（核反冲光产额与淬灭因子的测量）。
