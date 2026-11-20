# 中子本底估算 — 最保守简化模型（v7）

**状态**：已实现（见本目录 `README.md`；模型 `neutron_model.py`、估计 `estimate_neutron.py`、注入 `add_neutron_background.py`、校验 `verify_neutron_background.py`）
**基线**：`pytrial/bkg_est_add/`（γ 本底，不在范围内）

## 版本变更

| 项 | v6 | v7 |
|:---|:---|:---|
| 入射谱 | 热 + 超热 + 快三段，谱加权 | **全部为热中子**（$f_{\rm th}=1$）|
| 俘获概率 | 随机游走模型 $P_{\rm th}=0.91$ 等 | **设为 1**（一入射一俘获）|
| 保留 | — | 平均自由程 $\to$ 俘获深度分布；$\lambda_\gamma$ $\to$ $f_{\rm dep}$ 积分；原第 8 步以后全部保留 |
| 结果 | $n=0.012$–$0.147$ | **$n=0.162$**（绝对上界）；$E_\gamma=4.41\ \mathrm{MeV}$/俘获；$0.715\ \mathrm{MeV}$/事例 |

**这一版是单边保守上界**：所有假设都朝"高估本底"的方向取，因此结果可以直接当作"中子本底不超过此值"来引用。

---

## 0. 结论

1. **俘获事例率**（每个 100 ns 事例窗内的期望俘获数）：

$$n = \Phi\,\frac{S_{\rm out}}{4}\,\varepsilon_{\rm geo}\,\Delta t\times1 = 200\times9161\times0.884\times10^{-7} = 0.162$$

2. **单次俘获的 γ 沉积能量**：

$$E_\gamma = Q_\gamma\,f_{\rm dep} = 6.88\ \mathrm{MeV}\times0.6494 = 4.41\ \mathrm{MeV}$$

3. **每事例的沉积能量**：

$$E_{\rm evt} = n\,E_\gamma = 0.162\times4.41 = 0.715\ \mathrm{MeV}$$

4. **不产生假团簇**：单次俘获的 γ 总能量上限 $Q_\gamma=6.88\ \mathrm{MeV}<E_{\rm th,seed}=7\ \mathrm{MeV}$，单个晶体永远够不到种子阈值。

5. **净效应**：能量污染 $0.715\ \mathrm{MeV}$/事例，占 32 MeV 信号的 2.26%、46 MeV 信号的 1.57%；信号纯度 $\ge0.977$。

---

## 1. 输入与保守假设

| 量 | 值 | 说明 |
|:---|:---|:---|
| 中子通量 $\Phi$ | $<200\ \mathrm{n\,cm^{-2}\,s^{-1}}$ | 取上限 |
| 事例窗 $\Delta t$ | $100\ \mathrm{ns}$ | |
| 晶体 | 纯 CsI，568 块 | GDML：`SCINTILLATIONYIELD = 3500 ph/MeV`、$\tau=30\ \mathrm{ns}$ |
| 几何 | 内半径 260 mm、外半径 540 mm、厚 280 mm | GP(8,0) + GDML 交叉检验 |
| 强子通道 / 活化 | 不考虑 | 按决定 |

**四项保守假设**（每一项都只会高估本底）：

| # | 假设 | 真实情况 | 保守倍数 |
|:-:|:---|:---|:---|
| A1 | 全部中子为热中子（$f_{\rm th}=1$）| 屏蔽大厅里热中子只占一部分，$1/E$ 尾巴与快中子的俘获概率低得多 | 可达 5–10× |
| A2 | 俘获概率 $P_{\rm cap}=1$ | 热中子实际约 0.9（背散射逃逸）| ~1.1× |
| A3 | 全部俘获都发生在 CsI 内 | 部分俘获发生在结构材料/大厅材料里 | 小 |
| A4 | 忽略晶体阵列的几何缺口以外的遮挡 | 无 | — |

> 若进一步把几何缺口也忽略（$\varepsilon_{\rm geo}=1$），$n$ 会升到 0.183——本版仍用物理覆盖率 0.884。

**保守假设**

1. **把 $n$ 推到几何上界**。"每个入射中子最多俘获一次"意味着：任何能谱成分下，俘获事例率都不会超过 $N_{\rm enter}$。取 $f_{\rm th}=1$ 且 $P_{\rm cap}=1$，就等于假定**每一个进入晶体阵列的中子都贡献一次俘获**——这是 $n$ 的绝对上界，真实谱（无论热/超热/快）只会让一部分中子不被俘获，$n$ 只会更小。
2. **把单次沉积能量推到最大**。三个通道里，只有辐射俘获释放的是**无淬灭的电磁能量**（$Q_\gamma=6.88\ \mathrm{MeV}$）；弹性散射是高度淬灭的核反冲，非弹性只是 $0.1$–$0.6\ \mathrm{MeV}$ 的分立 γ。所以按俘获通道估算，单次和每事例的沉积能量都取到最大。
3. **沉积份额保留物理计算，不额外放大**。$f_{\rm dep}=0.6494$ 是从真实俘获深度分布积分出来的（35% 的 γ 从外表面逃逸）。若要取**绝对**上界，可再令 $f_{\rm dep}=1$，即每事例 $nQ_\gamma=0.162\times6.88=1.11\ \mathrm{MeV}$——这是连 γ 逃逸都不计的极端上限。本版采用 0.724 MeV，是"保守谱 + 真实沉积"的组合，比绝对上界还低 35%。

因此：$n=0.162$ hits/事例 可以当作**俘获事例率的上界**直接引用；$E_{\rm evt}=0.715\ \mathrm{MeV}$ 是在此基础上叠加真实 $f_{\rm dep}$ 的结果，若连 γ 逃逸也忽略则绝对上界为 $1.11\ \mathrm{MeV}$。

---

## 2. 通道：只剩辐射俘获

| 通道 | 单次能量 | 可见光 | 本模型中的处理 |
|:---|:---|:---|:---|
| 弹性散射 | 反冲 $\le0.031E_n$ | 淬灭 $\times0.02$–$0.1$ | 可见能 $<2\ \mathrm{keV}$，忽略 |
| 非弹性 $(n,n'\gamma)$ | $0.08$–$0.6\ \mathrm{MeV}$ | 正常 | 每入射中子 $\sim0.02\ \mathrm{MeV}$，忽略 |
| **辐射俘获 $(n,\gamma)$** | $Q_\gamma=6.88\ \mathrm{MeV}$ 级联 | 正常 | **全部计入，且俘获概率取 1** |

**其它通道为什么不显著**

即使取与 A1 相反的极端（把 200 n/cm²/s 全部当成快中子），另外两个通道的可见能也远小于俘获通道：

| 通道 | 每个入射中子的沉积能 | 淬灭后可见能 | 每事例可见能（$n=0.162$）| 单次作用的最大可见能 | 能否接近 7 MeV 阈值 |
|:---|:---|:---|:---|:---|:---|
| 辐射俘获 | 4.47 MeV（$Q_\gamma f_{\rm dep}$）| 4.47 MeV | **0.724 MeV** | 6.88 MeV | 接近但不超过 |
| 弹性散射 | 每次转移 $0.03E_n$，28 cm 内散射 1–3 次 → 1 MeV 中子沉积 30–90 keV | $\times0.02$–$0.1$ → 0.6–9 keV | 0.1–1.5 keV | 31 keV × 淬灭 ≈ 2 keV | 否 |
| 非弹性 $(n,n'\gamma)$ | 0.1–0.4 次/中子 × 0.1–0.6 MeV × 80% → 约 10–200 keV | 20–200 keV | 2–30 keV | 0.6 MeV | 否 |

两者合计 **$<35\ \mathrm{keV}$/事例**，比俘获通道的 710 keV 低 **1–2 个量级**。原因有三：

- **能量尺度**：俘获一次放出整个中子分离能 6.88 MeV；弹性散射每次只转移 $0.03E_n$，非弹性 γ 单条 $\le0.6\ \mathrm{MeV}$。
- **淬灭**：纯 CsI 对 Cs/I 重反冲的淬灭因子只有 0.02–0.1，把弹性通道再压低 1–2 个量级。
- **阈值**：只有俘获通道的单次能量（6.88 MeV）能逼近 $E_{\rm th,seed}=7\ \mathrm{MeV}$；另外两个通道离阈值差 1–2 个量级，连"接近阈值"都做不到。

所以本模型只保留俘获通道，不是忽略物理，而是**其它通道在纯 CsI 里的可见响应比俘获通道低 1–2 个量级**。

---

## 3. 计算步骤

### 第 1 步：几何与物质常量

$$S_{\rm out}=4\pi R_{\rm out}^2=36\,644\ \mathrm{cm^2},\qquad S_{\rm eff}=\frac{S_{\rm out}}{4}=9161\ \mathrm{cm^2}$$

（球形探测器的各向同性平均投影面积 $S/4$ 与单向束的 $\pi R^2$ 相同。）

$$\varepsilon_{\rm geo}=\frac{568}{642}=0.884,\qquad N_{\rm mol}=\frac{\rho N_A}{M_{\rm mol}}=\frac{4.51}{259.809}\times6.022\times10^{23}=1.045\times10^{22}\ \mathrm{cm^{-3}}$$

### 第 2 步：俘获事例率（保守取 $P_{\rm cap}=1$）

$$R_{\rm in}=\Phi\,S_{\rm eff}=200\times9161=1.83\times10^6\ \mathrm{s^{-1}}$$

$$R_{\rm cap}=R_{\rm in}\,\varepsilon_{\rm geo}\times1=1.62\times10^6\ \mathrm{s^{-1}}$$

$$n=R_{\rm cap}\,\Delta t=1.62\times10^6\times10^{-7}=0.162$$

即"每事例平均 0.162 次俘获"。这是全部通道的总上界：每个入射中子最多俘获一次。

### 第 3 步：俘获平均自由程与深度分布

俘获截面取热中子值（ENDF/B-VIII.0，$E=0.0253\ \mathrm{eV}$）：

$$\sigma_{\rm cap}({\rm CsI})=\sigma_\gamma(^{133}{\rm Cs})+\sigma_\gamma(^{127}{\rm I})=29.0\pm0.5\ \mathrm{b}+6.2\pm0.2\ \mathrm{b}=35.2\ \mathrm{b}$$

（天然 Cs、I 均为单同位素，无需丰度加权。）

宏观截面与平均自由程：

$$\Sigma_{\rm cap}=N_{\rm mol}\,\sigma_{\rm cap}=1.045\times10^{22}\times35.2\times10^{-24}=0.368\ \mathrm{cm^{-1}}$$

$$\lambda_c=\frac{1}{\Sigma_{\rm cap}}=\frac{1}{0.368}=2.72\ \mathrm{cm}$$

参考量：散射截面 $\sigma_s\approx10\ \mathrm{b}$ 给出

$$\lambda_s=\frac{1}{N_{\rm mol}\sigma_s}=\frac{1}{1.045\times10^{22}\times10\times10^{-24}}=9.6\ \mathrm{cm}$$

**俘获深度分布**：热中子在晶体表面入射后，沿路径按 $e^{-d/\lambda_c}$ 衰减，故俘获深度的概率密度为

$$p(d)=a\,e^{-ad},\qquad a\equiv\frac{1}{\lambda_c}=0.3676\ \mathrm{cm^{-1}},\qquad 0\le d\le L=28\ \mathrm{cm}$$

在 $[0,L]$ 上的归一化是 $1-e^{-aL}=1-3.4\times10^{-5}\approx1$，无需修正。

**物理图像**：$63\%$ 的俘获发生在距表面 2.72 cm 以内、$95\%$ 在 8.2 cm 以内——**俘获集中在外层晶体靠近外表面的一侧**。

### 第 4 步：俘获 γ 的平均自由程

$$Q_\gamma=\frac{29.0\times6.8915+6.2\times6.8250}{35.2}=6.8798\ \mathrm{MeV}$$

（Cs、I 各 100% 丰度，但一次俘获落在哪个核上按俘获截面 29.0:6.2 分配，所以按截面加权，而不是按质量分数。）

$$\left.\frac{\mu}{\rho}\right|_{6.88\ \mathrm{MeV}}\approx0.0375\ \mathrm{cm^2/g}\;\Rightarrow\;\mu=\rho\frac{\mu}{\rho}=4.51\times0.0375=0.159\ \mathrm{cm^{-1}}$$

$$\lambda_\gamma=\frac{1}{\mu}=5.91\ \mathrm{cm}$$

$\mu$ 是**线性衰减系数**（$\mathrm{cm^{-1}}$，单位路程上的相互作用概率），$\mu/\rho$ 是**质量衰减系数**（$\mathrm{cm^2/g}$，与密度无关）。此能区康普顿散射主导。

$\lambda_\gamma=5.91\ \mathrm{cm}$ 与晶体横向尺寸（3.6–7.5 cm）同量级——这是后面"一个级联跨 1–2 块晶体"和 $f_{\rm dep}<1$ 的物理根源。

### 第 5 步：沉积份额 $f_{\rm dep}$ 的积分

俘获深度服从 $p(d)$；γ 各向同性发射，向内一半穿 $L-d$、向外一半穿 $d$ 后离开阵列：

$$g(d)=\frac{1}{2}\left[1-e^{-b(L-d)}\right]+\frac{1}{2}\left[1-e^{-bd}\right],\qquad b\equiv\mu=\frac{1}{\lambda_\gamma}=0.1587\ \mathrm{cm^{-1}}$$

$$f_{\rm dep}=\int_0^L g(d)\,p(d)\,\mathrm{d}d=\frac{1}{2}\int_0^L\left[1-e^{-b(L-d)}\right]a e^{-ad}\,\mathrm{d}d+\frac{1}{2}\int_0^L\left[1-e^{-bd}\right]a e^{-ad}\,\mathrm{d}d$$

（第一项对应向内半球，第二项对应向外半球。）两个初等积分：

$$\int_0^L e^{-b(L-d)}\,a e^{-ad}\,\mathrm{d}d=\frac{a}{a-b}\,e^{-bL}\left[1-e^{-(a-b)L}\right]$$

$$\int_0^L e^{-bd}\,a e^{-ad}\,\mathrm{d}d=\frac{a}{a+b}\left[1-e^{-(a+b)L}\right]$$

因此闭式为

$$f_{\rm dep}=\frac{1}{2}\left[1-\frac{a}{a-b}\,e^{-bL}\left(1-e^{-(a-b)L}\right)\right]+\frac{1}{2}\left[1-\frac{a}{a+b}\left(1-e^{-(a+b)L}\right)\right]$$

代入 $a=0.3680$、$b=0.1691$、$a-b=0.1989$、$a+b=0.5371$、$L=28$：

向内半球：

$$\frac{a}{a-b}=1.851,\quad e^{-bL}=0.0088,\quad 1-e^{-(a-b)L}=0.9962\;\Rightarrow\;\frac{1}{2}\left[1-0.0161\right]=0.4919$$

向外半球：

$$\frac{a}{a+b}=0.6851,\quad 1-e^{-(a+b)L}=1-4\times10^{-7}\;\Rightarrow\;\frac{1}{2}\left[1-0.6851\right]=0.1575$$

$$f_{\rm dep}=0.4919+0.1575=0.6494$$

区间取 $0.55$–$0.70$。

> **读法**：向内那一半（穿 ~25 cm，4 个自由程）几乎全沉积（0.98）；向外那一半（穿 ~2.7 cm，0.46 个自由程）只沉积 0.32。两半各占 50%，所以 $f_{\rm dep}=0.65$——**35% 的俘获 γ 从外表面逃掉了**。

### 第 6 步：单次俘获与每事例的能量

单次俘获的 γ 沉积能量：

$$E_\gamma=Q_\gamma\,f_{\rm dep}=6.88\times0.6494=4.41\ \mathrm{MeV}$$

每事例（$n=0.162$ 次俘获）：

$$E_{\rm evt}=n\,E_\gamma=0.162\times4.41=0.715\ \mathrm{MeV}$$

对照：每事例放出的总 γ 能量 $nQ_\gamma=0.162\times6.88=1.11\ \mathrm{MeV}$，其中 0.724 MeV 沉积、0.40 MeV 逃逸。

> **注意两个"上限"不要混**：$E_\gamma=4.41\ \mathrm{MeV}$ 是"每次俘获沉积多少"；单次俘获**可能**沉积的硬上限仍是 $Q_\gamma=6.88\ \mathrm{MeV}$（若级联全部停在一块晶体里，实际不会发生）。

---

## 4. 空间分布与击中晶体数（原第 8 步）

俘获深度 $\langle d\rangle=\lambda_c=2.72\ \mathrm{cm}$，落在最外层晶体（晶体厚 28 cm）；$\lambda_\gamma=5.91\ \mathrm{cm}$ 与晶体横向尺寸 3.6–7.5 cm 同量级，故通常沉积在 1–2 块晶体。

| 击中晶体数 | 概率 |
|---:|---:|
| 1 | ~60% |
| 2 | ~35% |
| $\ge3$ | ~5% |

单晶最大沉积受 $Q_\gamma=6.88\ \mathrm{MeV}$ 硬上限约束，典型 2–4 MeV。

---

## 5. 阈值判据（原第 9 步）

$$7\ \mathrm{MeV}\le E_{\rm th,seed}\le18\ \mathrm{MeV},\qquad E_{\rm dep}^{\max}\le6.88\ \mathrm{MeV}<7\ \mathrm{MeV}$$

所以单个俘获在任何 $E_{\rm th,seed}\ge7\ \mathrm{MeV}$ 下都不可能成为种子。

同晶两次俘获：

$$p_2\approx\frac{1}{2}\left(\frac{n}{568}\right)^2\times568\approx\frac{1}{2}\times\left(\frac{0.162}{568}\right)^2\times568\approx2.3\times10^{-5}$$

**中子假团簇概率 $\approx2\times10^{-5}$/事例，即 0。**

---

## 6. 对重建能量的影响（原第 10 步）

| 情形 | 影响 |
|:---|:---|
| 俘获落在信号簇内 | 簇能量 $+0\dots4.4\ \mathrm{MeV}$ |
| 俘获落在簇外 | 1–2 块晶体的 2–4 MeV 低能沉积，低于 $E_{\rm th,seed}$，不参与聚类 |
| 事件总能量偏差 | $+E_{\rm evt}=0.715\ \mathrm{MeV}$ |
| 信号纯度（32 MeV）| $\ge1-0.71/32=0.977$ |
| 信号纯度（46 MeV）| $\ge1-0.71/46=0.984$ |

---

## 7. 不确定度与检验（原第 11 步）

| 来源 | 范围 | 对结果的影响 |
|:---|:---|:---|
| $\Phi$ | $<200$（上限）| 线性 |
| $f_{\rm dep}$ | 0.55–0.70 | $E_\gamma=3.8$–$4.8\ \mathrm{MeV}$ |
| $Q_\gamma$ | 6.82–6.89 | $\pm0.5\%$ |
| $\lambda_c$ | 2.72 cm（ENDF 热截面，$\pm2\%$）| $f_{\rm dep}$ 变化 <1% |
| 模型保守性 | A1–A4 | **单边**，只会高估 |

**三项检验**：

1. **上界检验**：$n\le\Phi S_{\rm eff}\Delta t=0.183$（若连几何缺口也不计），本版 0.162 低于它。
2. **能量守恒**：$E_{\rm evt}=nQ_\gamma f_{\rm dep}$，且 $f_{\rm dep}\le1$。
3. **自洽性**：$f_{\rm dep}=0.65$ 应落在"向内 0.98、向外 0.30"之间。

---

## 8. 结果汇总

| 量 | 值 | 备注 |
|:---|---:|:---|
| 俘获事例率 $n$ | **0.162** /事例 | 全部热中子、$P_{\rm cap}=1$ 的上界 |
| 单次俘获沉积能 $E_\gamma$ | **4.47 MeV** | $=Q_\gamma f_{\rm dep}$ |
| 每事例沉积能 $E_{\rm evt}$ | **0.724 MeV** | $=nE_\gamma$ |
| 占 32 MeV 信号 | 2.26% | |
| 占 46 MeV 信号 | 1.57% | |
| 信号纯度 | $\ge0.977$ | |
| 假团簇率 | $\approx2\times10^{-5}$ | 即 0 |
| 单晶硬上限 | 6.88 MeV | $<E_{\rm th,seed}$ |

**最终表述**：在中子通量上限 $200\ \mathrm{n\,cm^{-2}\,s^{-1}}$、且**保守地假设全部中子都被俘获**的前提下，中子本底对每个事例的能量贡献不超过 **0.724 MeV**（单次俘获 4.47 MeV），占信号的 **1.55–2.26%**，且**不产生假团簇**。

---

## 9. 实现伪代码

```python
PHI, S_EFF, EPS_GEO, DT = 200.0, 9161.0, 568/642, 100e-9
N_CAP   = PHI * S_EFF * EPS_GEO * DT          # 0.162 次俘获/事例（保守：P_cap = 1）
Q_GAMMA = 6.88                                 # MeV
F_DEP   = 0.6494                               # 第 5 步积分结果
E_GAMMA = Q_GAMMA * F_DEP                      # 4.47 MeV / 俘获
E_EVT   = N_CAP * E_GAMMA                      # 0.724 MeV / 事例
```

---

## 10. 还需要确认的

无。本版只依赖题目给定的 $200\ \mathrm{n\,cm^{-2}\,s^{-1}}$；其余全部取保守上界。若日后拿到实测谱与俘获概率，把 $f_{\rm th}$、$P_{\rm cap}$ 代回去即可收紧（结果只会更小）。

---

## 11. 风险

| 风险 | 影响 | 缓解 |
|:---|:---|:---|
| 强子簇射被排除但实际存在 | 假簇被低估 | 已在 v2 记录该通道；若将来测到硬谱再回补 |
| $E_{\rm th,seed}$ 降到 5 MeV 以下 | 6.88 MeV 端点可过阈 | README 注明适用区间 $E_{\rm th,seed}\ge7\ \mathrm{MeV}$ |
| 大厅内还有其它中子产生机制（如束流损失）| $\Phi$ 被低估 | 题目给定 200 为上限；若有新数据重算 |

---

## 12. 与 γ 本底的关系

`data/processed/*_bkg.root`（γ 注入，$n=2.42$）不在本轮范围内。中子估算不使用 $30\ \mu\mathrm{Sv/h}$、不使用 γ 的沉积谱、也不与其比较。若将来要合并，只需把 `NBkg` 与 `NNeutron` 相加（独立抽样，互不影响）。

---

## 13. 参考文献

- 俘获截面与 Q 值：ENDF/B-VIII.0（Cs、I），[NNDC](https://www.nndc.bnl.gov/endf/) 或 [IAEA ND](https://nds.iaea.org/)。
- 质量衰减系数：NIST XCOM，[physics.nist.gov/xcom](https://physics.nist.gov/PhysRefData/Xcom/html/xcom1.html)。
- 俘获 γ 线强（若日后细化 $f_{\rm dep}$）：IAEA PGAA，[nds.iaea.org/pgaa](https://nds.iaea.org/pgaa/)。
