# 中子本底分析工作报告

## 中子源与能量估算

质子轰击石墨靶产生散裂中子，中子被混凝土/铁屏蔽与多次散射后到达实验大厅：

- 低能端：中子被含氢材料慢化，形成**热中子**（Maxwellian，最可几 0.0253 eV）与**超热中子**（谱形状近似 $1/E$）成分；
- 高能端：高能中子衰减长度更长，散裂中子谱相对**硬化**，但绝对通量大幅下降。

参考：PSI 590 MeV 回旋加速器大厅的测量（[*NIM A* **476** (2002) 58](https://pdf.sciencedirectassets.com/271580/1-s2.0-S0168900200X02471/1-s2.0-S0168900201013882/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEDIaCXVzLWVhc3QtMSJGMEQCIBhmhHZzwS49EQXIY72LU5RG%2Fc9%2FcvnUdNDr2ZqRpSBZAiBortcUjhpTm7LLFDawGa%2B4UEHrDdYg70abeFCM9FIbRCq8BQj6%2F%2F%2F%2F%2F%2F%2F%2F%2F%2F8BEAUaDDA1OTAwMzU0Njg2NSIMeXMJ02wHFY8pDVdFKpAFcVO84GnEpGtFAJoKZg56PG5FQ3KU3gfNCTIAQRPWGxuRLc1bGNpwEi89gbZviQzEnxG3pjBXfv8pLvqkSZ%2FtsXYG%2FaHoz%2B%2BFxPs0t57MCxb%2BLtonya0aImumo5fMpKHaKLJqK%2Bov96GwMYNT3RNID%2FONwV9E17QwuaF0TkLMA6J8o4GXD8bYDWJNXzCe71W3aH%2BtGkV0zC9dxnTEbYe0PmhYGo3i3MzljmIsu9%2FfUb6PhHr%2FEj7fL6Lj%2Fh8hrjy7TRbdkHSyppYI5UM6uuQqN97Q55e8zeULLXwjyQf4MtizxyQrmMZy9ksE%2BOH%2B8gByRCQxBTOdT8ZDZr0Y8r8f7chSAUczycBFQcoqHd3F%2B8zLb3tx7GnJi2YNmJ6YR77%2FnnTS1upvwSDRNdGqSOk8hyzw3Z%2BOgj7fGnoKH%2FkBFUlzt9ngwDQdz45vUShEoxBmf4Grxm%2FAlYz29tYbqJ0jHCXKqq%2F%2FhHuGEVf6SDn7wt08D8H8JHFuSOHToBrVSLVKTPQVVL%2BMgmDcDZT%2B765RLtJvT0ZEfuRg1G3GDYvEDszTsgHQ%2FbNvyxDUDE%2BjsiSFXsxWRhiW1W2QZrUhnfLWf7Qw3oDKlKxLUFzQ%2BjDMr2SHAvMhdaIhdnypggO4pZyQuWeIfc2nlAMxJmv%2FgPOd82l7u5gxFNbSG4%2FdkQKkH7sPDngUUXZKIsj7dUwKtX4OEL7URvq7MYn2t6lZpSmTsojGZiFBrblPG406MUAv7Eeem2HySPq%2FyROnjaUh7WIeL9qV9Iz2e%2FxBnvuDzoNYNPyOp%2By8zIH%2FrpH8aynvRqZI7boPsm8T3K5Eiox5Y6HjRiAND6dM6b2eDhX4cmOJCUx2GC9RPZftzpuC%2Fy%2Ff5W8wkNz11wY6sgHUEbmNRK7pKYgy3zps5N%2BzTR7XDx2pUCf1cewFjToqSiOfMoCXxPfiMD0tPC8WIzuEV3dNF%2FL5JgI%2Fxuk7EWkuFDTqyWGJerH0vbuv0f8Ay1jBf19eT34R%2FnzQlHtAMHq60M75iQcoCP2d%2B1VJoFg1my4pnyD%2FeIO1Rl9W3hu8NWN3KSLCD9oo4UTGyRrI1jGUs2HPTBvo8f7WUjaacMVXp147FCN0sR0rHevO4RBOPJOd&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20261118T094216Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTY3MRFQC2V%2F20261118%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=86d45f73a44075704596b5c5766700530e6014d7dfac5903292f6768bb32f3c1&hash=d2a0cdc1fc299c5a4924103d5cfe4050ceb28ff0e2cd6b52c61ba9ab5ffe00e9&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0168900201013882&tid=spdf-f4a0b257-5fb6-4774-8328-30da2b3ca805&sid=8fae189c51323743833ba8820109c07630aegxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&rh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=1908045504505e57535357&rr=a5c606321c5e3478&cc=cn&kca=eyJrZXkiOiJzblR4T3lYM0IxelNaa3Bodjd4MHh6S1JEaHZtVUZ3WXN4c3ppOUVaSTRyb3BBNUYrdkpoZi9ZSnpnMU9Ock5jL1AvNC9sS3RPTDRLbEhFQTdpeU5iMmVTUWZndmRJdnZlenpMOGhLdDg1clBKZFQ5WVhHaWhzM3RSL28wbWdwSjlmU2RXNXhBaFZjS2Z4QisyRkhrc1lFeHlkY1JuaUFIZDlRUEczR1RTcEVMcnkxMXZRPT0iLCJpdiI6IjQ0Y2E5NGZjOGFlMzIyMjhhNGY3YWY3MDM4YjQ0NTNjIn0=_1794994951857)）

**后续只使用总通量上限 200 n/cm²/s，不需要谱形。** 最保守假设是"全部中子都是热中子且都被俘获"。

| 能段 | 能量范围 | 典型来源 | 在纯 CsI 中的主要反应 |
|:---|:---|:---|:---|
| 热中子 | $\lesssim0.1\ \mathrm{eV}$ | 屏蔽体慢化后扩散回来 | **辐射俘获**（$\sigma$ 最大）|
| 超热中子 | $0.1\ \mathrm{eV}$–$100\ \mathrm{keV}$ | 慢化中的 $1/E$ 尾巴 | 俘获（谱形按照$1/v$中子速度下降 + 共振）、弹性散射 |
| 快中子 | $100\ \mathrm{keV}$–$600\ \mathrm{MeV}$ | 靶源直接泄漏 / 屏蔽穿透 | 弹性、非弹性、强子簇射 |

---

## 中子与 CsI 的反应通道

### 1. 弹性散射 $(n,n)$

中子与 Cs/I 核发生弹性碰撞，把一部分动能传给反冲核。两体运动学给出反冲能

$$E_R=\frac{2A}{(1+A)^2}E_n\left(1-\cos\theta_{\rm cm}\right),\qquad E_R^{\max}=\frac{4A}{(1+A)^2}E_n$$

对 $^{127}\mathrm{I}$（$A=127$）与 $^{133}\mathrm{Cs}$（$A=133$），$4A/(1+A)^2$ 分别是 0.031 与 0.030：**1 MeV 中子的最大反冲只有 31 keV**，各向同性散射的平均值是这个的一半（15.5 keV）。

反冲核会产生闪烁光，但**重带电粒子的光产额被强烈淬灭**：CsI 对 Cs/I 反冲的淬灭因子在反冲能 10–100 keV 区间测得为 0.06–0.18（arXiv:hep-ex/0005031）。因此 1 MeV 中子沉积的 ~30 keV 反冲能对应 **1.8–5.4 keV 可见能**。

### 2. 非弹性散射 $(n,n'\gamma)$

当中子能量超过靶核第一激发态时，可以把核激发到激发态；退激时放出 γ：

$$n+{}^{A}\mathrm{X}\;\to\;n+{}^{A}\mathrm{X}^{*}\;\to\;n+{}^{A}\mathrm{X}+\gamma$$

阈能由第一激发态退激发能量决定：$^{133}\mathrm{Cs}$ 是 80 keV ([ENSDF-133Cs](https://www.nndc.bnl.gov/ensnds/133/Cs/133cs_n_nPg.pdf))、$^{127}\mathrm{I}$ 是 58 keV（[ENSDF-127I](https://www.nndc.bnl.gov/ensnds/127/I/127i_n_nPg.pdf)）, 退激 γ 在 CsI 中的自由程 2–4 cm，沉积份额约 0.8。反应截面与中子速率成近似反比，1 MeV 处截面在 1 barn 的量级（[EDNF](https://www.nndc.bnl.gov/sigma/getInterpreted.jsp?evalid=15176&mf=3&mt=4)），约为弹性截面的 2–10%。

带电粒子通道 $(n,p)$、$(n,\alpha)$

阈能在 MeV 量级、截面在 mb 量级以下（ENDF），比弹性低 2–3 个量级；产物是重带电粒子，淬灭与弹性散射同量级。对 MeV 量级的中子可忽略。

### 3. 辐射俘获 $(n,\gamma)$

中子被核吸收、形成复合核，复合核退激放出一串 γ：

$$n+{}^{133}\mathrm{Cs}\to{}^{134}\mathrm{Cs}^{*}\to{}^{134}\mathrm{Cs}+\gamma$$

（退激放出 γ 级联。）

$$n+{}^{127}\mathrm{I}\to{}^{128}\mathrm{I}^{*}\to{}^{128}\mathrm{I}+\gamma$$

（退激放出 γ 级联。）

目前分析的**唯一能在纯 CsI 中释放 MeV 级、无淬灭电磁能量的通道**：一次俘获放出整个中子分离能（约 6.88 MeV，后有计算）。它的截面在热中子处最大（35.2 b/分子），并随中子速率按 $1/v$ 下降。

### 4. 强子簇射（暂时不考虑）

能量足够高的中子（$\gtrsim20\ \mathrm{MeV}$）可以在核内产生 π 介子、打碎核子，形成强子级联，在晶体里沉积 10–300 MeV。**唯一可能产生假团簇**的中子通道（单次能量可超过 7 MeV 阈值）。


---

## MeV量级的快中子不在 CsI 中热化而俘获

俘获截面只在热能附近才大（35.2 b），而MeV量级中子在CsI内发生散射消耗的动能很少，无法停止被俘获

**平均对数能降**：中子每次弹性碰撞平均损失的能量（对数）为

$$\xi=\left\langle\ln\frac{E_{\rm in}}{E_{\rm out}}\right\rangle=1+\frac{(A-1)^2}{2A}\ln\frac{A-1}{A+1}$$

对重核近似为 $\xi\approx2/(A+2/3)$。CsI 按丰度平均：

$$\xi(\mathrm{Cs})=0.0150,\qquad \xi(\mathrm{I})=0.0157,\qquad \xi(\mathrm{CsI})=0.5115\times0.0150+0.4885\times0.0157=0.0153$$

**热化所需碰撞数**：从 $E_0=1MeV$ 慢化到 $E_{\rm th}=0.0253eV$ 需要

$$N=\frac{\ln(E_0/E_{\rm th})}{\xi}=\frac{\ln(1\ \mathrm{MeV}/0.0253\ \mathrm{eV})}{0.0153}=\frac{17.5}{0.0153}\approx1140$$

约 1140 次碰撞。

**热化所需路径**：1 MeV 处 $\sigma_s\approx10$ b/分子（两个原子各在 5b 以内），

$$\lambda_s=\frac{1}{N_{\rm mol}\sigma_s}=\frac{1}{1.0454\times10^{22}\times10\times10^{-24}}=9.6\ \mathrm{cm}$$

$$L_{\rm th}=N\lambda_s\approx1140\times9.6\ \mathrm{cm}\approx110\ \mathrm{m}$$

晶体只有 28 cm 长，中子在里面平均只散射 $28/9.6\approx2.9$ 次，能量只损失约 4%。


---

## 最保守假设：全热中子 + 全俘获

只有热中子能在晶体内俘获。真实谱中热中子只占一部分，但本报告取**最保守的上界**：

- **假设1**：200 n/cm²/s 全部是热中子；
- **假设2**：每个进入晶体的热中子都被俘获（$P_{\rm cap}=1$，真实值约 0.9，因为会有背散射逃逸）。

---

## 俘获截面与深度分布

### 1 热俘获截面

天然 Cs 是 100% $^{133}\mathrm{Cs}$、天然 I 是 100% $^{127}\mathrm{I}$（两者都是单同位素元素），所以分子截面就是两个核截面直接相加。在热中子参考点 $E_{\rm th}=0.0253\ \mathrm{eV}$（2200 m/s）：

$$\sigma_{\rm cap}(\mathrm{CsI})=\sigma_\gamma(^{133}\mathrm{Cs})+\sigma_\gamma(^{127}\mathrm{I})=29.0\pm0.5\ \mathrm{b}+6.2\pm0.2\ \mathrm{b}=35.2\ \mathrm{b/分子}$$

随着能量增加俘获截面按中子速率反比下降（不考虑共振）。
**0.0253 eV**： $kT$ 在 293.6 K（20 °C）的取值。


### 2 俘获深度分布

$$\Sigma_{\rm cap}=N_{\rm mol}\,\sigma_{\rm cap}=1.0454\times10^{22}\times35.2\times10^{-24}=0.3680\ \mathrm{cm^{-1}}$$

$$\lambda_\mathrm{cap}=\frac{1}{\Sigma_{\rm cap}}=\frac{1}{0.3680}=2.718\ \mathrm{cm}$$

也就是热中子在 CsI 中平均走 2.72 cm 就被俘获


一束热中子进入晶体后，强度按 $e^{-x/\lambda_\mathrm{cap}}$ 衰减；俘获深度按照指数分布：

$$p(d)=\frac{1}{\lambda_\mathrm{cap}}e^{-d/\lambda_\mathrm{cap}},\qquad 0\le d\le L=28\ \mathrm{cm}$$

由它得到：

$$P(d<\lambda_\mathrm{cap})=1-e^{-1}=63.2\%,\qquad P(d<3\lambda_\mathrm{cap})=1-e^{-3}=95.0\%$$

即 **63% 的俘获发生在距表面 2.72 cm 以内、95% 在 8.15 cm 以内**——俘获集中在外层晶体靠近外表面（或内表面）的薄层里。


---

## 俘获 γ 的能量与衰减

### 1 一次俘获放出 γ 能量：$Q_\gamma$

[IAEA](nds.iaea.org/pgaa/data/gamdoc.pdf)数据查表：能量释放在约10keV-6.1MeV不等，其中主要是Cs原子有较高概率释放MeV量级的gamma

按照最高能gamma估计：若按照原子质量差估算所有gamma的能量总和：$$Q=\left[M(A)+M(n)-M(A+1)\right]c^2$$

用 AME2020 的原子质量（原子质量已包含电子，两边电子数相同，电子质量抵消）：

| 反应 | 质量 (u) | $Q$ (MeV) |
|:---|:---|:---|
| $^{133}\mathrm{Cs}+n\to{}^{134}\mathrm{Cs}$ | 132.9054519610 + 1.00866491582 − 133.906718501 | **6.8915** |
| $^{127}\mathrm{I}+n\to{}^{128}\mathrm{I}$ | 126.9044719 + 1.00866491582 − 127.9058099 | **6.8250** |

一次俘获落在 Cs 上还是 I 上，概率正比于两者的俘获截面（29.0 : 6.2），所以平均每次俘获的 γ 能量是

$$Q_\gamma=\frac{29.0\times6.8915+6.2\times6.8250}{35.2}=6.8798\ \mathrm{MeV}$$


### 2 gamma的质量衰减系数 $\mu/\rho$ 

γ 在介质中的衰减服从指数律：

$$I(x)=I_0e^{-\mu x}$$

$\mu$ 是**线性衰减系数**（$\mathrm{cm^{-1}}$），正比于介质密度。

$$\frac{\mu}{\rho}\qquad[\mathrm{cm^2/g}]$$

$$\mu=\rho\cdot\frac{\mu}{\rho}$$

$\mu/\rho$ 的微观含义：$\mu=n_{\rm atom}\sigma_{\rm tot}$，其中 $\sigma_{\rm tot}$ 是 γ 与一个原子的总相互作用截面。所以

$$\sigma_{\rm tot}=\frac{\mu}{n_{\rm atom}}$$

NIST 的 CsI 质量衰减系数表（<http://physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/cesium.html>）在 6.88 MeV 的gamma处给出 $\mu/\rho\approx0.0375\ \mathrm{cm^2/g}$

$$\mu=4.51\times0.0375=0.1691\ \mathrm{cm^{-1}},\qquad \lambda_\gamma=\frac{1}{\mu}=5.913\ \mathrm{cm}$$

交叉检验：$\mu/n_{\rm atom}=0.1691/2.091\times10^{22}=8.1\times10^{-24}\ \mathrm{cm^2}=8.1\ \mathrm{b/原子}$，即 16.2 b/分子。CsI 有 108 个电子，摊到每个电子 0.15 b——正是 7 MeV 处康普顿散射的单电子截面量级。

$\lambda_\gamma=5.91\ \mathrm{cm}$ 与晶体横向尺寸（3.6–7.5 cm）同量级，也许可以解释一个级联通常跨 1–2 块晶体。

---

## 沉积份额 $f_{\rm dep}$ 的积分

俘获发生在深度 $d$（分布 $p(d)$）。γ 各向同性发射：一半向内，穿过 $L-d$ 的 CsI；一半向外，只穿过 $d$ 就离开晶体阵列。单次俘获的沉积份额是

$$g(d)=\frac{1}{2}\left[1-e^{-b(L-d)}\right]+\frac{1}{2}\left[1-e^{-bd}\right],\qquad b\equiv\mu=0.1691\ \mathrm{cm^{-1}}$$

对深度分布求平均：

$$f_{\rm dep}=\int_0^L g(d)\,p(d)\,\mathrm{d}d=\frac{1}{2}\int_0^L\left[1-e^{-b(L-d)}\right]a e^{-ad}\,\mathrm{d}d+\frac{1}{2}\int_0^L\left[1-e^{-bd}\right]a e^{-ad}\,\mathrm{d}d$$

（第一项对应向内半球，第二项对应向外半球。）

其中 $a=1/\lambda_\mathrm{cap}=0.3680\ \mathrm{cm^{-1}}$。

$$\int_0^L e^{-b(L-d)}\,a e^{-ad}\,\mathrm{d}d=\frac{a}{a-b}e^{-bL}\left[1-e^{-(a-b)L}\right]$$

$$\int_0^L e^{-bd}\,a e^{-ad}\,\mathrm{d}d=\frac{a}{a+b}\left[1-e^{-(a+b)L}\right]$$

$$f_{\rm dep}=\frac{1}{2}\left[1-\frac{a}{a-b}e^{-bL}\left(1-e^{-(a-b)L}\right)\right]+\frac{1}{2}\left[1-\frac{a}{a+b}\left(1-e^{-(a+b)L}\right)\right]$$

代入 $a=0.3680$、$b=0.1691$、$a-b=0.1989$、$a+b=0.5371$、$L=28$：

| 项 | 计算 | 结果 |
|:---|:---|:---|
| 向内半球 | $\frac{a}{a-b}=1.851$，$e^{-bL}=0.0088$，$1-e^{-(a-b)L}=0.9962$ | 0.4919 |
| 向外半球 | $\frac{a}{a+b}=0.6851$，$1-e^{-(a+b)L}=1-4\times10^{-7}$ | 0.1575 |
| **合计** | | **0.6494** |


**读法**：向内那一半要穿过约 25 cm（4 个自由程），几乎全沉积（0.98）；向外那一半穿过约 2.7 cm（0.46 个自由程），沉积 0.32。两半各占 50%，所以 $f_{\rm dep}=0.65$

---

## 俘获事例率与能量沉积

### 1 入射率

球壳探测器对**各向同性**通量的有效面积是 $S/4$（对球体它恰好等于单向束的 $\pi R^2$）：

$$S_{\rm out}=4\pi R_{\rm out}^2=4\pi\times54^2=36\,643.5\ \mathrm{cm^2}$$

$$\frac{S_{\rm out}}{4}=9160.9\ \mathrm{cm^2}$$

$$r_{\rm in}=\Phi\frac{S_{\rm out}}{4}=200\times9160.9=1.8322\times10^6\ \mathrm{s^{-1}}$$

### 2 俘获率与每事例俘获数

几何覆盖率：理想 GP(8,0) 有 642 个晶胞，探测器上下游各一个缺口、实际 568 块：

$$\varepsilon=\frac{568}{642}=0.8847$$

$$r_{\rm cap}=r_{\rm in}\,\varepsilon\times1=1.8322\times10^6\times0.8847=1.6210\times10^6\ \mathrm{s^{-1}}$$

100 ns 事例窗（连续束，稳态，$n=r\Delta t$）：

$$n=r_{\rm cap}\,\Delta t=1.6210\times10^6\times100\times10^{-9}=0.16210$$

即**每事例平均 0.162 次俘获**，约 **15% 的事例含至少一次俘获**。泊松分布：$P(0)=0.8504$、$P(1)=0.1378$、$P(2)=0.0112$。

### 3 单次俘获与每事例的沉积能量

平均每次俘获的沉积能量（分布呈指数）
$$E_\gamma=Q_\gamma f_{\rm dep}=6.8798\times0.6494=4.4674\ \mathrm{MeV}$$

乘以俘获率，平均每事例沉积能量
$$E_{\rm evt}=nE_\gamma=0.16210\times4.4674=0.7242\ \mathrm{MeV}$$


单次俘获的 γ 总能量上限是 $Q_\gamma=6.8798\ \mathrm{MeV}$（若级联全部停在一块晶体里；实际不可能），低于目前理想的种子阈值范围 $E_{\rm th,seed}=7$–18 MeV。因此：

- **单次俘获永远不能成为种子**，中子不产生假团簇；
- 同一晶体在一个 100 ns 窗内被两次俘获的概率 $\approx2\times10^{-5}$。

---


## 本底注入

### 注入方法

对每个事例：

1. 抽俘获数 $k\sim\mathrm{Poisson}(n)$（平均 0.162，85% 的事例 $k=0$）；
2. 在 568 块晶体中**均匀、有放回**地抽 $k$ 个编号（通量各向同性、阵列是球壳，每块晶体被从外击中的概率相同）；
3. 对每次俘获：从俘获深度分布抽 $d\sim\mathrm{Exp}(\lambda_c)$，按 $E=Q_\gamma\,g(d)$ 得到沉积能（均值 4.4674 MeV，范围 3.41–6.24 MeV）；
4. `Edeps[mod] += E`

**单晶注入**：一次俘获的全部沉积能放在**一块**晶体里。这比物理情况（级联通常跨 1–2 块）更保守

### 注入结果

| 数据集 | 事例数 | 注入俘获数 | 注入总能量 (MeV) | 平均每次俘获 (MeV) |
|:---|---:|---:|---:|---:|
| `ep_iso` | 1,000,000 | 161,747 | 723 307 | 4.4718 |
| `gamma_iso` | 1,000,000 | 162,692 | 726 947 | 4.4682 |
| `mu_rest` | 1,000,000 | 161,513 | 721 735 | 4.4686 |

（平均每次俘获与模型 4.4674 MeV 一致；总能量取自各文件的 `bkg_model_neutron` 元数据。）

### 校验

`verify_neutron_background.py` 用输出文件里保存的种子**重新生成**全部随机数，再逐事例比较原始文件与注入文件：

- `NNeutron == 0` 的事例必须与原始文件逐位相同；
- `NNeutron > 0` 的事例的差值必须与重新生成的沉积能完全一致；
- `NNeutron` 分布与 `Poisson(n)` 的一致性（χ² 检验）；
- 晶体编号使用频率的均匀性（χ² 检验）；
- 平均每次俘获能量与 $Q_\gamma f_{\rm dep}$ 的一致性；
- **最大单次沉积 ≤ $Q_\gamma$**

| 数据集 | 注入俘获 | 含俘获事例 | Poisson p | 均匀性 p | 平均能量 (MeV) | 最大单次 (MeV) | 结论 |
|:---|---:|---:|---:|---:|---:|---:|:---:|
| `ep_iso` | 161,747 | 149,389 (14.94%) | 0.719 | 0.444 | 4.4718 ± 0.0018 | 6.235 | PASS |
| `gamma_iso` | 162,692 | 150,016 (15.00%) | 0.213 | 0.322 | 4.4682 ± 0.0018 | 6.235 | PASS |
| `mu_rest` | 161,513 | 149,269 (14.93%) | 0.029 | 0.911 | 4.4686 ± 0.0018 | 6.235 | PASS |


### 对聚类的影响（冒烟测试）

在 `mu_rest` 的前 2000 个事例上，用[`Eth_seed=7`、`Eth_wave=0.05`、`Nth=0`]分别对纯信号文件与加本底文件各跑一遍：

| 文件 | 团簇总数 | 每事例团簇数分布 |
|:---|---:|:---|
| 信号| 1756 | {1: 1719, 2: 17, 3: 1} |
| 信号+中子本底 | 1756 | {1: 1719, 2: 17, 3: 1} |

**重建团簇数完全一致**（0 个事例改变），单次俘获沉积能 ≤6.88 MeV < $E_{\rm th,seed}=7$ MeV，几乎无法独立起簇，只可能中子本底混入事例团簇污染能量。

## 假设与口径

1. 弹性/非弹性：淬灭后的可见能只有 keV–几十 keV 量级，不考虑，比俘获通道低 1–3 个量级。
2. 强子簇射：暂时没考虑。唯一能产生假团簇的通道，若将来测到硬谱需要回补。
3. 谱：全部为热中子，保守估计
4. 俘获概率：1，保守估计
5. 单晶注入：保守，最大化了单晶能量
