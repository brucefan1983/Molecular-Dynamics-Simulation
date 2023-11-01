
# 《分子动力学模拟》第五章：输运性质的分子动力学模拟


# Table of contents
- [时间关联函数](#时间关联函数)
- [平衡态分子动力学模拟](#平衡态分子动力学模拟)
  - [扩散系数](#扩散系数)
  - [粘滞系数](#粘滞系数)
  - [热导率](#热导率)
- [齐性非平衡分子动力学模拟](#齐性非平衡分子动力学模拟)
- [非齐性非平衡分子动力学模拟](#非齐性非平衡分子动力学模拟)
- [习题](#习题)


根据热力学第二定律，一个系统在不受外场作用时，若其内部有热力学性质的不均匀性，则它一定处于非平衡的状态，并有向平衡态靠近的趋势。这种由热力学性质的不均匀性导致的热力学过程叫做输运过程，相应的现象叫做输运现象。例如，温度的不均匀性导致能量的输运（热传导现象）；粒子数密度的不均匀性导致粒子的输运（扩散现象）。将一个系统置于两个温度不同的热源之间，最终会在系统内建立一个稳定的~(不随时间变化的) 温度分布。我们说这样的系统处于一个稳态~(steady state)，但不处于一个平衡态~(equilibrium state)。稳态和平衡态都是不依赖与时间的，但前者属于非平衡态。
  
## 时间关联函数

在统计力学中，设有两个依赖于时间的物理量 $A(t)$ 和 $B(t)$，我们定义这两个量之间的时间关联函数 (time correlation function) $C(t)$ 为：

$$
C(t) = \left\langle A(t_0) B(t_0+t) \right\rangle.
$$

对这个公式的说明如下：
* 如果两个物理量相同， $A=B$，那么上式代表物理量 $A(t)$ 的'''自关联函数''' (auto-correlation function)。
* 上式中 $t$ 是某个时间间隔，叫做'''关联时间''' (correlation time)，而时间关联函数 $C(t)$ 是关联时间 $t$ 的函数。
* 尖括号在统计物理中代表系综平均，但在分子动力学模拟中一般代表“时间”平均，其中“时间”指的不是上述“关联时间” $t$，而是“时间原点” $t_0$。一般将 $t_0$ 写作 0，即常写

$$
C(t) = \left\langle A(0) B(t) \right\rangle。
$$

* 我们这里考虑的是平衡系统，即在每一个时间原点 $t_0$，系统都处于平衡态。所以，不同的时间原点在物理上是等价的，从而可以对时间原点求平均。在分子动力学模拟中，我们只能得到一条离散的相轨迹，故在实际计算中，与尖括号对应的时间平均将由求和方式表示。

首先，我们假设先在控制温度的情况下让系统达到平衡。然后在不控制温度的情况下（即在微正则系综）模拟了 $N_p$ 步（下标 $p$ 是 production 的意思），步长为 $\Delta t$。与计算静态热力学量的情形类似，我们不需要将每一步的数据都保存（因为相邻步的数据有关联性，保存得过频无益）。我们假设每 $N_s$ (下标 $s$ 是 sampling 的意思) 步保存一次数据，并称 $N_s$ 为取样间隔 (sampling interval)。我们假设 $N_p$ 是 $N_s$ 的整数倍，记

$$
\frac{N_p}{N_s} = N_d
$$

代表记录数据的总步数（下标 $d$ 是 data 的意思），并记

$$
N_s \Delta t = \Delta \tau.
$$

我们假设保存了如下 $N_d$ 个物理量 $A$ 和 $B$ 的值：

$$
A(\Delta \tau), A(2 \Delta \tau), A(3 \Delta \tau), \cdots, A(N_d \Delta \tau).
$$

$$
B(\Delta \tau), B(2 \Delta \tau), B(3 \Delta \tau), \cdots, B(N_d \Delta \tau).
$$

根据上面的讨论与记号，我们可以将关联时间为 $n \Delta \tau$ 的关联函数表达为：

$$
C(n \Delta \tau) = \frac{1}{N_d - n}\sum_{m=1}^{N_d - n} A(m \Delta \tau) B(m\Delta \tau + n \Delta \tau).
$$

其中， $N_d - n$ 是求时间平均时用的时间原点数目。

## 平衡态分子动力学模拟

这里要讲解更多理论基础（待写）

格林-久保公式实际上是一类公式，它们将非平衡过程的输运系数与平衡态中相应物理量的涨落相联系。格林-久保公式是说，输运系数等于自关联函数对关联时间的积分。例如，扩散系数是速度自关联函数的积分；粘性系数是压力自关联函数的积分；热导率是热流自关联函数的积分、等等。

### 扩散系数

#### 扩散系数的格林久保公式和爱因斯坦公式

速度自关联是单粒子关联函数。这就是说，我们可以为单个粒子定义速度自关联函数。对于粒子 $i$, 我们定义朝 $x$ 方向的速度自关联函数为：

$$
\langle v_{xi}(0) v_{xi}(t) \rangle.
$$

可以定义整体体系的 VAC 的平均值：


$$
\text{VAC} _{xx}(t) = \frac{1}{N} \sum _{i=1}^{N} \langle v _{xi}(0) v _{xi}(t) \rangle.
$$

上式中的时间平均 ($\langle \rangle$) 和空间平均 (对粒子的平均) 是可以交换的：

$$
\text{VAC} _{xx}(t) = \left\langle \frac{1}{N} \sum _{i=1}^{N}  v _{xi}(0) v _{xi}(t) \right\rangle.
$$

得到 VAC 之后，对关联时间积分就得到跑动扩散系数

$$
D _{xx}(t) = \int_0^{t} dt' ~\text{VAC} _{xx}(t')
$$

可以证明，格林-久保公式等价于爱因斯坦公式：

$$
D_{xx}(t) = \frac{1}{2} \frac{d}{dt} \Delta x^2(t),
$$

其中 $\Delta x^2(t)$ 是均方位移 （MSD），定义为

$$
\Delta x^2(t) =
\left\langle
\frac{1}{N} \sum _{i=1}^{N}  \left[ x _i(t) - x _i(0) \right]^2
\right\rangle =
\frac{1}{N} \sum _{i=1}^{N}
\left\langle \left[ x i(t) - x _i(0) \right]^2 \right\rangle.
$$

证明留作习题。

#### 例子：计算液态硅的扩散系数

#### 例子：计算水的扩散系数

#### 振动态密度

VAC 还可以用来计算振动态密度 （VDOS）[cite Dickey and Paskin]。VDOS 是归一化的 VAC 的傅里叶变换：

$$
\rho_x(\omega) = \int_{-\infty}^{\infty} dt e^{i\omega t}~\text{VAC}_{xx}(t).
$$

$\text{VAC} _{xx}(t)$ 应该被理解为 

$$
\text{VAC} _{xx}(t)/\text{VAC} _{xx}(0)
$$

利用自关联函数的对称性：

$$\text{VAC} _{xx}(-t) = \text{VAC} _{xx}(t)$$ 

我们有

$$
\rho _x(\omega) = \int _{-\infty}^{\infty} dt \cos (\omega t)~\text{VAC} _{xx}(t).
$$

因为我们只能在  $N_c$ 个离散的时间点计算 VAC，所以上述积分可以写成如下离散余弦变换：

$$
\rho_x(\omega) \approx \sum_{n_c=0}^{N_c-1} (2-\delta_{n_c0}) \Delta \tau \cos (\omega n_c \Delta \tau)~\text{VAC}_{xx}(n_c \Delta \tau).
$$

$\delta_{n_c0}$ 是克罗内克符号， 而因子 $(2-\delta_{n_c0})$ 表示时刻  $t = 0$ 只有一个关联数据，而其它时刻都有两个等价的数据。

为了消除吉布斯振荡，得到较为光滑的曲线，一般需要施加一个所谓的窗口函数。我们可以施加如下Hann窗口函数： 

$$
\rho_x(\omega) \approx \sum_{n_c=0}^{N_c-1}
(2-\delta_{n_c0}) \Delta \tau
\cos (\omega n_c \Delta \tau)~\text{VAC}_{xx}(n_c \Delta \tau) H(n_c);
$$

$$
H(n_c) = \frac{1}{2}
\left[ \cos \left( \frac{\pi n_c}{N_c} \right) + 1 \right].
$$

根据逆傅里叶变换：

$$
\text{VAC} _{xx}(t) = \int _{-\infty}^{\infty} \frac{d\omega}{2\pi} e^{-i\omega t}\rho_x(\omega).
$$

我们有

$$
1 = \text{VAC} _{xx}(0) = \int _{-\infty}^{\infty} \frac{d\omega}{2\pi}\rho_x(\omega).
$$

因为 $\rho_x(-\omega)=\rho_x(\omega)$ 所以有

$$
\int_{0}^{\infty}
 \frac{d\omega}{2\pi}\rho_x(\omega) = \frac{1}{2}.
$$

待写：多组分体系的公式需要推广。要使用质量加权的VAC。

#### 例子：计算石墨烯的振动态密度并讨论热熔的量子修正

#### 例子：计算水的振动态密度并讨论质量加权的意义

### 粘滞系数
 
粘滞系数可以表达为压强自关联函数的积分。记 $\alpha\beta$ 方向的压强与体积的乘积为 $S_{\alpha\beta}$ ，可以定义如下格林-久保积分公式：

$$
\eta _{\alpha\beta}(t) = \frac{1}{k _{\rm B} TV} \int_0^t \langle S _{\alpha\beta}(0) S _{\alpha\beta}(t') \rangle dt'
$$

剪切粘滞系数为

$$
\eta _{\rm S}(t) =  \frac{\eta _{xy}(t) +  \eta _{yz}(t) + \eta _{zx}(t) }{3}
$$

longitudinal 粘滞系数为

$$
\eta _{\rm L}(t) =  \frac{\eta _{xx}(t) +  \eta _{yy}(t) + \eta _{zz}(t) }{3}
$$

Bulk 粘滞系数为

$$
\eta _{\rm B}(t) =  \eta _{\rm L}(t) - \frac{4}{3} \eta _{\rm S}(t)
$$


### 热导率

热传导现象的宏观规律由傅里叶定律描述。傅里叶定律是说热流密度（热通量） $\vec{J}/V$，即单位时间穿过单位面积的热量，正比于温度梯度 $\vec{\nabla} T$：

$$
J _{\mu} = - \sum _{\nu} \kappa _{\mu\nu} \frac{\partial T}{\partial r^{\nu}}.
$$


这里的 $\kappa$ 就反映了热量输运的难易程度： $\kappa$ 越大代表热量越容易被输运。这样的物理量被称为输运系数~(transport coefficient)。具体到热传导，输运系数 $\kappa$ 叫做热导率~(thermal conductivity)。注意等式右边有个负号，它表示热量的传导方向与温度梯度的方向相反，指向温度降低的方向~(一个物理量的梯度的方向指向它增加的方向)。在国际单位制中，温度梯度的单位为~K/m，热流密度的单位是 W/m $^2$ ，故热导率的单位是 $\text{W}\text{m}^{-1}\text{K}^{-1}$ 。

对热导率的计算有如下格林-久保公式：

$$
\kappa_{\mu\nu}(t) = \frac{1}{k_B T^2 V} \int_0^{t} dt' C_{\mu\nu} (t').
$$

其中， $\kappa_{\mu\nu}(t)$ ($\mu, \nu = x, y, z$) 是热导率张量， $t'$ 是关联时间， $k_B$ 是 Boltzmann 常数,  $T$ 是温度,  $V$ 是体积， $C_{\mu\nu}(t)$ 是热流自关联函数（heat current autocorrelation function，常简称为 HCACF）。上式计算的跑动热导率（running thermal conductivity）。

热流自关联函数的表达式如下：

$$
 C_{\mu\nu}(t) = \langle J_{\mu}(0) J_{\nu}(t) \rangle.
$$

其中，尖括号表示统计平均，在分子动力学模拟中指对时间原点的平均， $J_{\mu}$ ($\mu = x, y, z$)  是热流。热流的表达式在第三章有详细讨论。

如果研究的是三维的各向同性（isotropic）的系统，则热导率张量的非对角分量一定为零，并可将最终计算的热导率取为对角分量的平均值：

$$
\kappa = \frac{\kappa_{xx} + \kappa_{yy} + \kappa_{zz}}{3}.
$$

使用格林-久保方法时要注意边界条件的选取：
* 如果模拟的是三维块体（bulk）系统，则每个方向都要使用周期边界条件。
* 如果要研究准两维系统（如薄膜或者两维材料），则在垂直于薄膜的方向用自由边界条件，在平行于薄膜的方向用周期边界条件。而且，此时垂直方向热导率的计算结果的意义比较特殊（我们后面会讨论）。
* 如果要研究的是准一维系统（如纳米线或者纳米管），则在垂直于线或者管的方向都要用自由边界条件，在平行于线或者管的方向用周期边界条件。


## 齐性非平衡分子动力学模拟

Consider a system of $N$ particles described by the general Hamiltonian


$$
H(\{\vec{r}_i,\vec{p}_i\}) = \sum _{i} \frac{\vec{p}_i^2}{2m_i} + U(\{\vec{r}_i\}),
$$

with the equations of motion 

$$
d\vec{r}_i/dt = \vec{p}_i/m_i
$$


$$
d\vec{p}_i/dt = \vec{F}_i
$$

$\vec{r}_i$ 是坐标

$m_i$ 是质量

$\vec{p}_i$ 是动量

$\vec{F}_i$ 是力


引入驱动力，使得运动方程变为


$$
\frac{d\vec{r}_i}{dt} = \frac{\vec{p}_i}{m_i} + \mathbf{C}_i(\{\vec{r}_i,\vec{p}_i\}) \cdot \vec{F} _{\rm e};
$$

$$
\frac{d\vec{p}_i}{dt} = \vec{F}_i + \mathbf{D}_i(\{\vec{r}_i,\vec{p}_i\}) \cdot \vec{F} _{\rm e}.
$$

Here 

$\mathbf{C} _i(\{\vec{r} _i,\vec{p} _i\})$ 和 $\mathbf{D} _i(\{\vec{r} _i,\vec{p} _i\})$  都是二阶张量。

$\vec{F}_{\rm e}$ 是矢量。

考察哈密顿量的时间导数：

$$
\frac{dH(\{\vec{r}_i,\vec{p}_i\})}{dt} = \vec{J} _{\rm d} \cdot \vec{F} _{\rm e},
$$

%
where $\vec{J}_{\rm d}=\vec{J}_{\rm d}(\{\vec{r}_i,\vec{p}_i\})$ is called the dissipative flux vecor. In terms of the dissipative flux, the nonequilibrium ensemble average $\langle \rangle_{\rm ne}$ of a general vecor physical quantity $\vec{A}(\{\vec{r}_i,\vec{p}_i\})$ at time $t$ after switching on the external driving force can be written as \cite{evans1990book,tuckerman2010book} ($k_{\rm B}T$ is the thermal energy)
%
\begin{align}
\langle \vec{A}(t)\rangle_{\rm ne}=\langle \vec{A}(0)\rangle 
+\left( \int_0^tdt'\frac{\langle \vec{A}(t')\otimes \vec{J}_{\rm d}(0)\rangle}{k_{\rm B}T} \right)
\cdot \vec{F}_{\rm e}.
\label{equation:A(t)}
\end{align}
%
Here, $\langle \vec{A}(0)\rangle$ is the usual equilibrium ensemble average of $\vec{A}$ and $\langle \vec{A}(t')\otimes \vec{J}_{\rm d}(0)\rangle$ is the equilibrium time correlation function between $\vec{A}$ and $\vec{J}_{\rm d}$.

The central idea of the HNEMD method by Evans \cite{evans1982pla} is to set both $\vec{A}$ and $\vec{J}_{\rm d}$ in Eq. (\ref{equation:A(t)}) to the heat current operator $\vec{J}_{\rm q}$, giving (note that $\langle \vec{J}_{\rm q}(0)\rangle=0$)
%
\begin{equation}
\langle \vec{J}_{\rm q}(t)\rangle_{\rm ne}
=
\left(
\frac{1}{k_{\rm B}T}\int_0^tdt'\langle \vec{J}_{\rm q}(t')\otimes \vec{J}_{\rm q}(0)\rangle
\right)
\cdot \vec{F}_{\rm e}.
\label{equation:J(t)}
\end{equation}
%
where $\langle \vec{J}_{\rm q}(t')\otimes \vec{J}_{\rm q}(0)\rangle$ is the equilibrium heat current autocorrelation function. Setting $\vec{J}_{\rm d}$ to $\vec{J}_{\rm q}$ fixes the equations of motion, as we will discuss soon. According to the Green-Kubo relation \cite{green1954jcp,kubo1957jpsj,mcquarrie2000book}, the quantity in the parentheses is related to the (running) thermal conductivity tensor, 
%
\begin{equation}
   \kappa^{\mu\nu}(t) = \frac{1}{k_{\rm B}T^2V} \int_0^t dt'\langle J_{\rm q}^{\mu}(t') J_{\rm q}^{\nu}(0)\rangle,
\end{equation}
%
$V$ being the system volume. Therefore, Eq. (\ref{equation:J(t)}) can be interpreted as 
%
\begin{equation}
\frac{\langle J^{\mu}_{\rm q}(t)\rangle_{\rm ne}}{TV} = \sum_{\nu} \kappa^{\mu\nu}(t) F_{\rm e}^{\nu}.    
\end{equation}
%
Working with principal axes \cite{nye1957book}, the thermal conductivity tensor is diagonal and the thermal conductivity $\kappa$ in a given direction is given by 
%
\begin{equation}
\label{equation:kappa}
\kappa(t) = \frac{\langle J_{\rm q}(t)\rangle_{\rm ne}}{TV F_{\rm e}}.    
\end{equation}
%
The running thermal conductivity $\kappa(t)$ calculated using this equation will show large fluctuations and it is not easy to judge when $\kappa(t)$ has converged. One can circumvent this difficulty by redefining $\kappa(t)$ as the following cumulative average:
%
\begin{equation}
\label{equation:kappa_prime}
\kappa(t)= \frac{1}{t} \int_0^{t} ds \frac{\langle J_{\rm q}(s)\rangle_{\rm ne}}{TV F_{\rm e}}.
\end{equation}
%
A similar definition has been implicitly used in previous works \cite{mandadapu2009jcp,dongre2017msmse} on the HNEMD method. 

To complete the derivation of the generalized HNEMD method, we need to determine the equations of motion, which are the foundation of the MD simulations. They are closely related to the heat current $\vec{J}_{\rm q}$ when the dissipative flux $\vec{J}_{\rm d}$ defined in Eq. (\ref{equation:dHdt}) is chosen to be the same as $\vec{J}_{\rm q}$. We discuss the heat current and the equations of motion next.

The general heat current formulae in MD simulations have been discussed in Ref. \cite{fan2015prb} in great detail. For a general many-body potential with the total potential energy $U=\sum_i U_i(\{\vec{r}_{ij}\}_{j\neq i})$, the heat current can be written as \cite{fan2015prb}
%
\begin{equation}
\label{equation:J_many-body}
\vec{J}_{\rm q} = \vec{J}_{\rm q}^{\rm kin}  + \vec{J}_{\rm q}^{\rm pot} 
= \sum_i \frac{\vec{p}_i}{m_i} E_i
+ \sum_{i,j\neq i} \frac{\vec{p}_i}{m_i} \cdot \left(\frac{\partial U_j}{\partial \vec{r}_{ji}} \otimes \vec{r}_{ij} \right),
\end{equation}
%
where $E_i=\vec{p}_i^2/2m_i+U_i$ is the total energy of particle $i$ and $U_i$ is the potential energy. The position difference is defined as $\vec{r}_{ij} \equiv \vec{r}_j - \vec{r}_i$. The equations of motion are constructed to make the dissipative flux $\vec{J}_{\rm d}$ identical to the heat current $\vec{J}_{\rm q}$. Evans chose the term $\mathbf{C}_i(\{\vec{r}_i,\vec{p}_i\})=0$. Then, the time derivative of the Hamiltonian (\ref{equation:H}) can be derived from the equations of motion (\ref{equation:eom-r}) and (\ref{equation:eom-p}) to be
%
\begin{equation}
\frac{dH}{dt} =\sum_i \frac{\vec{p}_i}{m_i} \cdot \left(\mathbf{D}_i\cdot \vec{F}_{\rm e} \right).    
\end{equation}
%
Comparing this with Eqs. (\ref{equation:dHdt}) and (\ref{equation:J_many-body}) and setting $\vec{J}_{\rm d}=\vec{J}_{\rm q}$, we have
%
\begin{equation}
\label{equation:driving_force}
\mathbf{D}_i \cdot \vec{F}_{\rm e} = E_i \vec{F}_{\rm e} +  \sum_{j \neq i} \left(\frac{\partial U_j}{\partial \vec{r}_{ji}} \otimes \vec{r}_{ij}\right) \cdot \vec{F}_{\rm e}.    
\end{equation}
%
This driving force will be added to the total force for particle $i$. Because the summation $\sum_i \mathbf{D}_i \cdot \vec{F}_{\rm e} \neq 0$, the total momentum of the system will not be conserved under this driving force. To restore momentum conservation, one needs to subtract the mean force of the total system from the force on each particle. Formally, this is equivalent to modifying the driving force to 
%
\begin{align}
\label{equation:DF-many-body}
\mathbf{D}_i \cdot \vec{F}_{\rm e}
&= E_i\vec{F}_{\rm e} -\frac{1}{N}\sum_jE_j \vec{F}_{\rm e} \nonumber \\
&+ \sum_{j \neq i} \left(\frac{\partial U_j}{\partial \vec{r}_{ji}} \otimes \vec{r}_{ij}\right) \cdot \vec{F}_{\rm e} \nonumber \\
&- \frac{1}{N}\sum_j \sum_{k \neq j} \left(\frac{\partial U_k}{\partial \vec{r}_{kj}} \otimes \vec{r}_{jk}\right) \cdot \vec{F}_{\rm e}.
\end{align}
%
One can easily verify that for two-body potentials, Eq. (\ref{equation:DF-many-body}) reduces to that by Evans \cite{evans1982pla}. However, we emphasize that the heat current formula for two-body potentials does not apply to many-body potentials \cite{fan2015prb}. One also needs to apply a thermostat to keep the temperature of the system at the target. To this end, we use the Nos\'{e}-Hoover chain thermostat \cite{tuckerman2010book} here. 

## 非齐性非平衡分子动力学模拟

讲热导率模拟的 NEMD 方法

## 习题

### 证明扩散系数的格林久保公式和爱因斯坦公式的等价性。

证明如下：

我们从坐标与速度的关系出发：

$$
x_i(t) - x_i(0) = \int_{0}^{t}dt' v_{xi}(t'),
$$

可以得到：

$$
[x_i(t) - x_i(0)]^2 =
\int_{0}^{t} dt' v_{xi}(t') \int_{0}^{t}dt''  v_{xi}(t'')=
\int_{0}^{t} dt' \int_{0}^{t}dt'' v_{xi}(t') v_{xi}(t'').
$$

那么，MSD可以表达为：

$$
\Delta x^2(t) =
\frac{1}{N}\sum_{i=1}^{N}
\int_{0}^{t} dt' \int_{0}^{t}dt''
\left\langle v_{xi}(t') v_{xi}(t'') \right\rangle.
$$

利用求导的莱布尼茨规则，可得：


$$
D_{xx}(t) = \frac{1}{2} \frac{d}{dt} \Delta x^2(t) =
\frac{1}{N}\sum_{i=1}^{N}
\int_{0}^{t} dt'
\left\langle v_{xi}(t) v_{xi}(t') \right\rangle,
$$

上式也可以写为：

$$
D_{xx}(t) =
\frac{1}{N}\sum_{i=1}^{N}
\int_{0}^{t} dt'
\left\langle v_{xi}(0) v_{xi}(t'-t) \right\rangle.
$$

令 $\tau=t'-t$ 可得

$$
D_{xx}(t) =
\frac{1}{N}\sum_{i=1}^{N}
\int_{-t}^{0} d\tau
\left\langle v_{xi}(0) v_{xi}(\tau) \right\rangle,
$$

进一步推到可得

$$
D_{xx}(t) =
\frac{1}{N}\sum_{i=1}^{N}
\int_{-t}^{0} d\tau
\left\langle v_{xi}(-\tau) v_{xi}(0) \right\rangle.
$$

再令 $t'=-\tau$ 可得

$$
D_{xx}(t) =
\frac{1}{N}\sum_{i=1}^{N}
\int_{0}^{t} dt'
\left\langle v_{xi}(t') v_{xi}(0) \right\rangle
=\int_0^t dt' ~\text{VAC}_{xx}(t').
$$

这就从基于MSD的爱因斯坦公式推导出了基于VAC的格林-久保公式。所以，它们是等价的。


