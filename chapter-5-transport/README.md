
# 《分子动力学模拟入门》第五章：输运性质的分子动力学模拟


# Table of contents
- [时间关联函数](#时间关联函数)
- [线性响应理论](#线性响应理论)
  - [扩散系数](#扩散系数)
  - [粘滞系数](#粘滞系数)
  - [热导率](#热导率)
  
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

* 我们这里考虑的是平衡系统，即在每一个时间原点 $t_0$，系统都处于平衡态。所以，不同的时间原点在物理上是等价的，从而可以对时间原点求平均。

在分子动力学模拟中，我们只能得到一条离散的相轨迹，故在实际计算中，与尖括号对应的时间平均将由求和方式表示。

首先，我们假设先在控制温度的情况下让系统达到平衡。然后在不控制温度的情况下（即在微正则系综）模拟了 $N_p$ 步（下标 $p$ 是 production 的意思），步长为 $\Delta t$。与计算静态热力学量的情形类似，我们不需要将每一步的数据都保存（因为相邻步的数据有关联性，保存得过频无益）。我们假设每 $N_s$ (下标 $s$ 是 sampling 的意思) 步保存一次数据，并称 $N_s$ 为取样间隔 (sampling interval)。我们假设 $N_p$ 是 $N_s$ 的整数倍，记

$$
\frac{N_p}{N_s} = N_d
$$

代表记录数据的总步数（下标 $d$ 是 data 的意思），并记

$$
N_s \Delta t = \Delta \tau.
$$

根据上面的讨论与记号，我们可以将关联时间为 $n \Delta \tau$ 的关联函数表达为：

$$
C(n \Delta \tau) = \frac{1}{M}\sum_{m=1}^{N_d - n} A(m \Delta \tau) B(m\Delta \tau+t).
$$

其中，$N_d - n$ 是求平均时用的时间原点数目。

## 线性响应理论

这里要讲解更多理论基础（代写）

格林-久保公式实际上是一类公式，它们将非平衡过程的输运系数与平衡态中相应物理量的涨落相联系。格林-久保公式是说，输运系数等于自关联函数对关联时间的积分。例如，扩散系数是速度自关联函数的积分；粘性系数是压力自关联函数的积分；热导率是热流自关联函数的积分、等等。

### 扩散系数

Velocity autocorrelation (VAC) is an important quantity in MD simulations. On the one hand, its integral with respect to the correlation time gives the running diffusion constant, which is equivalent to that obtained by a time derivative of the mean square displacement (MSD). On the other hand, its Fourier transform is the phonon density of states (PDOS).

The VAC is a single-particle correlation function. This means that we can define the VAC for individual particles. For particle $i$, the VAC along the $x$ direction is defined as
\begin{equation}
\langle v_{xi}(0) v_{xi}(t) \rangle.
\end{equation}
Then, one can define the mean VAC for any number of particles. In the current version of GPUMD, it is assumed that one wants to calculate the mean VAC in the whole simulated system:
\begin{equation}
\boxed{
\text{VAC}_{xx}(t) =
\frac{1}{N} \sum_{i=1}^{N} \langle v_{xi}(0) v_{xi}(t) \rangle
}.
\end{equation}
The order between the time-average (denoted by $\langle \rangle$) and the space-average (the average over the particles) can be changed:
\begin{equation}
\boxed{
\text{VAC}_{xx}(t) =
\left\langle \frac{1}{N} \sum_{i=1}^{N}  v_{xi}(0) v_{xi}(t) \right\rangle
}.
\end{equation}
Using the same conventions as in the case of HAC calculations, we have
the following explicit expression for the VAC:
\begin{equation}
\label{equation:VAC}
\text{VAC}_{xx}(n_c\Delta \tau) = \frac{1}{(N_d-N_c)N}
\sum_{m=0}^{N_d-N_c-1} \sum_{i=1}^{N}
v_{xi}(m\Delta \tau) v_{xi}((m+n_c)\Delta \tau),
\end{equation}
where $n_c = 0, 1, 2, \cdots, N_c-1$.
The algorithm for calculating the VAC is quite similar to that for calculating the HAC and it thus omitted.

After obtaining the VAC, we can calculate the running diffusion constant $D_{xx}(t)$ as
\begin{equation}
\boxed{
D_{xx}(t) = \int_0^{t} dt' ~\text{VAC}_{xx}(t')
}.
\end{equation}
One can prove that this is equivalent to the time-derivative of the MSD, i.e., the Einstein formula:
\begin{equation}
\boxed{
D_{xx}(t) = \frac{1}{2} \frac{d}{dt} \Delta x^2(t)
},
\end{equation}
where the MSD $\Delta x^2(t)$ is defined as
\begin{equation}
\boxed{
\Delta x^2(t) =
\left\langle
\frac{1}{N} \sum_{i=1}^{N}  \left[ x_i(t) - x_i(0) \right]^2
\right\rangle =
\frac{1}{N} \sum_{i=1}^{N}
\left\langle
  \left[ x_i(t) - x_i(0) \right]^2
\right\rangle
}.
\end{equation}

Here is the proof. Starting from the relation between position and velocity,
\begin{equation}
x_i(t) - x_i(0) = \int_{0}^{t}dt' v_{xi}(t'),
\end{equation}
we have
\begin{equation}
[x_i(t) - x_i(0)]^2 =
\int_{0}^{t} dt' v_{xi}(t') \int_{0}^{t}dt''  v_{xi}(t'')=
\int_{0}^{t} dt' \int_{0}^{t}dt'' v_{xi}(t') v_{xi}(t'').
\end{equation}
Then, the MSD can be expressed as
\begin{equation}
\Delta x^2(t) =
\frac{1}{N}\sum_{i=1}^{N}
\int_{0}^{t} dt' \int_{0}^{t}dt''
\left\langle v_{xi}(t') v_{xi}(t'') \right\rangle.
\end{equation}
Using Lebniz's rule, we have
\begin{equation}
D_{xx}(t) = \frac{1}{2} \frac{d}{dt} \Delta x^2(t) =
\frac{1}{N}\sum_{i=1}^{N}
\int_{0}^{t} dt'
\left\langle v_{xi}(t) v_{xi}(t') \right\rangle,
\end{equation}
which can be rewritten as
\begin{equation}
D_{xx}(t) =
\frac{1}{N}\sum_{i=1}^{N}
\int_{0}^{t} dt'
\left\langle v_{xi}(0) v_{xi}(t'-t) \right\rangle.
\end{equation}
Letting $\tau=t'-t$, we get (note that here $t$ is considered as a constant)
\begin{equation}
D_{xx}(t) =
\frac{1}{N}\sum_{i=1}^{N}
\int_{-t}^{0} d\tau
\left\langle v_{xi}(0) v_{xi}(\tau) \right\rangle,
\end{equation}
which can be rewritten as
\begin{equation}
D_{xx}(t) =
\frac{1}{N}\sum_{i=1}^{N}
\int_{-t}^{0} d\tau
\left\langle v_{xi}(-\tau) v_{xi}(0) \right\rangle.
\end{equation}
Letting $t'=-\tau$, we finally get
\begin{equation}
D_{xx}(t) =
\frac{1}{N}\sum_{i=1}^{N}
\int_{0}^{t} dt'
\left\langle v_{xi}(t') v_{xi}(0) \right\rangle
=\int_0^t dt' ~\text{VAC}_{xx}(t').
\end{equation}
We thus have derived the Green-Kubo formula from the Einstein formula.

In summary,
\begin{itemize}
\item The derivative of half of the MSD gives the running diffusion coefficient.
\item The integral of the VAC gives the running diffusion coefficient.
\item One can obtain the MSD by integrating the VAC twice (numerically).
\end{itemize}


It is interesting that the same VAC can be used to compute the PDOS, as first demonstrated by Dickey and Paskin. The PDOS is simply the Fourier transform of the normalized VAC:
\begin{equation}
\rho_x(\omega) = \int_{-\infty}^{\infty} dt e^{i\omega t}~\text{VAC}_{xx}(t).
\end{equation}
Here, $\text{VAC}_{xx}(t)$ should be understood as the normalized function $\text{VAC}_{xx}(t)/\text{VAC}_{xx}(0)$. Although it looks simple, it does not mean that one can get the correct PDOS by a naive fast Fourier transform (FFT) routine. Actually, this computation is very cheap and we do not need FFT at all. What we need is a discrete cosine transform. To see this, we first note that, by definition, $\text{VAC}_{xx}(-t) = \text{VAC}_{xx}(t)$. Using this, we have
\begin{equation}
\rho_x(\omega) = \int_{-\infty}^{\infty} dt \cos (\omega t)~\text{VAC}_{xx}(t).
\end{equation}
Because we only have the VAC data at the $N_c$ discrete time points, the above integral is approximated by the following discrete cosine transform:
\begin{equation}
\rho_x(\omega) \approx \sum_{n_c=0}^{N_c-1}
(2-\delta_{n_c0}) \Delta \tau
\cos (\omega n_c \Delta \tau)~\text{VAC}_{xx}(n_c \Delta \tau).
\end{equation}
Here, $\delta_{n_c0}$ is the Kronecker $\delta$ function and the factor $(2-\delta_{n_c0})$ accounts for the fact that there is only one point for $t = 0$ and there are two equivalent points for $t \neq 0$. Last, we note that a window function is needed to suppress the unwanted Gibbs oscillation in the calculated PDOS. In GPUMD, the Hann window $H(n_c)$ is applied:
\begin{equation}
\rho_x(\omega) \approx \sum_{n_c=0}^{N_c-1}
(2-\delta_{n_c0}) \Delta \tau
\cos (\omega n_c \Delta \tau)~\text{VAC}_{xx}(n_c \Delta \tau) H(n_c);
\end{equation}
\begin{equation}
H(n_c) = \frac{1}{2}
\left[ \cos \left( \frac{\pi n_c}{N_c} \right) + 1 \right].
\end{equation}

Here are some comments on the normalization of the PDOS. In the literature, one usually uses an arbitrary unit for the PDOS, but it actually has a dimension of [time], and an appropriate unit for it can be 1/THz or ps. The normalization of $\rho_x(\omega)$ can be determined by the inverse Fourier transform:
\begin{equation}
\text{VAC}_{xx}(t)
 = \int_{-\infty}^{\infty} \frac{d\omega}{2\pi} e^{-i\omega t}\rho_x(\omega).
\end{equation}
As we have normalized the VAC, we have
\begin{equation}
1 = \text{VAC}_{xx}(0)
 = \int_{-\infty}^{\infty}
 \frac{d\omega}{2\pi}\rho_x(\omega).
\end{equation}
Because $\rho_x(-\omega)=\rho_x(\omega)$, we have
\begin{equation}
\int_{0}^{\infty}
 \frac{d\omega}{2\pi}\rho_x(\omega) = \frac{1}{2}.
\end{equation}
The calculated PDOS should meet this normalization condition (approximately).

### 粘滞系数

### 热导率

根据热力学第二定律，一个系统在不受外场作用时，若其内部有热力学性质的不均匀性，则它一定处于非平衡的状态，并有向平衡态靠近的趋势。这种由热力学性质的不均匀性导致的热力学过程叫做输运过程~(transport process)，相应的现象叫做输运现象。例如，温度的不均匀性导致能量的输运（热传导现象）；粒子数密度的不均匀性导致粒子的输运（扩散现象）。将一个系统置于两个温度不同的热源之间，最终会在系统内建立一个稳定的~(不随时间变化的) 温度分布。我们说这样的系统处于一个稳态~(steady state)，但不处于一个平衡态~(equilibrium state)。稳态和平衡态都是不依赖与时间的，但前者属于非平衡态。

上述不均匀性都是由相应的不均匀的物理量的梯度来量化的。这里，我们只研究热输运，而且假设输运方向沿着一个特定方向（假设是~$x$ 方向）的情形。热传导现象的宏观规律由傅里叶定律描述。
傅里叶定律是说热流密度~(heat flux，或者~heat current) $J$，即单位时间穿过单位面积的热量，在数量上正比于温度梯度~$\frac{dT}{dx}$：
\begin{equation}
J = - \kappa \frac{dT}{dx}.
\end{equation}
这里的~$\kappa$ 就反映了热量输运的难易程度：$\kappa$ 越大代表热量越容易被输运。这样的物理量被称为输运系数~(
transport coefficient)。具体到热传导，输运系数~$\kappa$ 叫做热导率~(thermal conductivity)。注意等式右边有个负号，它表示热量的传导方向与温度梯度的方向相反，指向温度降低的方向~(一个物理量的梯度的方向指向它增加的方向)。在国际单位制中，温度梯度的单位为~K/m，热流密度的单位是~W/m$^2$，故热导率的单位是~$\text{W}\text{m}^{-1}\text{K}^{-1}$。

对热导率的计算有如下格林-久保公式：
\begin{equation}
\kappa_{\mu\nu}(t) = \frac{V}{k_B T^2} \int_0^{t} dt' C_{\mu\nu} (t').
\end{equation}
其中，$\kappa_{\mu\nu}(t)$ ($\mu, \nu = x, y, z$) 是热导率张量，$t'$ 是关联时间，$k_B$ 是~Boltzmann 常数, $T$ 是温度, $V$ 是体积，$C_{\mu\nu}(t)$ 是热流自关联函数（heat current autocorrelation function，常简称为~HCACF）。上式计算的跑动热导率（running thermal conductivity）。关于跑动热导率的技术细节，将在具体的范例中讨论。

热流自关联函数的表达式如下：
\begin{equation}
 C_{\mu\nu}(t) = \langle J_{\mu}(0) J_{\nu}(t) \rangle.
\end{equation}
其中，尖括号表示统计平均，在分子动力学模拟中指对时间原点的平均，$J_{\mu}$ ($\mu = x, y, z$) 是热流。下一节讨论关联函数；下下节讨论热流的具体表达式。

如果研究的是三维的各向同性（isotropic）的系统，则热导率张量的非对角分量一定为零，并可将最终计算的热导率取为对角分量的平均值：
\begin{equation}
\kappa = \frac{\kappa_{xx} + \kappa_{yy} + \kappa_{zz}}{3}.
\end{equation}
使用格林-久保方法时要注意边界条件的选取：
\begin{itemize}
\item 如果模拟的是三维块体（bulk）系统，则每个方向都要使用周期边界条件。
\item 如果要研究准两维系统（如薄膜或者两维材料），则在垂直于薄膜的方向用自由边界条件，在平行于薄膜的方向用周期边界条件。而且，此时垂直方向热导率的计算结果无意义。
\item 如果要研究的是准一维系统（如纳米线或者纳米管），则在垂直于线或者管的方向都要用自由边界条件，在平行于线或者管的方向用周期边界条件。而且，此时只有平行方向的热导率结果才有意义。
\end{itemize}

热流定义为能量密度矩的时间导数：
\begin{equation}
\vec{J} = \frac{1}{V} \frac{d}{d t} \sum_i \vec{r}_i E_i.
\end{equation}
其中，
\begin{equation}
E_i = \frac{1}{2} m_i \vec{v}_i^2 + U_i
\end{equation}
是第~$i$ 个粒子的总能量，$U_i$ 是势能部分，$\frac{1}{2} m_i \vec{v}_i^2$ 是动能部分。

首先，用求导的莱布尼茨法则可得
\begin{equation}
  \vec{J} = \sum_i \vec{v}_i E_i + \sum_i \vec{r}_i \frac{d}{d t} E_i
\end{equation}
通常将上式右边的两项分别称为动能项
\begin{equation}
\vec{J}_{\textmd{kin}} =  \sum_i \vec{v}_i E_i
\end{equation}
和势能项
\begin{equation}
\vec{J}_{\textmd{pot}} =  \sum_i \vec{r}_i \frac{d}{d t} E_i.
\end{equation}
合起来，有
\begin{equation}
\vec{J} =  \vec{J}_{\textmd{kin}} + \vec{J}_{\textmd{pot}}.
\end{equation}
动能项不需要再推导了。利用动能定理
\begin{equation}
\frac{d}{dt}\left(\frac{1}{2}m_i\vec{v}_i^2\right)=
\vec{F}_i \cdot \vec{v}_i,
\end{equation}
可以将势能项写为
\begin{equation}
  \vec{J}_{\textmd{pot}} = \sum_i \vec{r}_i (\vec{F}_i \cdot \vec{v}_i)
  + \sum_i \vec{r}_i \frac{d U_{i}}{d t}.
\end{equation}


在两体势（two-body potentials）情形中，系统总势能可表达为
\begin{equation}
  U =  \frac{1}{2} \sum_{i}  \sum_{j \neq i} U_{ij}(r_{ij}).
\end{equation}
可以推导如下力的表达式：
\begin{equation}
  \vec{F}_i = \sum_{j \neq i} \vec{F}_{ij},
\end{equation}
\begin{equation}
  \vec{F}_{ij}
  = \frac{\partial U_{ij}}{\partial \vec{r}_{ij}}
  = - \vec{F}_{ji}.
\end{equation}
其中，$\vec{F}_{ij}$ 是第~$i$ 个粒子受到的来自于第~$j$ 个粒子的力，
\begin{equation}
 \vec{r}_{ij} \equiv \vec{r}_{j} - \vec{r}_{i},
\end{equation}
是从第~$i$ 个粒子指向第~$j$ 个粒子的位置差矢量。

\vspace{0.5cm}
\hrule
\hrule
练习。定义
\begin{equation}
U_i = \frac{1}{2}\sum_{j\neq i} U_{ij},
\end{equation}
推导如下公式
\begin{equation}
  \vec{J}_{\textmd{pot}}
  = \frac{1}{2}\sum_i \sum_{j \neq i}
    \vec{r}_{i} [\vec{F}_{ij} \cdot (\vec{v}_i + \vec{v}_j)].
\end{equation}
在分子动力学模拟中，如果采用了周期边界条件，一个宏观物理量就不可能依赖于绝对坐标~$\vec{r}_i$, 而应该依赖于相对坐标。证明上面的热流表达式等价于
\begin{equation}
  \vec{J}_{\textmd{pot}}
  = - \frac{1}{4}\sum_i \sum_{j \neq i}
    \vec{r}_{ij} [\vec{F}_{ij} \cdot (\vec{v}_i + \vec{v}_j)].
\end{equation}
再证明上式等价于下面的不怎么对称的形式：
\begin{equation}
  \vec{J}_{\textmd{pot}}
  = - \frac{1}{2}\sum_i \sum_{j \neq i}
    \vec{r}_{ij} [\vec{F}_{ij} \cdot \vec{v}_i].
\end{equation}
\hrule
\hrule
\vspace{0.5cm}

在分子动力学模拟中，位力~(virial) 张量定义为
\begin{equation}
 \textbf{W} = \sum_i \textbf{W}_{i},
\end{equation}
其中，$\textbf{W}_i$ 是单粒子位力：
\begin{equation}
 \label{equation:per_atom_virial}
 \textbf{W}_i
 = -\frac{1}{2} \sum_{j \neq i}\vec{r}_{ij} \otimes \vec{F}_{ij}.
\end{equation}
于是有：
\begin{equation}
\label{equation:j_pot_pair_stress}
  \vec{J}_{\textmd{pot}} = \sum_{i} \textbf{W}_{i} \cdot \vec{v}_i.
\end{equation}
这个公式就是~LAMMPS 中用的热流公式，但它只适用于两体势，不适用于多体势（many-body potentials）。对多体势的讨论，请参看我的论文。



