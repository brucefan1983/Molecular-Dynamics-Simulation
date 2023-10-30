
# 《分子动力学模拟》第六章：分子动力学模拟的高级课题


本章介绍分子动力学模拟中的几个高级课题，包括蒙特卡洛-分子动力学混合模拟（一种以分子动力学模拟为主线，以蒙特卡洛模拟为辅助的模拟）、路径积分分子动力学、自由能计算、紧束缚分子动力学以及相关电子结构和输运的计算等。

# Table of contents
- [蒙特卡洛与分子动力学混合模拟](#蒙特卡洛与分子动力学混合模拟)
  - [正则系综](#正则系综)
  - [半巨正则系综](#半巨正则系综)
  - [方差约束的半巨正则系综](#方差约束的半巨正则系综)
- [路径积分分子动力学](#路径积分分子动力学)
  - [基本概念](#基本概念)
  - [路径积分分子动力学的算法](#路径积分分子动力学的算法)
  - [基本物理量的计算](#基本物理量的计算)
  - [路径积分分子动力学的应用](#路径积分分子动力学的应用)
- [自由能计算](#自由能计算)
  - [自由能微扰理论](#自由能微扰理论)
  - [热力学积分方法](#热力学积分方法)
  - [基于非平衡模拟的热力学积分方法](#基于非平衡模拟的热力学积分方法)
- [紧束缚分子动力学](#紧束缚分子动力学)
  - [紧束缚模型](#紧束缚模型)
  - [基于紧束缚模型的电子输运性质](#基于紧束缚模型的电子输运性质)

## 蒙特卡洛与分子动力学混合模拟

### 正则系综

Canonical ensemble

### 半巨正则系综

Semi-grand canonical ensemble

### 方差约束的半巨正则系综

Variance-constrained semi-grand canonical ensemble


## 路径积分分子动力学

### 基本概念

量子力学基础

待写。快速过渡到 Feynman 路径积分量子力学。

量子-经典对应

讨论 Chandler 和 Wolynes 1981 年 的工作。

 路径积分分子动力学的概念

讨论 Parrinello 和 Rahman 1984 年的工作。

RPMD 的概念

[Craig 和 Manolopoulos](https://doi.org/10.1063/1.1777575) 于2004 提出了 ring-polymer MD （RPMD).

CMD 的概念

[Jianshu Cao 和 Gregory A. Voth](https://doi.org/10.1063/1.467175) 于 1994 提出 centroid MD (CMD).

（但本书可能不打算针对 CMD 编程）

### 路径积分分子动力学的算法

首先讲 Ceriotti 等人针对PIMD的 PILE (path integral Langevin equation)。

[Efficient stochastic thermostatting of path integral molecular dynamics](https://doi.org/10.1063/1.3489925)

然后讲RPMD 以及 TRPMD (thermostatted RPMD) 的实现。

RPMD: [Quantum statistics and classical mechanics: Real time correlation functions from ring polymer molecular dynamics](https://doi.org/10.1063/1.1777575)

Lecture notes:
https://www.tugraz.at/fileadmin/user_upload/Institute/PTC/WTC/WTC_2019/2_Manolopoulos_19.L2.pdf

TRPMD: [How to remove the spurious resonances from ring polymer molecular dynamics](https://doi.org/10.1063/1.4883861)

这里要重点介绍 Korol 等人的 Cayley 变换。[Cayley modification for strongly stable path-integral and ring-polymer molecular dynamics](https://doi.org/10.1063/1.5120282)

用简谐振子为例编程实现。MATLAB 即可。

### 基本物理量的计算

势能的计算
动能的计算

位力的计算
热流的计算

### 路径积分分子动力学的应用

晶体的热熔

简写振子的能量。

![energy_ho](src/pimd/energy.png)

晶体的热膨胀

水的结构性质

用公开的 [NEP 势](https://gitlab.com/brucefan1983/nep-data) 计算水的径向分布函数，结果如下图所示。

![water_rdf](src/pimd/water_nep/water_rdf.jpg)

## 自由能计算

### 自由能微扰理论

### 热力学积分方法

### 基于非平衡模拟的热力学积分方法


## 紧束缚分子动力学

### 紧束缚模型

我们以sp3轨道为例介绍紧束缚哈密顿量的构造。

每个原子有四个轨道

$$
    \{ |s\rangle, |p_x\rangle, |p_y\rangle, |p_z\rangle \}.
$$

如果有 $N$ 个原子，那么哈密顿矩阵的维度为  $M = 4 N$ 。


For each atom $i$, we have a diagonal matrix $H_{i}$,

$$
    H^{ii} =
    \left(
        \begin{array}{cccc}
            E_s & 0   & 0   & 0   \\
            0   & E_p & 0   & 0   \\
            0   & 0   & E_p & 0   \\
            0   & 0   & 0   & E_p \\
        \end{array}
    \right).
$$

For each neighbor pair of atoms $i$ and $j$, there are two hopping matrix
$H_{ij}$ and $H_{ji}$ which are transpose of each other. The matrix $H_{ij}$
can be written as

$$
    H^{ij} =
    \left(
        \begin{array}{cccc}
            H^{ij}_{ss} & H^{ij}_{sx} & H^{ij}_{sy} & H^{ij}_{sz} \\
            H^{ij}_{xs} & H^{ij}_{xx} & H^{ij}_{xy} & H^{ij}_{xz} \\
            H^{ij}_{ys} & H^{ij}_{yx} & H^{ij}_{yy} & H^{ij}_{yz} \\
            H^{ij}_{zs} & H^{ij}_{zx} & H^{ij}_{zy} & H^{ij}_{zz} \\
        \end{array}
    \right).
$$

To calculate the matrix elements of $H_{ij}$, we need first define some angle
variables (not to be treated as usual triangular functions),

$$
    \cos(x) &=& x_{ij} / r_{ij},                          \\
    \cos(y) &=& y_{ij} / r_{ij},                          \\
    \cos(z) &=& z_{ij} / r_{ij},                          \\
    \sin^2(x) &=& 1 - \cos^2(x) = \cos^2(y) + \cos^2(z),  \\
    \sin^2(y) &=& 1 - \cos^2(y) = \cos^2(z) + \cos^2(x),  \\
    \sin^2(z) &=& 1 - \cos^2(z) = \cos^2(x) + \cos^2(y),
$$

where $r_{ij}$ is the distance of the two atoms, $x_{ij}$, $y_{ij}$, and
$z_{ij}$ are components of the position difference
$\textbf{r}_{ij} = \textbf{r}_{j} - \textbf{r}_{i}$ of the two atoms.

Using these angle variables and the hopping parameters $V_{ss\sigma}$,
$V_{sp\sigma}$, $V_{pp\sigma}$, and $V_{pp\pi}$, the matrix elements of
$H_{ij}$ can be written as


$$
    H^{ij}_{ss} &=& V_{ss\sigma},                                  \\
    H^{ij}_{xx} &=& V_{pp\sigma} \cos^2(x) + V_{pp\pi} \sin^2(x),  \\
    H^{ij}_{yy} &=& V_{pp\sigma} \cos^2(y) + V_{pp\pi} \sin^2(y),  \\
    H^{ij}_{zz} &=& V_{pp\sigma} \cos^2(z) + V_{pp\pi} \sin^2(z),  \\
    H^{ij}_{sx} &=& V_{sp\sigma} \cos(x),                          \\
    H^{ij}_{sy} &=& V_{sp\sigma} \cos(y),                          \\
    H^{ij}_{sz} &=& V_{sp\sigma} \cos(z),                          \\
    H^{ij}_{xy} &=& (V_{pp\sigma} - V_{pp\pi}) \cos(x) \cos(y),    \\
    H^{ij}_{yz} &=& (V_{pp\sigma} - V_{pp\pi}) \cos(y) \cos(z),    \\
    H^{ij}_{zx} &=& (V_{pp\sigma} - V_{pp\pi}) \cos(z) \cos(x),    \\
    H^{ij}_{xs} &=& - H^{ij}_{sx},                                 \\
    H^{ij}_{ys} &=& - H^{ij}_{sy},                                 \\
    H^{ij}_{zs} &=& - H^{ij}_{sz},                                 \\
    H^{ij}_{yx} &=& H^{ij}_{xy},                                   \\
    H^{ij}_{zy} &=& H^{ij}_{yz},                                   \\
    H^{ij}_{xz} &=& H^{ij}_{zx}.
$$




Band energy and the Feynman-Hellmann force using the direct
         Diagonalization method


Direct diagonalization can be easily performed by using linear algebra packages
such as Lapack. However, the complexity of the method is $O(N_{atom}^3)$, which
is very costly for large systems. It will takes several minutes for a single
MD step for a system of 1000 atoms.


The band energy $E_b$ of the system in a state can be written as

$$
    E_b &=& 2 \sum_{n = 1}^{N_{level}} f_n \langle n|H|n \rangle \\ \nonumber
        &=& 2 \sum_{n = 1}^{N_{level}} f_n \epsilon_n,
$$

where

$$
    f_n
    = \frac{1}{\exp\left(\frac{\epsilon_n - \mu}{k_B T}\right) + 1}
$$

is the Fermi-Dirac distribution function. $T$ is temperature and $k_B$ is
the Boltzmann's constant. The factor of 2 in the above equation accounts
for spin degeneracy.

To determine the Fermi-Dirac function, we should calculate the chemical
potential $\mu$, which is related to the number of valence electrons
$N_{electron}$ (which is $4 N_{atom}$ for Carbon) as follows,

$$
    N_{electron} = 2 \sum_{n = 1}^{N_{level}} f_n
$$


The Feynman-Hellmann force is the negative of the gradient of the band energy.
For atom $i$, we have,

$$
\textbf{f}_i = - \frac{\partial}{\partial \textbf{r}_i} E_b = - 2 \frac{\partial}{\partial \textbf{r}_i} \sum_n f_n \langle n|H|n \rangle            
= - 2 \sum_{j \alpha} \sum_{k \beta} \sum_n f_n C^n_{j \alpha} C^n_{k \beta} \frac{\partial}{\partial \textbf{r}_i} H^{jk}_{\alpha \beta}
$$


By defining the density matrix
\begin{equation}
    \rho^{jk}_{\alpha \beta} = \sum_n f_n C^n_{j \alpha} C^n_{k \beta},
\end{equation}
we can express the force on atom $i$ as
\begin{equation}
    \textbf{f}_i
    = - 2 \sum_{j \alpha} \sum_{k \beta} \rho^{jk}_{\alpha \beta}
    \frac{\partial}{\partial \textbf{r}_i} H^{jk}_{\alpha \beta}
    = -2 \textbf{Tr}(\rho \frac{\partial}{\partial \textbf{r}_i} H).
\end{equation}
Since $\frac{\partial}{\partial \textbf{r}_i} H^{jk}_{\alpha \beta} \neq 0$
only when $i = j$ or $i = k$, the above equation can be simplified to be
\begin{equation}
    \textbf{f}_i=
    - 2 \sum_k \sum_{\alpha \beta} \rho^{ik}_{\alpha \beta}
    \frac{\partial}{\partial \textbf{r}_i} H^{ik}_{\alpha \beta}
    - 2 \sum_j \sum_{\alpha \beta} \rho^{ji}_{\alpha \beta}
    \frac{\partial}{\partial \textbf{r}_i} H^{ji}_{\alpha \beta}
    = - 4 \sum_j \sum_{\alpha \beta} \rho^{ij}_{\alpha \beta}
    \frac{\partial}{\partial \textbf{r}_i} H^{ij}_{\alpha \beta}.
\end{equation}

The calculation of the Feynman-Hellmann matrix
\begin{equation}
    \textbf{K}^{ij}_{\alpha \beta} \equiv
    \frac{\partial}{\partial \textbf{r}_i} H^{ij}_{\alpha \beta}
\end{equation}
is straightforward. Since
\begin{equation}
    H^{ij}_{\alpha \beta} = H^{ij}_{\alpha \beta}(r_{ij} = r_0) s(r_{ij}),
\end{equation}
we have
\begin{equation}
    \textbf{K}^{ij}_{\alpha \beta} =  s(r_{ij})
    \frac{\partial}{\partial \textbf{r}_i} H^{ij}_{\alpha \beta}(r_{ij} = r_0)
    + H^{ij}_{\alpha \beta}(r_{ij} = r_0)
    \frac{\partial}{\partial \textbf{r}_i} s(r_{ij}).
\end{equation}
The remaining task is to calculate the matrix
\begin{equation}
    \textbf{G}^{ij}_{\alpha \beta} =
    \frac{\partial}{\partial \textbf{r}_i} H^{ij}_{\alpha \beta}(r_{ij} = r_0).
\end{equation}
By introducing some vectors
\begin{equation}
    \textbf{e}_{sx} = \left(\sin^2(x), -\cos(x)\cos(y), -\cos(x)\cos(z)\right),
\end{equation}
\begin{equation}
    \textbf{e}_{sy} = \left(-\cos(y)\cos(x), \sin^2(y), -\cos(y)\cos(z)\right),
\end{equation}
\begin{equation}
    \textbf{e}_{sz} = \left(-\cos(z)\cos(x), -\cos(z)\cos(y), \sin^2(z)\right),
\end{equation}
\begin{equation}
    \textbf{e}_{xx} = 2 \cos(x) \textbf{e}_{sx}
\end{equation}
\begin{equation}
    \textbf{e}_{yy} = 2 \cos(y) \textbf{e}_{sy}
\end{equation}
\begin{equation}
    \textbf{e}_{zz} = 2 \cos(z) \textbf{e}_{sz}
\end{equation}
\begin{equation}
    \textbf{e}_{xy} =
    \left(
        (\sin^2(x)-\cos^2(x)) \cos(y),
        (\sin^2(y)-\cos^2(y)) \cos(x),
        - 2 \cos(x) \cos(y) \cos(z)
    \right),
\end{equation}
\begin{equation}
    \textbf{e}_{yz} =
    \left(
        - 2 \cos(x) \cos(y) \cos(z),
        (\sin^2(y)-\cos^2(y)) \cos(z),
        (\sin^2(z)-\cos^2(z)) \cos(y)
    \right),
\end{equation}
\begin{equation}
    \textbf{e}_{zx} =
    \left(
        (\sin^2(x)-\cos^2(x)) \cos(z),
        - 2 \cos(x) \cos(y) \cos(z),
        (\sin^2(z)-\cos^2(z)) \cos(x)
    \right),
\end{equation}
we can write the matrix elements of this matrix to be
\begin{eqnarray}
    \textbf{G}^{ij}_{ss} &=& 0,                                           \\
    \textbf{G}^{ij}_{xx} &=& \frac{1}{r_{ij}}
                             (V_{pp\sigma} - V_{pp\pi})\textbf{e}_{xx},   \\
    \textbf{G}^{ij}_{yy} &=& \frac{1}{r_{ij}}
                             (V_{pp\sigma} - V_{pp\pi})\textbf{e}_{yy},   \\
    \textbf{G}^{ij}_{zz} &=& \frac{1}{r_{ij}}
                             (V_{pp\sigma} - V_{pp\pi})\textbf{e}_{zz},   \\
    \textbf{G}^{ij}_{sx} &=& \frac{1}{r_{ij}}
                             V_{sp\sigma} \textbf{e}_{sx},                \\
    \textbf{G}^{ij}_{sy} &=& \frac{1}{r_{ij}}
                             V_{sp\sigma} \textbf{e}_{sy},                \\
    \textbf{G}^{ij}_{sz} &=& \frac{1}{r_{ij}}
                             V_{sp\sigma} \textbf{e}_{sz},                \\
    \textbf{G}^{ij}_{xy} &=& \frac{1}{r_{ij}}
                             (V_{pp\sigma} - V_{pp\pi}) \textbf{e}_{xy},  \\
    \textbf{G}^{ij}_{yz} &=& \frac{1}{r_{ij}}
                             (V_{pp\sigma} - V_{pp\pi}) \textbf{e}_{yz},  \\
    \textbf{G}^{ij}_{zx} &=& \frac{1}{r_{ij}}
                             (V_{pp\sigma} - V_{pp\pi}) \textbf{e}_{zx},  \\
    \textbf{G}^{ij}_{xs} &=& - \textbf{F}^{ij}_{sx},                      \\
    \textbf{G}^{ij}_{ys} &=& - \textbf{F}^{ij}_{sy},                      \\
    \textbf{G}^{ij}_{zs} &=& - \textbf{F}^{ij}_{sz},                      \\
    \textbf{G}^{ij}_{yx} &=& \textbf{F}^{ij}_{xy},                        \\
    \textbf{G}^{ij}_{zy} &=& \textbf{F}^{ij}_{yz},                        \\
    \textbf{G}^{ij}_{xz} &=& \textbf{F}^{ij}_{zx}.
\end{eqnarray}


### 紧束缚分子动力学

要让一个紧束缚模型用于分子动力学模拟，必须再在band能量的基础上加入一个所谓的排斥势：


$$
U_{\rm tot} = U_{\rm bs} + U_{\rm rep}
$$

我们这里只考虑一个具体的碳材料的紧束缚模型 (C H Xu, C Z Wang, C T Chan and K M Ho, A transferable tight-binding potential for carbon, https://iopscience.iop.org/article/10.1088/0953-8984/4/28/006)

在这个模型中，排斥势被定义为类似于 EAM 势的多体势的形式：

$$
U_{\rm rep} = \sum_i f \left(  \sum_j \phi(r_{ij}) \right)
$$

函数 $\phi$ 的形式为：

$$
\phi(r) = \phi_0 (d_0/r)^m e^{ m [ -(r/d_c)^{m_c} + (d_0/d_c)^{m_c}]}
$$

函数 $s$ 的形式为：

$$
s(r) = (d_0/r)^n e^{ n [ -(r/r_c)^{n_c} + (r_0/r_c)^{n_c}]}
$$



### 基于紧束缚模型的电子输运性质

#### 与分子动力学模拟耦合的线性标度量子输运

这部分介绍最近在GPUMD实现的LSQT方法，详见 [Z Fan, Y Xiao, Y Wang, P Ying, S Chen, H Dong, Combining linear-scaling quantum transport and machine-learning molecular dynamics to study thermal and electronic transports in complex materials](https://arxiv.org/abs/2310.15314)。

类似于热导率的格林-久保公式，电导率可以表达为电流自关联的积分。因为电流密度等于电子电量乘以速度，所以我们也可以用速度自关联进行讨论。下面是用速度自关联表达的电导率公式：

$$
    \Sigma(E,t)=\frac{2e^2}{\Omega} \int_0^{t} \mathrm{Tr} \left[\delta (E-\hat{H}) \mathrm{Re} (\hat{V}\hat{V}(\tau)) \right] d\tau,
$$

该式表明，电导率是积分上限 $t$ 和能量 $E$ 的函数。这里 $e$ 就是基本电荷（即一个电子的电荷），  $\Omega$ 是体系的体积， $\hat{H}$ 是电子的紧束缚哈密顿量算符， $\hat{V}$ 是电子的速度算符，  $\delta(E-\hat{H})$  是能量分辨算符。

$$\hat{V}(\tau)=e^{i\hat{H}\tau} \hat{V} e^{-i\hat{H}\tau}$$

是时间演化后的速度算符。将原子体系的时间演化和电子体系的时间演化耦合起来，就能有效地描述电声耦合。

需要进一步讨论 LSQT 的算法才能讲清楚。


电子态密度的计算是类似的，所以我们先来讨论它的计算

$$
    \rho(E)=\frac{2}{\Omega}  \mathrm{Tr} \left[\delta (E-\hat{H})  \right].
$$

#### 应用实例：石墨烯纳米结构的热电输运

作为例子，我们研究石墨烯反点格，也叫做石墨烯纳米网，其结构如下图所示：

 
该模型的坐标文件 model.xyz 可在本章程序库找到。该结构有 187200 个原子，在 xy 平面内的尺寸为 88.5 纳米乘以 76.7 纳米。对于两维材料，通常认为地取一个厚度，我们这里按照文献中的约定将厚度取为 0.335 纳米。 

对于该结构，我们使用一个非常简单的单  $p_z$ 轨道的紧束缚模型。 在该模型中，只有最近邻的碳原子之间才有如下跃迁矩阵元：

$$
    H_{ij} = t_0 \left( \frac{r_0}{r_{ij}}\right)^2,
$$

这里 $t_0=-2.7$ eV,  $r_0=0.142$ 纳米。 $r_{ij}$ 是原子 $i$ 和 $j$ 之间的距离。

该体系的实空间哈密顿算符可写为：

$$
    \hat{H} = \sum_{i,j} H_{ij} |i \rangle\langle j|;
$$

如果假设输运方向为 x， 那么速度算符可写为：

$$
    \hat{V} = \frac{i}{\hbar}\sum_{i,j} (x_j-x_i)H_{ij} |i \rangle\langle j|,
$$

其中 $x_i$ 是 $i$ 原子的 x 坐标。
