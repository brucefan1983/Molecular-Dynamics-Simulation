
# 《分子动力学模拟入门》第三章：多体势函数

本章介绍分子动力学模拟中的势函数，重点讨论几个典型的多体势函数，包括 Embedded-atom-method (EAM) 势、Tersoff 势以及一个机器学习势。

# Table of contents
- [经典势函数的一般性质](#经典势函数的一般性质)
  - [两体势与多体势的定义](#两体势与多体势的定义)
  - [多体势中力的表达式](#多体势中力的表达式)
  - [多体势中位力和热流的表达式](#多体势中位力和热流的表达式)
- [两个典型的经验多体势](#两个典型的经验多体势)
  - [EAM势](#EAM势)
    - [EAM势的通用表达式](#EAM势的通用表达式)
    - [一个解析EAM势的编程实现](#一个解析EAM势的编程实现)
  - [Tesoff势](#Tersoff势)
    - [势函数形式](#势函数形式)
    - [编程实现](#编程实现)
- [NEP机器学习势](#NEP机器学习势)
  - [NEP机器学习势的人工神经网络模型](#NEP机器学习势的人工神经网络模型)
  - [NEP机器学习势的描述符](#NEP机器学习势的描述符)
  - [NEP机器学习势的训练](#NEP机器学习势的训练)
    - [损失函数](#损失函数)
    - [自然演化策略](#自然演化策略)
  - [NEP机器学习势的编程实现](#NEP机器学习势的编程实现)
- [习题](#习题)


## 经典势函数的一般性质

经典势函数，简称经典势（classical potentials)，是相对于量子力学势（quantum-mechanical potentials) 来说的。一个经典多粒子系统的势能可以写为各个粒子坐标的函数。文献中也将经典势函数称为力场（force field）或者势能面（potential energy surface）。一般来说，物理学家更乐意用势函数的称呼，而化学家更乐意用力场的称呼。笔者不是化学家，故后面基本上都用势函数的说法。

相对于量子力学的计算，经典势要高效得多。大多数的量子力学计算的复杂度都不低于 $\mathcal{O}(N^3)$, 而大部分经典势的计算复杂度都是线性标度的，即 $\mathcal{O}(N)$。这里的 $N$ 是体系的原子数。正是计算的高效性让经典势成为分子动力学模拟中重要的研究对象。

从某种程度上说，经典势可以分为经验势（empirical potentials）和机器学习势（machine-leared potentials）两大类。其中，经验势一般是用具有物理或者化学意义的数学函数表达的，而机器学习势则涉及某种机器学习模型（比如线性回归和人工神经网络等）。

无论是经验势，还是机器学习势，都包含一些待定参数。经验势中的参数个数一般较少（几个到几十个不等），且一般通过拟合实验结果来确定。机器学习势中的参数个数一般较多（几百个到上百万个不等），且一般通过拟合量子力学计算结果来确定。

### 两体势与多体势的定义

如果可以将体系的势能写成

$$
U= \frac{1}{2}\sum_{i}\sum_{j \neq i} U_{ij} \left(r_{ij} \right),
$$

其中，

$$
U_{ij} \left(r_{ij} \right) = U_{ji} \left(r_{ji} \right),
$$

那么我们称该系统的相互作用势能为两体势。其中， $U_{ij} \left(r_{ij} \right)$ 代表粒子 $i$ 和 $j$ 之间的相互作用势能，仅仅依赖于两粒子的相对距离 $r_{ij}$。前两章讨论的 LJ 势就是一个典型的两体势。两体势系统的势能也可以写成如下等价的形式：

$$
U= \sum_{i}\sum_{j > i} U_{ij} \left(r_{ij} \right).
$$

如果一个体系的势能无法写成以上形式，那么我们称该势能为多体势。相对而言，多体势比两体势更加接近量子力学计算的结果，故在各种材料体系中应用得较为成功。

迄今已有众多的经验多体势，其中在材料模拟领域最为广泛使用的当属 EAM 势和 Tersoff 势。在有机物模拟领域，所谓的“拓扑力场”更为常用。这些拓扑力场往往是多种经验势的组合(包括两体的和多体的），包括比如键长项、键角项、二面角项、范德瓦尔斯作用项、静电作用项等。这里“拓扑”的含义指的就是键长、键角和二面角等相互作用是针对具有固定拓扑结构的体系定义的。比如，某三个原子之间的键角相互作用总是存在，而不管它们的位置如何。由于笔者对拓扑力场没有研究经验，故本书不讨论拓扑力场。

### 多体势中力的表达式

为了推导多体势的一系列表达式，我们假设一个多体势系统的总能量可以写为各个粒子的能量之和

$$
U = \sum_i U_i.
$$

其中， $U_i$ 称为粒子 $i$ 的能量，它依赖于各个从 $i$ 指向其它粒子的位置矢量差 $\vec{r}_{ij}$：

$$
U_i = U_i\left(\vec{r}_{ij}\right).
$$

该表达式显然满足空间平移不变性，但我们还没有对其施加空间转动不变性。后面我们会看到，EAM 势和 Tersoff 势都满足这个假设，但它们都有额外的限制。以后我们还会看到，最近发展迅猛的机器学习势也满足这个假设。

从以上假设出发，可以推导出如下力的表达式：

$$
\vec{F}_{i} = \sum_{j \neq i} \vec{F}_{ij};
$$

$$
\vec{F}_{ij} = - \vec{F}_{ji} =
\frac{\partial U_{i}}{\partial \vec{r}_{ij}} -
\frac{\partial U_{j}}{\partial \vec{r}_{ji}} =
\frac{\partial \left(U_{i} + U_{j}\right) }{\partial \vec{r}_{ij}}.
$$

这里，


$$
\partial U_{i}/\partial \vec{r}_{ij} =
\partial U_{i}/\partial x_{ij} \vec{e}_x +
\partial U_{i}/\partial y_{ij} \vec{e}_y +
\partial U_{i}/\partial z_{ij} \vec{e}_z
$$

以上结果由[笔者和合作者](https://doi.org/10.1103/PhysRevB.92.094301) 于 2015 年推导出来，详细证明如下。

我们从保守力的定义出发。粒子 $i$ 的力等于体系总势能对粒子坐标的梯度的负值：

$$
\vec{F}_{i} = - \frac{\partial U}{\partial \vec{r}_{i}}
$$

代入总能量表达式，得

$$
\vec{F}_{i} = - \frac{\partial \sum_j U_j}{\partial \vec{r}_{i}}
$$

注意，为了避免混淆指标，上式中的求和不能写成原先的 $\sum_i U_i$，这是在推导公式时要特别注意的。接下来的任务就是推导上式中的偏导数了。为此，我们注意到 $U_{j}$ 是所有 $\vec{r}_{jk}$ 的函数，于是有

$$
    \frac{\partial U_j}{\partial \vec{r}_{i}} = \sum_k \frac{\partial U_j}{\partial \vec{r}_{jk}} \frac{\partial \vec{r}_{jk}}{\partial \vec{r}_{i}}
$$

因为

$$
\frac{\partial \vec{r}_{jk}}{\partial \vec{r}_{i}} = \frac{\partial (\vec{r}_{k} - \vec{r}_{j})}{\partial \vec{r}_{i}} = \frac{\partial \vec{r}_{k}}{\partial \vec{r}_{i}} - \frac{\partial \vec{r}_{j}}{\partial \vec{r}_{i}} = \delta_{ki}-
\delta_{ji},
$$

所以有

$$
    \frac{\partial U_j}{\partial \vec{r}_{i}} =  \frac{\partial U_j}{\partial \vec{r}_{ji}} - \sum_k \frac{\partial U_j}{\partial \vec{r}_{jk}} \delta_{ji}
$$

$$
\vec{F}_{i} = - \sum_j \left(\frac{\partial U_j}{\partial \vec{r}_{ji}} - \sum_k \frac{\partial U_j}{\partial \vec{r}_{jk}} \delta_{ji}\right) = \sum_k \frac{\partial U_i}{\partial \vec{r}_{ik}} - \sum_j \frac{\partial U_j}{\partial \vec{r}_{ji}}
= \sum_j \left(\frac{\partial U_i}{\partial \vec{r}_{ij}} - \frac{\partial U_j}{\partial \vec{r}_{ji}} \right).
$$

根据以上公式，我们以说多体势的力满足牛顿第三定律的弱形式，但不一定满足牛顿第三定律的强形式。也就是说，我们可以定义两个粒子之间的力  $\vec{F}_{ij}$ 和 $\vec{F}_{ji}$, 他们大小相等、方向相反，但不一定作用在两个粒子所在直线上。

### 多体势中位力和热流的表达式

从力的表达式出发，可以推导出如下位力的表达式：

$$
\mathbf{W} = \sum_i \mathbf{W}_{i}.
$$

$$
\mathbf{W}_{i} = \vec{r}_{ij} \otimes \frac{\partial U_j}{\partial \vec{r}_{ji}} .
$$

我们还可以推导出如下和相互作用有关的热流的表达式：

$$
\vec{J} = \sum_i \vec{J}_{i}.
$$

$$
\vec{J}_i = \mathbf{W}_{i} \cdot \vec{v}_i.
$$

这里的 $\vec{v}_i$ 是粒子 $i$ 的速度。

关于位力和热流的定义以及推导，涉及到很多的统计物理知识，笔者有空后会补上。

## 两个典型的经验多体势

迄今已有众多的经验多体势，其中在材料模拟领域最为广泛使用的当属 EAM 势和 Tersoff 势。

### EAM势

#### EAM势的通用表达式

EAM 势由若干人同时提出，包括 [Daw & Baskes](https://doi.org/10.1103/PhysRevLett.50.1285) 以及 [Finnis & Sinclair](https://doi.org/10.1080/01418618408244210)。

在EAM势中，原子 $i$ 的势能为

$$
U_i = \frac{1}{2} \sum_{j\neq i} \phi(r_{ij}) + F (\rho_i).
$$

这里，含有 $\phi(r_{ij})$ 的部分是两体势， $F(\rho_i)$ 即为嵌入势。嵌入势是 $i$ 粒子处电子密度 $\rho_i$ 的函数。粒子 $i$ 所在点的电子密度是由它的邻居贡献的：

$$
\rho_i = \sum_{j\neq i} f(r_{ij}).
$$

所以，对单元素体系来说，一个EAM势完全由如下三个函数确定： $\phi(r_{ij})$, $F(\rho_i)$ 和 $f(r_{ij})$。他们都是一元函数，可以用解析表达式确定，但为了提高计算速度和通用性，目前最为广泛使用的方式是用样条插值表示他们。

可以推导如下表达式：

$$
\frac{\partial U_i}{\partial \vec{r}_{ij}}
= \frac{1}{2}  \phi'(r_{ij})  \frac{\partial r_{ij}} {\partial \vec{r}_{ij}} +
F'(\rho_i)  f'(r_{ij}) \frac{\partial r_{ij}} {\partial \vec{r}_{ij}}.
$$

#### 一个解析EAM势的编程实现

### Tersoff势

#### 势函数形式

Tersoff 势有几个稍有不同的变体。我这里介绍 [Tersoff 在 1989 年发表的一篇文章中使用的形式](https://doi.org/10.1103/PhysRevB.39.5566)。为简单起见，我们考虑单种元素的势函数。本书不会涉及多元素体系的 Tersoff 势。

粒子 $i$ 的势能可以写为：

$$
U_i =  \frac{1}{2} \sum_{j \neq i} f_C(r_{ij}) \left[ f_R(r_{ij}) - b_{ij} f_A(r_{ij}) \right].
$$

其中， $f_{C}$ 是一个截断函数，当 $r_{ij} < R$ 时取值为 1，当 $r_{ij} > S$ 时取值为 0。在这之间，该函数为

$$
f_{C}(r_{ij}) = \frac{1}{2}
\left[
1 + \cos \left( \pi \frac{r_{ij} - R_{ij}}{S_{ij} - R_{ij}} \right)
\right].
$$

排斥函数 $f_{R}$ 和吸引函数 $f_{A}$ 为

$$
f_{R}(r) = A e^{-\lambda r_{ij}};
$$

$$
f_{A}(r) = B e^{-\mu r_{ij}}.
$$

键序为

$$
\label{equation:bij}
b_{ij} = \left(1 + \beta^{n} \zeta^{n}_{ij}\right)^{-\frac{1}{2n}},
$$

其中，

$$
\zeta_{ij} = \sum_{k\neq i, j}f_C(r_{ik}) g_{ijk},
$$

$$
g_{ijk} = 1 + \frac{c^2}{d^2} - \frac{c^2}{d^2+(h-\cos\theta_{ijk})^2}.
$$

在以上表达式中，有如下参数： $A$, $B$, $\lambda$, $\mu$, $\beta$, $n$, $c$, $d$, $h$, $R$, $S$。

#### 编程实现

给出一个 C++ 程序。

## NEP机器学习势

机器学习势近年来得到了迅猛的发展。当前主流的机器学习势框架由 [Behler 和 Parrinello](https://doi.org/10.1103/PhysRevLett.98.146401) 建立。该机器学习势采用人工神经网络作为机器学习模型，这也是当前主流的机器学习模型之一。在该机器学习势框架中，体系的总势能依然表达为各个原子的势能的和。其中，每个原子的势能则表达为同一个神经网络模型的输出。作为神经网络模型的输入，则采用若干所谓的对称函数，也称为描述符。所有的描述符函数组成一个矢量，称为描述符矢量。该矢量不是三维空间的，而是一个抽象空间的。描述符矢量中的每一个分量都是原子坐标的某种函数，而且都满足空间平移和转动不变性，也满足同类原子置换的不变性。每一个原子的描述符都反映了该原子的局部化学环境，而人工神经网络模型就是将不同的化学环境与不同的原子势能联系起来。

当前已有众多的机器学习势方法，它们在各种细节上有些差别，但总的来说都很相似。笔者于2021年开发了一个名为 [NEP (Neuroevolution potential)](https://doi.org/10.1103/PhysRevB.104.104309) 的机器学习势方法，本节将以其为例介绍更多机器学习势的细节。

### NEP机器学习势的人工神经网络模型

人工神经网络是一个比较简单但应用广泛的机器学习模型。从数学函数的角度来看，一个神经网络模型是一个多元函数：

$$
U = U(q_0, q_1, q_2, ...).
$$

这里的自变量 $q_i$ 组成一个高维空间的矢量，构成了神经网络的输入层，而 $U$ 就是神经网络的输出层。在我们讨论的机器学习势中，输入层就是前面提到的描述符矢量，而输出层是一个标量，就是某个原子的势能。从输入层到输出层的映射，即函数 $U$，一般来说是一个非线性函数，它的具体形式取决于神经网络中有多少个隐藏层（hidden layers）。在NEP机器学习势中，所用神经网络仅有一个隐藏层。用更多的隐藏层有可能提高势函数的精度，但也会增加势函数的计算量。在经过大量测试与权衡后，笔者决定仅用一个隐藏层。

在仅有一个隐藏层的情况下，NEP势函数的神经网络模型可表达为如下复合函数


$$
U = \sum_{\mu=1}^{N_\mathrm{neu}}w^{(1)}_{\mu} x_{\mu} - b^{(1)},
$$

$$
x_{\mu} = \tanh\left(\sum_{\nu=1}^{N_\mathrm{des}} w^{(0)}_{\mu\nu} q_{\nu} - b^{(0)}_{\mu}\right),
$$

这里的 $N_\mathrm{des}$ 是输入层的维度，也就是描述符的维度，而 $N_\mathrm{neu}$ 是隐藏层的维度，即隐藏层的神经元个数。该复合函数可以形象地由下图表示：

![nep](fig/nep.png)

从以上表达式可以看出，从输入层到隐藏层，首先经过了一个线性变换

$$
y_{\mu} = \sum_{\nu=1}^{N_\mathrm{des}} w^{(0)}_{\mu\nu} q_{\nu} - b^{(0)}_{\mu},
$$

然后经过了一个非线性变换

$$
x_{\mu} = \tanh (y_{\mu}).
$$

这里的非线性变换函数 $\tanh$ 是神经网络中的激活函数 (activation function)。笔者尝试过其它的激活函数，没有发现有明显更好的。从隐藏层到输出层，仅有一个线性变换。以上公式中, $w^{(0)}$ 是连接输入层和隐藏层的权重, $w^{(1)}$ 是连接隐藏层和输出层的权重, $b^{(0)}$ 是隐藏层的偏置, $b^{(1)}$ 是输出层的偏置。这些参数都是需要通过训练来确定的，稍后将对此进行进一步讨论。


### NEP机器学习势的描述符

机器学习势中的描述符是用来描述原子周围的环境的，即一个原子的近邻原子的分布情况。有一种描述符只依赖于邻居原子与所考虑的中心原子的距离，称为径向描述符。另外一种描述符与邻居原子相对于中心原子的方向也有关，称为角度描述符。下面我们一一介绍。

#### 径向描述符

对于某个中心原子 $i$，NEP 中的径向描述符分量定义为

$$
q^i_n = \sum_i g_n(r_{ij})
$$

也就是说，我们一共有若干（由下标 $n$ 标记）径向描述符分量。每一个径向描述符分量是某个径向函数 $g_n(r_{ij})$ 对近邻的求和。对比前面的 EAM 势，会发现一个径向描述符分量就类似于 EAM 势中的电子密度。之所以需要用若干径向描述符分量，是为了对近邻原子相对于中心原子的距离进行细致的刻画。

那么径向函数 $g_n(r_{ij})$ 应该有怎样的形式呢？在 [Behler 和 Parrinello](https://doi.org/10.1103/PhysRevLett.98.146401) 提出的方案中，径向函数是一些高斯分布函数，它们具有不同的分布中心和宽度。在 NEP 中，径向函数 $g_n(r_{ij})$ 是另一套径向基函数的线性叠加：

$$
   g_n(r_{ij}) = \sum_{k=0}^{N_\mathrm{bas}^\mathrm{R}} c^{ij}_{nk} f_k(r_{ij}),
$$

其中的径向基函数通过 Chebyshev 多项式定义：

$$
   f_k(r_{ij}) = \frac{1}{2}
   \left[T_k\left(2\left(r_{ij}/r_\mathrm{c}^\mathrm{R}-1\right)^2-1\right)+1\right]
   f_\mathrm{c}(r_{ij}),
$$

这里的 $T_k(x)$ 是 $k$-阶第一类 Chebyshev 多项式。上式中的 $f_\mathrm{c}(r_{ij})$ 是一个光滑截断函数，类似于 Tersoff 势的截断函数：

$$
   f_\mathrm{c}(r_{ij}) 
   = \begin{cases}
   \frac{1}{2}\left[
   1 + \cos\left( \pi \frac{r_{ij}}{r_\mathrm{c}^\mathrm{R}} \right) 
   \right],& r_{ij}\leq r_\mathrm{c}^\mathrm{R}; \\
   0, & r_{ij} > r_\mathrm{c}^\mathrm{R}.
   \end{cases}
$$

下图展示了前 5 个径向基函数，其中的 $f_0$ 实际上就是截断函数，因为 $T_0(x) = 1$。图中的横坐标是约化距离 $r_{ij}/r_\mathrm{c}^\mathrm{R}$，其中 $r_\mathrm{c}^\mathrm{R}$ 是径向描述符的截断距离。

![radial_basis_functions](fig/basis.png)

#### 角度描述符

NEP 中的角度描述符包含所谓的三体、四体和五体描述符。我们这里仅介绍三体描述符，关于NEP中四体和五体描述符的细节，请参考如下文章 [GPUMD: A package for constructing accurate machine-learned potentials and performing highly efficient atomistic simulations](https://doi.org/10.1063/5.0106617).

角度描述符需要有键角的信息，所以表达式中一定要有至少三个原子之间的键角。类似于径向描述符，我们可以构造出如下三体角度描述符：

$$
   q^i_{nl} = \sum_{j \neq i} \sum_{k \neq i} g_n(r_{ij}) g_n(r_{ij}) P_l(\theta_{ijk}),
$$

可以看到，三体角度描述符相比两体的径向描述符多了一个下标 $l$。这个下标是勒让德多项式 $P_l(\theta_{ijk})$ 的阶数。勒让德多项式是键角 $\theta_{ijk}$ 的函数。该键角是以原子 $i$ 为中心，以键 $ij$ 和 $ik$ 为两边的夹角。该式中的函数 $g_n(r_{ij})$ 和径向描述符中的对应函数一致，只不过可能有不同的截断距离 $r_\mathrm{c}^\mathrm{A}$。

上述表达式中的双重求和保证了同类原子置换的不变性，但导致计算量正比于近邻个数的平方。有一个办法可以在不改变结果的前提下降低计算复杂度，那就是利用球谐函数的加法定理。利用该定理（练习题），可以将上述三体角度描述符写成如下等价的形式：

$$
   q^i_{nl} = \sum_{m=-l}^l (-1)^m A^i_{nlm} A^i_{nl(-m)},
$$

$$
   A^i_{nlm} = \sum_{j\neq i} g_n(r_{ij}) Y_{lm}(\theta_{ij},\phi_{ij}),
$$

其中， $Y_{lm}(\theta_{ij},\phi_{ij})$ 是球谐函数。 $\theta_{ij}$ 是极角， $\phi_{ij}$ 是方位角。







### NEP机器学习势的训练

#### 损失函数

#### 自然演化策略

### NEP机器学习势的编程实现

NEP 机器学习势已经在 [GPUMD 程序包](https://github.com/brucefan1983/GPUMD) 中实现，可用该程序包的 `nep` 可执行文件训练，并由该程序包的 `gpumd` 可执行文件进行分子动力学模拟。本章暂不讨论 GPUMD 程序包的使用。

NEP 机器学习势目前也有一个独立的 C++ 编程实现，见 [NEP_CPU 程序包](https://github.com/brucefan1983/NEP_CPU)。该程序包给出了一个名为`NEP3` 的 C++ 类，见程序包中 `src` 文件夹内的 `nep.cpp` 和`nep.h`。我们下面的测试将使用该 C++ 实现。

## 习题

1. 将本章的 tersoff_md.cpp 从单元素体系推广到多元素体系，使其能模拟比如 SiC 和 SiGe。

2. 将本章的 eam_md.cpp 从单元素体系推广到多元素体系，使其能模拟高熵合金。

3. 利用球谐函数的加法定理，将三体角度描述符

$$
   q^i_{nl} = \sum_{j \neq i} \sum_{k \neq i} g_n(r_{ij}) g_n(r_{ij}) P_l(\theta_{ijk}),
$$

推导为如下形式：

$$
   q^i_{nl} = \sum_{m=-l}^l (-1)^m A^i_{nlm} A^i_{nl(-m)},
$$

$$
   A^i_{nlm} = \sum_{j\neq i} g_n(r_{ij}) Y_{lm}(\theta_{ij},\phi_{ij}).
$$

