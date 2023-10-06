
# 《分子动力学模拟》第三章：多体势函数

本章介绍分子动力学模拟中的势函数，重点讨论几个典型的多体势函数，包括 Embedded-atom-method (EAM) 势、Tersoff 势以及一个机器学习势。

# Table of contents
- [经典势函数的一般性质](#经典势函数的一般性质)
  - [两体势与多体势的定义](#两体势与多体势的定义)
  - [多体势中力的表达式](#多体势中力的表达式)
  - [多体势中位力和热流的表达式](#多体势中位力和热流的表达式)
- [两个典型的经验多体势](#两个典型的经验多体势)
  - [EAM势](#EAM势)
  - [Tesoff势](#Tersoff势)
- [NEP机器学习势](#NEP机器学习势)
  - [NEP机器学习势的人工神经网络模型](#NEP机器学习势的人工神经网络模型)
  - [NEP机器学习势的描述符](#NEP机器学习势的描述符)
  - [NEP机器学习势的训练](#NEP机器学习势的训练)
  - [NEP机器学习势的编程实现](#NEP机器学习势的编程实现)
- [习题](#习题)


## 经典势函数的一般性质

经典势函数，简称经典势（classical potentials)，是相对于量子力学势（quantum-mechanical potentials) 来说的。一个经典多粒子系统的势能可以写为各个粒子坐标的函数。文献中也将经典势函数称为力场（force field）或者势能面（potential energy surface）。本书将根据上下文的需要随意使用这些类似的术语。

相对于量子力学的计算，经典势要高效得多。大多数的量子力学计算的复杂度都不低于 $\mathcal{O}(N^3)$, 而大部分经典势的计算复杂度都是线性标度的，即 $\mathcal{O}(N)$。这里的 $N$ 是体系的原子数。正是计算的高效性让经典势成为分子动力学模拟中重要的研究对象。

从某种程度上说，经典势可以分为经验势（empirical potentials）和机器学习势（machine-learned potentials）两大类。其中，经验势一般是用具有物理或化学意义的数学函数表达的，而机器学习势则涉及某种机器学习模型（如线性回归和人工神经网络等）。

无论是经验势，还是机器学习势，都包含一些待定参数。经验势中的参数个数一般较少（几个到几十个不等），且一般通过拟合实验结果来确定。机器学习势中的参数个数一般较多（几百个到上百万个不等），且一般通过拟合量子力学计算结果来确定。

### 两体势与多体势的定义

如果可以将体系的势能写成

$$
U= \frac{1}{2}\sum _{i}\sum _{j \neq i} U _{ij} \left(r _{ij} \right),
$$

其中，

$$
U _{ij} \left(r _{ij} \right) = U _{ji} \left(r _{ji} \right),
$$

那么我们称该系统的相互作用势能为两体势。其中， $U _{ij} \left(r _{ij} \right)$ 代表粒子 $i$ 和 $j$ 之间的相互作用势能，仅仅依赖于两粒子的相对距离 $r _{ij}$。前一章讨论的 LJ 势就是一个典型的两体势。两体势系统的势能也可以写成如下等价的形式：

$$
U= \sum _{i}\sum _{j > i} U _{ij} \left(r _{ij} \right).
$$

如果一个体系的势能无法写成以上形式，那么我们称该势能为多体势。相对而言，多体势比两体势更加接近量子力学计算的结果，故在各种材料体系中应用得较为成功。

迄今已有众多的经验多体势，其中在材料模拟领域最为广泛使用的是EAM（embedded atom method）势和 Tersoff 势。在有机物模拟领域，目前最为常用的是不考虑化学键断裂与重组以及电荷转移的分子力场。这些分子力场往往是多种经验势的组合(包括两体的和多体的），包括比如键长项、键角项、二面角项、范德瓦尔斯作用项、静电作用项等。这里的键长、键角和二面角相互作用都是针对固定的化学键定义的（由所谓的“拓扑”结构定义）。比如，某三个原子之间的键角相互作用总是存在，而不管它们的位置如何。由于笔者对该类型的分子力场没有研究经验，故本书不讨它们。

### 多体势中力的表达式

为了推导多体势的一系列表达式，我们假设一个多体势系统的总能量可以写为各个粒子的能量之和

$$
U = \sum _i U _i.
$$

其中， $U _i$ 称为粒子 $i$ 的能量，它依赖于各个从 $i$ 指向其它粒子的位置矢量差 $\vec{r} _{ij} \equiv \vec{r} _{j} - \vec{r} _{i}$：

$$
U _i = U _i\left(\lbrace \vec{r} _{ij} \rbrace\right).
$$

该表达式显然满足空间平移不变性，但我们还没有对其施加空间转动不变性。后面我们会看到，EAM 势和 Tersoff 势都满足空间转动不变性，但它们都有额外的限制。以后我们还会看到，最近发展迅猛的机器学习势也满足空间转动不变性。

从以上假设出发，可以推导出如下力的表达式：

$$
\vec{F} _{i} = \sum _{j \neq i} \vec{F} _{ij};
$$

$$
\vec{F} _{ij} = - \vec{F} _{ji} =
\frac{\partial U _{i}}{\partial \vec{r} _{ij}} -
\frac{\partial U _{j}}{\partial \vec{r} _{ji}} =
\frac{\partial \left(U _{i} + U _{j}\right) }{\partial \vec{r} _{ij}}.
$$

这里，


$$
\partial U _{i}/\partial \vec{r} _{ij} =
\partial U _{i}/\partial x _{ij} \vec{e} _x +
\partial U _{i}/\partial y _{ij} \vec{e} _y +
\partial U _{i}/\partial z _{ij} \vec{e} _z
$$

以上结果由[笔者和合作者](https://doi.org/10.1103/PhysRevB.92.094301) 于 2015 年推导出来，详细证明如下。

我们从保守力的定义出发。粒子 $i$ 的力等于体系总势能对粒子坐标的梯度的负值：

$$
\vec{F} _{i} = - \frac{\partial U}{\partial \vec{r} _{i}}
$$

代入总能量表达式，得

$$
\vec{F} _{i} = - \frac{\partial \sum _j U _j}{\partial \vec{r} _{i}}
$$

注意，为了避免混淆指标，上式中的求和不能写成原先的 $\sum _i U _i$，这是在推导公式时要特别注意的。接下来的任务就是推导上式中的偏导数了。为此，我们注意到 $U _{j}$ 是所有 $\vec{r} _{jk}$ 的函数，于是有

$$
    \frac{\partial U _j}{\partial \vec{r} _{i}} = \sum _k \frac{\partial U _j}{\partial \vec{r} _{jk}} \frac{\partial \vec{r} _{jk}}{\partial \vec{r} _{i}}
$$

因为

$$
\frac{\partial \vec{r} _{jk}}{\partial \vec{r} _{i}} = \frac{\partial (\vec{r} _{k} - \vec{r} _{j})}{\partial \vec{r} _{i}} = \frac{\partial \vec{r} _{k}}{\partial \vec{r} _{i}} - \frac{\partial \vec{r} _{j}}{\partial \vec{r} _{i}} = \delta _{ki}-
\delta _{ji},
$$

所以有

$$
    \frac{\partial U _j}{\partial \vec{r} _{i}} =  \frac{\partial U _j}{\partial \vec{r} _{ji}} - \sum _k \frac{\partial U _j}{\partial \vec{r} _{jk}} \delta _{ji}
$$

$$
\vec{F} _{i} = - \sum _j \left(\frac{\partial U _j}{\partial \vec{r} _{ji}} - \sum _k \frac{\partial U _j}{\partial \vec{r} _{jk}} \delta _{ji}\right) = \sum _k \frac{\partial U _i}{\partial \vec{r} _{ik}} - \sum _j \frac{\partial U _j}{\partial \vec{r} _{ji}}
= \sum _j \left(\frac{\partial U _i}{\partial \vec{r} _{ij}} - \frac{\partial U _j}{\partial \vec{r} _{ji}} \right).
$$

定义

$$
\vec{F} _{ij} = - \vec{F} _{ji} =
\frac{\partial U _{i}}{\partial \vec{r} _{ij}} -
\frac{\partial U _{j}}{\partial \vec{r} _{ji}} =
\frac{\partial \left(U _{i} + U _{j}\right) }{\partial \vec{r} _{ij}}.
$$

则有

$$
\vec{F} _{i} = \sum _{j \neq i} \vec{F} _{ij}.
$$

证毕。

根据以上公式，我们可以说多体势的力满足牛顿第三定律的弱形式，但不一定满足牛顿第三定律的强形式。也就是说，我们可以定义两个粒子之间的力  $\vec{F} _{ij}$ 和 $\vec{F} _{ji}$, 他们大小相等、方向相反，但不一定作用在两个粒子所在直线上。

### 多体势中位力和热流的表达式

从力的表达式出发，可以推导出如下位力的表达式：

$$
\mathbf{W} = \sum _i \mathbf{W} _{i}.
$$

$$
\mathbf{W} _{i} = \sum _j \vec{r} _{ij} \otimes \frac{\partial U _j}{\partial \vec{r} _{ji}} .
$$

我们还可以推导出如下和相互作用有关的热流的表达式：

$$
\vec{J} = \sum _i \vec{J} _{i}.
$$

$$
\vec{J} _i = \mathbf{W} _{i} \cdot \vec{v} _i.
$$

这里的 $\vec{v} _i$ 是粒子 $i$ 的速度。

关于位力，我们从其定义出发：

$$
\mathbf{W} = \sum _i \vec{r}_i \otimes \vec{F}_i.
$$

该定义看起来很简单，但是它并不适合编程实现，因为它依赖于原子的绝对坐标。我们需要由此出发推导一个仅依赖于原子相对坐标的表达式。

将力的表达式代入上式，得

$$
\mathbf{W} = \sum _i  \sum _{j \neq i} \vec{r}_i \otimes \vec{F} _{ij}.
$$

交换两个求和哑指标，可得

$$
\mathbf{W} = \sum _i  \sum _{j \neq i} \vec{r}_j \otimes \vec{F} _{ji}.
$$

将以上两式相加再除以二可得

$$
\mathbf{W} = -\frac{1}{2} \sum _i  \sum _{j \neq i} \vec{r} _{ij} \otimes \vec{F} _{ij}.
$$

该式已经可以很方便地编程实现了。但是，为了方便地从位力计算热流，我们后面还需要推导出另一个等价的位力表达式。我们从如下定义出发推导热流的表达式：

$$
\vec{J} = \frac{d}{dt} \sum _i \vec{r} _{i} E _{i}.
$$

其中, $E _{i}$ 是粒子 $i$ 的总能量（动能与势能之和）。 上述求和所代表的物理量可称为能量矩 （energy moment)，而热流就是能量矩的时间变化率。将上式的求导展开得

$$
\vec{J} = \sum _i \vec{v} _{i} E _{i} + \sum _i \vec{r} _{i} \frac{d}{dt} E _{i}.
$$

上式右边第一项叫做对流项，而第二项叫做势能项。我们可以记为

$$
\vec{J} = \vec{J}^{\rm conv} + \vec{J}^{\rm pot}
$$

对流项无需进一步推导，而利用动能定理

$$
\frac{d}{dt} K _{i} = \vec{F} _i \cdot \vec{v} _i
$$

可将势能项的热流写为：

$$
\vec{J}^{\rm pot} = \sum _i \vec{r} _{i} \vec{F} _i \cdot \vec{v} _i + \vec{r} _{i} \frac{d}{dt} U _{i} .
$$

根据力的表达式，我们有

$$
\sum _i \vec{r} _{i} \vec{F} _i \cdot \vec{v} _i = \sum _i \sum _{j \neq i} \vec{r} _{i} \vec{F} _{ij} \cdot \vec{v} _i  
= \sum _i \sum _{j \neq i} \vec{r} _{i} \left(\frac{\partial U _{i}}{\partial \vec{r} _{ij}} - \frac{\partial U _{j}}{\partial \vec{r} _{ji}}\right) \cdot \vec{v} _i .
$$

根据势能 $U _{i}$ 的表达式，我们还有

$$
\sum _i \vec{r} _{i} \frac{d}{dt} U _{i} = \sum _i \sum _{j \neq i} \vec{r} _{i} \frac{\partial U _{i} }{\partial \vec{r} _{ij}} \cdot ( \vec{v} _{j} - \vec{v} _{i}).
$$

将以上两式相加，得到如下势能项的热流：


$$
\vec{J}^{\rm pot} = \sum _i \sum _{j \neq i} \vec{r} _{i} \left(\frac{\partial U _{i}}{\partial \vec{r} _{ij}} \cdot \vec{v} _j - \frac{\partial U _{j}}{\partial \vec{r} _{ji}} \cdot \vec{v} _i \right)
$$

类似于位力的推导，我们可以将上式用相对坐标表达：

$$
\vec{J}^{\rm pot} = -\frac{1}{2} \sum _i \sum _{j \neq i} \vec{r} _{ij} \left(\frac{\partial U _{i}}{\partial \vec{r} _{ij}} \cdot \vec{v} _j - \frac{\partial U _{j}}{\partial \vec{r} _{ji}} \cdot \vec{v} _i \right)
$$

该表达式涉及到一个粒子及其邻居的速度，不利于编程实现。可以通过交换哑指标的方式将上式改写为如下等价的形式：

$$
\vec{J}^{\rm pot} = \sum _i \sum _{j \neq i} \vec{r} _{ij} \left(\frac{\partial U _{j}}{\partial \vec{r} _{ji}} \cdot \vec{v} _i \right)
$$

注意到上式等价于

$$
\vec{J}^{\rm pot} = \sum _i \sum _{j \neq i} \left( \vec{r} _{ij} \otimes \frac{\partial U _{j}}{\partial \vec{r} _{ji}} \right) \cdot \vec{v} _i 
$$

还注意到位力可以写为如下等价的形式：

$$
\mathbf{W} = \sum _i \sum _{j \neq i} \vec{r} _{ij} \otimes \frac{\partial U _j}{\partial \vec{r} _{ji}}.
$$

若定义如下单粒子的位力

$$
\mathbf{W}_i = \sum _{j \neq i} \vec{r} _{ij} \otimes \frac{\partial U _j}{\partial \vec{r} _{ji}}.
$$

则有

$$
\vec{J}^{\rm pot} = \sum _i \mathbf{W}_i \cdot \vec{v} _i 
$$

也就是说，势能项的热流是与单粒子位力紧密相关的。但是，要注意的是，上述表达式中的单粒子位力必须是我们定义的形式。实际上，文献中有不少使用错误的热流公式的例子，都是因为没有仔细推导，或者盲目地信任已有的编程实现。到目前为止，流行的 LAMMPS 程序中的热流计算对大部分的多体势函数来说都是错误的。相比之下，GPUMD 程序中的热流计算都使用了上述正确的公式。


## 两个典型的经验多体势

迄今已有众多的经验多体势，其中在材料模拟领域最为广泛使用的当属 EAM 势和 Tersoff 势。

### EAM势

EAM 势由若干人同时提出，包括 [Daw & Baskes](https://doi.org/10.1103/PhysRevLett.50.1285) 以及 [Finnis & Sinclair](https://doi.org/10.1080/01418618408244210)。

在EAM势中，原子 $i$ 的势能为

$$
U _i = \frac{1}{2} \sum _{j\neq i} \phi(r _{ij}) + F (\rho _i).
$$

这里，含有 $\phi(r _{ij})$ 的部分是两体势， $F(\rho _i)$ 即为嵌入势。嵌入势是 $i$ 粒子处电子密度 $\rho _i$ 的函数。粒子 $i$ 所在点的电子密度是由它的邻居贡献的：

$$
\rho _i = \sum _{j\neq i} f(r _{ij}).
$$

所以，对单元素体系来说，一个EAM势完全由如下三个函数确定： $\phi(r _{ij})$, $F(\rho _i)$ 和 $f(r _{ij})$。他们都是一元函数，可以用解析表达式确定，但为了提高计算速度和通用性，常用样条插值表示他们。GPUMD 程序中实现了 [Xiaowang Zhou 等人的解析版本](https://doi.org/10.1103/PhysRevB.69.144113)。

可以推导如下表达式：

$$
\frac{\partial U _i}{\partial \vec{r} _{ij}}
= \frac{1}{2}  \phi'(r _{ij})  \frac{\partial r _{ij}} {\partial \vec{r} _{ij}} +
F'(\rho _i)  f'(r _{ij}) \frac{\partial r _{ij}} {\partial \vec{r} _{ij}}.
$$

有了上述表达式，EAM 势的编程实现就很直接了。

### Tersoff势

Tersoff 势有几个稍有不同的变体。我这里介绍 [Tersoff 在 1989 年发表的一篇文章中使用的形式](https://doi.org/10.1103/PhysRevB.39.5566)。为简单起见，我们考虑单种元素的势函数。本书不会涉及多元素体系的 Tersoff 势。

粒子 $i$ 的势能可以写为：

$$
U _i =  \frac{1}{2} \sum _{j \neq i} f _C(r _{ij}) \left[ f _R(r _{ij}) - b _{ij} f _A(r _{ij}) \right].
$$

其中， $f _{C}$ 是一个截断函数，当 $r _{ij} < R$ 时取值为 1，当 $r _{ij} > S$ 时取值为 0。在这之间，该函数为

$$
f _{C}(r _{ij}) = \frac{1}{2}
\left[
1 + \cos \left( \pi \frac{r _{ij} - R}{S - R} \right)
\right].
$$

排斥函数 $f _{R}$ 和吸引函数 $f _{A}$ 为

$$
f _{R}(r _{ij}) = A e ^{-\lambda r _{ij}};
$$

$$
f _{A}(r _{ij}) = B e ^{-\mu r _{ij}}.
$$

键序为

$$
b _{ij} = \left(1 + \beta ^{n} \zeta ^{n} _{ij} \right) ^{-\frac{1}{2n}},
$$

其中，

$$
\zeta _{ij} = \sum _{k\neq i, j}f _C(r _{ik}) g _{ijk},
$$

$$
g _{ijk} = 1 + \frac{c^2}{d^2} - \frac{c^2}{d^2+(h-\cos\theta _{ijk})^2}.
$$

在以上表达式中，有如下参数： $A$, $B$, $\lambda$, $\mu$, $\beta$, $n$, $c$, $d$, $h$, $R$, $S$。

在 Tersoff 势的表达式中，有一个角 $\theta _{ijk}$ ，它指的是键 $ij$ 和 $ik$ 的夹角。显然有

$$
\cos\theta _{ijk} = \frac{\vec{r} _{ij} \cdot \vec{r} _ {ik}}{r _{ij}  r _ {ik} }.
$$

根据我们对多体势的推导，在编程实现时需要使用 $\frac{\partial U _i}{\partial \vec{r} _{ij}}$ 的表达式，它可以在[作者一篇文章的附录](https://doi.org/10.1103/PhysRevB.92.094301)找到。

## NEP机器学习势

机器学习势近年来得到了迅猛的发展。当前主流的机器学习势框架由 [Behler 和 Parrinello](https://doi.org/10.1103/PhysRevLett.98.146401) 建立。该机器学习势采用人工神经网络作为机器学习模型，这也是当前主流的机器学习模型之一。在该机器学习势框架中，体系的总势能依然表达为各个原子的势能的和。其中，每个原子的势能则表达为同一个神经网络模型的输出。作为神经网络模型的输入，则采用若干所谓的对称函数，也称为描述符。所有的描述符函数组成一个矢量，称为描述符矢量。该矢量不是三维空间的，而是一个抽象空间的。描述符矢量中的每一个分量都是原子坐标的某种函数，而且都满足空间平移和转动不变性，也满足同类原子置换的不变性。每一个原子的描述符都反映了该原子的局部化学环境，而人工神经网络模型就是将不同的化学环境与不同的原子势能联系起来。

当前已有众多的机器学习势方法，它们在各种细节上有些差别，但总的来说都很相似。笔者于2021年开发了一个名为 [NEP (Neuroevolution potential)](https://doi.org/10.1103/PhysRevB.104.104309) 的机器学习势方法，本节将以其为例介绍更多机器学习势的细节。

### NEP机器学习势的人工神经网络模型

人工神经网络是一个比较简单但应用广泛的机器学习模型。从数学函数的角度来看，一个神经网络模型是一个多元函数：

$$
U = U(q _0, q _1, q _2, ...).
$$

这里的自变量 $q _i$ 组成一个高维空间的矢量，构成了神经网络的输入层，而 $U$ 就是神经网络的输出层。在我们讨论的机器学习势中，输入层就是前面提到的描述符矢量，而输出层是一个标量，就是某个原子的势能。从输入层到输出层的映射，即函数 $U$，一般来说是一个非线性函数，它的具体形式取决于神经网络中有多少个隐藏层（hidden layers）。在NEP机器学习势中，所用神经网络仅有一个隐藏层。用更多的隐藏层有可能提高势函数的精度，但也会增加势函数的计算量。在经过大量测试与权衡后，笔者决定仅用一个隐藏层。

在仅有一个隐藏层的情况下，NEP势函数的神经网络模型可表达为如下复合函数


$$
U = \sum _{\mu=1} ^{N _\mathrm{neu}}w ^{(1)} _{\mu} x _{\mu} - b ^{(1)},
$$

$$
x _{\mu} = \tanh\left(\sum _{\nu=1} ^{N _\mathrm{des}} w ^{(0)} _{\mu\nu} q _{\nu} - b ^{(0)} _{\mu}\right),
$$

这里的 $N _\mathrm{des}$ 是输入层的维度，也就是描述符的维度，而 $N _\mathrm{neu}$ 是隐藏层的维度，即隐藏层的神经元个数。该复合函数可以形象地由下图表示：

![nep](fig/nep.png)

从以上表达式可以看出，从输入层到隐藏层，首先经过了一个线性变换

$$
y _{\mu} = \sum _{\nu=1} ^{N _\mathrm{des}} w ^{(0)} _{\mu\nu} q _{\nu} - b ^{(0)} _{\mu},
$$

然后经过了一个非线性变换

$$
x _{\mu} = \tanh (y _{\mu}).
$$

这里的非线性变换函数 $\tanh$ 是神经网络中的激活函数 (activation function)。笔者尝试过其它的激活函数，没有发现有明显更好的。从隐藏层到输出层，仅有一个线性变换。以上公式中, $w ^{(0)}$ 是连接输入层和隐藏层的权重, $w ^{(1)}$ 是连接隐藏层和输出层的权重, $b ^{(0)}$ 是隐藏层的偏置, $b ^{(1)}$ 是输出层的偏置。这些参数都是需要通过训练来确定的，稍后将对此进行进一步讨论。


### NEP机器学习势的描述符

机器学习势中的描述符是用来描述原子周围的环境的，即一个原子的近邻原子的分布情况。有一种描述符只依赖于邻居原子与所考虑的中心原子的距离，称为径向描述符。另外一种描述符与邻居原子相对于中心原子的方向也有关，称为角度描述符。下面我们一一介绍。

#### 径向描述符

对于某个中心原子 $i$，NEP 中的径向描述符分量定义为

$$
q^i _n = \sum _i g _n(r _{ij})
$$

也就是说，我们一共有若干（由下标 $n$ 标记）径向描述符分量。每一个径向描述符分量是某个径向函数 $g _n(r _{ij})$ 对近邻的求和。对比前面的 EAM 势，会发现一个径向描述符分量就类似于 EAM 势中的电子密度。之所以需要用若干径向描述符分量，是为了对近邻原子相对于中心原子的距离进行细致的刻画。

那么径向函数 $g _n(r _{ij})$ 应该有怎样的形式呢？在 [Behler 和 Parrinello](https://doi.org/10.1103/PhysRevLett.98.146401) 提出的方案中，径向函数是一些高斯分布函数，它们具有不同的分布中心和宽度。在 NEP 中，径向函数 $g _n(r _{ij})$ 是另一套径向基函数的线性叠加：

$$
   g _n(r _{ij}) = \sum _{k=0} ^{N _\mathrm{bas}^\mathrm{R}} c ^{ij} _{nk} f _k(r _{ij}),
$$

其中的径向基函数通过 Chebyshev 多项式定义：

$$
   f _k(r _{ij}) = \frac{1}{2}
   \left[T _k\left(2\left(r _{ij}/r _\mathrm{c}^\mathrm{R}-1\right)^2-1\right)+1\right]
   f _\mathrm{c}(r _{ij}),
$$

这里的 $T _k(x)$ 是 $k$-阶第一类 Chebyshev 多项式。上式中的 $f _\mathrm{c}(r _{ij})$ 是一个光滑截断函数，类似于 Tersoff 势的截断函数：

$$
   f _\mathrm{c}(r _{ij}) 
   = \begin{cases}
   \frac{1}{2}\left[
   1 + \cos\left( \pi \frac{r _{ij}}{r _\mathrm{c}^\mathrm{R}} \right) 
   \right],& r _{ij}\leq r _\mathrm{c}^\mathrm{R}; \\
   0, & r _{ij} > r _\mathrm{c}^\mathrm{R}.
   \end{cases}
$$

下图展示了前 5 个径向基函数，其中的 $f _0$ 实际上就是截断函数，因为 $T _0(x) = 1$。图中的横坐标是约化距离 $r _{ij}/r _\mathrm{c}^\mathrm{R}$，其中 $r _\mathrm{c}^\mathrm{R}$ 是径向描述符的截断距离。

![radial _basis _functions](fig/basis.png)

#### 角度描述符

NEP 中的角度描述符包含所谓的三体、四体和五体描述符。我们这里仅介绍三体描述符，关于NEP中四体和五体描述符的细节，请参考如下文章 [GPUMD: A package for constructing accurate machine-learned potentials and performing highly efficient atomistic simulations](https://doi.org/10.1063/5.0106617).

角度描述符需要有键角的信息，所以表达式中一定要有至少三个原子之间的键角。类似于径向描述符，我们可以构造出如下三体角度描述符：

$$
   q^i _{nl} = \sum _{j \neq i} \sum _{k \neq i} g _n(r _{ij}) g _n(r _{ij}) P _l(\theta _{ijk}),
$$

可以看到，三体角度描述符相比两体的径向描述符多了一个下标 $l$。这个下标是勒让德多项式 $P _l(\theta _{ijk})$ 的阶数。勒让德多项式是键角 $\theta _{ijk}$ 的函数。该键角是以原子 $i$ 为中心，以键 $ij$ 和 $ik$ 为两边的夹角。该式中的函数 $g _n(r _{ij})$ 和径向描述符中的对应函数一致，只不过可能有不同的截断距离 $r _\mathrm{c}^\mathrm{A}$。

上述表达式中的双重求和保证了同类原子置换的不变性，但导致计算量正比于近邻个数的平方。有一个办法可以在不改变结果的前提下降低计算复杂度，那就是利用球谐函数的加法定理。利用该定理（练习题），可以将上述三体角度描述符写成如下等价的形式：

$$
   q^i _{nl} = \sum _{m=-l}^l (-1)^m A^i _{nlm} A^i _{nl(-m)},
$$

$$
   A^i _{nlm} = \sum _{j\neq i} g _n(r _{ij}) Y _{lm}(\theta _{ij},\phi _{ij}),
$$

其中， $Y _{lm}(\theta _{ij},\phi _{ij})$ 是球谐函数。 $\theta _{ij}$ 是极角， $\phi _{ij}$ 是方位角。


### NEP机器学习势的训练

机器学习势的训练也称为拟合，目的是根据一些物理量的参考值优化NEP机器学习势模型中的参数，包括人工神经网络中的权重和偏置参数，以及描述符中的 $g _n$ 函数所涉及的一些线性叠加参数。一般来说，一个元素一共有好几千个参数，对多元素体系则有更多参数。

#### 损失函数

为了优化模型中的参数，我们需要定义一个损失函数，优化的方向是得到尽可能小的损失函数值。也就是说，大的损失函数值对应精度差的模型，而小的损失函数值对应精度高的模型。很自然地，损失函数应该定义为由NEP计算的某些物理量相对参考值的“距离”。通常用方均根误差（root-mean-square error) 表示。

那么，我们要考虑哪些物理量呢？我们知道，通过势函数，我们可以计算每个原子的能量、受力（矢量）以及位力（张量）。在NEP中我们就是使用这三个物理量进行训练。一般来说，参考值是通过量子力学密度泛函理论（DFT) 的方法计算的，其中能量和位力的一般来说只是针对整个构型来说的，而不是针对一个构型中的单个原子来说的。那么，我们需要针对 NEP 计算出整个构型的能量和位力，再和参考值对比。相反，力的对比可针对单个原子实现，因为 DFT 也可以计算单个原子的受力。

根据以上考虑，我们在NEP中提出了如下损失函数：

$$
L(\boldsymbol{z}) = L _{\rm e}(\boldsymbol{z}) + L _{\rm f}(\boldsymbol{z}) + L _{\rm v}(\boldsymbol{z}) + L _1(\boldsymbol{z}) + L _2(\boldsymbol{z})
$$

$$
L _{\rm e}(\boldsymbol{z}) 
= \lambda _\mathrm{e} 
\left( \frac{1}{N _\mathrm{str}}\sum _{n=1} ^{N _\mathrm{str}} \left( U^\mathrm{NEP}(n,\boldsymbol{z}) - U^\mathrm{tar}(n)\right)^2
\right) ^{1/2} 
$$

$$
L _{\rm f}(\boldsymbol{z}) = \lambda _\mathrm{f} \left( \frac{1}{3N} \sum _{i=1} ^{N} \left( \boldsymbol{F} _i^\mathrm{NEP}(\boldsymbol{z}) - \boldsymbol{F} _i^\mathrm{tar}\right)^2 \right) ^{1/2} 
$$

$$
L _{\rm v}(\boldsymbol{z}) = \lambda _\mathrm{v} \left( 
   \frac{1}{6N _\mathrm{str}}
   \sum _{n=1} ^{N _\mathrm{str}} \sum _{\mu\nu} \left( W _{\mu\nu}^\mathrm{NEP}(n,\boldsymbol{z}) - W _{\mu\nu}^\mathrm{tar}(n)\right)^2
   \right) ^{1/2} 
$$

$$
L _1(\boldsymbol{z}) = \lambda _1 \frac{1}{N _\mathrm{par}} \sum _{n=1} ^{N _\mathrm{par}} |z _n| 
$$

$$
L _2(\boldsymbol{z}) = \lambda _2 \left(\frac{1}{N _\mathrm{par}} \sum _{n=1} ^{N _\mathrm{par}} z _n^2\right) ^{1/2}.
$$

在公式中有如下符号：
- $\boldsymbol{z}$ 是所有被优化参数组成的抽象矢量。
- $N _\mathrm{str}$ 是训练集中结构的数目。
- $N$ 是训练集中原子的数目。
- $U^\mathrm{NEP}(n,\boldsymbol{z})$ 和 $W _{\mu\nu}^\mathrm{NEP}(n,\boldsymbol{z})$ 是通过 NEP 计算出的结构 $n$ 的能量和位力， $U^\mathrm{tar}(n)$ 和 $W _{\mu\nu}^\mathrm{tar}(n)$ 是对应的参考值。
- $\boldsymbol{F} _i^\mathrm{NEP}(\boldsymbol{z})$ 是通过 NEP 计算出的原子 $i$ 的力矢量， $\boldsymbol{F} _i^\mathrm{tar}$ 是对应的参考值。
- $L _1(\boldsymbol{z})$ 和 $L _1(\boldsymbol{z})$ 代表 $\mathcal{L} _1$ 和 $\mathcal{L} _2$ 正则化。

#### 可分离的自然演化策略

NEP 机器学习势得名于将演化算法用于训练神经网络势函数。我们所用的演化算法称为可分离的自然演化策略 (separable natural evolution strategy, 简称为 SNES)，由 Tom Schaul 等人在文章 [High dimensions and heavy tails for natural evolution strategies](https://doi.org/10.1145/2001576.2001692) 中提出。这是一种不需要使用任何解析梯度的优化算法，编程实现非常简单。在 GPUMD 中实现的 NEP 势函数没有使用任何第三方机器学习程序包，就是得益于该演化算法的简单性。下面列出使用 SNES 训练 NEP 势函数的算法流程：

- 初始化。

- 对演化代数进行迭代，直到完成特定代数 $N_{\rm gen}$
  - 产生一个种群
  - 计算种群中各个个体的损失函数
  - 计算自然梯度
  - 根据自然梯度更新解的平均值和方差
 
该算法的参数选取如下： 



### NEP机器学习势的编程实现

NEP 机器学习势已经在 [GPUMD 程序包](https://github.com/brucefan1983/GPUMD) 中实现，可用该程序包的 `nep` 可执行文件训练，并由该程序包的 `gpumd` 可执行文件进行分子动力学模拟。

NEP 机器学习势目前也有一个独立的 C++ 编程实现，见 [NEP_CPU 程序包](https://github.com/brucefan1983/NEP_CPU)。该程序包给出了一个名为`NEP3` 的 C++ 类，见程序包中 `src` 文件夹内的 `nep.cpp` 和`nep.h`。

## 习题

1. 推导 Tersoff 势 和 NEP 势的表达势  $\partial U_i / \partial \vec{r} _{ij}$ 。结果可分别参考 [2015年](https://doi.org/10.1103/PhysRevB.92.094301) 和 [2021年](https://doi.org/10.1103/PhysRevB.104.104309) 的文章。

2. 利用球谐函数的加法定理，将三体角度描述符

$$
   q^i _{nl} = \sum _{j \neq i} \sum _{k \neq i} g _n(r _{ij}) g _n(r _{ij}) P _l(\theta _{ijk}),
$$

推导为如下形式：

$$
   q^i _{nl} = \sum _{m=-l}^l (-1)^m A^i _{nlm} A^i _{nl(-m)},
$$

$$
   A^i _{nlm} = \sum _{j\neq i} g _n(r _{ij}) Y _{lm}(\theta _{ij},\phi _{ij}).
$$

