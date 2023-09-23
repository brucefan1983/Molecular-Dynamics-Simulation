# 《分子动力学模拟入门》第三章：一个简单的分子动力学模拟程序

# Table of contents
- [简单分子动力学模拟的基本要素](#简单分子动力学模拟的基本要素)
 	- [分子动力学模拟的定义](#分子动力学模拟的定义)
 	- [初始条件](#初始条件)
 	- [边界条件](#边界条件)
 	- [相互作用](#相互作用)
 	- [运动方程的数值积分](#运动方程的数值积分)
- [开发一个简单的分子动力学模拟程序](#开发一个简单的分子动力学模拟程序)
  - [程序中使用的单位制](#程序中使用的单位制)
  - [本章程序的源代码解析](#本章程序的源代码解析)
    - [主函数](#主函数)
    - [内存分配](#内存分配)
    - [输入处理](#输入处理)
    - [速度初始化](#速度初始化)
    - [运动方程的数值积分](#运动方程的数值积分)
    - [求势能与力](#求势能与力)
    - [程序的编译与运行](#程序的编译与运行)
  - [能量守恒的测试](#能量守恒的测试)
- [三斜盒子](#三斜盒子)
 	- [三斜盒子的定义](#三斜盒子的定义)
 	- [三斜盒子情况下的周期边界条件](#三斜盒子情况下的周期边界条件)
- [近邻列表的创建](#近邻列表的创建)
  - [为什么要用近邻列表](#为什么要用近邻列表)
  - [自动判断何时更新近邻列表](#自动判断何时更新近邻列表)
  - [构建近邻列表的平方标度算法](#构建近邻列表的平方标度算法)
  - [构建近邻列表的线性标度算法](#构建近邻列表的线性标度算法)
  - [程序速度测试](#程序速度测试)
- [GPUMD程序简介](#GPUMD程序简介)

本章将从一个简单的分子动力学模拟程序开始，带领读者走进分子动力学模拟的世界。在以后的章节，我们将逐步深入地探讨分子动力学模拟的若干重要课题。

## 简单分子动力学模拟的基本要素

### 分子动力学模拟的定义

作者曾经在一篇博文给分子动力学模拟下过一个定义：

分子动力学模拟是一种数值计算方法，在这种方法中，我们对一个具有一定初始条件和边界条件且具有相互作用的多粒子系统的运动方程进行数值积分，得到系统在相空间中的一条离散的轨迹，并用统计力学的方法从这条相轨迹中提取出有用的物理结果。

在本章余下的部分，我们会一一考察上述定义中的重要概念，如初始条件、边界条件、相互作用、运动方程、数值积分等。这里，我们首先讨论上述定义中的“统计力学的方法”。

在一个特定的平衡态统计系综中，一个物理量 $A$ 的统计平均值可表达为

$$
\langle A \rangle _{\rm ensemble} = \int f(p,q) A(p,q) dpdq
$$

其中， $f(p,q)$ 是该系综的分布函数， $q$ 和 $p$ 代表广义坐标和广义动量的集合。这样的统计平均称为系综平均。然而，根据我们的定义，在分子动力学模拟中并没有使用“系综”（即系统的集合），而是仅对一个系统的时间演化过程进行分析。那么，上述定义中的“统计力学的方法”指的是什么呢？

这里的“统计力学的方法”指的是时间平均，即

$$
\langle A \rangle _{\rm time} = \lim _{t\to\infty}\frac{1}{t}\int_0^t A(t') dt'
$$

我们知道，随着时间的推移，系统将历经一条相轨迹。因为这条相轨迹永不与自身相交，如果它不代表一个周期性运动的话，那么随着时间的增加，这条相轨迹应该历经越来越多的相点。一个自然的假设是，当时间趋近于无穷大时，这条相轨迹将遍历系统的所有相点。这就是各态历经假设。在此假设下，系综平均与时间平均等价。

从实验的角度来说，对热力学体系的实验测量采用了时间平均而非系综平均。所以，分子动力学模拟实际上比基于系综的统计力学更加贴近实验。从计算机模拟的角度来看，除了采用时间平均的分子动力学模拟，还有直接采用系综平均的蒙特卡罗模拟，但本书不讨论蒙特卡罗模拟。

根据上述定义，我们可以设想一个典型的、简单的分子动力学模拟有如下大致的计算流程：
- 初始化。设置系统的初始条件，具体包括各个粒子的位置矢量和速度矢量。
- 时间演化。根据系统中的粒子所满足的相互作用规律，由牛顿定律确定所有粒子的运动方程（二阶常微分方程组），并对运动方程进行数值积分，即不断地更新每个粒子的坐标和速度。最终, 我们将得到一系列离散的时刻系统在相空间中的位置，即一条离散的相轨迹。
- 测量。用统计力学的方法分析相轨迹所蕴含的物理规律。

###初始条件

初始化指的是确定一个初始的相空间点，包括各个粒子初始的坐标和速度。在分子动力学模拟中，我们需要对 $3N$（ $N$ 是粒子数目）个二阶常微分方程进行数值积分。因为每一个二阶常微分方程的求解都需要有两个初始条件，所以我们需要确定 $6N$ 个初始条件： $3N$ 个初始坐标分量和 $3N$ 个初始速度分量。

#### 坐标初始化

坐标的初始化指的是为系统中的每个粒子选定一个初始的位置坐标。分子动力学模拟中如何初始化位置主要取决于所要模拟的体系。例如，如果要模拟固态氩，就得让各个氩原子的位置按面心立方结构排列。如果要模拟的是液态或者气态物质，那么初始坐标的选取就可以比较随意了。重要的是，在构造的初始结构中，任何两个粒子的距离都不能太小，因为这可能导致有些粒子受到非常大的力，以至于让后面的数值积分变得非常不稳定。坐标的初始化也常被称为建模，往往需要用到一些专业的知识，例如固体物理学（晶体学）中的知识。本章将通过一个程序介绍固态氩的建模。

#### 速度初始化

我们知道，任何经典热力学系统在平衡时各个粒子的速度满足麦克斯韦分布。然而，作为初始条件，我们并不一定要求粒子的速度满足麦克斯韦分布。最简单的速度初始化方法是产生 $3N$ 个在某个区间均匀分布的随机速度分量，再通过如下几个基本条件对速度分量进行修正。

第一个条件是让系统的总动量为零。也就是说，我们不希望系统的质心在模拟的过程中跑动。分子间作用力是所谓的内力，不会改变系统的整体动量，即系统的整体动量是守恒的。只要初始的整体动量为零，在分子动力学模拟的时间演化过程中整体动量将保持为零。如果整体动量明显偏离零（相对于所用浮点数精度来说），则说明模拟出了问题。这正是判断程序是否有误的标准之一。

第二个条件是系统的总动能应该与所选定的初始温度对应。我们知道，在经典统计力学中，能量均分定理成立，即粒子的哈密顿量中每一个具有平方形式的能量项的统计平均值都等于 $k _{\rm B} T/2$。其中， $k _{\rm B}$ 是玻尔兹曼常数， $T$ 是系统的绝对温度。所以，在将质心的动量取为零之后就可以对每个粒子的速度进行一个标度变换，使得系统的初始温度与所设定的温度一致。假设我们设置的目标温度是 $T_0$，那么对各个粒子的速度做如下变换即可让系统的温度从 $T$ 变成 $T_0$：

$$
\vec{v}_i \rightarrow \vec{v}_i'= \vec{v}_i\sqrt{\frac{T_0}{T}}    
$$

容易验证（习题），在做上式中的变换之前，如果系统的总动量已经为零，那么在做这个变换之后，系统的总动量也为零。

第三个可选的条件是角动量的初始化。我们将在一个习题中研究这个问题。

### 边界条件

在我们对分子动力学模拟的定义中，除了初始条件，还提到了边界条件。边界条件对常微分方程的求解并不是必要的，但在分子动力学模拟中通常会根据所模拟的物理体系选取合适的边界条件，以期得到更合理的结果。边界条件的选取对粒子间作用力的计算也是有影响的。常用的边界条件有好几种，但我们这里只先讨论其中的一种：周期边界条件，也称为 Born-von Karman 边界条件。在计算机模拟中，模拟的系统尺寸一定是有限的，通常比实验中对应的体系的尺寸小很多。选取周期边界条件通常可以让模拟的体系更加接近于实际的情形，因为原本有边界的系统在应用了周期边界条件之后，“似乎”没有边界了。当然，并不能说应用了周期边界条件的系统就等价于无限大的系统，只能说周期边界条件的应用可以部分地消除边界效应，让所模拟系统的性质更加接近于无限大系统的性质。通常，在这种情况下，我们要模拟几个不同大小的系统，分析所得结果对模拟尺寸的依赖关系。

在计算两个粒子，如粒子 $i$ 和粒子 $j$ 的距离时，就要考虑周期边界条件带来的影响。举个一维的例子。假设模拟在一个长度为 $L_x$ 的模拟盒子中进行，采用周期边界条件时，可以将该一维的盒子想象为一个圆圈。假设 $L_x=10$（任意单位），第 $i$ 个粒子的坐标 $x_i=1$，第 $j$ 个粒子的坐标 $x_j=8$，则这两个粒子的距离是多少呢？如果忽略周期边界条件，那么答案是 $|x_j-x_i|=7$，而且 $j$ 粒子在 $i$ 粒子的右边（坐标值大的一边）。但是，在采取周期边界条件时，也可认为 $j$ 粒子在 $i$ 粒子的左边，且坐标值可以平移至 $8-10=-2$。这样， $j$ 与 $i$ 的距离是 $|x_j-x_i|=3$，比平移 $j$ 粒子之前两个粒子之间的距离要小。在我们的模拟中，总是采用最小镜像约定 [Equation of state calculations by fast computing machines](https://doi.org/10.1063/1.1699114)：在计算两个粒子的距离时，总是取最小的可能值。定义

$$
x_j-x_i \equiv x _{ij}    
$$

则这个约定等价于如下规则：如果 $x _{ij}<-L_x/2$，则将 $x _{ij}$ 换为 $x _{ij}+L_x$；如果 $x _{ij}>+L_x/2$，则将 $x _{ij}$换为 $x _{ij}-L_x$。最终效果就是让变换后的 $x _{ij}$ 的绝对值不大于 $L_x/2$。

很容易将上述讨论推广到二维和三维的情形。例如，在二维的情形中，可以将一个周期的模拟盒子想象为一个环面，就像一个甜甜圈或一个充了气的轮胎的表面。在三维的情形中，可以将一个周期的模拟盒子想象为一个三维环面，而最小镜像约定可以表达为：
- 如果 $x _{ij}<-L_x/2$，则将 $x _{ij}$ 换为 $x _{ij}+L_x$；如果 $x _{ij}>+L_x/2$，则将 $x _{ij}$ 换为 $x _{ij}-L_x$。
- 如果 $y _{ij}<-L_y/2$，则将 $y _{ij}$ 换为 $y _{ij}+L_y$；如果 $y _{ij}>+L_y/2$，则将 $y _{ij}$ 换为 $y _{ij}-L_y$。
- 如果 $z _{ij}<-L_z/2$，则将 $z _{ij}$ 换为 $z _{ij}+L_z$；如果 $z _{ij}>+L_z/2$，则将 $z _{ij}$ 换为 $z _{ij}-L_z$

这里，我们假设了三维模拟盒子中 3 个共点的边的长度分别为 $L_x$、 $L_y$ 和 $L_z$，且两两相互垂直（所谓的正交模拟盒子）。如果有任意两个共点的边不是相互垂直的，情况就要复杂一些。本章仅讨论正交盒子的情形，以后再讨论非正交盒子的情形。

### 相互作用

宏观物质的性质在很大程度上是由微观粒子之间的相互作用力决定的。所以，对粒子间相互作用力的计算在分子动力学模拟中是至关重要的。粒子间有何种相互作用不是分子动力学模拟本身所能回答的；它本质上是一个量子力学的问题。在经典分子动力学模拟中，粒子间的相互作用力常常由一个或多个经验势函数描述。经验势函数能够在某种程度上反映出某些物质的某些性质。近年来，机器学习也广泛地用于构造更加准确的势函数。在本章，我们只介绍一个称为 Lennard-Jones 势的简单势函数（简称为 LJ 势）[On the determination
of molecular fields. I. From the variation of the viscosity of a gas with temperature](https://doi.org/10.1098/rspa.1924.0081) and [On the determination of molecular fields. —II. From the equation of state of a gas](https://doi.org/10.1098/rspa.1924.0082)。在以后的章节中，我们将介绍更多的经验势函数以及机器学习势函数。

考虑系统中的任意粒子对 $i$ 和 $j$，它们之间的相互作用势能可以写为

$$
U _{ij}(r _{ij})=4\epsilon
\left(
\frac{\sigma^{12}}{r _{ij}^{12}}-\frac{\sigma^{6}}{r _{ij}^{6}}
\right).    
$$

其中， $\epsilon$ 和 $\sigma$ 是势函数中的参数，分别具有能量和长度的量纲； $r _{ij}=|\vec{r}_j-\vec{r}_i|$ 是两个粒子间的距离。

LJ 势比较适合描述惰性元素组成的物质。它是最早提出的两体势函数之一。所谓两体势，指的是两个粒子 $i$ 和 $j$ 之间的相互作用势仅依赖于它们之间的距离 $r _{ij}$，不依赖于系统中其他粒子的存在与否及具体位置。本章只讨论两体势，后续的章节会讨论一些多体势，即非两体势。对于两体势函数，我们可以将整个系统的总势能 $U$ 写为

$$
U=\sum _{i=1}^N U_i;
$$
$$
U_i= \frac{1}{2} \sum _{j \neq i} U _{ij}(r _{ij}).
$$

将以上两式合起来，可以写成

$$
U=\frac{1}{2}\sum _{i=1}^N  \sum _{j \neq i} U _{ij}(r _{ij})
$$

上面的 $U_i$ 可以称为粒子 $i$ 的势能。也可以将总势能写为如下形式：

$$
U=\sum _{i=1}^N \sum _{j > i} U _{ij}(r _{ij})
$$

根据力的定义可得（习题），粒子 $i$ 所受总的力为

$$
\vec{F} _{i} = \sum _{j \neq i} \vec{F} _{ij}
$$

$$
\vec{F} _{ij} =
\frac{\partial U _{ij}(r _{ij})}{\partial r _{ij}}
\frac{\vec{r} _{ij} }{r _{ij}}
$$

其中，我们定义了一个表示粒子间相对位置的符号

$$
\vec{r} _{ij} \equiv \vec{r}_j - \vec{r}_i
$$

显然，牛顿第三定律成立：

$$
\vec{F} _{ij} = - \vec{F} _{ji}.
$$

进一步推导可得 LJ 势中力的明确表达式：

$$
\vec{F} _{ij} = \vec{r} _{ij}
\left( 
\frac{24 \epsilon \sigma^6} {r _{ij}^8} - \frac{48 \epsilon \sigma^{12}} {r _{ij}^{14}} 
\right).
$$

通常，为了节约计算，我们会对势函数进行一个截断，即认为当两个原子之间的距离大于某个截断距离 $R _{\rm c}$ 时，它们之间的相互作用势能和力都是零：

$$
U _{ij}(r _{ij}) = 0  \quad (r _{ij} > R _{\rm c})
$$

$$
\vec{F} _{ij} = \vec{0}  \quad (r _{ij} > R _{\rm c})
$$

### 运动方程的数值积分

我们知道在经典力学中，粒子的运动方程可以用牛顿第二定律表达。例如，对于第 $i$ 个粒子，其运动方程为

$$
m_i \frac{d^2\vec{r}_i}{dt^2} = \vec{F}_i
$$

这是一个二阶常微分方程，我们可以把它改写为两个一阶常微分方程：

$$
\frac{d\vec{r}_i}{dt} = \vec{v}_i
$$

$$
\frac{d\vec{v}_i}{dt} = \frac{\vec{F}_i}{m_i}
$$

对运动方程进行数值积分的目的就是在给定的初始条件下找到各个粒子在一系列离散的时间点的坐标和速度值。我们假设每两个离散的时间点之间的间隔是固定的，记为 $\Delta t$，称为时间步长。在分子动力学模拟中使用的数值积分方法有很多种，本章只介绍所谓的“速度-Verlet”积分方法，其推导过程见习题。在该方法中，第 $i$ 个粒子在时刻 $t+\Delta t$ 的速度 $\vec{v}_i(t+\Delta t)$ 和位置 $\vec{r}_i(t+\Delta t)$ 分别由以下两式给出：

$$
\vec{v}_i(t+\Delta t)=\vec{v}_i(t)+\frac{1}{2}\frac{\vec{F}_i(t)+\vec{F}_i(t+\Delta t)}{m_i}\Delta t
$$

$$
\vec{r}_i(t+\Delta t)
=\vec{r}_i(t)
+\vec{v}_i(t)\Delta t
+\frac{1}{2}\frac{\vec{F}_i(t)}{m_i}(\Delta t)^2
$$

由以上两式可以看出， $t+\Delta t$ 时刻的坐标仅依赖于 $t$ 时刻的坐标、速度和力，但 $t+\Delta t$ 时刻的速度依赖于 $t$ 时刻的速度、力及 $t+\Delta t$ 时刻的力。所以，从算法的角度来说，以上两式应该对应如下的计算流程：

第一步：部分地更新速度并完全地更新坐标（注意，我们引入了一个中间的速度变量 $\vec{v}_i(t+\Delta t/2)$）：

$$
\vec{v}_i(t) \rightarrow \vec{v}_i(t+\Delta t/2)=\vec{v}_i(t)+\frac{1}{2}\frac{\vec{F}_i(t)}{m_i}\Delta t
$$
$$
\vec{r}_i(t)\rightarrow \vec{r}_i(t+\Delta t)
=\vec{r}_i(t)
+\vec{v}_i(t+\Delta t/2)\Delta t
$$

第二步：用更新后的坐标计算新的力

$$
\vec{F}_i(t)\rightarrow \vec{F}_i(t+\Delta t)
$$

第三步：用更新后的力完成速度的更新：

$$
\vec{v}_i(t+\Delta t/2) \rightarrow \vec{v}_i(t+\Delta t)=\vec{v}_i(t+\Delta t/2)+\frac{1}{2}\frac{\vec{F}_i(t+\Delta t)}{m_i}\Delta t
$$

完成上述计算之后，粒子的坐标、速度、和力都从 $t$ 时刻的更新为 $t+\Delta t$ 时刻的。这就是一个时间步的计算。反复执行这样的计算流程，系统的微观状态就会不断地随时间变化，从而得到一条相空间的轨迹。系统所有的宏观性质都包含在相轨迹中。

## 开发一个简单的分子动力学模拟程序

本节给出一个简单的用 C++ 开发的分子动力学模拟程序。我们选择用 C++ 语言开发程序，因为这是作者最熟悉的编程语言。因为分子动力学模拟一般较为耗时，所以用高效的编译型语言开发较为合适。除了 C++ 语言，C 语言和 Fortran 语言也很高效。以后我们还将介绍分子动力学模的 CUDA 编程实现。另外，在完成分子动力学模拟的计算后，我们往往需要对得到的数据进行可视化，此时使用一门解释性语言更为方便。作者使用 Matlab 语言进行数据的后处理和可视化。

### 程序中使用的单位制

我们的分子动力学模拟程序只涉及经典力学和热力学，故只需要用到 4 个基本物理量的单位。我们选择如下 4 个基本单位来确定各个物理量的数值：
- 能量：电子伏特（记号为 eV），约为 $1.6\times 10^{-19}$ J。
- 长度：埃（angstrom，记号为 A），即 $10^{-10}$ m。
- 质量：原子质量单位（atomic mass unit，记号为 amu），约为 $1.66 \times 10^{-27}$ kg。
- 温度：开尔文（记号为 K）。

用这样的基本单位，可使程序中大部分物理量的数值都接近 1。我们称这样的单位为该程序的“自然单位”。

从以上基本单位可以推导出程序中其他相关物理量的单位：
- 力。因为力乘以距离等于功（能量），故力的单位是能量单位除以长度单位，即 eV/A。
- 速度。因为动能正比于质量乘以速度的平方，故速度的单位是能量单位除以质量单位再开根号，即 eV $^{1/2}$ amu $^{-1/2}$。
- 时间。因为长度等于速度乘以时间，故时间的单位是长度单位除以速度单位，即A amu $^{1/2}$ eV $^{-1/2}$，约为 $1.018051 \times 10^{1}$ fs（fs 指飞秒，即 $10^{-15}$ s）。
- 玻尔兹曼常数 $k _{\rm B}$。这是一个很重要的常数，它在国际单位制中约为 $1.38\times 10^{-23}$ J/K，对应于程序自然单位制的 $8.617343 \times 10^{-5}$ eV/K。

### 本章程序的源代码解析

本章的程序很简单，一共只有 200 行代码，故将所有代码写在一个源文件 `ljmd.cpp`。下面，我们详细地讲解该程序。

#### 主函数

我们从该文件的主函数 `main()` 看起。下面是主函数的全部代码：

```C++
int main(int argc, char** argv)
{
  if (argc != 5) {
    printf("usage: %s numCells numSteps temperature timeStep\n", argv[0]);
    exit(1);
  }
  const int numCells = atoi(argv[1]);
  const int numSteps = atoi(argv[2]);
  const double temperature = atof(argv[3]);
  double timeStep = atof(argv[4]);
  timeStep /= TIME_UNIT_CONVERSION; // from fs to natural unit

  Atom atom;
  allocateMemory(numCells, atom);
  initializePosition(numCells, atom);
  initializeVelocity(temperature, atom);

  const clock_t tStart = clock();
  std::ofstream ofile("energy.txt");
  for (int step = 0; step < numSteps; ++step) {
    integrate(true, timeStep, atom);  // step 1 in the book
    findForce(atom);                  // step 2 in the book
    integrate(false, timeStep, atom); // step 3 in the book
    if (step % Ns == 0) {
      ofile << findKineticEnergy(atom) << " "
            << std::accumulate(atom.pe.begin(), atom.pe.end(), 0.0)
            << std::endl;
    }
  }
  ofile.close();
  const clock_t tStop = clock();
  const float tElapsed = float(tStop - tStart) / CLOCKS_PER_SEC;
  std::cout << "Time used = " << tElapsed << " s" << std::endl;

  return 0;
}
```

首先，程序从命令行读入四个参数，分别是面心立方 固态氩单胞在每个方向的个数（`numCells`）、整个分子动力学模拟的步数（`numSteps`）、体系的目标温度（`temperature`）和数值积分的时间步长（`timeStep`）。在读入时间步长后，立刻将其单位从输入的 fs 转换至程序的自然单位。这种单位转换只需要在处理输入和输出时实施，在程序的其他地方任何物理量的单位都将是我们定义的自然单位。

接着，程序定义了一个结构体 `Atom` 的变量 `atom`。该结构体类型中定义了程序中用到的大部分数据。我们这里用了 C++ 标准模板库中的 `std::vector` 来表示一些数组。

```C++
struct Atom {
  int number; // 总的粒子数（原子数）
  double box[6]; // 前三个数是三个盒子长度；后三个数是盒子长度的一半
  // 质量、坐标、速度、力、势能：
  std::vector<double> mass, x, y, z, vx, vy, vz, fx, fy, fz, pe;
};
```


在定义 `atom` 后，我们调用函数 `allocateMemory()` 对一些数组分配内存，然后调用函数 `initializePosition()` 初始化模拟体系的原子坐标，并调用函数 `initializeVelocity()` 初始化体系的速度。

接下来，做一个次数为 `numSteps` 的循环进行时间演化。在该循环中，我们实现前面讲述的速度-Verlet积分算法，一共三个步骤：
- 语句 \verb"integrate(true, timeStep, atom);" 实现速度的部分更新和坐标的完全更新。
- 语句 \verb"findForce(atom);" 实现力的更新。
- 语句 \verb"integrate(false, timeStep, atom);" 实现剩下的速度更新

在循环过程中，每 `Ns=100` 步计算一次体系的总动能和总势能并输出到文件 `energy.txt`。程序将对演化过程计时，在结束程序之前报道演化过程所花的总时间。

#### 内存分配

本书代码基本上都用 C++ 的 `std::vector` 容器来处理动态的内存分配与释放，而不是用 `malloc` 和 `free` 或者 `new` 和 `delete` 来处理。在函数 `allocateMemory()` 中，我们用 `std::vector` 的成员函数 `resize()` 分配内存，并同时将每个数组元素初始化。除了原子质量初始化为 40 amu 之外，其它物理量都初始化为零。

```C++
void allocateMemory(const int numCells, Atom& atom)
{
  const int numAtomsPerCell = 4;
  atom.number = numCells * numCells * numCells * numAtomsPerCell;
  atom.mass.resize(atom.number, 40.0); // argon mass
  atom.x.resize(atom.number, 0.0);
  atom.y.resize(atom.number, 0.0);
  atom.z.resize(atom.number, 0.0);
  atom.vx.resize(atom.number, 0.0);
  atom.vy.resize(atom.number, 0.0);
  atom.vz.resize(atom.number, 0.0);
  atom.fx.resize(atom.number, 0.0);
  atom.fy.resize(atom.number, 0.0);
  atom.fz.resize(atom.number, 0.0);
  atom.pe.resize(atom.number, 0.0);
}
```

#### 输入处理

函数 `initializePosition()` 根据每个方向扩包的次数 `numCells` 对体系的坐标进行初始化。同时，与盒子大小有关的数据 `box` 也在该函数中得到初始化。

```C++
void initializePosition(const int numCells, Atom& atom)
{
  const int numAtomsPerCell = 4;
  const double latticeConstant = 5.385;
  atom.box[0] = atom.box[1] = atom.box[2] = latticeConstant * numCells;
  atom.box[3] = atom.box[0] * 0.5;
  atom.box[4] = atom.box[1] * 0.5;
  atom.box[5] = atom.box[2] * 0.5;
  const double x0[numAtomsPerCell] = {0.0, 0.0, 0.5, 0.5};
  const double y0[numAtomsPerCell] = {0.0, 0.5, 0.0, 0.5};
  const double z0[numAtomsPerCell] = {0.0, 0.5, 0.5, 0.0};
  int n = 0;
  for (int ix = 0; ix < numCells; ++ix) {
    for (int iy = 0; iy < numCells; ++iy) {
      for (int iz = 0; iz < numCells; ++iz) {
        for (int i = 0; i < numAtomsPerCell; ++i) {
          atom.x[n] = (ix + x0[i]) * latticeConstant;
          atom.y[n] = (iy + y0[i]) * latticeConstant;
          atom.z[n] = (iz + z0[i]) * latticeConstant;
          ++n;
        }
      }
    }
  }
}
```

#### 速度初始化

下面是速度初始化的函数 `initializeVelocity()`。在该函数中，首先利用随机数获得分布在 -1 到 1 之间的随机速度分量，同时计算体系的质心速度 `centerOfMassVelocity`。注意，函数 `rand()` 返回一个从零到 `RAND_MAX` 之间的整数。 接着，对速度进行修正，使得体系的整体动量为零。最后，调用 `scaleVelocity()` 函数对速度进行标度变换，使得体系的温度达到目标值 `T0`。

```C++
void initializeVelocity(const double T0, Atom& atom)
{
#ifndef DEBUG
  srand(time(NULL));
#endif
  double centerOfMassVelocity[3] = {0.0, 0.0, 0.0};
  double totalMass = 0.0;
  for (int n = 0; n < atom.number; ++n) {
    totalMass += atom.mass[n];
    atom.vx[n] = -1.0 + (rand() * 2.0) / RAND_MAX;
    atom.vy[n] = -1.0 + (rand() * 2.0) / RAND_MAX;
    atom.vz[n] = -1.0 + (rand() * 2.0) / RAND_MAX;
    centerOfMassVelocity[0] += atom.mass[n] * atom.vx[n];
    centerOfMassVelocity[1] += atom.mass[n] * atom.vy[n];
    centerOfMassVelocity[2] += atom.mass[n] * atom.vz[n];
  }
  centerOfMassVelocity[0] /= totalMass;
  centerOfMassVelocity[1] /= totalMass;
  centerOfMassVelocity[2] /= totalMass;
  for (int n = 0; n < atom.number; ++n) {
    atom.vx[n] -= centerOfMassVelocity[0];
    atom.vy[n] -= centerOfMassVelocity[1];
    atom.vz[n] -= centerOfMassVelocity[2];
  }
  scaleVelocity(T0, atom);
}
```

下面是函数 `scaleVelocity()` 的定义。在该函数中，首先调用函数 `findKineticEnergy()` 根据当前的速度计算当前的动能，进而得到当前的温度 `temperature`，然后根据公式计算速度标度变换的因子 `scaleFactor`，对速度进行标度变换。

```C++
void scaleVelocity(const double T0, Atom& atom)
{
  const double temperature =
    findKineticEnergy(atom) * 2.0 / (3.0 * K_B * atom.number);
  double scaleFactor = sqrt(T0 / temperature);
  for (int n = 0; n < atom.number; ++n) {
    atom.vx[n] *= scaleFactor;
    atom.vy[n] *= scaleFactor;
    atom.vz[n] *= scaleFactor;
  }
}
```

下面是通过速度计算动能的函数 `findKineticEnergy()`。注意，该函数计算的结果将通过一个值返回。

```C++
double findKineticEnergy(const Atom& atom)
{
  double kineticEnergy = 0.0;
  for (int n = 0; n < atom.number; ++n) {
    double v2 = atom.vx[n] * atom.vx[n] + atom.vy[n] * atom.vy[n] +
                atom.vz[n] * atom.vz[n];
    kineticEnergy += atom.mass[n] * v2;
  }
  return kineticEnergy * 0.5;
}
```

#### 运动方程的数值积分

我们用函数 `integrate()` 来实现速度-Verlet 积分算法的两个步骤。该函数的第一个输入参数 `isStepOne` 是一个布尔型变量。当该变量为真时，就实行速度-Verlet 积分算法的第一个步骤，即部分地将速度更新，并完全地将坐标更新。当该变量为假时，就实行速度-Verlet 积分算法的第二个步骤，只更新速度，不再更新坐标。

```C++
void integrate(const bool isStepOne, const double timeStep, Atom& atom)
{
  const double timeStepHalf = timeStep * 0.5;
  for (int n = 0; n < atom.number; ++n) {
    const double mass_inv = 1.0 / atom.mass[n];
    const double ax = atom.fx[n] * mass_inv;
    const double ay = atom.fy[n] * mass_inv;
    const double az = atom.fz[n] * mass_inv;
    atom.vx[n] += ax * timeStepHalf;
    atom.vy[n] += ay * timeStepHalf;
    atom.vz[n] += az * timeStepHalf;
    if (isStepOne) {
      atom.x[n] += atom.vx[n] * timeStep;
      atom.y[n] += atom.vy[n] * timeStep;
      atom.z[n] += atom.vz[n] * timeStep;
    }
  }
}
```

#### 求势能与力

函数 `findForce()` 负责求体系中各个粒子的势能和受到的力。在该函数的开头，我们定义了若干常量。这种常量的计算将在编译期间就完成。在循环之前尽可能多地计算常量可以省去很多不必要的计算。我们用了固态氩的 LJ 参数 $\epsilon = 0.01032e$ eV， $\sigma = 3.405$ A。并将截断距离取为 $R_c = 10$ A。

接着，我们将每个原子的势能和受力初始化为零，因为在后面的循环中我们将对势能和力进行累加。

接下来，是一个两重循环，因为们要计算每一对粒子之间的相互作用力。不过这里的两重循环有些特殊，排除了 `i >= j` 的可能性，这就是利用牛顿第三定律节约一半的计算量。

在循环体中，首先计算相对位置 $\vec{r} _{ij}$  并对其实施最小镜像约定。紧接着计算两个粒子距离的平方并忽略截断距离之外的粒子对。最后，通过非常节约的方式计算两个粒子之间的势能和相互作用力，并存储在对应的数组中。这部分的计算要避免使用耗时的 `sqrt()` 函数 和 `pow()` 函数。在 LJ 势的编程中，虽然从公式来看好像需要使用，但是仔细思考后会发现这些都是可以避免的。还有一点值得注意，那就是除法运算大概是乘法运算的几倍耗时，所以在编写程序时，要将除法运算的个数最小化。

```C++
void findForce(Atom& atom)
{
  const double epsilon = 1.032e-2;
  const double sigma = 3.405;
  const double cutoff = 10.0;
  const double cutoffSquare = cutoff * cutoff;
  const double sigma3 = sigma * sigma * sigma;
  const double sigma6 = sigma3 * sigma3;
  const double sigma12 = sigma6 * sigma6;
  const double e24s6 = 24.0 * epsilon * sigma6;
  const double e48s12 = 48.0 * epsilon * sigma12;
  const double e4s6 = 4.0 * epsilon * sigma6;
  const double e4s12 = 4.0 * epsilon * sigma12;
  for (int n = 0; n < atom.number; ++n)
    atom.fx[n] = atom.fy[n] = atom.fz[n] = atom.pe[n] = 0.0;

  for (int i = 0; i < atom.number - 1; ++i) {
    for (int j = i + 1; j < atom.number; ++j) {
      double xij = atom.x[j] - atom.x[i];
      double yij = atom.y[j] - atom.y[i];
      double zij = atom.z[j] - atom.z[i];
      applyMic(atom.box, xij, yij, zij);
      const double r2 = xij * xij + yij * yij + zij * zij;
      if (r2 > cutoffSquare)
        continue;

      const double r2inv = 1.0 / r2;
      const double r4inv = r2inv * r2inv;
      const double r6inv = r2inv * r4inv;
      const double r8inv = r4inv * r4inv;
      const double r12inv = r4inv * r8inv;
      const double r14inv = r6inv * r8inv;
      const double f_ij = e24s6 * r8inv - e48s12 * r14inv;
      atom.pe[i] += e4s12 * r12inv - e4s6 * r6inv;
      atom.fx[i] += f_ij * xij;
      atom.fx[j] -= f_ij * xij;
      atom.fy[i] += f_ij * yij;
      atom.fy[j] -= f_ij * yij;
      atom.fz[i] += f_ij * zij;
      atom.fz[j] -= f_ij * zij;
    }
  }
}
```

最小镜像约定的实施由如下两个函数实现。注意，这里的函数参数用了 C++ 里面的引用（Reference）。

```C++
void applyMicOne(const double length, const double halfLength, double& x12)
{
  if (x12 < -halfLength)
    x12 += length;
  else if (x12 > +halfLength)
    x12 -= length;
}

void applyMic(const double box[6], double& x12, double& y12, double& z12)
{
  applyMicOne(box[0], box[3], x12);
  applyMicOne(box[1], box[4], y12);
  applyMicOne(box[2], box[5], z12);
}
```

#### 程序的编译与运行

本书所开发的 C++ 程序都可以在 Linux 和 Windows 操作系统使用。我们推荐使用 GCC 工具。在命令行可以用如下方式编译本章的程序：
```shell
    $ g++ -O3 ljmd.cpp -o ljmd
```
其中，`-O3` 选项表示优化等级。编译完成后，将生成名为 `ljmd` 的可执行文件（在 Windows 中为 `ljmd.exe`)。

然后，就可以在命令行使用该程序：
```shell
    $ ljmd  numCells numSteps temperature timeStep
```
如果使用时忘了给命令行参数，程序会提示正确的用法。

具体地，笔者用如下命令运行程序：
```shell
    $ ljmd 4 20000 60 5
```
也就是说，考虑原子数为 $4^3 \times 4 = 256$ 的固态氩体系，运行 $2\times 10^4$ 步，初始温度为 60 K，积分步长为 5 fs。该模拟在笔者的计算机中运行的时间约为 10 秒钟。

### 能量守恒的测试

程序每隔100步输出系统的总动能 $K(t)$ 和总势能 $U(t)$，它们都是时间 $t$ 的函数。对于大小有限的体系，它们都是随时间 $t$ 涨落的。然而，根据能量守恒定律，系统动能和势能的和，即总能量 $E(t)=K(t)+V(t)$，应该是不随时间变化的。当然，我们的模拟中使用了具有一定误差的数值积分方法，故总能量也会有一定大小的涨落。这个涨落主要与积分的时间步长有关系。一般来说，积分的时间步长越大，总能量的涨落越大。

![energy](fig/energy.png)

   （a) 体系总动能随时间的变化；（b) 体系总势能随时间的变化；（c) 体系总能随时间的变化；（d) 相对总能量随时间的变化。


上图（a-c）给出了系统的总动能、总势能和总能量随时间变化的情况。可以看出动能是正的，势能是负的，涨落相对较大；总能是负的，涨落较小。细心的读者可以注意到，体系的总动能在很短的时间内突然降低了大约一半。这是因为，我们的模拟体系一开始是完美的面心立方晶格，每个原子都处于受力为零的平衡状态，没有振动产生的额外势能。我们用某个温度（在我们的例子中是 60 K）初始化了原子的速度，故该体系一开始是有一定的动能的。假设每个原子都会在其平衡位置做简谐振动（这就是所谓的简谐近似，它对很多问题的研究来说是一个很好的出发点），故根据能量均分定理，在达到热力学平衡后体系的动能会基本上等于体系的振动势能。也就是说，随着时间的推移，体系的动能平均值会减半，减少的部分变成了原子的振动势能。从图（a-b）可以看到，这个过程是很快的。该过程实际上就是一个从非平衡态跑向平衡态的过程。在这个例子中，这个过程只需要 1 ps 的量级。

图 （d）给出了 $(E(t)-\langle E\rangle)/|\langle E\rangle|$，即总能的相对涨落值。因为总能量在模拟的初期也有一个突然的变化，我们在计算总能的相对涨落时去掉了第一组输出的能量值（这相当于去掉了一个远离平衡态的时刻）。总能量相对涨落值在 $10^{-4}$ 量级。对于很小的体系来说，这是一个合理的值。

## 三斜盒子

### 三斜盒子的定义

我们首先将上一章的正交盒子进行推广。在正交盒子中，我们只有三个关于盒子的自由度，
分别是三个盒子边长 $L_x$、 $L_y$ 和 $L_z$。对于这样的盒子，我们其实假设了 $L_x$ 朝着 $x$ 的方向，
$L_y$ 朝着 $y$ 的方向， $L_z$ 朝着 $z$ 的方向。其实，一个盒子中共点的三条有向线段可以表示成矢量，
可记为 $\vec{a}$、 $\vec{b}$、 $\vec{c}$。它们可以用分量表示为

$$
   \vec{a} = a_x \vec{e}_x + a_y \vec{e}_y + a_z \vec{e}_z;
$$

$$
   \vec{b} = b_x \vec{e}_x + b_y \vec{e}_y + b_z \vec{e}_z;
$$

$$
   \vec{c} = c_x \vec{e}_x + c_y \vec{e}_y + c_z \vec{e}_z.
$$

这样的具有9个自由度的盒子，叫做三斜盒子（triclinic box）。
对于正交盒子，我们有

$$
   \vec{a} = L_x \vec{e}_x + 0 \vec{e}_y + 0 \vec{e}_z;
$$

$$
   \vec{b} = 0 \vec{e}_x + L_y \vec{e}_y + 0 \vec{e}_z;
$$

$$
   \vec{c} = 0 \vec{e}_x + 0 \vec{e}_y + L_z \vec{e}_z.
$$

我们可以将三个矢量的一共 9 个分量组合成一个矩阵，称为盒子矩阵，记为

$$
    H = \left(
    \begin{array}{ccc}
    a_x & b_x & c_x \\
    a_y & b_y & c_y \\
    a_z & b_z & c_z 
    \end{array}
    \right).
$$

对于正交盒子，该矩阵为对角矩阵：

$$
    \left(
    \begin{array}{ccc}
    L_x & 0 & 0 \\
    0 & L_y & 0 \\
    0 & 0 & L_z 
    \end{array}
    \right).
$$

读者也许要问，盒子矩阵为什么不定义为上述定义的转置？实际上，这只是一个约定，也许只是在后面的计算中比较方便而已。读者可以思考如果在定义中加一个转置，后面的讨论会如何改变。

根据高等数学知识，我们知道三斜盒子的体积 $V$ 等于上述矩阵行列式 $\det(H)$ 的绝对值：

$$
V = |\det(H)|.
$$

对于 $3 \times 3$ 的矩阵的行列式，有如下计算公式：

$$
\det(H) = a_x  (b_y  c_z - c_y  b_z) +
          b_x  (c_y  a_z - a_y  c_z) +
          c_x  (a_y  b_z - b_y  a_z).
$$

据此可写出如下 C++ 函数：

```C++
double getDet(const double* box)
{
  return box[0] * (box[4] * box[8] - box[5] * box[7]) +
         box[1] * (box[5] * box[6] - box[3] * box[8]) +
         box[2] * (box[3] * box[7] - box[4] * box[6]);
}
```

这里我们假设矩阵 $H$ 中的 9 个元素是按行主序（row-major order）的次序存储在数组 `double box[18]` 的前 9 个元素中的。
该数组有 18 个元素，这里只用了 9 个，剩下的 9 个元素将用于存放矩阵 $H$ 的逆:

$$
G = H^{-1}.
$$

矩阵求逆的 C++ 函数为：

```C++
void getInverseBox(double* box)
{
  box[9] = box[4] * box[8] - box[5] * box[7];
  box[10] = box[2] * box[7] - box[1] * box[8];
  box[11] = box[1] * box[5] - box[2] * box[4];
  box[12] = box[5] * box[6] - box[3] * box[8];
  box[13] = box[0] * box[8] - box[2] * box[6];
  box[14] = box[2] * box[3] - box[0] * box[5];
  box[15] = box[3] * box[7] - box[4] * box[6];
  box[16] = box[1] * box[6] - box[0] * box[7];
  box[17] = box[0] * box[4] - box[1] * box[3];
  double det = getDet(box);
  for (int n = 9; n < 18; ++n) {
    box[n] /= det;
  }
}
```


### 三斜盒子情况下的周期边界条件

为了理解盒子矩阵的作用，我们假设将一个坐标 $\vec{r}=(x,y,z)$ 表示为：

$$
\vec{r} = s_a \vec{a} + s_b \vec{b} + s_c \vec{c}.
$$

可以看出， 

$$
(s_a,s_b,s_c)=(0,0,0),(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),(1,1,1)
$$ 

分别代表盒子的 8 个顶点。那么，如果要求 $0\leq s_a\leq 1$,  $0\leq s_b\leq 1$,  $0\leq s_c\leq 1$，上述坐标就完全在盒子内（或者表面）。这里的坐标分量 $(s_a,s_b,s_c)$ 称为“分数坐标”。上式可用矩阵表示如下

$$
\left(
    \begin{array}{c}
    x  \\
    y  \\
    z 
    \end{array}
    \right)=
    \left(
    \begin{array}{ccc}
    a_x & b_x & c_x \\
    a_y & b_y & c_y \\
    a_z & b_z & c_z 
    \end{array}
    \right) 
    \left(
    \begin{array}{c}
    s_a  \\
    s_b  \\
    s_c 
    \end{array}
    \right).
$$

简写为：

$$
\left(
    \begin{array}{c}
    x  \\
    y  \\
    z 
    \end{array}
    \right)=
    H 
    \left(
    \begin{array}{c}
    s_a  \\
    s_b  \\
    s_c 
    \end{array}
    \right).
$$

反过来，我们有

$$
\left(
    \begin{array}{c}
    s_a  \\
    s_b  \\
    s_z 
    \end{array}
    \right)=
    G 
    \left(
    \begin{array}{c}
    x  \\
    y  \\
    z 
    \end{array}
    \right).
$$

也就是说，盒子矩阵的逆矩阵 $G$ 可以将一个盒子内的坐标变换为相对于盒子的分数坐标。所有的分数坐标形成一个边长为 1 的立方体。
结合上一章关于正交盒子的最小镜像约定，我们可以得到如下的关于三斜盒子的最小镜像约定的算法：

- 对于粒子 $i$ 和 $j$ 的相对坐标 

$$ 
\vec{r} _{ij} = (x _{ij}, y _{ij}, z _{ij}) 
$$
  
  首先用盒子逆矩阵 $G$ 将其变换为分数相对坐标 
  
$$ 
\vec{s} _{ij} = (\xi _{ij}, \eta _{ij}, \zeta _{ij})
$$
  
  变换关系为：
  
  $$
\left(
    \begin{array}{c}
    \xi _{ij}  \\
    \eta _{ij}  \\
    \zeta _{ij} 
    \end{array}
    \right)=
    G
    \left(
    \begin{array}{c}
    x _{ij}  \\
    y _{ij}  \\
    z _{ij}
    \end{array}
    \right).
$$
  
- 对分数相对坐标实施如下最小镜像约定操作。例如，当 $\xi _{ij}<-1/2$ 时，将其换为 $\xi _{ij}+1$；当 $\xi _{ij}>1/2$ 时，将其换为 $\xi _{ij}-1$。
- 将实施了最小镜像约定操作的分数相对坐标变换到普通相对坐标：

$$
\left(
    \begin{array}{c}
    x _{ij}  \\
    y _{ij}  \\
    z _{ij} 
    \end{array}
    \right)=
    H 
    \left(
    \begin{array}{c}
    \xi _{ij}  \\
    \eta _{ij}  \\
    \zeta _{ij}
    \end{array}
    \right).
$$

该算法由如下两个函数实现：

```C++
void applyMicOne(double& x12)
{
  if (x12 < -0.5)
    x12 += 1.0;
  else if (x12 > +0.5)
    x12 -= 1.0;
}

void applyMic(const double* box, double& x12, double& y12, double& z12)
{
  double sx12 = box[9] * x12 + box[10] * y12 + box[11] * z12;
  double sy12 = box[12] * x12 + box[13] * y12 + box[14] * z12;
  double sz12 = box[15] * x12 + box[16] * y12 + box[17] * z12;
  applyMicOne(sx12);
  applyMicOne(sy12);
  applyMicOne(sz12);
  x12 = box[0] * sx12 + box[1] * sy12 + box[2] * sz12;
  y12 = box[3] * sx12 + box[4] * sy12 + box[5] * sz12;
  z12 = box[6] * sx12 + box[7] * sy12 + box[8] * sz12;
}
```


## 近邻列表的创建

### 为什么要用近邻列表？

根据上一章的程序，我们知道，如果用一个自变量为原子数 $N$ 的函数表示程序的计算量 $y$，那么该函数一定是：

$$
    y = c_0 + c_1 N + c_2 N^2.
$$

其中， $c_2N^2$ 是求能量和力的部分的， $c_1N$ 主要对应运动方程积分的部分，而 $c_0$ 对应其它和 $N$ 无关的部分。当 $N$ 很大时， $c_2N^2$ 占主导，我们说该程序的计算复杂度是 $\mathcal{O}(N^2)$，或者说具有平方标度。对于这样的算法，粒子数很大时计算量会变得很大，以至于难以承受。

我们注意到，在上一章的 `findForce` 函数中，即使两个原子离得很远，我们也要对它们的距离进行判断，决定是否考虑它们之间的相互作用力。如果我们事先将在某一距离范围内的原子对记录下来，并在 `findForce` 函数中查看记录，就可以在 `findForce` 函数中仅考虑在某一距离范围内的原子对。这个“记录”就是近邻列表（neighbor list）。这个距离叫做近邻列表的截断距离。它应该比势函数的截断距离大一些，否则我们在每一步都需要先建立近邻列表，然后在 `findForce` 中使用。记势函数的截断距离为 $r _{\rm c}$，近邻列表的截断距离为 $R _{\rm c}$，一般来说用 $R _{\rm c}-r _{\rm c}=1$ Angstrom 是个不错的选择。当 $R _{\rm c}>r _{\rm c}$ 时，一个构建好的近邻列表将在若干步内都是“安全”的，不需要每一步都更新。如果用下面将要介绍的平方标度算法，整个程序的计算复杂度依然是 $\mathcal{O}(N^2)$，但平方项的系数 $c_2$ 将变小，从而有效地提升效率。如果用下面将要介绍的线性标度算法，整个程序的计算复杂度将变为 $\mathcal{O}(N)$，对于原子数很多的情形将具有很大的优势。我们稍后将给出具体的测试结果。

### 自动判断何时更新近邻列表

记 $R _{\rm c}-r _{\rm c}=\delta$，我们可以构造如下的算法，实现近邻列表更新的自动判断：
- 在程序的开头定义一套额外的坐标 $\{\vec{r}^0_i\}$，初始化为 0。
  
- 在积分过程的每一步，对每个原子 $i$ 计算距离 $d_i = |\vec{r}_i - \vec{r}^0_i|$，然后计算出这些距离中的最大值 $d _{\rm max}$。若 $2d _{\rm max}>\delta$，则更新近邻列表，同时更新 $\{\vec{r}^0_i\}$：

$$
\vec{r}^0_i = \vec{r}_i
$$
    
这是一个非常保守的更新判据。读者可以思考更加好的方案（使得更新频率降低）。

### 构建近邻列表的平方标度算法

我们先讨论一个简单的平方标度算法。首先，我们定义近邻列表。一个近邻列表指定了研究体系中每个原子的近邻个数，即与某个中心原子距离小于  $R _{\rm c}$ 的原子的个数。我们记原子 $i$ 的近邻个数为 $NN _{i}$。除此以外，确定一个近邻列表还需要知道原子 $i$ 的所有这 $NN _{i}$ 个近邻的指标。我们记原子 $i$ 的第 $k$ 个邻居的指标为 $NL _{ik}$。因为我们在求力的时候将利用牛顿第三定律，所以在构建近邻列表时也要求一个原子的近邻的指标大于原子本身的指标，即 $i < NL _{ik}$。这样定义的近邻列表最早由 Verlet 提出 [Computer "Experiments" on Classical Fluids. I. Thermodynamical Properties of Lennard-Jones Molecules](https://doi.org/10.1103/PhysRev.159.98)，所以称为 Verlet 近邻列表。

一个很自然的构建近邻列表的算法是检验所有的粒子对的距离。这显然是一个
$O(N^2)$ 复杂度的算法。我们的 C++ 实现如下：

```C++
void findNeighborON2(Atom& atom)
{
  const double cutoffSquare = atom.cutoffNeighbor * atom.cutoffNeighbor;
  std::fill(atom.NN.begin(), atom.NN.end(), 0);

  for (int i = 0; i < atom.number - 1; ++i) {
    const double x1 = atom.x[i];
    const double y1 = atom.y[i];
    const double z1 = atom.z[i];
    for (int j = i + 1; j < atom.number; ++j) {
      double xij = atom.x[j] - x1;
      double yij = atom.y[j] - y1;
      double zij = atom.z[j] - z1;
      applyMic(atom.box, xij, yij, zij);
      const double distanceSquare = xij * xij + yij * yij + zij * zij;
      if (distanceSquare < cutoffSquare) {
        atom.NL[i * atom.MN + atom.NN[i]++] = j;
        if (atom.NN[i] > atom.MN) {
          std::cout << "Error: number of neighbors for atom " << i
                    << " exceeds " << atom.MN << std::endl;
          exit(1);
        }
      }
    }
  }
}
```

下面是对该函数的一些说明：
- 我们调用了 C++ 标准库的 `std::fill()` 函数将每个原子的近邻个数置零，为后面的累加做准备。
- 语句

```C++
atom.NL[i * atom.MN + atom.NN[i]++] = j;
```

等价于如下两句：

```C++
atom.NL[i * atom.MN + atom.NN[i]] = j;
atom.NN[i]++;  
```

- 该函数对近邻列表的存储空间进行判断：如果任何原子的近邻个数超出了设定的最大值 `atom.MN`，就报告一个错误消息并退出程序。该最大值在该程序中设置为 1000，但读者可以适当地改动。这个数值设置得过大，就会浪费内存，设置得过小，就容易引发此处的错误。读者可以思考如下问题：如何修改程序，使得程序能够自动地调整相关数组的内存分配量，并总是让程序顺利地执行，而不会因为 `atom.MN` 设置得过小而退出程序，也不会因为 `atom.MN` 设置得过大而浪费太多的内存？

### 构建近邻列表的线性标度算法

上述简单的近邻列表算法已经能够加速程序了，但它还是一个 $\mathcal{O}(N^2)$ 复杂度的算法。
本节介绍一个 $\mathcal{O}(N)$ 复杂度的算法，即所谓的线性标度算法。
该算法的主要思想有以下两点：（1）整个模拟盒子被划分为一系列小胞（cell），每个胞的任何厚度不小于近邻列表的截断距离；
（2）对于每个原子，只需要在 27 个小胞（一个是原子所在的胞，另外 26 个是与该胞紧挨着的）中寻找近邻。 
在分子动力学模拟中最早提出此类方法的可能是 [Quentrec 和 Brot](https://doi.org/10.1016/0021-9991(73)90046-6) 。

先看第一点。首先，我们需要确定将整个体系划分为多少个小胞。我们针对一般的三斜盒子来讨论。对于三斜盒子，我们需要计算三个厚度，然后除以近邻截断，并向下取整，得到每个盒子矢量方向的小胞个数：

$$
N _{a} = \lfloor (V / A _{bc}) / R_c \rfloor;
$$

$$
N _{b} = \lfloor (V / A _{ca}) / R_c \rfloor;
$$

$$
N _{c} = \lfloor (V / A _{ab}) / R_c \rfloor.
$$

总的小盒子个数为 $N _{\rm cell} = N_a N_b N_c$.

计算盒子厚度的 C++ 函数如下：

```C++
float getArea(const double* a, const double* b)
{
  const double s1 = a[1] * b[2] - a[2] * b[1];
  const double s2 = a[2] * b[0] - a[0] * b[2];
  const double s3 = a[0] * b[1] - a[1] * b[0];
  return sqrt(s1 * s1 + s2 * s2 + s3 * s3);
}

void getThickness(const Atom& atom, double* thickness)
{
  double volume = abs(getDet(atom.box));
  const double a[3] = {atom.box[0], atom.box[3], atom.box[6]};
  const double b[3] = {atom.box[1], atom.box[4], atom.box[7]};
  const double c[3] = {atom.box[2], atom.box[5], atom.box[8]};
  thickness[0] = volume / getArea(b, c);
  thickness[1] = volume / getArea(c, a);
  thickness[2] = volume / getArea(a, b);
}
```

在确定小盒子的个数之后，我们就可以确定每个原子处于哪个小胞了。对应的 C++ 函数如下：

```C++
void findCell(
  const double* box,
  const double* thickness,
  const double* r,
  double cutoffInverse,
  const int* numCells,
  int* cell)
{
  double s[3];
  s[0] = box[9] * r[0] + box[10] * r[1] + box[11] * r[2];
  s[1] = box[12] * r[0] + box[13] * r[1] + box[14] * r[2];
  s[2] = box[15] * r[0] + box[16] * r[1] + box[17] * r[2];
  for (int d = 0; d < 3; ++d) {
    cell[d] = floor(s[d] * thickness[d] * cutoffInverse);
    if (cell[d] < 0)
      cell[d] += numCells[d];
    if (cell[d] >= numCells[d])
      cell[d] -= numCells[d];
  }
  cell[3] = cell[0] + numCells[0] * (cell[1] + numCells[1] * cell[2]);
}
```

第 9-12 行计算与原子坐标 $\vec{r}$ 对应的分数坐标 $\vec{s}$。
第 14 行计算原子在某个盒子矢量方向的小盒子指标。
第 15-18 行保证计算的小盒子指标不超出界限。
最后，第 20 行根据三维的小盒子指标计算一个一维指标。

有了以上准备之后，我们考虑第二点，即在 27 个小盒子内寻找一个原子的近邻。
实现该部分的 C++ 函数如下。第 5-6 行计算三个盒子厚度。
第 8-12 行计算每个盒子矢量方向的小盒子个数。
第 17 行定义一个数组，代表每个小盒子中的原子个数 $C_n (0 \leq n \leq N _{\rm cell} - 1)$，
其计算过程在第 20-24 行。
第 18 行定义一个数组 $S_n (0 \leq n \leq N _{\rm cell} - 1)$，其元素定义为：

$$
S_0 = 0;
$$

$$
S_n = \sum _{m=0}^{n-1} C_m \quad (1 \leq n \leq N _{\rm cell}-1).
$$

用计算机的术语来说，数组 $S_n$ 是数组 $C_n$ 的前缀和（prefix sum），
其计算过程在第 26-28 行。第 29 行将数组 $C_n$ 的元素置零，因为后面又要对它进行累加。
第 32-39 行计算根据小盒子指标排列的原子指标，保存为一个长度为原子数 $N$ 的数组 $I$。
其中，从 $I _{S_n}$ 到 $I _{S_n + C_n - 1}$ 的元素就是处于小盒子 $n$ 的原子的指标。
第 41-85 行根据以上几个辅助数组的信息构建 Verlet 近邻列表，其算法和前面的 $\mathcal{O}(N)$ 算法差不多，
只不过我们只需要考虑 27 个小盒子，而不是整个大盒子。

```C++
void findNeighborON1(Atom& atom)
{
  const double cutoffInverse = 1.0 / atom.cutoffNeighbor;
  double cutoffSquare = atom.cutoffNeighbor * atom.cutoffNeighbor;
  double thickness[3];
  getThickness(atom, thickness);

  int numCells[4];

  for (int d = 0; d < 3; ++d) {
    numCells[d] = floor(thickness[d] * cutoffInverse);
  }

  numCells[3] = numCells[0] * numCells[1] * numCells[2];
  int cell[4];

  std::vector<int> cellCount(numCells[3], 0);
  std::vector<int> cellCountSum(numCells[3], 0);

  for (int n = 0; n < atom.number; ++n) {
    const double r[3] = {atom.x[n], atom.y[n], atom.z[n]};
    findCell(atom.box, thickness, r, cutoffInverse, numCells, cell);
    ++cellCount[cell[3]];
  }

  for (int i = 1; i < numCells[3]; ++i) {
    cellCountSum[i] = cellCountSum[i - 1] + cellCount[i - 1];
  }

  std::fill(cellCount.begin(), cellCount.end(), 0);

  std::vector<int> cellContents(atom.number, 0);

  for (int n = 0; n < atom.number; ++n) {
    const double r[3] = {atom.x[n], atom.y[n], atom.z[n]};
    findCell(atom.box, thickness, r, cutoffInverse, numCells, cell);
    cellContents[cellCountSum[cell[3]] + cellCount[cell[3]]] = n;
    ++cellCount[cell[3]];
  }

  std::fill(atom.NN.begin(), atom.NN.end(), 0);

  for (int n1 = 0; n1 < atom.number; ++n1) {
    const double r1[3] = {atom.x[n1], atom.y[n1], atom.z[n1]};
    findCell(atom.box, thickness, r1, cutoffInverse, numCells, cell);
    for (int k = -1; k <= 1; ++k) {
      for (int j = -1; j <= 1; ++j) {
        for (int i = -1; i <= 1; ++i) {
          int neighborCell = cell[3] + (k * numCells[1] + j) * numCells[0] + i;
          if (cell[0] + i < 0)
            neighborCell += numCells[0];
          if (cell[0] + i >= numCells[0])
            neighborCell -= numCells[0];
          if (cell[1] + j < 0)
            neighborCell += numCells[1] * numCells[0];
          if (cell[1] + j >= numCells[1])
            neighborCell -= numCells[1] * numCells[0];
          if (cell[2] + k < 0)
            neighborCell += numCells[3];
          if (cell[2] + k >= numCells[2])
            neighborCell -= numCells[3];

          for (int m = 0; m < cellCount[neighborCell]; ++m) {
            const int n2 = cellContents[cellCountSum[neighborCell] + m];
            if (n1 < n2) {
              double x12 = atom.x[n2] - r1[0];
              double y12 = atom.y[n2] - r1[1];
              double z12 = atom.z[n2] - r1[2];
              applyMic(atom.box, x12, y12, z12);
              const double d2 = x12 * x12 + y12 * y12 + z12 * z12;
              if (d2 < cutoffSquare) {
                atom.NL[n1 * atom.MN + atom.NN[n1]++] = n2;
                if (atom.NN[n1] > atom.MN) {
                  std::cout << "Error: number of neighbors for atom " << n1
                            << " exceeds " << atom.MN << std::endl;
                  exit(1);
                }
              }
            }
          }
        }
      }
    }
  }
}
```

### 程序速度测试
  
![energy](examples/neighbor.png)

用本章的程序 `md2.cpp` 进行测试，得到如图 1 的结果。
此图对比了三种近邻列表方案下程序跑 1000 步所花的时间和体系原子数的关系。
这里的三种近邻列表方案分别是：
- 不使用近邻列表（类似于上一章的程序，但注意我们从本章开始使用了三斜的盒子）；
- 使用 $\mathcal{O}(N^2)$ 计算复杂度的近邻列表构建方法；
- 使用 $\mathcal{O}(N)$ 计算复杂度的近邻列表构建方法。

该图的结果是符合预期了，它展示了如下特征：
- 不使用近邻列表时，程序的计算量正比于原子数的平方，拥有典型的 $\mathcal{O}(N^2)$ 计算复杂度。
- 使用 $\mathcal{O}(N^2)$ 计算复杂度的近邻列表构建方法时，程序的计算量在小体系的极限下正比于原子数，但在大体系的极限下正比于原子数平方。
- 使用 $\mathcal{O}(N)$ 计算复杂度的近邻列表构建方法时，程序的计算量始终正比于原子数，具有 $\mathcal{O}(N)$ 复杂度。
所以，在研究较大体系时，一定要采用 $\mathcal{O}(N)$ 计算复杂度的近邻列表构建方法。

## GPUMD程序简介

本章介绍的近邻列表技术只是提高 MD 程序效率的一个方面。在实际的 MD 程序中，往往还需要使用并行编程的技术进一步提高效率。并行编程是和计算机硬件紧密相关的。对于 CPU 计算来说，MPI (message passing interface) 是最常用的选择，它可以在一个节点内并行，也可以在不同的节点之间并行。另一个并行技术，OpenMP，就只能在节点内并行。当前，图形处理器（GPU) 相对于 CPU 有很大的性能优势。Nvidia 的GPU可使用 CUDA (Compute Unified Device Architecture) 进行编程，而 AMD 的GPU可使用HIP (Heterogeneous Interface for Portability) 进行编程。

当前学术界使用最为广泛的两个 MD 程序是 [LAMMPS](https://www.lammps.org) 和 [Gromacs](https://www.gromacs.org)，它们都广泛使用了 MPI、OpenMP 和 CUDA 等并行技术加速。笔者主导开发了另一个 MD 程序 [GPUMD](https://www.gpumd.org)，大量使用了 CUDA 编程，获得了很高的计算性能。本书除了使用自编 C++ 程序之外，还会适当使用 GPUMD。所以，本节简要介绍 GPUMD。

GPUMD 是开源程序，在Github托管，地址为：[https://www.gpumd.org](https://github.com/brucefan1983/GPUMD)

从 Github 下载程序包之后，解压缩，从终端进入 `src` 文件夹，敲 `make` 即可完成程序的编译。这要求读者的计算机有 CUDA 编程环境和 Nvidia 的 GPU。关于 CUDA 编程的基础知识，请参考笔者的另一本书 《CUDA编程：基础与实践》。编译好之后，会在 `src` 文件夹内产生 `gpumd` 和 `nep` 两个可执行文件。其中，`gpumd` 可执行文件就类似于我们目前使用过的自编 MD 程序，需要`model.xyz`和`run.in`两个输入文件，而`nep`可执行文件是用来训练机器学习势函数的（将在下一章介绍）。我们将在后面的章节通过具体的例子逐步介绍 GPUMD 程序的使用。感兴趣的读者可以浏览该程序的使用手册：[https://www.gpumd.org](https://gpumd.org/)

## 习题

### 自动更新近邻判据的改进

在正文中，我们仅仅计算了体系中移动距离最大的粒子。其实，计算体系中移动距离最大和次大的粒子，然后判断何时更新近邻列表，有可能提升效率。请据此修改本章程序，并将测试结果和正文的结果进行比较。

### 小盒子边长为截断距离一半的情况

在正文中，我们让小盒子的边长等于或大于近邻列表的截断距离。其实，将小盒子的边长取得小一点，有可能提升效率。请修改本章程序，将小盒子的边长设置为近邻列表截断距离的一半。将测试结果和正文的结果进行比较。

## $\delta$ 参数的优化

在正文中，我们假设近邻 $\delta$ 为 1 Å。请测试该参数对计算时间的影响，找出若干情况下最优的 $\delta$ 值。

