
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

## 基于紧束缚模型的电子输运性质

### 与分子动力学模拟耦合的线性标度量子输运

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
