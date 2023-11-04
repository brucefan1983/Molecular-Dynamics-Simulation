# 《分子动力学模拟》第一章：经典力学

本章程序 harmonic_oscillator.m 的结果。

## 简谐振子运动的数值求解程序


## 源代码

见 src/harmonic_oscillator.m

我们用一个简谐振子模型来展示速度-Verlet算法的实现。简谐振子偏离平衡位置的坐标为 $x$ ， 受力为 $-kx$ ， $k$ 是弹簧的劲度系数。


## 结果分析

下面是坐标和速度随时间变化的图：

![position_and_momentum](src/position_and_momentum.png)


下面是简谐振子的相空间图：

![phase_space](src/phase_space.png)


下面的结果展示了总能量守恒：

![phase_space](src/fig-chapter1-energy_conservation.png)
