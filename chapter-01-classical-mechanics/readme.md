# 《分子动力学模拟》第一章：经典力学基础

在学习分子动力学模拟之前，需要熟悉经典力学、热力学和统计物理。本章将简要地回顾经典力学的知识并通过简谐振子开始编写一个最简单的类似于分子动力学模拟的小程序。

# Table of contents
- [牛顿力学](#牛顿力学)
 	- [质点力学](#质点力学)
 	- [粒子系力学](#粒子系力学)
 	- [牛顿运动方程的数值积分](#牛顿运动方程的数值积分)
 	- [简谐振子运动的数值求解](#简谐振子运动的数值求解)
- [分析力学](#分析力学)
  - [拉格朗日方程](#拉格朗日方程)
  - [哈密顿方程](#哈密顿方程)
  - [相空间](#相空间)
  - [哈密顿体系运动方程的数值积分](#哈密顿体系运动方程的数值积分)

## 牛顿力学

### 质点力学

理论上，牛顿力学里面最基本的研究对象是质点，它是指一个有特定的质量但其大小对所研究的问题不重要的物体。有人也称它为一个粒子。这也将是分子动力学模拟中最重要的研究对象之一。 一个粒子在三维空间的运动由一个位置函数 $\vec{r}(t)$ 完全地描述。也就是说，在一个时刻 $t$ ，粒子的位置由一个有三个分量的矢量 

$$
\vec{r}(t)=x(t)\vec{e}_x + y(t) \vec{e}_y + z(t) \vec{e}_z
$$ 

给出。当然，如果不指定一个参考系，我们将无法确定这三个分量的值。参考系的数学本质是坐标系。我们说粒子的位置函数完全地描述了粒子的运动性质，是因为其它的运动性质，比如速度和加速度，都可以由位置函数导出。

速度函数定义为位置函数对时间的一阶导数：

$$
\vec{v}(t) = \frac{d \vec{r}(t)}{d t} = \dot{\vec{x}}(t)
=\frac{dx}{dt} \vec{e}_x + \frac{dy}{dt} \vec{e}_y + \frac{dz}{dt} \vec{e}_z.
$$

加速度函数定义为位置函数对时间的二阶导数，或者速度函数对时间的一阶导数：

$$
\vec{a}(t) = \frac{d ^2\vec{r}(t)}{d t^2} = \frac{d \vec{v}(t)}{d t}
= \ddot{\vec{r}}(t) = \dot{\vec{v}}(t).
$$

以上对粒子运动的描述就是所谓的运动学。下面我们讨论质点的动力学。

动力学研究物体运动及其改变的原因。牛顿运动定律让我们可以在一定的条件下确定一个粒子的位置函数 $\vec{r}(t)$。

牛顿第一定律。一个粒子在施加于它的合外力 $\vec{F}$ 为零时保持静止或者一个恒定的速度 $\vec{v}$ :

$$
\vec{F} = \vec{0} \Rightarrow  \vec{v} = \text{const}.
$$

该定律也叫惯性定律，或者惰性定律，意思是在合外力为零的情况下，一个粒子懒得改变它原来的运动状态：原来静止的话就继续静止；原来运动的话就继续原来的运动。

牛顿第二定律。一个粒子的动量 $\vec{p}$ 的时间变化率 $d\vec{p}/dt$ 正比于作用于它的合外力 $\vec{F}$：

$$
\vec{F} = \frac{d \vec{p}}{d t}.
$$

一个粒子的动量定义为其惯性质量 $m$ 与速度 $\vec{v}$ 的乘积：

$$
\vec{p}=m\vec{v}
$$

如果粒子的惯性质量是个常数的话，那么上式就变成

$$
\vec{F} = m\frac{d \vec{v}}{d t} = m \vec{a}.
$$

我们强调了这里的 $m$ 是惯性质量，原因是牛顿力学里面还有另外一个质量，叫做引力质量，虽然它们是等价的。

牛顿第三定律。每一个作用力都有一个与之大小相等、方向相反的反作用力。例如，如果有一个由粒子 2 作用在粒子 1 上的作用力  $\vec{F}_1{}_2$ ，那么一定同时存在一个由粒子 1 作用在粒子 2 上的反作用力  $\vec{F}_2{}_1$ ，它们大小相等，方向相反：

```math
\vec{F}_1{}_2 = -\vec{F}_2{}_1.
```

当然，作用力和反作用力是相对的；把一对作用力中的哪一个叫做作用力，哪一个叫做反作用力是随意的。上式表达的是牛顿第三定律的弱形式。还有一个强形式的牛顿第三定律：

$$
\vec{F}_1{}_2 \propto \vec{r}_1{}_2 \equiv \vec{r}_2 - \vec{r}_1
$$

显然，满足了牛顿第三定律的强形式，就一定满足了牛顿第三定律的弱形式，但反之不然。

我们在陈述牛顿第二定律时已经定义了动量。根据牛顿第二定律，如果作用在一个粒子上的合力为零，那么该粒子的动量将保持为一个常数：

$$
\vec{F} = \vec{0} \Rightarrow \vec{p} = \text{常数}.
$$

这就叫做动量守恒定律。

一个粒子的角动量被定义为：

$$
\vec{L} = \vec{r} \times \vec{p}.
$$

因为位置 $\vec{r}$ 依赖于坐标系原点的选取，角动量也依赖于坐标原点的选择。类似地，我们定义作用在粒子上的力矩 $\vec{\tau}$ 为：

$$
\vec{\tau} = \vec{r} \times \vec{F}.
$$

其中， $\vec{F}$ 是作用在粒子上的力。让我们来考察角动量的时间变化率：

$$
\frac{d \vec{L}}{dt} = \frac{d (\vec{r} \times \vec{p})}{dt}
= \vec{v} \times \vec{p} + \vec{r} \times \vec{F}
= \vec{r} \times \vec{F} = \vec{\tau}.
$$

这就是说，角动量的变化率等于力矩。显然，如果作用在粒子上的力矩为零，则粒子的角动量是常数：

$$
\vec{\tau} = \vec{0} \Rightarrow \vec{L} = \text{常数}.
$$

这叫做角动量守恒定律。

除了动量，还可以定义一个仅仅与质量和速度有关的物理量，叫做动能：

$$
T = \frac{1}{2}m \vec{v}^2 =
\frac{1}{2}m \vec{v} \cdot \vec{v}
= \frac{1}{2}m (v_x^2 + v_y^2 + v_z^2).
$$

在一个力 $\vec{F}$ 的作用下，如果粒子运动了一个微分位移

$$
d\vec{x} = dx \vec{e}_1 + dy \vec{e}_2 + dz \vec{e}_3
$$

那么我们定义该力对该粒子做的微功为：

$$
d W = \vec{F} \cdot d\vec{x} = F_x dx + F_y dy + F_z dz.
$$

进一步推导可得

$$
d W
= \vec{F} \cdot \vec{v} dt
= (\vec{F} dt) \cdot \vec{v}
= d\vec{p} \cdot \vec{v}
= m d\vec{v} \cdot \vec{v}
= d \left( \frac{1}{2}  m \vec{v}^2 \right)
= dT.
$$

所以，在一个过程中，外力对一个粒子做的功等于粒子动能的改变量。

外力 $\vec{F}$ 做的功一般来说依赖于路径的选取。然而，如果外力做的功与具体的路径无关，而只与起点和终点有关，那么由矢量分析的定理可知，该外力可以写成一个标量场 $U(\vec{r})$ 的梯度的负值：

$$
\vec{F} = -\nabla U(\vec{r}).
$$

该标量场 $U(\vec{r})$ 叫做粒子的势能场，简称为势能，或者进一步简称为势。与此对应的力称为保守力。该保守力沿任意微小路径对粒子做的总功为：

$$
dW = - \nabla U(\vec{r}) \cdot d\vec{r} = - dU.
$$

将该式与前面得到的 $dW = dT$ 比较可得：

$$
dT + dU = 0.
$$

这就是说，在一个保守力的作用下，任意过程中粒子的动能与势能的和都不改变。这个不改变的量称为机械能。所以，在保守力的作用下，粒子的机械能是守恒的。

### 粒子系力学

在具体讨论粒子系的动力学行为之前必须先弄清楚内力和外力的区别。内力是系统中某个粒子作用于另一个粒子的，而外力来自于系统之外。内力满足牛顿第三定律。虽然外界施与粒子 $i$ 一个力 $\vec{F}^{\text{ext}}_{i}$ 的同时，粒子 $i$ 同时也施与外界一个大小相等、方向相反的反作用力，但由于我们的系统不包含外界，我们通常不会对外力利用第三定律。以上就是内力与外力的区别。

对任意一个粒子 $i$，我们可以写下它的动力学方程：

$$
m_i \ddot{\vec{r}}_i = \sum _{{ j \neq i }} \vec{F}_i{}_j + \vec{F}^{\text{ext}}_i.
$$

将上式左右两边都对指标 $i$ 求和，得

$$
\sum_i m_i \ddot{\vec{r}}_i = \sum_i \sum _{{j \neq i}} \vec{F}_i{}_j + \sum_i \vec{F} ^{\text{ext}} _i.
$$

根据牛顿第三定律，等号右边的第一项等于零。系统的总质量定义为 $m = \sum_i m_i$。如果定义一个平均坐标

$$
\vec{r} = \frac{\sum_i m_i \vec{r}_i}{m},
$$

那么我们有

$$
m \ddot{\vec{r}} = \sum_i \vec{F}^{\text{ext}}_{i}.
$$

这个式子看上去很像一个质量为系统总质量 $m$，坐标为平均坐标 $\vec{x}$ 的粒子的动力学方程；该粒子所受的合外力为整个粒子系所受的合外力。我们称这个等效的粒子为粒子系的质心。质心的质量就是整个系统的质量；质心坐标就是上述平均坐标；质心的运动满足上述等效的牛顿第二定律。

由质心坐标可以定义质心速度

$$
\dot{\vec{r}} = \frac{\sum_i m_i \dot{\vec{r}}_i}{m}
$$

和质心动量

$$
\vec{p} = m \dot{\vec{r}} = \sum_i m_i \dot{\vec{r}}_i.
$$

于是，质心的牛顿第二定律可用质心动量表达为：

$$
\frac{d \vec{p}}{dt} = \sum_i \vec{F}^{\text{ext}}_{i}.
$$

如果系统受到的合外力为零，那么系统的质心动量（即系统的总动量）是守恒的。这就是质点系的动量守恒定律。

类似地，由牛顿第三定律可以证明内力对系统的总力矩的贡献也是零，即系统所受的总力矩等于外力的总力矩 $\vec{\tau}^{\text{ext}}$：

$$
\vec{\tau}^{\text{ext}} = \sum_i^N \vec{r}_i \times \vec{F}_i^{\text{ext}}.
$$

如果定义系统的总角动量为

$$
\vec{L} = \sum_i^N \vec{L}_i = \sum_i^N \vec{r}_i \times \vec{p}_i,
$$

那么可以证明如下的角动量定理：

$$
\frac{d \vec{L}}{dt} = \vec{\tau}^{\text{ext}}.
$$

这就是说，外力产生的总力矩等于系统总角动量的时间变化率。如果外力产生的总力矩等于零，则有系统的总角动量守恒。这就是质点系的角动量守恒定理。

对于有相互作用的多粒子系统，如果内力和外力都是保守力，那么系统的总势能可以表达为

$$
U = \sum_i U_i + \frac{1}{2} \sum_i^N \sum_{j\neq i}^N U_i{}_j,
$$

其中， $U_i$ 是第 $i$ 个粒子在外力场中的势能， $U_i{}_j$ 是系统中由 $i$ 与 $j$ 的相互作用导致的势能。


### 牛顿运动方程的数值积分

给定一个多粒子体系的初始状态（坐标和速度），根据各个粒子之间的相互作用力就可以预测该体系的运动状态，即任意时刻各个粒子的坐标和速度。该预测过程本质上就是对运动方程的数值积分。

我们对粒子 $i$ 在 $t+\Delta t$ 时刻的坐标做泰勒级数展开：

$$
\vec{r}_i(t+\Delta t) \approx \vec{r}_i(t) + \vec{v}_i(t) \Delta t + \frac{1}{2} \frac{\vec{F}_i(t)}{m_i} \Delta t^2.
$$

我们也可考虑一个过去的时刻 $t-\Delta t$ 并做类似的展开

$$
\vec{r}_i(t-\Delta t) \approx \vec{r}_i(t) - \vec{v}_i(t) \Delta t + \frac{1}{2} \frac{\vec{F}_i(t)}{m_i} \Delta t^2.
$$

由以上两式可得

$$
\vec{r}_i(t+\Delta t) \approx 2\vec{r}_i(t) - \vec{r}_i(t-\Delta t) + \frac{\vec{F}_i(t)}{m_i} \Delta t^2.
$$

这就是所谓的 [Verlet积分算法](https://doi.org/10.1103/PhysRev.159.98)，它只涉及坐标，不涉及速度。如果要获得速度，需要通过如下差分求得：

$$
\vec{v}_i(t) \approx \frac{ \vec{r}_i(t+\Delta t) - \vec{r}_i(t-\Delta t) }{2\Delta t}.
$$

Verlet积分算法中的速度计算涉及到时间上相差 $2\Delta t$ 的坐标，不是很方便。为了得到更加方便的速度计算方式，我们考虑如下两个展开：

$$
\vec{r}_i(t+\Delta t) \approx \vec{r}_i(t) + \vec{v}_i(t) \Delta t + \frac{1}{2} \frac{\vec{F}_i(t)}{m_i} \Delta t^2.
$$

$$
\vec{r}_i(t) \approx \vec{r}_i(t+\Delta t) - \vec{v}_i(t+\Delta t) \Delta t + \frac{1}{2} \frac{\vec{F}_i(t+\Delta t)}{m_i} \Delta t^2.
$$

由此可得到如下速度计算公式：

$$
\vec{v}_i(t+\Delta t) \approx \vec{v}_i(t) + \Delta t \frac{ \vec{F}_i(t) + \vec{F}_i(t+\Delta t) }{2m_i}.
$$

以上就是 [速度-Verlet 积分算法](https://doi.org/10.1063/1.442716)。

可以看出， $t+\Delta t$ 时刻的坐标仅依赖于 $t$ 时刻的坐标、速度和力，但 $t+\Delta t$ 时刻的速度依赖于 $t$ 时刻的速度、力及 $t+\Delta t$ 时刻的力。所以，从算法的角度来说，速度-Verlet 积分算法对应如下的计算流程：

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

### 简谐振子运动的数值求解

我们用一个简谐振子模型来展示速度-Verlet算法的实现。简谐振子偏离平衡位置的坐标为 $x$ ， 受力为 $-kx$ ， $k$ 是弹簧的劲度系数。

Matlab 代码如下：
```matlab
m=1; k=1; dt=0.01; n_step=1000;
v=0; x=1; 
v_vector=zeros(n_step,1); x_vector=zeros(n_step,1);
for step=1:n_step
    v=v-(dt/2)*k*x; 
    x=x+dt*v;
    v=v-(dt/2)*k*x; 
    v_vector(step,:)=v; x_vector(step,:)=x;
end
```

下面是坐标和速度随时间变化的图：


![position_and_momentum](src/position_and_momentum.png)

## 分析力学

### 拉格朗日方程

我们知道，一个由 $N$ 个质点构成的力学系统需要用 $3N$ 个坐标分量及其相应的速度分量来描述其运动状态。我们说，该系统的力学自由度为 $3N$。 但如果是一个由这 $N$ 个质点构成的刚体，则自由度只有 6：3 个平动自由度和 3 个转动自由度。刚体的自由度之所以小于自由质点系统的自由度，是因为刚体系统中有约束。

约束的分类很复杂，但我们只关注最简单的一种：稳定的几何约束。一个原本有 $M$ 个自由度的体系，若受到 $m$ 个几何约束，则其自由度为 $s=M-m$。要描述该体系，我们不一定需要使用原来的 $M$ 个坐标和速度，而可以构造 $s$ 个新的坐标和速度。这 $s$ 个新的坐标不一定限于在原来的 $M$ 个老坐标中挑选，而完全可以另行选取。这样选取的新坐标叫做广义坐标。广义坐标不一定能三个一组地构成矢量，而且某个广义坐标也不一定具有长度的量纲。相应地，某个广义速度也不一定具有速度的量纲。

将 $s$ 个广义坐标 $q_1, q_2, \cdots q_s$ 的集合简记为 $q$。对于保守力体系，可以从牛顿力学出发推导出如下拉格朗日方程：

$$
\frac{d}{dt}
\left(\frac{\partial (T-U)}{\partial \dot{q}_{\alpha}}\right)
-\frac{\partial (T-U)}{\partial q _{\alpha} } = 0.
\quad (\alpha = 1, 2, \cdots, s)
$$

这里出现的动能与势能的差是一个很重要的量，叫做拉格朗日量，记为

$$
L(q, \dot{q}) = T(\dot{q}) - U(q).
$$

用拉格朗日量可以将拉格朗日方程写成更加简洁的形式：

$$
\boxed{
\frac{d}{dt}
\left(\frac{\partial L}{\partial \dot{q} _{\alpha} }\right)
-\frac{\partial L}{\partial q _{\alpha} } = 0.
\quad (\alpha = 1, 2, \cdots, s)
}
$$

例子：简谐振子的拉格朗日量为

$$
L = \frac{1}{2} m \dot{x}^2 - \frac{1}{2} k x^2
$$

由此可以推导出简谐振子的运动方程

$$
m \ddot{x} + k x = 0.
$$


### 哈密顿方程

由于拉格朗日量具有能量的量纲，故当广义坐标具有长度的量纲时，偏导数 $\frac{\partial L}{\partial \dot{q}_{\alpha}}$ 具有动量的量纲。我们称这个偏导数为广义动量，记为

$$
\boxed{p_{\alpha} = \frac{\partial L}{\partial \dot{q}_{\alpha}}}.
$$

注意，当广义坐标不具有长度的量纲时，相应的广义动量也不具有动量的量纲。运用广义动量可以将拉格朗日方程写成更加简洁的形式：

$$
\boxed{
\dot{p} _{\alpha}
-\frac{\partial L}{\partial q _{\alpha} } = 0.
\quad (\alpha = 1, 2, \cdots, s)
}
$$

虽然广义动量出现在拉格朗日方程中，但我们要清楚的是拉格朗日量是广义坐标和广义速度的函数，而不是广义动量的函数。用勒让德变换可以改变一个多元函数的独立变量。为此，我们先写下拉格朗日量的全微分：

$$
d L(q, \dot{q}) = \frac{\partial L}{\partial q _{\alpha} } d q _{\alpha}  + \frac{\partial L}{\partial \dot{q} _{\alpha} } d \dot{q} _{\alpha} = \dot{p} _{\alpha} d q _{\alpha} + p _{\alpha} d \dot{q} _{\alpha}.
$$

勒让德变换能将独立变量从广义速度变成广义动量。勒让德变换是指将函数减去它的某个独立变量与它对该独立变量的偏导数的乘积。选取广义速度为该独立变量，并将变换后的函数取一个相反数，则有如下勒让德变换

$$
L \rightarrow p_{\alpha} \dot{q}_{\alpha} - L.
$$

变换后的函数

$$
H = p_{\alpha} \dot{q}_{\alpha} - L.
$$

称为哈密顿量，它的全微分为：

$$
d H = p _{\alpha} d\dot{q} _{\alpha} + \dot{q} _{\alpha} d p _{\alpha} - (\dot{p} _{\alpha} d q _{\alpha} + p _{\alpha} d \dot{q} _{\alpha}) =\dot{q} _{\alpha} d p _{\alpha} - \dot{p} _{\alpha} d q _{\alpha}.
$$

看来哈密顿量确实是广义动量以及广义坐标的函数，而不是广义速度的函数了。

既然哈密顿量是广义动量和广义坐标的函数，那么根据全微分的定义，我们又有

$$
d H = \frac{\partial H}{\partial q _{\alpha}} d q _{\alpha} + \frac{\partial H}{\partial p _{\alpha} } d p _{\alpha}.
$$

对比以上两个式子，我们得到如下两组重要的方程：

$$
\dot{q} _{\alpha} = \frac{\partial H}{\partial p _{\alpha}}, \quad (\alpha = 1 , 2, \cdots, s)
$$

$$
\dot{p} _{\alpha} = - \frac{\partial H}{\partial q _{\alpha}}. \quad (\alpha = 1 , 2, \cdots, s)
$$

这 $2s$ 个一阶微分方程组称为哈密顿正则方程。这里的“正则”意为“简单且对称”。

对于由 $N$ 个粒子组成的体系，根据上述哈密顿量的定义可知：

$$
H = \sum_i \frac{\vec{p}_i^2}{2m_i} + U(\vec{r}_1, \vec{r}_2, \cdots, \vec{r}_N)
$$

由此可见，哈密顿就是体系的总能量。

例子：简谐振子的哈密顿量为

$$
H = \frac{p^2}{2m} + \frac{1}{2} k x^2.
$$

由此可得运动方程为：

$$
\dot{x} = \frac{p}{m}
$$

$$
\dot{p} = - k x.
$$

这都和牛顿力学的结果是一致的。

### 相空间

广义坐标和广义动量是相互独立的变量。我们可以将 $s$ 个广义坐标看成一个 $2s$ 维“空间”的“坐标”。这个抽象的“空间”叫做相空间。一组给定的广义坐标和广义动量叫做相空间的一个相点。另外，根据哈密顿正则方程，只要给定一个初始条件，即初始时刻的各个广义坐标和广义动量，就可以唯一地确定任意时刻的各个广义坐标和广义动量，即如下 $2s$ 个函数：

$$
q _{\alpha} = q _{\alpha}(t), \quad p _{\alpha} = p _{\alpha}(t). \quad (\alpha = 1, 2, \cdots, s)
$$

这 $2s$ 个函数将在 $2s$ 维的相空间给出一条轨迹，叫做相轨迹。随着时间的推移，一条相轨迹会在相空间跑动，历经很多不同的相点。一个系统中不同的初始条件会给出不同的相轨迹，而两条不同的相轨迹绝对不会相交于某一个相点。

例子：简谐振子的相空间。

![phase_space](src/phase_space.png)

对于哈密顿体系，我们还可以证明，相空间是不可压缩的。为此，我们先将广义坐标和动量简写为如下的 $2s$ 分量的矢量：

$$
x \equiv (q_1, q_2, \cdots, q_s, p_1, p_2, \cdots, p_s).
$$

这个矢量就代表一个相空间点，它代表相空间的“坐标”。 根据哈密顿正则方程，该相空间“坐标”的速度为

$$
\dot{x} = \left(\frac{\partial H}{\partial p_1}, \frac{\partial H}{\partial p_2}, \cdots, 
\frac{\partial H}{\partial p_s}, -\frac{\partial H}{\partial q_1}, -\frac{\partial H}{\partial q_2}, \cdots, -\frac{\partial H}{\partial q_s}\right).
$$

在流体力学中，流体的不可压缩性指的是其中的流速场的散度为零（即没有源和汇）。将相空间类比为流体，那么相空间的不可压缩性指的是：

$$
\nabla_{x} \cdot \dot{x} = 0.
$$

这是很显然的：

$$
\nabla_{x} \cdot \dot{x} =
\sum_{\alpha} \frac{\partial^2 H}{\partial p_{\alpha} \partial q_{\alpha}} -
\frac{\partial^2 H}{\partial q_{\alpha} \partial p_{\alpha}} = 0.
$$

所以，哈密顿体系的相空间是不可压缩的。


哈密顿体系相空间的不可压缩性也意味着刘伟尔定理。首先定义相空间的体积元：

$$
d x = dx_1  dx_2 \cdots dx_{2s}
$$

记某个初始时刻 $t=0$ 的相空间点为 $x_0$，时刻 $t$ 的相空间点为 $x_t$， 刘伟尔定理是说

$$
d x_t = dx_0
$$

即相空间的体积元的体积是守恒的。 

为了证明该等式，我们首先注意到，两个体积元可由一个雅可比行列式联系：

$$
d x_t = \det(J) dx_0
$$

为了确定雅可比行列式，我们先计算它的时间导数：

$$
\frac{d }{dt} \det(J) = \frac{d}{dt} e^{\mathrm{Tr}[\ln (J)]}
 = \det(J)  \mathrm{Tr}[ \frac{dJ}{dt} J^{-1}]
$$

上面利用了等式

$$
\det(J) = e^{\mathrm{Tr}[\ln (J)]}
$$

将求迹展开得

$$
\frac{d \det(J)}{dt} 
= \det(J)  \sum_{k,l} \frac{\partial \dot{x}_t^k}{\partial x_0^l} \frac{\partial x_0^l}{\partial x_t^k}
$$

利用求导的链式法则得

$$
\frac{d \det(J)}{dt} 
= \det(J)  \sum_{k} \frac{\partial \dot{x}_t^k}{\partial x_t^k} 
$$

根据相空间体积的不可压缩性

$$
\sum_{k} \frac{\partial \dot{x}_t^k}{\partial x_t^k} = 0 
$$ 

可知

$$
\frac{d \det(J)}{dt} = 0.
$$

也就是说，雅可比行列式不随时间变化。将时刻 $t$ 取 0 可知它是个常数：

$$
\det(J) = 1
$$

于是，我们就证明了刘伟尔定理：

$$
d x_t = dx_0
$$

### 哈密顿体系运动方程的数值积分

一个一般的物理量可以表达为相空间坐标的函数

$$
A=A(q_{\alpha}, p_{\alpha})
$$

我们求它的时间导数

$$
\frac{dA}{dt}= \sum_{\alpha} \left( \frac{\partial A}{\partial q _{\alpha} } \dot{q} _{\alpha} 
+\frac{\partial A}{\partial p _{\alpha} } \dot{p} _{\alpha} \right)
$$

利用哈密顿正则方程，可得

$$
\frac{dA}{dt}= \sum_{\alpha} \left( \frac{\partial A}{\partial q _{\alpha} } 
\frac{\partial H}{\partial p _{\alpha} }
-\frac{\partial A}{\partial p _{\alpha} } 
\frac{\partial H}{\partial q _{\alpha} } \right)
$$

定义任意两个物理量之间的泊松括号

$$
\lbrace A,B \rbrace = \sum_{\alpha} \left( \frac{\partial A}{\partial q _{\alpha} } 
\frac{\partial B}{\partial p _{\alpha} }
-\frac{\partial A}{\partial p _{\alpha} } 
\frac{\partial B}{\partial q _{\alpha} } \right)
$$

我们有

$$
\frac{dA}{dt}= \lbrace A,H \rbrace.
$$

如果将物理量 $A$ 取为哈密顿量本身，则有

$$
\frac{dH}{dt} = \lbrace H,H \rbrace = 0.
$$

由此可见，由哈密顿描述的力学体系的能量是守恒的。

一个物理量 $A$ 和哈密顿量之间的泊松括号运算也常用刘维尔算符表示

$$
\frac{dA}{dt} = \lbrace A,H \rbrace \equiv iLA
$$

该方程的形式解为

$$
A(t)=e^{iL}A(t=0)
$$

可以证明，刘维尔算符是厄米算符，故上述指数算符是幺正算符，称为经典演化算符（学过量子力学的读者会发现该指数算符和量子演化算符类似）。

将物理量 $A$ 取为相空间坐标，则有

$$
x_t = e^{iL} x_0
$$

该式只是形式解；我们通常是无法解析求解经典演化算符的。这也就是要研究近似的数值计算方法的根本原因所在。

为了推导近似的计算方法，我们首先注意到，可将刘维尔算符写成：


$$
iL = \sum_{\alpha} \left( \frac{\partial }{\partial q _{\alpha} } 
\frac{\partial H}{\partial p _{\alpha} }
-\frac{\partial }{\partial p _{\alpha} } 
\frac{\partial H}{\partial q _{\alpha} } \right)
$$

我们先将其写成两部分的和：

$$
iL = iL_1 + iL_2
$$

$$
iL_1 = \sum_{\alpha} 
\left( 
\frac{\partial }{\partial q _{\alpha} } 
\frac{\partial H}{\partial p _{\alpha} } 
\right)
$$

$$
iL = \sum_{\alpha} 
\left( 
-\frac{\partial }{\partial p _{\alpha} } 
\frac{\partial H}{\partial q _{\alpha} } 
\right)
$$

以上两个刘维尔算符是非对易的，即

$$
iL_1 iL_2 - iL_2 iL_1 \equiv [iL_1, iL_2] \neq 0
$$

这可由一维简谐振子很容易地验证。

因为以上两个刘维尔算符不对易，我们就不能将经典演化算符分开，即：

$$
e^{iL_1 + iL_2} \neq e^{iL_1} e^{iL_2}
$$

实际上，两个部分的时间演化算符 $e^{iL_1}$ 和 $e^{iL_2}$ 都是可精确求解的。那么我就希望找到一种方法，使得可以用 $e^{iL_1}$ 和 $e^{iL_2}$ 近似地表达 $e^{iL_1 + iL_2}$ 。对称 Trotter 定理（引用）提供了一种近似方法：

$$
e^{A + B} = \lim_{P \to \infty} \left[ e^{B/2P} e^{A/P} e^{B/2P} \right]^P
$$


未完待续。


The above velocity-Verlet integrator can be derived by finite-difference method (Taylor series expansion), but a more general method, which can be generalized to more sophisticated situations, is the classical time-evolution operator approach, or the Liouville operator approach \cite{tuckerman2010}. In this approach, the time-evolution of a classical system by one step can be formally expressed as
\begin{equation}
\left(
\begin{array}{c}
\vect{r}_i(t+\Delta t) \\
\vect{p}_i(t+\Delta t)
\end{array}
\right) =
e^{iL\Delta t}
\left(
\begin{array}{c}
\vect{r}_i(t) \\
\vect{p}_i(t)
\end{array}
\right),
\end{equation}
where $\vect{p}_i$ is the momentum of particle $i$ and
$e^{iL\Delta t}$ is called the classical evolution operator, which is the classical counterpart of the quantum evolution operator. The operator $iL$ in the exponent of the evolution operator is called the Liouville operator and is defined by
\begin{equation}
iL (\text{anything}) = \{\text{anything}, H\} \equiv
\sum_{i=1}^N
\left(
\frac{\partial H}{\partial \vect{p}_i} \cdot
\frac{\partial  }{\partial \vect{r}_i}  -
\frac{\partial H}{\partial \vect{r}_i} \cdot
\frac{\partial  }{\partial \vect{p}_i}
\right) (\text{anything}).
\end{equation}
Here, $H$ is the Hamiltonian of the system. Because
\begin{equation}
\frac{\partial H}{\partial \vect{p}_i} =
\frac{\vect{p}_i}{m_i} ~ \text{and} ~
-\frac{\partial H}{\partial \vect{r}_i} =
\vect{F}_i,
\end{equation}
we have
\begin{equation}
iL = iL_1 + iL_2,
\end{equation}
\begin{equation}
iL_1 = \sum_{i=1}^N
\frac{\vect{p}_i}{m_i} \cdot
\frac{\partial }{\partial \vect{r}_i},
\end{equation}
\begin{equation}
iL_2 = \sum_{i=1}^N
\vect{F}_i \cdot
\frac{\partial }{\partial \vect{p}_i}.
\end{equation}
Here, we have divided the Liouville operator into two parts. In general, $iL_1$ and $iL_2$ do not commute, and therefore
$e^{iL \Delta t} \neq e^{iL_1 \Delta t}e^{iL_2 \Delta t}$. However, there is an important theorem called the Trotter theorem, which can be used to derive the following approximation:
\begin{equation}
e^{iL \Delta t} \approx
e^{iL_2 \Delta t/2} e^{iL_1 \Delta t}e^{iL_2 \Delta t/2}.
\end{equation}
Now, we can express the one-step integration as
\begin{equation}
\left(
\begin{array}{c}
\vect{r}_i(t+\Delta t) \\
\vect{p}_i(t+\Delta t)
\end{array}
\right) \approx
e^{iL_2 \Delta t/2} e^{iL_1 \Delta t}e^{iL_2 \Delta t/2}
\left(
\begin{array}{c}
\vect{r}_i(t) \\
\vect{p}_i(t)
\end{array}
\right).
\end{equation}
To make further derivations, we note that for an arbitrary constant $c$, we have
\begin{equation}
e^{c \frac{\partial}{\partial x}} x = x+c.
\end{equation}
Applying this identity to the right most operator in the above equation, we have
\begin{equation}
\left(
\begin{array}{c}
\vect{r}_i(t+\Delta t) \\
\vect{p}_i(t+\Delta t)
\end{array}
\right) \approx
e^{iL_2 \Delta t/2} e^{iL_1 \Delta t}
\left(
\begin{array}{c}
\vect{r}_i(t) \\
\vect{p}_i(t) + \frac{\Delta t}{2} \vect{F}_i(t)
\end{array}
\right).
\end{equation}
Then, applying the operator $e^{iL_1 \Delta t}$, we have
\begin{equation}
\left(
\begin{array}{c}
\vect{r}_i(t+\Delta t) \\
\vect{p}_i(t+\Delta t)
\end{array}
\right) \approx
e^{iL_2 \Delta t/2}
\left(
\begin{array}{c}
\vect{r}_i(t) + \Delta t \frac{\vect{p}_i(t) + \frac{\Delta t}{2} \vect{F}_i(t)}{m_i} \\
\vect{p}_i(t) + \frac{\Delta t}{2} \vect{F}_i(t)
\end{array}
\right).
\end{equation}
Last, applying the remaining operator $e^{iL_2 \Delta t/2}$, we have
\begin{equation}
\left(
\begin{array}{c}
\vect{r}_i(t+\Delta t) \\
\vect{p}_i(t+\Delta t)
\end{array}
\right) \approx
\left(
\begin{array}{c}
\vect{r}_i(t) + \Delta t \frac{\vect{p}_i(t) + \frac{\Delta t}{2} \vect{F}_i(t)}{m_i} \\
\vect{p}_i(t) +
\frac{\Delta t}{2} \vect{F}_i(t) +
\frac{\Delta t}{2} \vect{F}_i(t+\Delta t)
\end{array}
\right).
\end{equation}
It is clear that this equation is equivalent to Eqs. (\ref{equation:velocity-Verlet-velocity}) and (\ref{equation:velocity-Verlet-position}).

\begin{algorithm}[htb]
\caption{The whole time-stepping in the $NVE$ ensemble. }
\label{algorithm:integration_NVE}
\begin{algorithmic}[1]
\State update the velocities partially
\begin{equation}
\vect{v}_i \leftarrow \vect{v}_i + \frac{1}{2} \frac{\vect{F}_i}{m_i} \Delta t
\end{equation}
\State update the positions completely
\begin{equation}
\vect{r}_i \leftarrow \vect{r}_i + \vect{v}_i \Delta t
\end{equation}
\State update the forces
\begin{equation}
\vect{F}_i \leftarrow \vect{F}_i(\{\vect{r}_i\})
\end{equation}
\State complete updating the velocities
\begin{equation}
\vect{v}_i \leftarrow \vect{v}_i + \frac{1}{2} \frac{\vect{F}_i}{m_i} \Delta t
\end{equation}
 \end{algorithmic}
\end{algorithm}

We can see that in the velocity-Verlet integrator, the position updating can be done in one step, but the velocity updating can only be done by two steps, one before force updating and the other after it. Algorithm \ref{algorithm:integration_NVE} gives the pseudo code for the complete time-stepping in the $NVE$ ensemble, including force updating.
