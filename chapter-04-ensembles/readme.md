
# 《分子动力学模拟》第四章：控温与控压算法

本章介绍几个常用的控温与控压算法。

# Table of contents
- [控温算法](#控温算法)
  - [Berendsen控温算法](#Berendsen控温算法)
  - [BDP控温算法](#BDP控温算法)
  - [Nose-Hoover-chain控温算法](#Nose-Hoover-chain控温算法)
  - [郎之万控温算法](#郎之万控温算法)
- [控压算法](#控压算法)
  - [Berendsen控压算法](#Berendsen控压算法)
  - [SCR控压算法](#SCR控压算法)
  - [MTTK控压算法](#MTTK控压算法)

## 控温算法

### Berendsen控温算法

在 Berendsen 控温算法中，用如下方式对原子速度进行变换：

$$
\vec{v} _i^{\text{scaled}} = \vec{v} _i \sqrt{1 + \frac{\Delta t}{\tau _T}  \left(\frac{T_0}{T} - 1\right)}.
$$

此处, $\tau_T$ 是一个时间参数, $T_0$ 是目标温度,  $T$ 当前根据速速计算出来的瞬时温度。一般来说 $\tau_T$ 取 100 倍步长比较合适。

如果 $\tau_T=\Delta t$, 那么上述公式变为简单的速度重标

$$
\vec{v} _i^{\text{scaled}} = \vec{v} _i \sqrt{\frac{T_0}{T}}.
$$

结合文献讲一讲 Berendsen 控温算法的优缺点。优点是编程实现简单、一般情况下很稳定，适合从非平衡态到平衡态过渡期间的模拟，缺点是该算法并不给出真正的正则系综，而且可能会出现所谓的“飞冰”现象。

### Nose-Hoover-chain 控温算法

运动方程：

$$
\frac{d}{dt} \vec{r}_i = \frac{\vec{p}_i}{m_i},
$$

$$
\frac{d}{dt} \vec{p}_i = \vec{F}_i - \frac{\pi_0}{Q_0} \vec{p}_i,
$$

$$
\frac{d}{dt} \eta_k = \frac{\pi_k}{Q_k} ~(k = 0, 1, \cdots, M-1),
$$

$$
\frac{d}{dt} \pi_0 =
2\left(
\sum_i \frac{\vec{p}_i^2}{2m_i} - dN\frac{k_BT}{2}
\right)
- \frac{\pi_1}{Q_1} \pi_0,
$$

$$
\frac{d}{dt} \pi_k =
2\left( \frac{\pi_{k-1}^2}{2Q_{k-1}} - \frac{k_BT}{2} \right)
- \frac{\pi_{k+1}}{Q_{k+1}} \pi_{k} ~(k = 1, 2, \cdots, M-2),
$$

$$
\frac{d}{dt} \pi_{M-1} =
2\left( \frac{\pi_{M-2}^2}{2Q_{M-2}} - \frac{k_BT}{2} \right).
$$

The optimal choice  for the thermostat masses is

$$
Q_0 = dNk_BT\tau^2,
$$

$$
Q_k = k_BT\tau^2 ~(k = 1, 2, \cdots, M-1),
$$

where $\tau$ is a time parameter, whose value is usually chosen by try and error in practice. A good choice is $\tau = 100 \Delta t$, where $\Delta t$ is the time step for integration.

总刘维尔算符为

$$
iL = iL_1 + iL_2 + iL_T,
$$

$$
iL_1 = \sum_{i=1}^N
\frac{\vec{p}_i}{m_i} \cdot
\frac{\partial }{\partial \vec{r}_i},
$$

$$
iL_2 = \sum_{i=1}^N
\vec{F}_i \cdot
\frac{\partial }{\partial \vec{p}_i}.
$$

$$
iL_T =
\sum_{k=0}^{M-1} \frac{\pi_k}{Q_k}\frac{\partial}{\partial \eta_k} +
\sum_{k=0}^{M-2} \left(G_k - \frac{\pi_{k+1}}{Q_{k+1}}\pi_k\right)
\frac{\partial}{\partial \pi_k}                               +
G_{M-1} \frac{\partial}{\partial \pi_{M-1}}                        -
\sum_{i=0}^{N-1}
\frac{\pi_0}{Q_0} \vec{p}_i \cdot \frac{\partial}{\partial \vec{p}_i}.
$$

也就是说，和 NVE 系综相比 （第一章）， NVT 系综的刘维尔算符多了和热浴变量有关的一项 $iL_T$ 。

The total time-evolution operator $e^{iL\Delta t} $ for one step can be factorized using the Trotter theorem as in the case of the $NVE$ ensemble:


$$
e^{iL} \approx
e^{iL_T\Delta t/2}
e^{iL_2\Delta t/2}
e^{iL_1\Delta t}
e^{iL_2\Delta t/2}
e^{iL_T\Delta t/2}.
$$

Comparing this with the factorization in the $NVE$ ensemble, we see that we only need to apply the operator $e^{iL_T\Delta t/2}$ before and after applying the usual velocity-Verlet integrator in the $NVE$ ensemble.

The operator $e^{iL_T\Delta t/2}$ can be further factorized into some elementary factors using the Trotter theorem. First, we define the following decomposition of the operator $iL_T$:

$$
iL_T =  iL_{T1} + iL_{T2} + iL_{T3},
$$

$$
iL_{T1} =
\sum_{k=0}^{M-1} \frac{\pi_k}{Q_k}\frac{\partial}{\partial \eta_k},
$$

$$
iL_{T2} =
\sum_{k=0}^{M-2} \left(G_k - \frac{\pi_{k+1}}{Q_{k+1}}\pi_k\right)
\frac{\partial}{\partial \pi_k}                               +
G_{M-1} \frac{\partial}{\partial \pi_{M-1}},
$$

$$
iL_{T3} = -\sum_{i=0}^{N-1}
\frac{\pi_0}{Q_0} \vec{p}_i \cdot \frac{\partial}{\partial \vec{p}_i}.
$$

We can then make the following factorization:

$$
e^{iL_T\Delta t/2} \approx
e^{iL_{T2}\Delta t/4}
e^{iL_{T3}\Delta t/2}
e^{iL_{T1}\Delta t/2}
e^{iL_{T2}\Delta t/4}.
$$

There are still a few terms in $iL_{T2}$ and we need to factorize $e^{iL_{T2}\Delta t/4}$ further. We can factorize the $e^{iL_{T2}\Delta t/4}$ term on the right of the above equation as


$$
e^{iL_{T2}\Delta t/4} \approx
\prod_{k=0}^{M-2}
\left(
e^{-\frac{\Delta t}{8} \frac{\pi_{k+1}}{Q_{k+1}}\pi_k \frac{\partial}{\partial \pi_k} }
e^{\frac{\Delta t}{4} G_k \frac{\partial}{\partial \pi_k} }
e^{-\frac{\Delta t}{8} \frac{\pi_{k+1}}{Q_{k+1}}\pi_k \frac{\partial}{\partial \pi_k} }
\right)
e^{\frac{\Delta t}{4} G_{M-1} \frac{\partial}{\partial \pi_{M-1}} }
$$

and correspondingly factorize that on the left as

$$
e^{iL_{T2}\Delta t/4} \approx
e^{\frac{\Delta t}{4} G_{M-1} \frac{\partial}{\partial \pi_{M-1}} }
\prod_{k=M-2}^{0}
\left(
e^{-\frac{\Delta t}{8} \frac{\pi_{k+1}}{Q_{k+1}}\pi_k \frac{\partial}{\partial \pi_k} }
e^{\frac{\Delta t}{4} G_k \frac{\partial}{\partial \pi_k} }
e^{-\frac{\Delta t}{8} \frac{\pi_{k+1}}{Q_{k+1}}\pi_k \frac{\partial}{\partial \pi_k} }
\right).
$$

It can be shown that the effect of the operator $e^{cx\frac{\partial}{\partial x}}$ on $x$ is to scale it by a factor of $e^c$:

$$
e^{cx\frac{\partial}{\partial x}} x = e^c x.
$$

Therefore, the effect of the operator $e^{iL_{T3}\Delta t/2}$ is to scale the momenta of all the particles in the system by a uniform factor $e^{-(\pi_0/Q_0)\Delta t/2}$. Although this operator appears in the factorization of $e^{iL_{T}\Delta t/2}$, its does not affect the thermostat variables. Therefore, when applying the operator $e^{iL_{T}\Delta t/2}$, we only need to update the variables related to the thermostats and save this factor for later use when we update the variables for the particles. In this way, the update for the thermostat variables and that for the particle variables are separated.

### BDP控温算法

### 郎之万控温算法

## 控压算法


控压算法是实现 NPT系综的关键。在 NPT 系综中，我们同时对温度和压强进行控制，
使得体系在平衡之后具有确定的温度平均值和压强平均值。
注意，我们不能说在 NPT 系综具有恒定的温度和压强，因为它们总是有涨落的。

### Berendsen 控压算法

The Berendsen barostat [Berendsen1984]_ is used with the Berendsen thermostat discussed above.
The barostat scales the box and positions as follows:

.. math::

   \left(
   \begin{array}{ccc}
   a_x^{\rm scaled} & b_x^{\rm scaled} & c_x^{\rm scaled} \\
   a_y^{\rm scaled} & b_y^{\rm scaled} & c_y^{\rm scaled} \\
   a_z^{\rm scaled} & b_z^{\rm scaled} & c_z^{\rm scaled} 
   \end{array}
   \right)
   =
   \left(
   \begin{array}{ccc}
   \mu_{xx} & \mu_{xy} & \mu_{xz} \\
   \mu_{yx} & \mu_{yy} & \mu_{yz} \\
   \mu_{zx} & \mu_{zy} & \mu_{zz} \\
   \end{array}
   \right)
   \left(
   \begin{array}{ccc}
   a_x & b_x & c_x \\
   a_y & b_y & c_y \\
   a_z & b_z & c_z 
   \end{array}
   \right)

and

.. math::

   \left(
   \begin{array}{c}
   x^{\rm scaled}_i \\
   y^{\rm scaled}_i \\
   z^{\rm scaled}_i
   \end{array}
   \right)
   =
   \left(
   \begin{array}{ccc}
   \mu_{xx} & \mu_{xy} & \mu_{xz} \\
   \mu_{yx} & \mu_{yy} & \mu_{yz} \\
   \mu_{zx} & \mu_{zy} & \mu_{zz} \\
   \end{array}
   \right)
   \left(
   \begin{array}{c}
   x_i \\
   y_i \\
   z_i
   \end{array}
   \right).

We consider the following three pressure-controlling conditions:

* *Condition 1*:
  The simulation box is *orthogonal* and only the hydrostatic pressure (trace of the pressure tensor) is controlled.
  The simulation box must be periodic in all three directions.
  The scaling matrix only has nonzero diagonal components and the diagonal components can be written as:

  .. math::

     \mu_{xx}=\mu_{yy}=\mu_{zz}= 1-\frac{\beta_{\rm hydro} \Delta t}{3 \tau_p} (p^{\rm target}_{\rm hydro} - p^{\rm instant}_{\rm hydro}).

* *Condition 2*:
  The simulation box is *orthogonal* and the three diagonal pressure components are controlled independently.
  The simulation box can be periodic or non-periodic in any of the three directions.
  Pressure is only controlled for periodic directions.
  The diagonal components of the scaling matrix can be written as:

  .. math::

     \mu_{xx}= 1-\frac{\beta_{xx} \Delta t}{3 \tau_p} (p^{\rm target}_{xx} - p^{\rm instant}_{xx}) \\
     \mu_{yy}= 1-\frac{\beta_{yy} \Delta t}{3 \tau_p} (p^{\rm target}_{yy} - p^{\rm instant}_{yy}) \\
     \mu_{zz}= 1-\frac{\beta_{zz} \Delta t}{3 \tau_p} (p^{\rm target}_{zz} - p^{\rm instant}_{zz}).

* *Condition 3*:
  The simulation box is *triclinic* and the 6 nonequivalent pressure components are controlled independently. The
  simulation box must be periodic in all three directions.
  The scaling matrix components are:

  .. math::

     \mu_{\alpha\beta}= 1-\frac{\beta_{\alpha\beta} \Delta t}{3 \tau_p} (p^{\rm target}_{\alpha\beta} - p^{\rm instant}_{\alpha\beta}).

The parameter :math:`\beta_{\alpha\beta}` is the isothermal compressibility, which is the inverse of the elastic modulus.
:math:`\Delta t` is the time step and :math:`\tau_p` is the pressure coupling time (relaxation time).

### SCR控压算法

Berendsen 方法并不给出真正的 NPT 系综。2020年，Bernetti 和 Bussi 对Berendsen系综进行了推广，提出了随机盒子重标 （stochastic cell rescaling ）控压算法，详见：[Mattia Bernetti and Giovanni Bussi, Pressure control using stochastic cell rescalin](https://doi.org/10.1063/5.0020514)。该控压算法于 BDP 控温算法结合可以实现 NPT 系综。

SCR 控压算法中的盒子变换矩阵是Berendsen控压算法中的盒子变换矩阵加上一个随机项。该随机项为：

$$
\mu^{\rm stochastic} _{\alpha\beta} 
= \sqrt{
\frac{1}{D _{\rm couple}}}
\sqrt{ 
\frac{\beta _{\alpha\beta} \Delta t}{3\tau_p} 
\frac{2k _{\rm B} T^{\rm target}}{V} 
} R _{\alpha\beta}.
$$

类似于 BDP 控温算法， $R_{\alpha\beta}$ 是平均值为零、方差为一的正则分布的随机数。

$D_{\rm couple}$ 是一个和盒子“耦合”有关的自由度。

### MTTK控压算法

从 Andersen 到 Parrinello-Rahman， 到 MTTK.
