
# 《分子动力学模拟入门》第四章：控温与控压算法

本章介绍几个常用的控温与控压算法。

# Table of contents
- [控温算法](#控温算法)
  - [Berendsen控温算法](#Berendsen控温算法)
  - [Nose-Hoover-chain控温算法](#Nose-Hoover-chain控温算法)
  - [BDP控温算法](#BDP控温算法)
  - [郎之万控温算法](#郎之万控温算法)
- [控压算法](#控压算法)
  - [Berendsen控压算法](#Berendsen控压算法)
  - [SCR控压算法](#SCR控压算法)


## 控温算法

### Berendsen控温算法

Using the Berendsen thermostat, the integration algorithm in the $NVT$ ensemble only requires an extra scaling of all the velocity components. For the $NPT$ ensemble, the Berendsen barostat requires an extra scaling of positions and box lengths. The Berendsen thermostat and barosat are very suitable for equilibrating the system to a target temperature and pressure.

The velocities are scaled in the Berendsen thermostat in the following way:
\begin{equation}
\vec{v}_i^{\text{scaled}}
= \vec{v}_i
\sqrt{1 + \alpha_T \left(\frac{T_0}{T} - 1\right)}.
\end{equation}
Here, $\alpha_T$ is a dimensionless parameter, $T_0$ is the target temperature, and $T$ is the instant temperature calculated from the current velocities $\{ \vec{v}_i \}$. The parameter $\alpha_T$ should be positive and not larger than 1. When $\alpha_T=1$, the above formula reduces to the simple velocity-scaling formula:
\begin{equation}
\vec{v}_i^{\text{scaled}}
= \vec{v}_i \sqrt{\frac{T_0}{T}}.
\end{equation}
A smaller $\alpha_T$ represents a weaker coupling between the system and the thermostat. Practically, any value of $\alpha_T$ in the range of $0.001 \sim 1$ can be used.

### Nose-Hoover-chain控温算法

The Nos\'e-Hoover chain method  is more suitable for calculating equilibrium properties in a specific ensemble. In the current version of GPUMD, only the Nos\'e-Hoover chain thermostat is implemented. We hope to implement the Nos\'e-Hoover chain barostat in a future version.

The equations of motion in the Nos\'{e}-Hoover chain method are
\begin{equation}
\frac{d}{dt} \vec{r}_i = \frac{\vec{p}_i}{m_i},
\end{equation}
\begin{equation}
\frac{d}{dt} \vec{p}_i = \vec{F}_i - \frac{\pi_0}{Q_0} \vec{p}_i,
\end{equation}
\begin{equation}
\frac{d}{dt} \eta_k = \frac{\pi_k}{Q_k} ~(k = 0, 1, \cdots, M-1),
\end{equation}
\begin{equation}
\frac{d}{dt} \pi_0 =
2\left(
\sum_i \frac{\vec{p}_i^2}{2m_i} - dN\frac{k_BT}{2}
\right)
- \frac{\pi_1}{Q_1} \pi_0,
\end{equation}
\begin{equation}
\frac{d}{dt} \pi_k =
2\left( \frac{\pi_{k-1}^2}{2Q_{k-1}} - \frac{k_BT}{2} \right)
- \frac{\pi_{k+1}}{Q_{k+1}} \pi_{k} ~(k = 1, 2, \cdots, M-2),
\end{equation}
\begin{equation}
\frac{d}{dt} \pi_{M-1} =
2\left( \frac{\pi_{M-2}^2}{2Q_{M-2}} - \frac{k_BT}{2} \right).
\end{equation}
The optimal choice  for the thermostat masses is
\begin{equation}
Q_0 = dNk_BT\tau^2,
\end{equation}
\begin{equation}
Q_k = k_BT\tau^2 ~(k = 1, 2, \cdots, M-1),
\end{equation}
where $\tau$ is a time parameter, whose value is usually chosen by try and error in practice. A good choice is $\tau = 100 \Delta t$, where $\Delta t$ is the time step for integration.

An integration scheme for the $NVT$ ensemble using the Nos\'{e}-Hoover chain can also be formulated using the approach of the time-evolution operator. The total Liouville operator for the equations of motion in the Nos\'{e}-Hoover chain method is 
\begin{equation}
iL = iL_1 + iL_2 + iL_T,
\end{equation}
\begin{equation}
iL_1 = \sum_{i=1}^N
\frac{\vec{p}_i}{m_i} \cdot
\frac{\partial }{\partial \vec{r}_i},
\end{equation}
\begin{equation}
iL_2 = \sum_{i=1}^N
\vec{F}_i \cdot
\frac{\partial }{\partial \vec{p}_i}.
\end{equation}
\begin{equation}
iL_T =
\sum_{k=0}^{M-1} \frac{\pi_k}{Q_k}\frac{\partial}{\partial \eta_k} +
\sum_{k=0}^{M-2} \left(G_k - \frac{\pi_{k+1}}{Q_{k+1}}\pi_k\right)
     \frac{\partial}{\partial \pi_k}                               +
G_{M-1} \frac{\partial}{\partial \pi_{M-1}}                        -
\sum_{i=0}^{N-1}
\frac{\pi_0}{Q_0} \vec{p}_i \cdot \frac{\partial}{\partial \vec{p}_i}.
\end{equation}
That is, the Liouville operator for the $NVT$ ensemble contains an extra term $iL_T$ related to the thermostat variables, which is absent from that for the $NVE$ ensemble.

The total time-evolution operator $e^{iL\Delta t} $ for one step can be factorized using the Trotter theorem as in the case of the $NVE$ ensemble:
\begin{equation}
e^{iL} \approx
e^{iL_T\Delta t/2}
e^{iL_2\Delta t/2}
e^{iL_1\Delta t}
e^{iL_2\Delta t/2}
e^{iL_T\Delta t/2}.
\end{equation}
Comparing this with the factorization in the $NVE$ ensemble, we see that we only need to apply the operator $e^{iL_T\Delta t/2}$ before and after applying the usual velocity-Verlet integrator in the $NVE$ ensemble.

The operator $e^{iL_T\Delta t/2}$ can be further factorized into some elementary factors using the Trotter theorem. First, we define the following decomposition of the operator $iL_T$:
\begin{equation}
iL_T =  iL_{T1} + iL_{T2} + iL_{T3},
\end{equation}
\begin{equation}
iL_{T1} =
\sum_{k=0}^{M-1} \frac{\pi_k}{Q_k}\frac{\partial}{\partial \eta_k},
\end{equation}
\begin{equation}
iL_{T2} =
\sum_{k=0}^{M-2} \left(G_k - \frac{\pi_{k+1}}{Q_{k+1}}\pi_k\right)
     \frac{\partial}{\partial \pi_k}                               +
G_{M-1} \frac{\partial}{\partial \pi_{M-1}},
\end{equation}
\begin{equation}
iL_{T3} = -\sum_{i=0}^{N-1}
\frac{\pi_0}{Q_0} \vec{p}_i \cdot \frac{\partial}{\partial \vec{p}_i}.
\end{equation}
We can then make the following factorization:
\begin{equation}
e^{iL_T\Delta t/2} \approx
e^{iL_{T2}\Delta t/4}
e^{iL_{T3}\Delta t/2}
e^{iL_{T1}\Delta t/2}
e^{iL_{T2}\Delta t/4}.
\end{equation}
There are still a few terms in $iL_{T2}$ and we need to factorize $e^{iL_{T2}\Delta t/4}$ further. We can factorize the $e^{iL_{T2}\Delta t/4}$ term on the right of the above equation as
\begin{equation}
e^{iL_{T2}\Delta t/4} \approx
\prod_{k=0}^{M-2}
\left(
e^{-\frac{\Delta t}{8} \frac{\pi_{k+1}}{Q_{k+1}}\pi_k \frac{\partial}{\partial \pi_k} }
e^{\frac{\Delta t}{4} G_k \frac{\partial}{\partial \pi_k} }
e^{-\frac{\Delta t}{8} \frac{\pi_{k+1}}{Q_{k+1}}\pi_k \frac{\partial}{\partial \pi_k} }
\right)
e^{\frac{\Delta t}{4} G_{M-1} \frac{\partial}{\partial \pi_{M-1}} }
\end{equation}
and correspondingly factorize that on the left as
\begin{equation}
e^{iL_{T2}\Delta t/4} \approx
e^{\frac{\Delta t}{4} G_{M-1} \frac{\partial}{\partial \pi_{M-1}} }
\prod_{k=M-2}^{0}
\left(
e^{-\frac{\Delta t}{8} \frac{\pi_{k+1}}{Q_{k+1}}\pi_k \frac{\partial}{\partial \pi_k} }
e^{\frac{\Delta t}{4} G_k \frac{\partial}{\partial \pi_k} }
e^{-\frac{\Delta t}{8} \frac{\pi_{k+1}}{Q_{k+1}}\pi_k \frac{\partial}{\partial \pi_k} }
\right).
\end{equation}

It can be shown that the effect of the operator $e^{cx\frac{\partial}{\partial x}}$ on $x$ is to scale it by a factor of $e^c$:
\begin{equation}
e^{cx\frac{\partial}{\partial x}} x = e^c x.
\end{equation}
Therefore, the effect of the operator $e^{iL_{T3}\Delta t/2}$ is to scale the momenta of all the particles in the system by a uniform factor $e^{-(\pi_0/Q_0)\Delta t/2}$. Although this operator appears in the factorization of $e^{iL_{T}\Delta t/2}$, its does not affect the thermostat variables. Therefore, when applying the operator $e^{iL_{T}\Delta t/2}$, we only need to update the variables related to the thermostats and save this factor for later use when we update the variables for the particles. In this way, the update for the thermostat variables and that for the particle variables are separated.

### BDP控温算法

### 郎之万控温算法

## 控压算法


控压算法是实现 NPT系综的关键。在 NPT 系综中，我们同时对温度和压强进行控制，
使得体系在平衡之后具有确定的温度平均值和压强平均值。
注意，我们不能说在NPT系综具有恒定的温度和压强，因为它们总是有涨落的。

### Berendsen控压算法

In the Berendsen barostat algorithm, the particle positions and box length in a given direction are scaled if periodic boundary conditions are applied to that direction. The scaling of the positions reads
\begin{equation}
\vec{r}_i^{\text{scaled}}
= \vec{r}_i \left[ 1-\alpha_p (\vec{p}_0 - \vec{p}) \right].
\end{equation}
Here, $\alpha_p$ is a parameter and $\vec{p}_0$ ($\vec{p}$) is the target (instant) pressure in the three directions. The parameter $\alpha_p$ is not dimensionless, and it requires some try-and-error to find a good value of it for a given system. A harder/softer system requires a smaller/larger value of $\alpha_p$. In the unit system adopted by GPUMD, it is recommended that $\alpha_p = 10^{-4} \sim 10^{-2}$.
Only directions with periodic boundary conditions will be affected by the barostat.

### SCR控压算法
