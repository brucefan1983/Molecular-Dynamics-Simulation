# 第2章 热力学与统计力学基础回顾




## 经典力学基础

## 热力学基础

### 基本概念

一个**热力学系统**(thermodynamic system)是大量粒子的集合，常简称为\textbf{系统}(system)。系统之外的所有物质称为\textbf{环境}(environment)。

如果系统与环境既不能交换能量也不能交换物质，则称其为一个\textbf{孤立系统}(isolated system)；若可交换能量但不能交换物质，则称为\textbf{闭合系统}(closed system)；若即可交换能量也可交换物质，则称为\textbf{开放系统}(open system)。交换的能量可以是\textbf{热能}(thermal energy)，称为\textbf{热量}(heat)，也可以是\textbf{机械能}(mechanical energy)，称为\textbf{功}
(work)。

如果我们说一个系统达到了\textbf{热力学平衡}(thermodynamic equilibrium)，那么它的各个部分(叫做子系统)之间一定同时达到了(1)\textbf{热平衡}(thermal equilibrium)、(2)\textbf{力学平衡}(mechanical equilibrium) 以及(3)\textbf{扩散平衡}(diffusive equilibrium)。


热力学中一共有四条基本经验定律，其中在时间上最后提出但在逻辑上排在第一个的叫\textbf{热力学第零定律}(the zeroth law of thermodynamics)。 该定律说如果两个系统分别与第三个系统能达到热平衡，那么这两个系统也一定能达到热平衡。这说明处于热平衡的两个系统之间有一个相等的量。这个量就是热力学中最重要的物理量之一：
\textbf{温度}(temperature)。

通过规定两个特定系统的温度值，可以建立\textbf{温标}(temperature scale)。一般选择一个大气压下水的冰点和沸点作为这两个特定系统。Celsius 温标规定冰点的温度为~0 $^{\circ}$C，沸点的温度为~100 $^{\circ}$C，而且将这中间的温度均分为~100 份。Kelvin 温标规定冰点的温度为~273.15 K，沸点的温度为~373.15 K，而且将这中间的温度均分为~100 份。显然: $T (^{\circ}\text{C}) = T (\text{K}) -  273.15$.

虑孤立和闭合系统，达到热力学平衡时，可由三个物理量描述。它们是：压强~$p$，体积~$V$, 和温度~$T$。给定这三个量，就确定了系统的一个\textbf{状态}(state)。事实上，这三个量并非完全独立，而是由一个叫做\textbf{状态方程}(equation of state)的方程联系着的。该方程可写为
\begin{equation}
\label{equation:equation_of_state_general}
f(p, V, T) = 0.
\end{equation}
其中， $f$ 是一个特定的(三元)函数。当然我们可以将其中任何一个量写成另外两个量的(两元)函数，比如~$V=V(p,T)$。这样的两元函数在由~$p$、$V$、$T$ 这三个参数构成的三维空间中可以表示为一个曲面。系统的某一个状态就对应于该曲面上的一个点。

如果系统的状态发生了变化，我们称系统经历了一个\textbf{过程}(process)。如果系统在其状态发生变化时始终无限接近平衡态，那么系统经历的过程为\textbf{准静态过程}(quasi-static process)。显然，准静态过程对应于状态方程曲面上的一条曲线。

\textbf{理想气体}(ideal gas)是严格满足~\textbf{Boyle 定律}(即等温时~$p \propto 1/V $)，\textbf{Charles 定律}(即等容时 ~$p \propto T$) 以及~\textbf{Gay-Lussac 定律}(即等压时~$V \propto T$)的气体。通过这些定律以及标准状况(压强为一个大气压，温度为 0 摄氏度)下的摩尔体积($22.4\times10^{-3}$ m$^3$/mol)可以定出\textbf{理想气体的状态方程}：
\begin{equation}
pV = \nu R T.
\end{equation}
这里，$\nu$ 是\textbf{物质的量}(amount of substance)，$R=8.31$ J/mol/K 是\textbf{理想气体普适常量}(universal constant of ideal gas)。它与~Boltzmann 常量 $k_B = 1.38$ J/K 由下式联系：
\begin{equation}
R \equiv N_A k_B.
\end{equation}
这里，$N_A = 6.02\times 10^{23}$ mol$^{-1}$ 是~Avogadro 常量。注意到~$\nu N_A =N$ 就是系统的\textbf{粒子数}，我们可将理想气体的状态方程改写成
\begin{equation}
p V = N k_B T.
\end{equation}
再定义
\begin{equation}
n=N/V
\end{equation}
为系统的\textbf{数密度}(number density)，我们又可以将理想气体的状态方程改写成
\begin{equation}
p = n k_B T.
\end{equation}
对于单原子理想气体，可以证明：
\begin{equation}
p = \frac{1}{3}n m \overline{\vec{v}^2}.
\end{equation}
其中，$m$是一个原子的质量，$\overline{\vec{v}^2}$ 是原子速度平方的平均值。于是，我们有
\begin{equation}
\frac{3}{2}k_BT = \frac{1}{2} m \overline{\vec{v}^2}.
\end{equation}
这就是温度的微观意义之一：温度是分子平均平动动能的量度。


\textbf{热力学第一定律}(the first law of thermodynamics)是说，在一个过程中，如果环境对系统做了功，也传了热，那么系统的\textbf{内能}(internal energy)的增加量$\Delta E$就等于环境对系统做的功$W$ 和传给系统的热$Q$的和：
\begin{equation}
\Delta E = Q + W.
\end{equation}
当然，如果系统对环境做功，则约定~$W<0$；如果系统传给环境热量，则约定~$Q<0$。这个定律其实就是\textbf{能量守恒定律}在热力学系统中的应用。如果考虑的是无限小过程，那么可以将热力学第一定律写成
\begin{equation}
d E = d Q + d W.
\end{equation}
功是一个过程量，即它依赖于具体的过程。或者说，功不能描述系统的状态。这点与压强和温度不同。

系统在吸热时一般温度会升高(有反例)。与功类似，热量也是与过程有关的。因此，将一定质量的系统的温度升高一定的值所需要的热量依赖于系统所经历的过程。指定一个过程，就可以定义一个叫\textbf{热容}(heat capacity)的物理量：
\begin{equation}
C = \frac{d Q}{d T}.
\end{equation}
常见的两个过程是\textbf{等体过程}和\textbf{等压过程}，对应于\textbf{等体热容}和\textbf{等压热容}。

如果体积固定，系统与外界互不做功，热力学第一定律告诉我们$d Q = d E$。从而，等体热容可表达为
\begin{equation}
C_V =\left(\frac{~\delta Q}{dT}\right)_V
= \left(\frac{\partial E}{\partial T}\right)_V.
\end{equation}

如果压强固定(但体积不固定)，系统要对环境做功~$pdV$，此时热力学第一定律告诉我们
$d Q = d E + p d V = d (E + p V)$。
从而，等压热容可表达为
\begin{equation}
C_p = \left(\frac{~\delta Q}{dT}\right)_p
       = \left(\frac{\partial H}{\partial T}\right)_p.
\end{equation}
其中，我定义了一个新的类似于内能的热力学函数，\textbf{焓}(enthalpy)：
\begin{equation}
H = E + p V.
\end{equation}

综上，我们知：\textbf{等体过程中系统吸收的热量等于其内能的增加量；等压过程中系统吸收的热量等于其焓的增加量。}


\section{统计力学基础}

待写。两三页的介绍足够。


\section{进一步阅读}

Tuckerman 的统计物理

\section{习题}
