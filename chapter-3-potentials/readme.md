\documentclass{ctexart}

\setCJKmainfont{BabelStone Han}

\setlength{\oddsidemargin}{6mm}
\setlength{\evensidemargin}{6mm}
\setlength{\topmargin}{12mm}
\textwidth 145mm
\textheight 195mm
\parindent=2\ccwd

\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{bm}
\usepackage{dcolumn}
\usepackage{glossaries}
\usepackage{graphicx}
\usepackage{multirow}
\usepackage{fancyvrb}
\usepackage[linesnumbered,ruled]{algorithm2e}
\usepackage{siunitx}

\usepackage[usenames,dvipsnames,svgnames]{xcolor}
\usepackage{hyperref}
\hypersetup{
    pdfnewwindow=true,      % links in new window
    colorlinks=true,        % false: boxed links; true: colored links
    linkcolor=Blue,         % color of internal links
    citecolor=Blue,         % color of links to bibliography
    filecolor=Blue,         % color of file links
    urlcolor=Blue           % color of external links
}

\usepackage{listings}	
\lstset{
    basicstyle=\ttfamily,
	language={[ISO]C++},			
	commentstyle=\color{orange},    
	keywordstyle=\color{purple},		
	stringstyle=\color{blue},		
	frame=single,			
	xleftmargin=1.5em,	
	xrightmargin=.5em,		
    aboveskip=2em,
    belowskip=2em,
	numbers=left,		
	numberstyle=\small,		
	showstringspaces=false, 
	breaklines=true,	
	columns=flexible,	
	tabsize=4				
}

\begin{document}
\tableofcontents

\begin{abstract}
《分子动力学模拟入门》第三章：多体势函数\\

本章介绍分子动力学模拟中势函数的一般性知识，重点讨论两个典型的多体势函数，包括 Embedded-atom-method (EAM) 势 和 Tersoff 势。其中，EAM 势广泛应用于金属材料，而 Tersoff 势广泛应用于半导体和绝缘体材料。
\end{abstract}


\section{经验势函数的一般性质}

一个经典多粒子系统的势能可以写为
\begin{equation}
U=U\left(\{\vec{r}_i\}_{i=1}^N\right).
\end{equation}
该势能应该在欧几里得群的操作下保持不变。据此可以证明系统势能仅仅依赖于两两粒子对之间的距离：
\begin{equation}
U=U\left(\{r_{ij}\}_{i<j}\right).
\end{equation}
如果可以进一步将体系的势能写成
\begin{equation}
\label{equation:U-two-body}
U= \frac{1}{2}\sum_{i}\sum_{j \neq i} U_{ij} \left(r_{ij} \right),
\end{equation}
其中，
\begin{equation}
U_{ij} \left(r_{ij} \right) = U_{ji} \left(r_{ji} \right),
\end{equation}
那么我们称该系统的相互作用势能为两体势。其中，$U_{ij} \left(r_{ij} \right)$ 代表粒子 $i$和 $j$ 之间的相互作用势能，仅仅依赖于两粒子的相对距离$r_{ij}$。两体势系统的势能也可以写成如下等价的形式：
\begin{equation}
U= \sum_{i}\sum_{j > i} U_{ij} \left(r_{ij} \right).
\end{equation}
如果一个体系的势能无法写成以上形式，那么我们称该势能为多体势。

原则上，一个材料体系中原子间的相互作用势是可以由量子力学计算出来的。相对而言，多体势比两体势更加接近量子力学计算的结果，故在各种材料体系中应用得较为成功。

\section{多体势中力的表达式}

为了推导多体势的一系列表达式，我们假设一个多体势系统的总能量可以写为各个粒子的能量之和
\begin{equation}
U = \sum_i U_i.
\end{equation}
其中，$U_i$ 称为粒子 $i$ 的能量，它依赖于各个从 $i$ 指向其它粒子的位置矢量差 $\{\vec{r}_{ij}\}_j$：
\begin{equation}
U_i = U_i\left(\{\vec{r}_{ij}\}_j\right).
\end{equation}
该表达式显然满足空间平移不变性，但我们还没有对其施加空间转动不变性。后面我们会看到，EAM 势和 Tersoff 势都满足这个假设，但它们都有额外的限制。以后我们还会看到，最近发展迅猛的机器学习势也满足这个假设。

从以上假设出发，可以推导出如下力的表达式：
\begin{equation}
\vec{F}_{i} = \sum_{j \neq i} \vec{F}_{ij};
\end{equation}
\begin{equation}
\vec{F}_{ij} = - \vec{F}_{ji} =
\frac{\partial U_{i}}{\partial \vec{r}_{ij}} -
\frac{\partial U_{j}}{\partial \vec{r}_{ji}} =
\frac{\partial \left(U_{i} + U_{j}\right) }{\partial \vec{r}_{ij}}.
\end{equation}
这里，
\begin{equation}
\partial U_{i}/\partial \vec{r}_{ij} =
\partial U_{i}/\partial x_{ij} \vec{e}_x +
\partial U_{i}/\partial y_{ij} \vec{e}_y +
\partial U_{i}/\partial z_{ij} \vec{e}_z
\end{equation}
以上结果由笔者于 2015 年推导出来 \cite{fan2015prb}，详细证明如下。

我们从保守力的定义出发。粒子 $i$ 的力等于体系总势能对粒子坐标的梯度的负值：
\begin{equation}
\vec{F}_{i} = - \frac{\partial U}{\partial \vec{r}_{i}}
\end{equation}
代入总能量表达式，得
\begin{equation}
\vec{F}_{i} = - \frac{\partial \sum_j U_j}{\partial \vec{r}_{i}}
\end{equation}
注意，为了避免混淆指标，上式中的求和不能写成原先的 $\sum_i U_i$，这是在推导公式时要特别注意的。接下来的任务就是推导 $\partial U_j/\partial \vec{r}_i$ 了。为此，我们注意到 $U_{j}$ 是所有 $\{\vec{r}_{jk}\}_k$ 的函数，于是有
\begin{equation}
    \frac{\partial U_j}{\partial \vec{r}_{i}} = \sum_k \frac{\partial U_j}{\partial \vec{r}_{jk}} \frac{\partial \vec{r}_{jk}}{\partial \vec{r}_{i}}
\end{equation}
因为
\begin{equation}
\frac{\partial \vec{r}_{jk}}{\partial \vec{r}_{i}} = \frac{\partial (\vec{r}_{k} - \vec{r}_{j})}{\partial \vec{r}_{i}} = \frac{\partial \vec{r}_{k}}{\partial \vec{r}_{i}} - \frac{\partial \vec{r}_{j}}{\partial \vec{r}_{i}} = \delta_{ki}-
\delta_{ji},
\end{equation}
所以有
\begin{equation}
    \frac{\partial U_j}{\partial \vec{r}_{i}} =  \frac{\partial U_j}{\partial \vec{r}_{ji}} - \sum_k \frac{\partial U_j}{\partial \vec{r}_{jk}} \delta_{ji}
\end{equation}
最终得
\begin{equation}
\vec{F}_{i} = - \sum_j \left(\frac{\partial U_j}{\partial \vec{r}_{ji}} - \sum_k \frac{\partial U_j}{\partial \vec{r}_{jk}} \delta_{ji}\right) = \sum_k \frac{\partial U_i}{\partial \vec{r}_{ik}} - \sum_j \frac{\partial U_j}{\partial \vec{r}_{ji}}
= \sum_j \left(\frac{\partial U_i}{\partial \vec{r}_{ij}} - \frac{\partial U_j}{\partial \vec{r}_{ji}} \right).
\end{equation}
证毕。

\section{EAM 势}

EAM 势由若干人同时提出 \cite{daw1984prb,finnis1984pma}。

原子 $i$ 的势能为
\begin{equation}
U_i = \frac{1}{2} \sum_{j\neq i} \phi(r_{ij}) + F (\rho_i).
\end{equation}
这里，含有 $\phi(r_{ij})$ 的部分是两体势，$F(\rho_i)$ 即为嵌入势。嵌入势是 $i$ 粒子处电子密度 $\rho_i$ 的函数。粒子 $i$ 所在点的电子密度是由它的邻居贡献的：
\begin{equation}
\rho_i = \sum_{j\neq i} f(r_{ij}).
\end{equation}


可以推导如下表达式：
\begin{equation}
\label{equation:eam_partial_force}
\frac{\partial U_i}{\partial \vec{r}_{ij}}
= \frac{1}{2}  \phi'(r_{ij})  \frac{\partial r_{ij}} {\partial \vec{r}_{ij}} +
F'(\rho_i)  f'(r_{ij}) \frac{\partial r_{ij}} {\partial \vec{r}_{ij}}.
\end{equation}

\section{Tersoff 势}

Tersoff 势有几个稍有不同的变体。我这里介绍 Tersoff 在 1989 年发表的一篇文章中使用的形式 \cite{tersoff1989prb}。为简单起见，我们考虑单种元素的势函数。

粒子$i$的势能可以写为：
\begin{equation}
U_i =  \frac{1}{2} \sum_{j \neq i} f_C(r_{ij}) \left[ f_R(r_{ij}) - b_{ij} f_A(r_{ij}) \right].
\end{equation}
其中，$f_{C}$是一个截断函数，当$r_{ij}<R$ 时取值为 1，当$r_{ij}>S$时取值为 0。在这之间，该函数为
\begin{equation}
f_{C}(r_{ij}) = \frac{1}{2}
\left[
1 + \cos \left( \pi \frac{r_{ij} - R_{ij}}{S_{ij} - R_{ij}} \right)
\right].
\end{equation}

排斥函数$f_{R}$和吸引函数$f_{A}$为
\begin{equation}
f_{R}(r) = A e^{-\lambda r_{ij}};
\end{equation}
\begin{equation}
f_{A}(r) = B e^{-\mu r_{ij}}.
\end{equation}

键序为
\begin{equation}
\label{equation:bij}
b_{ij} = \left(1 + \beta^{n} \zeta^{n}_{ij}\right)^{-\frac{1}{2n}},
\end{equation}
其中，
\begin{equation}
\zeta_{ij} = \sum_{k\neq i, j}f_C(r_{ik}) g_{ijk},
\end{equation}
\begin{equation}
g_{ijk} = 1 + \frac{c^2}{d^2} - \frac{c^2}{d^2+(h-\cos\theta_{ijk})^2}.
\end{equation}

在以上表达式中，有如下参数：$A$, $B$, $\lambda$, $\mu$, $\beta$, $n$, $c$, $d$, $h$, $R$, $S$。

未完待续。

\section{习题}

1. 证明公式 \ref{equation:eam_partial_force}.

\bibliographystyle{apalike}
\bibliography{ref}

\end{document}