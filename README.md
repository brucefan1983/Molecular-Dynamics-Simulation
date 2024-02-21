# Source code and examples for my on-going book《分子动力学模拟》 (Molecular dynamics simulation)

To be published in the summer of 2024.

## Chapter 1: Classical Physics

* [A Python code demonstrating the velocity-Verlet integrator using 1D harmonic oscillator](chapter-01-classical_physics/python-harmonic-oscillator)

## Chapter 2: Simple MD code

* [A C++ implemenation of simple MD without using neighbor list](chapter-02-simple_md/cpp-simpleMD)

## Chapter 3: Simulation box and neighbor list

* [A C++ implemenation of linear-scaling MD with neighbor list](chapter-03-linear_md/cpp-linearMD)

## Chapter 4: Empirical potentials

* [A C++ implemenation of Tersoff potential, validating force calculation using finite difference](chapter-04-empirical_potentials/cpp-tersoff-validation)

* [An example of using GPUMD to validate Tersoff force against the C++ code above](chapter-04-empirical_potentials/gpumd-tersoff)

## Chapter 5: Machine-learned potentials

* [An example of using GPUMD to train a NEP model for silicon crystal](chapter-05-machine_learned_potentials/gpumd-nep-training-Si)

## 讨论群
* 读者可加入作者的 QQ 群： 887975816。
* 该群主要讨论分子动力学模拟中的算法和作者开发的 GPUMD 程序 (https://gpumd.org) 的使用。

