# lsqt-matlab
A 200-line MATLAB code for a linear-scaling quantum transport method

* This code can be used to obtain intrinsic electronic transport properties of large systems described by a real-space tight-binding Hamiltonian. It differs from the GPUQT code (https://github.com/brucefan1983/gpuqt) in two ways:
    * GPUQT is more than two orders of magnitude faster.
    * GPUQT accepts inputs for a general Hamiltonian, but the current code is only for 
the square lattice model with Anderson disorder.

* Purposes of this code
    * Help the readers who are interested in this linear-scaling quantum transport method to better understand it. 
    * The students in my course are asked to use the code to reproduce the results in sections 5.1 of Ref. [2] below. 
This is an optional course project.

* The major references for the implementation are (check the references cited in the papers below for original works on this method):
    * [1] Z. Fan, A. Uppstu, T. Siro, and A. Harju, Efficient linear-scaling quantum transport calculations on graphics processing units and applications on electron transport in graphene, Comput. Phys. Commun. 185, 28 (2014).
    * [2] Z. Fan, V. Vierimaa, and Ari Harju, GPUQT: An efficient linear-scaling quantum transport code fully implemented on graphics processing units, arXiv:1705.01387 [physics.comp-ph], to be published in Comput. Phys. Commun.

## File organizations

* This package consists of the following MATLAB functions (the file names are the same as the function names):
    * lsqt (calls find_H, create_state, find_dos, find_vac, and find_msd)            
    * find_dos (calls find_moments and chebyshev_summation)                  
    * find_vac (calls evolve, find_moments and chebyshev_summation)              
    * find_msd (calls evolve, evolvex, find_moments and chebyshev_summation)           
    * find_H (calls nothing)
    * create_state (calls nothing)
    * evolve (calls nothing)     
    * evolvex (calls nothing)
    * chebyshev_summation (calls nothing)
    * find_moments (calls nothing)

* lsqt is the driver function which the user will call.

## Inputs and outputs of the driver function (Nt is the number of time points; Ne is the number of energy points)
* Inputs of the driver function 
    * Nx: number of lattice points in the x direction
    * Ny: number of lattice points in the y direction
    * W: Anderson disorder strength
    * E: energy points (1*Ne matrix)
    * E_max: energy scaling factor
    * dt: time steps (Nt*1 matrix)
    * M: order of Chebyshev polynomial expansion
    * flag_vac: if this is 1, calculate the VAC and related quantities
    * flag_msd: if this is 1, calculate the MSD and related quantities
* Outputs of the driver function 
    * dos: the density of states (DOS), 1*Ne matrix
    * vac: the velocity autocorrelation (VAC), Nt*Ne matrix
    * msd: the mean square displacement (MSD), Nt*Ne matrix
    * sigma_vac: conductivity from the VAC, Nt*Ne matrix
    * sigma_msd: conductivity from the MSD, Nt*Ne matrix
* Unit system
    * basic units:
        * reduced Planck constant hbar = 1
        * elementary charge e = 1
        * energy unit gamma is choosen by the user
        * length unit a is choosen by the user
    * derived units:
        * time: hbar/gamma
        * DOS: 1/gamma/a^2
        * VAC: a^2*gamma^2/hbar^2
        * MSD: a^2
        * electrical conductivity: e^2/hbar

## Contact

* Zheyong Fan: brucenju(at)gmail.com

