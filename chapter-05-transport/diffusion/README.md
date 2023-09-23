# phonon-dos-matlab
A 20-line MATLAB function for calculating phonon density of states from velocity data

* This code can be used to obtain the phonon density of states (PDOS) from velocity data obtained from molecular dynamics simulation:
    * First, the normalized velocity autocorrelation function (VACF) is calculated from the velocity data.
    * Then, the PDOS is calculated from the normalized VACF by discrete cosine tranform.

* The orignianl reference for this method is
    * [1] J. M. Dickey and A. Paskin, Computer Simulation of the Lattice Dynamics of Solids, Phys. Rev. 188, 1407 (1969).
* For more backgrounds, see section 3.8 of the manual of the GPUMD code (https://github.com/brucefan1983/GPUMD).

## File organizations

* There is a single Matlab function in the `find_pdos.m` file:
    * `[vacf_output,pdos]=find_pdos(v_all,Nc,dt,omega)`

* Inputs of the function 
    * `v_all`: all the velocity data, which is an `N*3*Nf` matlab array. Here, `N` is the number of atoms, `3` refers to the space dimension, and `Nf` is the number of velocity frames. The units of the velocity data are not important.
    * `Nc`: number of correlation steps.
    * `dt`: time interval between two velocity frames, in units of ps. The maximum correlation time is then `dt*Nc`.
    * `omega`: the `Nw` phonon angular frequency points you want to consider in units of THz, which can be a `1*Nw` or `Nw*1` array.

* Outputs of the function 
    * `vacf_output`: the normalized (dimensionless) VACF, which is a `1*Nc` array
    * `pdos`: the PDOS in units of 1/THz, which is a `Nw*1` array.

## Contact

* Zheyong Fan: brucenju(at)gmail.com; zheyong.fan(at)aalto.fi; zheyongfan(at)163.com
