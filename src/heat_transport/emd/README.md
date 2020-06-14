# heat-conductivity-emd
A 400-line C code for heat conductivity calculations using equilibrium molecular dynamics (EMD). 

## Features and Restrictions

* This is a small but working MD code, which is particularly useful to beginners 
  who are interested in thermal conductivity calculations using molecular dynamics simulations. 
  I am not aware of a shorter code which can do essentially the same calculations.

* A matlab script is provided to show how to properly analyze the computed data. 

* Without proper modifications, this code can only calculate the lattice thermal conductivity of solid argon. 
  
* This is not intended as a code for doing original research, but is useful for teaching.

## Theory
The method used is also called the Green-Kubo method, as the thermal conductivity 
is calculated according to a Green-Kubo relation. The simulation consists of three steps:
* The equilibrium stage: The system is equilibrated in an NVT ensemble
  (It does not matter much whether a true NVT ensemble is used or not here; so the simple velocity-rescaling method is used).
* The production stage: The total heat current of the system is calculated in equilibrium state (NVE ensemble) 
  with a given sampling interval.
* The post-processing stage: The heat current data are used to calculate the heat current 
  auto-correlation function, from which the running thermal conductivity is calculated 
  according to the Green-Kubo formula.
  
## File organizations

* The main code is a standalone C code:
  * kappa_emd.cpp

* The following Matlab script can be used to analyze the data:
  * plot_kappa.m

## Unit system

* I use the following basic units:
  * Length: Angstrom
  * Mass: amu (atomic mass unit)
  * Energy: eV
  * Temperature: K
  
* Other units are then derived.

## Compile and run

* Although the code was written in pure C, I used a suffix of .cpp 
  and one can therefore use g++ to compile it. In linux, it can be compiled using, e.g.
  * g++ -O3 kappa_emd.cpp
  
* Then, just run the executable:
  * ./a.out
  
* The simulation parameters are hard coded. Without modifying the code, one simulation takes about 20 seconds. 
  This method requires to do many independent simulations for a given set of parameters 
  to get statistically meaningful results. Different runs usually give different results 
  because the initial velocities are different from run to run.
  
* After running the C code one or more times, one can run the Matlab script to analyze the results. Two figures 
  will show up:
  * The first figure shows the normalized HCACF (heat current auto-correlation function) as a function of the 
    correlation time. It is conventional to show the normalized quantity here. I guess one reason is that this
    quantity does not have a nice dimension (unit).
  * The second figure shows the running thermal conductivity as a function of the correlation time. 
  * Note that it is not meaningful to do ensemble average for results obtained with different input parameters 
   (such as different temperatures).
  
* Beginners (including myself a few years ago) usually do not know how to report statistically meaningful results. 
  I recommend such readers check the last few lines in the Matlab script and 
  Section III.B of one of my recently published papers: [Phys. Rev. B 92, 094301 (2015)].

* The readers are encouraged to experiment with the many adjustable parameters defined in the beginning of the main function.

## Contact

* Zheyong Fan: brucenju(at)gmail.com; zheyong.fan(at)aalto.fi; zheyongfan(at)163.com
