# heat-conductivity-hnemd
A 350-line C code for heat conductivity calculations using homogeneous nonequilibrium molecular dynamics (HNEMD).

## Features and Restrictions

* This is a small but working MD code, which is particularly useful to beginners who are interested in thermal conductivity calculations using molecular dynamics simulations.

* A matlab script is provided to show how to properly analyze the computed data. 

* Without proper modifications, this code can only calculate the lattice thermal conductivity of solid argon. 
  
* This is not intended as a code for doing original research, but is useful for teaching.

## Theory
The method used is called the homogeneous nonequilibrium molecular dynamics (HNEMD) method and was proposed by Evans in 1982 [D. J. Evans, Physics Letters A 91, 457 (1982)]. The simulation consists of three steps:
* The equilibrium stage: The system is equilibrated in an NVT ensemble (It does not matter much whether a true NVT ensemble is used or not here; so the simple velocity-rescaling method is used).
* The production stage: An external driving force Fe is applied on each particle (the total force on the whole system is then set to zero) and the total heat current J of the system is calculated in an NVT ensemble (again the simple velocity-rescaling method suffices).
* The post-processing stage: The heat conductivity is calculated as J/T/V/Fe, where T is the temperature and V is the volume of the system.
  
## File organizations

* The main code is a standalone C code:
  * kappa_hnemd.cpp

* The following Matlab script in the "20" directory can be used to analyze the data:
  * plot_kappa.m

## Unit system

* I use the following basic units:
  * Length: Angstrom
  * Mass: amu (atomic mass unit)
  * Energy: eV
  * Temperature: K
  
* Other units are then derived.

## Compile and run

* Although the code was written in pure C, I used a suffix of .cpp and one can therefore use g++ to compile it. In linux, it can be compiled using, e.g.
  * g++ -O3 kappa_hnemd.cpp -o kappa_hnemd
  
* Then, go to the directory "20" and type 
  * ../kappa_hnemd < input 
  
  If you want to do a few independent runs using one command, you can execute the "run" script I prepared:
  * ./run
  
* The simulation parameters are read in from the file "input". There are only a few parameters to play with. Check the beginning of the main() function in the C code to understand the meanings (and units) of the input parameters.
  
* After running the C code, a file named kappa.txt will be generated and one can run the Matlab script to analyze the results. Two figures will show up:
  * The first figure shows the running average of the diagonal thermal conductivity component k_xx, which should be about 1.5 W/mK (LJ argon at 20 K).
  * The second figure shows the running average of the off-diagonal thermal conductivity component k_yx, which should be close to zero (the system is isotropic).
  
## References

* For the HNEMD method, see D. J. Evans, "Homogeneous NEMD algorithm for thermal conductivity--Application of non-canonical linear response theory", Physics Letters A 91, 457 (1982).

* For benchmark results for Lennard-Jones argon, see table 1 of A.J.H. McGaughey and M. Kaviany, "Thermal conductivity decomposition and analysis using molecular dynamics simulations. Part I. Lennard-Jones argon", International Journal of Heat and Mass Transfer 47, 1783 (2004).

## Contact

* Zheyong Fan: brucenju(at)gmail.com; zheyong.fan(at)aalto.fi; zheyongfan(at)163.com

