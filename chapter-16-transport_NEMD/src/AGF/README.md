# AGF-phonon-transport
A 200-line Matlab implementation of the Atomistic Green's function (AGF) method for phonon transport

## Features

* This is a simple but working code for calculating the phonon transmission in systems described by the Tersoff interatomic potential using the AGF (atomistic Green's function) method. This method is also called the NEGF (nonequilibrium Green's function) method.
* This code was developed based on the RGF-electron-matlab code: https://github.com/brucefan1983/RGF-electron-matlab. So you see that phonon transport is very similar to electron transport.
* This code was developed for the following paper: Zhen Li, Shiyun Xiong, Charles Sievers, Yue Hu, Zheyong Fan, Ning Wei, Hua Bao, Shunda Chen, Davide Donadio, and Tapio Ala-Nissila,
[Influence of Thermostatting on Nonequilibrium Molecular Dynamics Simulations of Heat Conduction in Solids](https://doi.org/10.1063/1.5132543),
J. Chem. Phys. **151**, 234105 (2019).

## Good references
* J.-S. Wang et al., "Quantum thermal transport in nanostructures", Eur. Phys. J. B 62, 381-404 (2008).
* S. Sadasivam, et al., "The atomistic Green's function method for interfacial phonon transport", Annual Review of Heat Transfer, 17, 89-145, (2014).
  
## File organizations

* Scripts for testing:
  * test_graphene.m: calculates and plots the phonon transmission in a graphene sheet

* The testing script calls
  * the function "find_r" to construct the simulation model
  * the function "find_T" to calculate the phonon transmission for the whole spectrum

* The function "find_T" calls
  * the "find_H" function to construct the dynamical matrix
  * the "find_Sigma" function to calculate the self-energies of the leads
  * the "find_T1" function to calculate the transmission for a single phonon frequency
  
* The "find_H" function calls 
  * the "find_neighbor" function to calculate the neighbor list
  * the "find_H1" function to calculate one element of the dynamical matrix
  
* The "find_H1" function calls 
  * the "find_E" function to calculate some energies

* The "find_Sigma" function calls 
  * the "find_g00" function to calculate the surface Green's function
  
## Unit system

* I use the following basic units in most parts of the code:
  * energy: eV
  * mass: amu
  * length: angstrom

## Uploaded results

* `graphene/`: results for graphene
* `SiNW/`: results for silicon nanowire (Obained by Shiyun Xiong using his code)

## Contact

* Zheyong Fan: brucenju(at)gmail.com
