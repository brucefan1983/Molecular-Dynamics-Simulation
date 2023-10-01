# MATLADY
A MATlab toolkit for LAttice DYnamics calculations based on empirical potentials

## Features

* This is a pedagogical code for teaching some basics of lattice dynamics.
  * Force constants are calculated by using finite difference from empirical potentials.
  * The dynamical matrices at a given set of k-points are constructed from the force constants.
  * From the eigenvalues and eigenvectors of the dynamics matrices, one can calculate many physical quantities. 
  * In the very first version, the code can only calculate the phonon dispersion for systems described by the Tersoff potential.
  * I plan to enrich the features little by little.
  
## File organizations

* Some scripts for testing:
  * test_Si.m: calculate the phonon dispersion of bulk silicon using the Tersoff-1989 potential;
  * test_graphene.m: calculate the phonon dispersion of 2D graphene using the Tersoff-2010 potential (with some modifications).

* The testing scripts call
  * the function "find_r" to build up the real space super cell;
  * the function "find_k" to build up the reciprocal space;
  * the driver function "matlady" to do the major calculations;
  * the "plot_dispersion" function to plot the phonon dispersion.

* The "matlady" function calls
  * the "find_neighbor" function to construct one neighbor list for force calculations and one for building up the dynamical matrix;
  * the "find_phi_all" function to calculate all the force constants that are needed;
     * the "find_phi_all" function calls the "find_phi_one" function
     * the "find_phi_one" function calls the "find_E" function
     * the "find_E" function calls the "find_E_tersoff" function
  * the "find_nu" function to construct the dynamical matrix and to calculate the phonon spectrum. 
  
## Unit system

* I use the following basic units:
  * Length: Angstrom
  * Mass: amu (atomic mass unit)
  * Energy: eV
  
* Other units are then derived.

## Running the examples

* Run the "test_Si.m" script to get the phonon dispersion for bulk silicon. This takes less than one second. 

* Run the "test_graphene.m" script to get the phonon dispersion for 2D graphene. This takes less than one second. 

## Contact

* Zheyong Fan: brucenju(at)gmail.com
