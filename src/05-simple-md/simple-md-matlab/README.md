# simple-md-matlab
A 100-line MATLAB implementation of molecular dynamics (MD) simulations with the Lennard-Jones (LJ) potential. 

## Features

* This is a small but working MD code, which is particularly useful to beginners.

* It contains two versions of the force-evaluation function, "find_force" and "find_force_vectorized".
  In the latter, the calculations of the pair forces between one atom and its neighboring atoms are vectorized.
  It turns out that the vectorized version is about 3 times as fast as the normal version. 
  
## File organizations

* There are two scripts:
  * test_md.m
  * test_md_vectorized.m

* The "test_md.m" script calls the "md" function in the "md.m" file.

* The "test_md_vectorized.m" script calls the "md_vectorized" function in the "md_vectorized.m" file.

* The "md" function calls the functions in the following files:
  * initialize_position.m
  * initialize_velocity.m
  * find_neighbor.m
  * find_force.m
  
* The "md_vectorized" function calls the functions in the following files:
  * initialize_position.m
  * initialize_velocity.m
  * find_neighbor.m
  * find_force_vectorized.m
  
## Unit system

* I use the following basic units:
  * Length: Angstrom
  * Mass: amu (atomic mass unit)
  * Energy: eV
  * Temperature: K
  
* Other units are then derived.

## Running the examples

* Run the "test_md.m" script or the "test_md_vectorized.m" sript in MATLAB, 
  which takes about 3 min and 1 min, respectively, in my laptop.
  
* Upon finished, two figures will show up. 
  The first figure shows the time evolution of the kinetic, potential, and total energies.
  The second figure shows the relative fluctuations of the total energy, 
  which are of the order of 1.0e-5, indicating good energy conservation.

## Contact

* Zheyong Fan: brucenju(at)gmail.com
