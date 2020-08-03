# ljmd.cpp
A `C++` implementation of simple molecular dynamics (MD) simulations with the Lennard-Jones (LJ) potential. 

## Compile

```
g++ -std=c++11 -O2 ljmd.cpp # using g++ from GCC
cl.exe /O2 ljmd.cpp # using cl.exe from MSVC
```

## Run

```
./ljmd nx Ne   # Linux
ljmd.exe nx Ne # Windows
# nx = number of cubic cells in x direction
# Ne = number of equilibrium steps = number of production steps
```
