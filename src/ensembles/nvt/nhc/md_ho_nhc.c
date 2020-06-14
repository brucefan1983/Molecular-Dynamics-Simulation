/*----------------------------------------------------------------------------80
    Copyright 2016 Zheyong Fan <brucenju@gmail.com>

    This is a simple C program demonstrating the molecular dynamics of a single
    harmonic oscillator in the canonical ensemble using the Nose-Hoover chain
    method.

    This is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You can find a copy of the GNU General Public License at
    <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------------*/


#include <stdio.h>
#include <math.h>


#define NUMBER_OF_STEPS 10000000
#define SAMPLE_INTERVAL 100


/*----------------------------------------------------------------------------80
    The Nose-Hoover chain algorithm

    See G. J. Martyna, M. E. Tuckerman, D. J. Tobias and M. L. Klein,
        "Explicit reversible integrators for extended systems dynamics",
        Molecular Physics, 87, 1117-1157 (1996).

    This paper has a lot of useful Fortran codes
------------------------------------------------------------------------------*/

double nhc
(
    int M, double* pos_eta, double *vel_eta, double *mas_eta,
    double Ek2, double kT, double dN, double dt2
)
{
    const double dt4 = dt2 * 0.5;
    const double dt8 = dt4 * 0.5;

    // update velocity of the last (M - 1) thermostat:
    double G = vel_eta[M - 2] * vel_eta[M - 2] / mas_eta[M - 2] - kT;
    vel_eta[M - 1] += dt4 * G;

    // update thermostat velocities from M - 2 to 0:
    for (int m = M - 2; m >= 0; m--)
    { 
        double tmp = exp(-dt8 * vel_eta[m + 1] / mas_eta[m + 1]);
        G = vel_eta[m - 1] * vel_eta[m - 1] / mas_eta[m - 1] - kT;
        if (m == 0) { G = Ek2 - dN  * kT; }
        vel_eta[m] = tmp * (tmp * vel_eta[m] + dt4 * G);   
    }

    // update thermostat positions from M - 1 to 0:
    for (int m = M - 1; m >= 0; m--)
    { 
        pos_eta[m] += dt2 * vel_eta[m] / mas_eta[m];  
    } 

    // compute the scale factor 
    double factor = exp(-dt2 * vel_eta[0] / mas_eta[0]); 

    // update thermostat velocities from 0 to M - 2:
    for (int m = 0; m < M - 1; m++)
    { 
        double tmp = exp(-dt8 * vel_eta[m + 1] / mas_eta[m + 1]);
        G = vel_eta[m - 1] * vel_eta[m - 1] / mas_eta[m - 1] - kT;
        if (m == 0) {G = Ek2 * factor * factor - dN * kT;}
        vel_eta[m] = tmp * (tmp * vel_eta[m] + dt4 * G);   
    }

    // update velocity of the last (M - 1) thermostat:
    G = vel_eta[M - 2] * vel_eta[M - 2] / mas_eta[M - 2] - kT;
    vel_eta[M - 1] += dt4 * G;

    return factor;
}


/*----------------------------------------------------------------------------80
    The main function
------------------------------------------------------------------------------*/


int main(void)
{
    // parameters
    double dN = 1.0; // one particle in one dimension here
    double dt = 0.001;
    double dt2 = dt * 0.5;
    double kT = 1.0; // temperature
    double x = 1.0; 
    double v = 0.0;  
    double mass = 1.0;
    double k_spring = 1.0;
    double f = - k_spring * x; // Harmonic oscillator
    const int M = 4;       // a good choice
    double tau = dt * 100; // a good choice (larger tau gives looser coupling)
    double pos_eta[M];
    double vel_eta[M];
    double mas_eta[M];

    vel_eta[0] = vel_eta[2] = +1.0; // good choice (but not important)
    vel_eta[1] = vel_eta[3] = -1.0; // good choice (but not important)

    for (int i = 0; i < M; i++)
    {
        pos_eta[i] = 0.0; // good choice (but not important)
        mas_eta[i] = kT * tau * tau;
        if (i == 0)
        {
            mas_eta[i] = dN * kT * tau * tau;
        }
    }
 
    double Ek2 = 0.0;
    double factor = 0.0; 

    for (int step = 0; step < NUMBER_OF_STEPS; ++step) 
    {
        // The invariant is not only the particle Hamiltonian
        double inv = mass * v * v * 0.5 + k_spring * x * x * 0.5;  
        inv += kT * dN * pos_eta[0];
        for (int m = 1; m < M; m++)
        {
            inv += kT * pos_eta[m];
        }
        for (int m = 0; m < M; m++)
        {
            inv += 0.5 * vel_eta[m] * vel_eta[m] / mas_eta[m];
        }

        if (step % SAMPLE_INTERVAL == 0)
        {
            printf("%25.15f%25.15f%25.15f\n", x, v, inv);
        } 

        // The first half of the thermostat integration
        Ek2 = v * v * mass;
        factor = nhc(M, pos_eta, vel_eta, mas_eta, Ek2, kT, dN, dt2);
        v *= factor;

        // The Velocity-Verlet integration
        v += dt2 * (f/mass);
        x += dt * v;
        f = - k_spring * x;
        v += dt2 * (f/mass);

        // The second half of the thermostat integration
        Ek2 = v * v * mass;
        factor = nhc(M, pos_eta, vel_eta, mas_eta, Ek2, kT, dN, dt2);
        v *= factor;
    }
}




