/*----------------------------------------------------------------------------80
    Copyright 2016 Zheyong Fan <brucenju@gmail.com>

    This is a simple C program demonstrating the molecular dynamics of a single
    harmonic oscillator in the micro-canonical ensemble using the
    velocity-Verlet method.

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


#define NUMBER_OF_STEPS 10000000
#define SAMPLE_INTERVAL 100


int main(void)
{
    // parameters
    double dt = 0.01;          // time step
    double dt2 = dt * 0.5;     // half of time step
    double x = 1.0;            // initial position
    double v = 0.0;            // initial velocity
    double mass = 1.0;         // mass of the particle
    double k_spring = 1.0;     // spring constant
    double f = - k_spring * x; // initial force
  
    for (int step = 0; step < NUMBER_OF_STEPS; ++step) 
    {
        double inv = mass * v * v * 0.5 + k_spring * x * x * 0.5;

        if (step % SAMPLE_INTERVAL == 0)
        {
            printf("%25.15e%25.15e%25.15e\n", x, v, inv);
        } 

        v += dt2 * (f / mass); // Velocity-Verlet: 1st step
        x += dt * v;           // Velocity-Verlet: 1st step
        f = - k_spring * x;    // force update
        v += dt2 * (f / mass); // Velocity-Verlet: 2nd step
    }
}


