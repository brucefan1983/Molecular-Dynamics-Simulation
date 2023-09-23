/*----------------------------------------------------------------------------80
    Copyright 2017 Zheyong Fan <brucenju@gmail.com>

    This is a simple C program demonstrating the NPT integrator based on the
	MTK equation.

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


#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#define NUMBER_OF_STEPS 10000000
#define SAMPLE_INTERVAL 1000
#define TWOPI           6.283185307179586


/*----------------------------------------------------------------------------80
    Reference:
    
    M. E. Tuckerman, Statistical Mechanics: Theory and Molecular Simulation
    Chapters 4 and 5.
    
    We try to reproduce Figure 5.3 in this book.
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


double sinh_factor(double x2)
{
    const double c2 = 1.0/(2.0*3.0);
    const double c4 = c2/(4.0*5.0);
    const double c6 = c4/(6.0*7.0);
    const double c8 = c6/(8.0*9.0);
    return 1.0 + x2*(c2 + x2*(c4 + x2*(c6 + x2*c8)));
}


double find_force(double x, double h) 
{
    double tmp = TWOPI / h;
    return - sin(x * tmp) / tmp;
}


double find_potential(double x, double h) 
{
    double tmp = TWOPI / h;
    return (1.0 - cos(x * tmp)) / tmp / tmp;
}


double find_pressure(double ek2, double x, double h) 
{
    double tmp = TWOPI / h;
    double potential = (1.0 - cos(x * tmp)) / tmp / tmp;
    return (ek2 - 2.0 * potential) / h;
}


int main(void)
{
	// box
    double h = 1.0;         // box length
    double pg = 0.0;        // box momentum
    double Wg = 18.0;       // box mass
    double dN_box = 1.0;    // # (degrees of freedom) of the box
    
    // particle
    double dN = 1.0;             // # (degrees of freedom) of the particle
    double dt = 0.005;           // time step
    double dt2 = dt * 0.5;       // half time step
    double kT = 1.0;             // target temperature
    double p0 = 1.0;             // target pressure
    double x = 0.2;              // position
    double v = 0.0;              // velocity
    double mass = 1.0;           // mass
    double f = find_force(x, h); // force

    // NHC for particle:
    const int M = 4;      
    double tau = 1; 
    double pos_eta[M];
    double vel_eta[M];
    double mas_eta[M];
    vel_eta[0] = vel_eta[2] = +1.0; 
    vel_eta[1] = vel_eta[3] = -1.0; 
    for (int i = 0; i < M; i++)
    {
        pos_eta[i] = 0.0; 
        mas_eta[i] = kT * tau * tau;
        if (i == 0) mas_eta[i] = dN * kT * tau * tau;
    }
    
    // NHC for box      
    const int M_baro = 4; 
    double tau_baro = 1;
    double pos_eta_baro[M_baro];
    double vel_eta_baro[M_baro];
    double mas_eta_baro[M_baro]; 
    vel_eta_baro[0] = vel_eta_baro[2] = +1.0; 
    vel_eta_baro[1] = vel_eta_baro[3] = -1.0;
    for (int i = 0; i < M_baro; i++)
    {
        pos_eta_baro[i] = 0.0;
        mas_eta_baro[i] = kT * tau_baro * tau_baro;
        if (i == 0) mas_eta_baro[i] = dN_box * kT * tau_baro * tau_baro;
    }

    // time evolution:
    FILE *fid = fopen("data.txt", "w");
    for (int step = 0; step < NUMBER_OF_STEPS; ++step) 
    {
        // output data:
        if (step % SAMPLE_INTERVAL == 0)
        {
            double inv = mass * v * v * 0.5 + find_potential(x, h); // particle
            inv += 0.5 * pg * pg / Wg + p0 * h; // box
            inv += kT * dN * pos_eta[0]; // NHC-particle
            for (int m = 1; m < M; m++) inv += kT * pos_eta[m];
            for (int m = 0; m < M; m++) 
			    inv += 0.5 * vel_eta[m] * vel_eta[m] / mas_eta[m];
            inv += kT * dN_box * pos_eta_baro[0]; // NHC-box
            for (int m = 1; m < M_baro; m++) inv += kT * pos_eta_baro[m];
            for (int m = 0; m < M_baro; m++)
                inv += 0.5*vel_eta_baro[m]*vel_eta_baro[m]/mas_eta_baro[m];
            fprintf(fid, "%25.15f%25.15f%25.15f%25.15f\n", x, h, inv, v);
        }
		
		double Ek2 = 0.0;
        double factor = 0.0; 
    
		// \exp^{iL_{T-b} dt2}
        Ek2 = pg * pg / Wg;
        factor = nhc
		(
		    M_baro, pos_eta_baro, vel_eta_baro, 
		    mas_eta_baro, Ek2, kT, dN_box, dt2
		);
        pg *= factor; 
        
        // \exp^{iL_{T-p} dt2}
        Ek2 = v * v * mass;
        factor = nhc(M, pos_eta, vel_eta, mas_eta, Ek2, kT, dN, dt2);
        v *= factor;
        
        // \exp^{iL_{g2} dt2}
        Ek2 = v * v * mass;
        pg += dt2 * (h * (find_pressure(Ek2, x, h) - p0) + Ek2);
        
        // \exp^{iL_{2} dt2}
        double x1 = -pg / Wg * dt2;
        double tmp = exp(x1);
        v = v * tmp * tmp + (dt2 / mass) * f * tmp * sinh_factor(x1*x1);
        
        //\exp^{iL_{1} dt}
        x1 = pg / Wg * dt2;
        tmp = exp(x1); 
        x =  x * tmp * tmp + dt * v * tmp * sinh_factor(x1*x1);
        if (x > h) x -= h;
        if (x < 0) x += h;
        
        //\exp^{iL_{1g} dt}
        h *= exp(pg / Wg * dt);

        // update force
        f = find_force(x, h); 
        
        // \exp^{iL_{2} dt2}
        x1 = -pg / Wg * dt2;
        tmp = exp(x1);
        v = v * tmp * tmp + (dt2 / mass) * f * tmp * sinh_factor(x1*x1);
        		 
        // \exp^{iL_{g2} dt2}
        Ek2 = v * v * mass;
        pg += dt2 * (h * (find_pressure(Ek2, x, h) - p0) + Ek2);

        // \exp^{iL_{T-p} dt2}
        Ek2 = v * v * mass;
        factor = nhc(M, pos_eta, vel_eta, mas_eta, Ek2, kT, dN, dt2);
        v *= factor;
        
        // \exp^{iL_{T-b} dt2}
        Ek2 = pg * pg / Wg;
        factor = nhc
		(
		    M_baro, pos_eta_baro, vel_eta_baro, 
		    mas_eta_baro, Ek2, kT, dN_box, dt2
		);
        pg *= factor;
    }
    
    fclose(fid);
    return 0;
}


