/*****************************************************************************80
    Copyright 2016 Zheyong Fan <brucenju@gmail.com>

    This is the source code of MD_DIFFUSION, which can output velocity
    and position data in a molecular dynamics simulation. 
    The other two related post-processing codes, namely, VAC and MSD, can
    calculate the velocity auto-correlation (VAC) using the velocity data and
    calculate the mean square displacement (MSD) using the position data.

    MD_DIFFUSION is a free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MD_DIFFUSION is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You can find a copy of the GNU General Public License at
    <http://www.gnu.org/licenses/>.
*******************************************************************************/


/*----------------------------------------------------------------------------80
    Author: Zheyong Fan <brucenju@gmail.com>
------------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------80
The units used in this program are as follows:
    (1) The basic units are chosen to be
        (A) Mechanics
            -- energy: eV (electron volt)
            -- length: A (Angstrom)
            -- mass: amu (atomic mass unit)
        (B) Thermodynamics
            -- temperature: K (Kelvin)
    (2) Some derived units are
        (A) Mechanics
            -- time: A amu^(1/2) eV^(-1/2) ~ 10 fs
            -- velocity: eV^(1/2) amu^(-1/2) ~ 0.1 A fs^(-1)
            -- acceleration: eV A^(-1) amu^(-1) ~ 0.01 A fs^(-2)
            -- force: eV A^(-1)
        (B) Thermodynamics
            Boltzmann constant: K_B ~ 0.8625 * 10^(-4) eV K^(-1)
------------------------------------------------------------------------------*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define K_B                   8.625e-5 // Boltzmann's constant in my unit
#define TIME_UNIT_CONVERSION  1.018e+1 // fs <-> my unit


// minimum image convention
void apply_mic
(
    double lx, double ly, double lz, double lxh, double lyh, 
    double lzh, double *x12, double *y12, double *z12
)
{
    if (*x12 < - lxh) {*x12 += lx;} else if (*x12 > + lxh) {*x12 -= lx;}
    if (*y12 < - lyh) {*y12 += ly;} else if (*y12 > + lyh) {*y12 -= ly;}
    if (*z12 < - lzh) {*z12 += lz;} else if (*z12 > + lzh) {*z12 -= lz;}
}


// pull the atoms back to the simulation box (needed for simulating fluids)
void apply_pbc
(
    int N, double lx, double ly, double lz, double *x, double *y, double *z
)
{
    for (int n = 0; n < N; n++)
    {
        if (x[n] < 0) {x[n] += lx;} else if (x[n] > + lx) {x[n] -= lx;}
        if (y[n] < 0) {y[n] += ly;} else if (y[n] > + ly) {y[n] -= ly;}
        if (z[n] < 0) {z[n] += lz;} else if (z[n] > + lz) {z[n] -= lz;}
    }
}


// neighbor list construction (slow for large system, but ok for teaching)
void find_neighbor
(
    int N, int *NN, int *NL, double *x, double *y, double *z, 
    double lx, double ly, double lz, int MN, double cutoff
)              
{
    double lxh = lx * 0.5;
    double lyh = ly * 0.5;
    double lzh = lz * 0.5; 
    double cutoff_square = cutoff * cutoff;
    for (int n = 0; n < N; n++) {NN[n] = 0;}
    for (int n1 = 0; n1 < N - 1; n1++)
    {  
        for (int n2 = n1 + 1; n2 < N; n2++)
        {   
            double x12 = x[n2] - x[n1];
            double y12 = y[n2] - y[n1];
            double z12 = z[n2] - z[n1];
            apply_mic(lx, ly, lz, lxh, lyh, lzh, &x12, &y12, &z12);
            double  distance_square = x12 * x12 + y12 * y12 + z12 * z12;
            if (distance_square < cutoff_square)
            {        
                NL[n1 * MN + NN[n1]] = n2;
                NN[n1]++;
                NL[n2 * MN + NN[n2]] = n1;
                NN[n2]++;
            }
            if (NN[n1] > MN)
            {
                printf("Error: cutoff for neighbor list is too large.\n");
                exit(1);
            }
        }
    } 
}


// position initialization (FCC crystal)
void initialize_position 
(
    int n0, int nx, int ny, int nz, double ax, double ay, double az, 
    double *x, double *y, double *z
)
{
    double x0[4] = {0.0, 0.0, 0.5, 0.5}; // only for simple FCC lattice
    double y0[4] = {0.0, 0.5, 0.0, 0.5}; 
    double z0[4] = {0.0, 0.5, 0.5, 0.0};
    int n = 0;
    for (int ix = 0; ix < nx; ++ix)
    {
        for (int iy = 0; iy < ny; ++iy)
        {
            for (int iz = 0; iz < nz; ++iz)
            {
                for (int i = 0; i < n0; ++i)
                {
                    x[n] = (ix + x0[i]) * ax;
                    y[n] = (iy + y0[i]) * ay;
                    z[n] = (iz + z0[i]) * az;
                    n++;
                }
            }
        }
    }
} 


// velocity scaling (the simplest way of controlling temperature)
void scale_velocity
(int N, double T_0, double *m, double *vx, double *vy, double *vz)
{  
    double temperature = 0.0;
    for (int n = 0; n < N; ++n) 
    {
        double v2 = vx[n] * vx[n] + vy[n] * vy[n] + vz[n] * vz[n];     
        temperature += m[n] * v2; 
    }
    temperature /= 3.0 * K_B * N;
    double scale_factor = sqrt(T_0 / temperature);
    for (int n = 0; n < N; ++n)
    { 
        vx[n] *= scale_factor;
        vy[n] *= scale_factor;
        vz[n] *= scale_factor;
    }
}  
   

// velocity initialization (with zero total momentum)  
void initialize_velocity
(int N, double T_0, double *m, double *vx, double *vy, double *vz)
{
    double momentum_average[3] = {0.0, 0.0, 0.0};
    for (int n = 0; n < N; ++n)
    { 
        vx[n] = -1.0 + (rand() * 2.0) / RAND_MAX; 
        vy[n] = -1.0 + (rand() * 2.0) / RAND_MAX; 
        vz[n] = -1.0 + (rand() * 2.0) / RAND_MAX;    
        
        momentum_average[0] += m[n] * vx[n] / N;
        momentum_average[1] += m[n] * vy[n] / N;
        momentum_average[2] += m[n] * vz[n] / N;
    } 
    for (int n = 0; n < N; ++n) 
    { 
        vx[n] -= momentum_average[0] / m[n];
        vy[n] -= momentum_average[1] / m[n];
        vz[n] -= momentum_average[2] / m[n]; 
    }
    scale_velocity(N, T_0, m, vx, vy, vz);
}


// force evaluation for Lennard-Jones potential
void find_force
(
    int N, int *NN, int *NL, int MN, double lx, double ly, double lz,
    double *x, double *y, double *z, double *fx, double *fy, double *fz,
    double cutoff
)
{
    // precompute something (a trick to make the code faster)
    const double epsilon = 1.032e-2;
    const double sigma = 3.405;
    const double sigma_3 = sigma * sigma * sigma;
    const double sigma_6 = sigma_3 * sigma_3;
    const double sigma_12 = sigma_6 * sigma_6;
    const double factor_1 = 24.0 * epsilon * sigma_6; 
    const double factor_2 = 48.0 * epsilon * sigma_12;

    // initialize force before accumulation
    for (int n = 0; n < N; ++n) { fx[n]=fy[n]=fz[n]=0.0; }

    // do somthing before the loop (a trick to make the code faster)
    double cutoff_square = cutoff * cutoff;
    double lxh = lx * 0.5;
    double lyh = ly * 0.5;
    double lzh = lz * 0.5;

    for (int i = 0; i < N; ++i)
    {
        for (int k = 0; k < NN[i]; k++)
        {
            int j = NL[i * MN + k];
            if (j < i) { continue; } // will use Newton's 3rd law

            double x_ij = x[j] - x[i];
            double y_ij = y[j] - y[i];
            double z_ij = z[j] - z[i];
            apply_mic(lx, ly, lz, lxh, lyh, lzh, &x_ij, &y_ij, &z_ij);

            double r_2 = x_ij * x_ij + y_ij * y_ij + z_ij * z_ij;
            if (r_2 > cutoff_square) { continue; }

            double r_4  = r_2 * r_2;
            double r_8  = r_4 * r_4;
            double r_14 = r_2 * r_4 * r_8;
            double f_ij = factor_1 / r_8 - factor_2 / r_14;
            fx[i] += f_ij * x_ij; fx[j] -= f_ij * x_ij; // use Newton's 3rd law
            fy[i] += f_ij * y_ij; fy[j] -= f_ij * y_ij;
            fz[i] += f_ij * z_ij; fz[j] -= f_ij * z_ij;  
        }
    }
}


// velocity-Verlet method of integration
void integrate
(
    int N, double time_step, double *m, double *fx, double *fy, double *fz, 
    double *vx, double *vy, double *vz, double *x, double *y, double *z, 
    double *x_msd, double *y_msd, double *z_msd, int flag
)
{
    double time_step_half = time_step * 0.5;
    for (int n = 0; n < N; ++n)
    {
        double mass_inv = 1.0 / m[n];
        double ax = fx[n] * mass_inv;
        double ay = fy[n] * mass_inv;
        double az = fz[n] * mass_inv;
        vx[n] += ax * time_step_half;
        vy[n] += ay * time_step_half;
        vz[n] += az * time_step_half;
        if (flag == 1) 
        { 
            x[n] += vx[n] * time_step; 
            y[n] += vy[n] * time_step; 
            z[n] += vz[n] * time_step; 
            x_msd[n] += vx[n] * time_step; 
            y_msd[n] += vy[n] * time_step; 
            z_msd[n] += vz[n] * time_step; 
        }
    }
}


// calculate MSD from position data
void find_msd
(int N, int Nd, int Nc, double delta_tau, double *x, double *y, double *z)
{  
    FILE *fid_out = fopen("msd.txt", "w");
    int M = Nd - Nc;
    for (int nc = 0; nc < Nc; nc++)
    {
        double msd_x = 0.0; double msd_y = 0.0; double msd_z = 0.0;
        for (int m = 0; m < M; m++)
        {
            double sum_x = 0.0; double sum_y = 0.0; double sum_z = 0.0;
            for (int i = 0; i < N; i++)
            {
                double dx = (x[m * N + i] - x[(nc + 1 + m) * N + i]);
                double dy = (y[m * N + i] - y[(nc + 1 + m) * N + i]);
                double dz = (z[m * N + i] - z[(nc + 1 + m) * N + i]);
                sum_x += dx * dx; sum_y += dy * dy; sum_z += dz * dz;
            }
            msd_x += sum_x; msd_y += sum_y; msd_z += sum_z;
        }
        msd_x /= N * M; msd_y /= N * M; msd_z /= N * M;
        fprintf
        (
            fid_out, "%25.15e%25.15e%25.15e%25.15e\n", 
            nc * delta_tau, msd_x, msd_y, msd_z
        );
    }
    fclose(fid_out);
}


// calculate VAC from velocity data
void find_vac
(int N, int Nd, int Nc, double delta_tau, double *vx, double *vy, double *vz)
{   
    FILE *fid_out = fopen("vac.txt", "w");
    int M = Nd - Nc;
    for (int nc = 0; nc < Nc; nc++)
    {
        double vac_x = 0.0; double vac_y = 0.0; double vac_z = 0.0;
        for (int m = 0; m < M; m++)
        {
            for (int i = 0; i < N; i++)
            {
                vac_x += vx[m * N + i] * vx[(m + nc) * N + i]; 
                vac_y += vy[m * N + i] * vy[(m + nc) * N + i];
                vac_z += vz[m * N + i] * vz[(m + nc) * N + i];
            }
        }
        vac_x /= N * M; vac_y /= N * M; vac_z /= N * M;
        fprintf
        (
            fid_out, "%25.15e%25.15e%25.15e%25.15e\n", 
            nc * delta_tau, vac_x, vac_y, vac_z
        );
    }
    fclose(fid_out);
}


// the main function
int main(int argc, char *argv[])
{
    //srand(time(NULL));

    int nx = 4; // number of unit cells in the x-direction
    int ny = 4; // number of unit cells in the y-direction
    int nz = 4; // number of unit cells in the z-direction
    int n0 = 4; // number of particles in the unit cell
    int N = n0 * nx * ny * nz; // total number of particles

    int Ne = 100000;  // number of steps in the equilibration stage
    int Np = 100000;  // number of steps in the production stage
    int Ns = 10;      // sampling interval
    int Nd = Np / Ns; // number of position/velocity data
    int Nc = Nd / 10; // number of correlation data
    double time_step = 10.0 / TIME_UNIT_CONVERSION; // time step
    
    double T_0 = 1.475 * 119.8; // temperature prescribed
    double ax = 11.6449;        // lattice constant in the x direction
    double ay = ax;             // lattice constant in the y direction
    double az = ax;             // lattice constant in the z direction
    double lx = ax * nx;        // box length in the x direction
    double ly = ay * ny;        // box length in the y direction
    double lz = az * nz;        // box length in the z direction

    double rcn = 15.0;     // cutoff distance for neighbor list
    double rcf = 10.0;     // cutoff distance for force
    int neighbor_update_interval = 10; // interval of neighbor list update
    int MN = 100;                      // amximal number of neighbors 
    
    // memory for neighbor list
    int *NN = (int*) malloc(N * sizeof(int));
    int *NL = (int*) malloc(N * MN * sizeof(int));

    // major data for the particles
    double *m  = (double*) malloc(N * sizeof(double)); // mass
    double *x  = (double*) malloc(N * sizeof(double)); // position
    double *y  = (double*) malloc(N * sizeof(double));
    double *z  = (double*) malloc(N * sizeof(double));
    double *vx = (double*) malloc(N * sizeof(double)); // velocity
    double *vy = (double*) malloc(N * sizeof(double));
    double *vz = (double*) malloc(N * sizeof(double));
    double *fx = (double*) malloc(N * sizeof(double)); // force
    double *fy = (double*) malloc(N * sizeof(double));
    double *fz = (double*) malloc(N * sizeof(double));

    // copy of position data for msd calculations (do not apply PBC to them)
    double *x_msd = (double*) malloc(N * sizeof(double)); 
    double *y_msd = (double*) malloc(N * sizeof(double)); 
    double *z_msd = (double*) malloc(N * sizeof(double)); 

    // data used for correlation function calculations
    double *x_all = (double*) malloc(Nd * N * sizeof(double)); 
    double *y_all = (double*) malloc(Nd * N * sizeof(double)); 
    double *z_all = (double*) malloc(Nd * N * sizeof(double)); 

    double *vx_all = (double*) malloc(Nd * N * sizeof(double)); 
    double *vy_all = (double*) malloc(Nd * N * sizeof(double)); 
    double *vz_all = (double*) malloc(Nd * N * sizeof(double)); 
    
    // initialize mass, position, and velocity
    for (int n = 0; n < N; ++n) { m[n] = 40.0; } // mass for argon atom
    initialize_position(n0, nx, ny, nz, ax, ay, az, x, y, z);
    for (int n = 0; n < N; n++) // make a copy
    { 
        x_msd[n] = x[n]; 
        y_msd[n] = y[n];
        z_msd[n] = z[n]; 
    }
    initialize_velocity(N, T_0, m, vx, vy, vz);

    // initialize neighbor list and force
    find_neighbor(N, NN, NL, x, y, z, lx, ly, lz, MN, rcn);
    find_force(N, NN, NL, MN, lx, ly, lz, x, y, z, fx, fy, fz, rcf);
 
    // equilibration
    clock_t time_begin = clock();
    for (int step = 0; step < Ne; ++step)
    { 
        if (0 == step % neighbor_update_interval)
        {
            find_neighbor(N, NN, NL, x, y, z, lx, ly, lz, MN, rcn);
        }

        integrate
        (
            N, time_step, m, fx, fy, fz, vx, vy, vz, x, y, z, 
            x_msd, y_msd, z_msd, 1
        );

        find_force(N, NN, NL, MN, lx, ly, lz, x, y, z, fx, fy, fz, rcf);

        integrate
        (
            N, time_step, m, fx, fy, fz, vx, vy, vz, x, y, z, 
            x_msd, y_msd, z_msd, 2
        );

        scale_velocity(N, T_0, m, vx, vy, vz); // control temperature

        apply_pbc(N, lx, ly, lz, x, y, z); // needed for simulating fluids
    } 
    clock_t time_finish = clock();
    double time_used = (time_finish - time_begin) / (double) CLOCKS_PER_SEC;
    fprintf(stderr, "time use for equilibration = %f s\n", time_used); 

    // production
    time_begin = clock();
    for (int step = 0; step < Np; ++step)
    {  
        if (0 == step % neighbor_update_interval)
        {
            find_neighbor(N, NN, NL, x, y, z, lx, ly, lz, MN, rcn);
        }

        integrate
        (
            N, time_step, m, fx, fy, fz, vx, vy, vz, x, y, z, 
            x_msd, y_msd, z_msd, 1
        );

        find_force(N, NN, NL, MN, lx, ly, lz, x, y, z, fx, fy, fz, rcf);

        integrate
        (
            N, time_step, m, fx, fy, fz, vx, vy, vz, x, y, z, 
            x_msd, y_msd, z_msd, 2
        );

        apply_pbc(N, lx, ly, lz, x, y, z); // needed for simulating fluids

        // record data
        if (0 == step % Ns) 
        {
            int offset = (step / Ns) * N;
            for (int n = 0; n < N; n++)
            { 
                x_all[n + offset]  = x_msd[n]; // record x_msd, not x
                y_all[n + offset]  = y_msd[n];
                z_all[n + offset]  = z_msd[n];
                vx_all[n + offset] = vx[n];
                vy_all[n + offset] = vy[n];
                vz_all[n + offset] = vz[n];
            }
        }
    } 
    time_finish = clock();
    time_used = (time_finish - time_begin) / (double) CLOCKS_PER_SEC;
    fprintf(stderr, "time use for production = %f s\n", time_used); 

    // free some memory
    free(NN); free(NL); free(m);  free(x);  free(y);  free(z);
    free(vx); free(vy); free(vz); free(fx); free(fy); free(fz);
    free(x_msd);  free(y_msd);  free(z_msd);

    // calculate MSD and VAC
    find_msd(N, Nd, Nc, Ns * time_step, x_all, y_all, z_all);
    find_vac(N, Nd, Nc, Ns * time_step, vx_all, vy_all, vz_all);

    // free some other memory
    free(x_all);  free(y_all);  free(z_all);
    free(vx_all); free(vy_all); free(vz_all);

    return 0;
}



