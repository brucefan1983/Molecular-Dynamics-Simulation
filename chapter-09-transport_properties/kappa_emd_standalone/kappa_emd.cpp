#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define K_B                   8.617343e-5 // Boltzmann's constant  
#define TIME_UNIT_CONVERSION  1.018051e+1 // fs     <-> my natural unit
#define KAPPA_UNIT_CONVERSION 1.573769e+5 // W/(mK) <-> my natural unit


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


void find_force
(
    int N, int *NN, int *NL, int MN, double lx, double ly, double lz,
    double *x, double *y, double *z, double *fx, double *fy, double *fz,
    double *vx, double *vy, double *vz, double *hc
)
{
    const double epsilon = 1.032e-2;
    const double sigma = 3.405;
    const double cutoff = sigma * 3.0;
    const double cutoff_square = cutoff * cutoff;
    const double sigma_3 = sigma * sigma * sigma;
    const double sigma_6 = sigma_3 * sigma_3;
    const double sigma_12 = sigma_6 * sigma_6;
    const double factor_1 = 24.0 * epsilon * sigma_6; 
    const double factor_2 = 48.0 * epsilon * sigma_12;

    // initialize heat current and force
    hc[0] = hc[1] = hc[2] = 0.0; 
    for (int n = 0; n < N; ++n) { fx[n]=fy[n]=fz[n]=0.0; }

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
            double f_dot_v
                = x_ij*(vx[i]+vx[j])+y_ij*(vy[i]+vy[j])+z_ij*(vz[i]+vz[j]);
            f_dot_v *= f_ij * 0.5;
            hc[0] -= x_ij * f_dot_v; // calculate heat current
            hc[1] -= y_ij * f_dot_v;
            hc[2] -= z_ij * f_dot_v;
        }
    }
}


void integrate
(
    int N, double time_step, double *m, double *fx, double *fy, double *fz, 
    double *vx, double *vy, double *vz, double *x, double *y, double *z, 
    int flag
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
        }
    }
}


void find_hac
(
    int Nc, int M, double *hx, double *hy, double *hz, double *hac_x, 
    double *hac_y, double *hac_z 
)
{
    for (int nc = 0; nc < Nc; nc++) // loop over the correlation time points
    {
        for (int m = 0; m < M; m++) // loop over the time origins
        {
            hac_x[nc] += hx[m] * hx[m + nc]; 
            hac_y[nc] += hy[m] * hy[m + nc]; 
            hac_z[nc] += hz[m] * hz[m + nc]; 
        }
        hac_x[nc] /= M; hac_y[nc] /= M; hac_z[nc] /= M;
    }
}


static void find_rtc
(
    int Nc, double factor, double *hac_x, double *hac_y, double *hac_z,
    double *rtc_x, double *rtc_y, double *rtc_z
)
{
    for (int nc = 1; nc < Nc; nc++)  
    {
        rtc_x[nc] = rtc_x[nc - 1] + (hac_x[nc - 1] + hac_x[nc]) * factor;
        rtc_y[nc] = rtc_y[nc - 1] + (hac_y[nc - 1] + hac_y[nc]) * factor;
        rtc_z[nc] = rtc_z[nc - 1] + (hac_z[nc - 1] + hac_z[nc]) * factor;
    }
}


void find_hac_kappa
(
    int Nd, int Nc, double dt, double T_0, double V, 
    double *hx, double *hy, double *hz
)
{
    double dt_in_ps = dt * TIME_UNIT_CONVERSION / 1000.0; // ps
    int M = Nd - Nc; // number of time origins

    double *hac_x = (double *)malloc(sizeof(double) * Nc);
    double *hac_y = (double *)malloc(sizeof(double) * Nc);
    double *hac_z = (double *)malloc(sizeof(double) * Nc);
    double *rtc_x = (double *)malloc(sizeof(double) * Nc);
    double *rtc_y = (double *)malloc(sizeof(double) * Nc);
    double *rtc_z = (double *)malloc(sizeof(double) * Nc);
    for (int nc = 0; nc < Nc; nc++) {hac_x[nc] = hac_y[nc] = hac_z[nc] = 0.0;}
    for (int nc = 0; nc < Nc; nc++) {rtc_x[nc] = rtc_y[nc] = rtc_z[nc] = 0.0;}

    find_hac(Nc, M, hx, hy, hz, hac_x, hac_y, hac_z);
    double factor = dt * 0.5 *  KAPPA_UNIT_CONVERSION / (K_B * T_0 * T_0 * V);
    find_rtc(Nc, factor, hac_x, hac_y, hac_z, rtc_x, rtc_y, rtc_z);

    FILE *fid = fopen("kappa.txt", "a");
    for (int nc = 0; nc < Nc; nc++) 
    {
        fprintf
        (
            fid, "%25.15e%25.15e%25.15e%25.15e%25.15e%25.15e%25.15e\n", 
            nc * dt_in_ps, 
            hac_x[nc], hac_y[nc], hac_z[nc],
            rtc_x[nc], rtc_y[nc], rtc_z[nc]
        );
    }
    fclose(fid);

    free(hac_x); free(hac_y); free(hac_z);
    free(rtc_x); free(rtc_y); free(rtc_z);
}


int main(int argc, char *argv[])
{
    srand(time(NULL));
    int nx = 4; // number of unit cells in the x-direction
    int ny = 4; // number of unit cells in the y-direction
    int nz = 4; // number of unit cells in the z-direction
    int n0 = 4; // number of particles in the unit cell
    int N = n0 * nx * ny * nz; // total number of particles
    int Ne = 20000;   // number of steps in the equilibration stage
    int Np = 20000;   // number of steps in the production stage
    int Ns = 10;      // sampling interval
    int Nd = Np / Ns; // number of heat current data
    int Nc = Nd / 10;   // number of correlation data
    int MN = 200;     // maximum number of neighbors for one particle

    // For LJ argon
    // Temperature (K)      20       30       40       50       60    
    // lattice constant (A) 5.284    5.305    5.329    5.356    5.385
    double T_0 = 60.0;    // temperature prescribed
    double ax = 5.385;    // lattice constant in the x direction
    double ay = ax;       // lattice constant in the y direction
    double az = ax;       // lattice constant in the z direction
    double lx = ax * nx;  // box length in the x direction
    double ly = ay * ny;  // box length in the y direction
    double lz = az * nz;  // box length in the z direction
    double cutoff = 12.0; // cutoff distance for neighbor list
    
    double time_step = 10.0 / TIME_UNIT_CONVERSION; // time step
    
    // fixed neighbor list
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
    double *hx = (double*) malloc(Nd * sizeof(double)); // heat current
    double *hy = (double*) malloc(Nd * sizeof(double));
    double *hz = (double*) malloc(Nd * sizeof(double));

    // initialize mass, position, and velocity
    for (int n = 0; n < N; ++n) { m[n] = 40.0; } // mass for argon atom
    initialize_position(n0, nx, ny, nz, ax, ay, az, x, y, z);
    initialize_velocity(N, T_0, m, vx, vy, vz);

    // initialize neighbor list and force
    find_neighbor(N, NN, NL, x, y, z, lx, ly, lz, MN, cutoff);
    double hc[3]; // heat current at a specific time point
    find_force(N, NN, NL, MN, lx, ly, lz, x, y, z, fx, fy, fz, vx, vy, vz, hc);
 
    // equilibration
    clock_t time_begin = clock();
    for (int step = 0; step < Ne; ++step)
    { 
        integrate(N, time_step, m, fx, fy, fz, vx, vy, vz, x, y, z, 1);
        find_force
        (N, NN, NL, MN, lx, ly, lz, x, y, z, fx, fy, fz, vx, vy, vz, hc);
        integrate(N, time_step, m, fx, fy, fz, vx, vy, vz, x, y, z, 2);
        scale_velocity(N, T_0, m, vx, vy, vz); // control temperature
    } 
    clock_t time_finish = clock();
    double time_used = (time_finish - time_begin) / (double) CLOCKS_PER_SEC;
    fprintf(stderr, "time use for equilibration = %f s\n", time_used); 

    // production
    time_begin = clock();
    int count = 0;
    for (int step = 0; step < Np; ++step)
    {  
        integrate(N, time_step, m, fx, fy, fz, vx, vy, vz, x, y, z, 1);
        find_force
        (N, NN, NL, MN, lx, ly, lz, x, y, z, fx, fy, fz, vx, vy, vz, hc);
        integrate(N, time_step, m, fx, fy, fz, vx, vy, vz, x, y, z, 2);
        if (0 == step % Ns) 
        { hx[count] = hc[0]; hy[count] = hc[1]; hz[count] = hc[2]; count++; }
    } 
    time_finish = clock();
    time_used = (time_finish - time_begin) / (double) CLOCKS_PER_SEC;
    fprintf(stderr, "time use for production = %f s\n", time_used); 

    // calculate heat current autocorrelation function and thermal conductivity
    find_hac_kappa(Nd, Nc, time_step * Ns, T_0, lx * ly * lz, hx, hy, hz);

    free(NN); free(NL); free(m);  free(x);  free(y);  free(z);
    free(vx); free(vy); free(vz); free(fx); free(fy); free(fz);
    free(hx); free(hy); free(hz);
    
    //system("PAUSE"); // For DEV-C++ in Windows; 

    return 0;
}
