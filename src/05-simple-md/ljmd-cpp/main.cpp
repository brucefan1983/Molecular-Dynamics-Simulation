#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <numeric>
#include <vector>

const double K_B = 8.617343e-5;
const double TIME_UNIT_CONVERSION = 1.018051e+1;

static void scale_velocity1(
  int N,
  double T_0,
  std::vector<double>& vx,
  std::vector<double>& vy,
  std::vector<double>& vz,
  std::vector<double>& ke)
{
  double temperature = std::accumulate(ke.begin(), ke.end(), 0.0);
  temperature /= 1.5 * K_B * N;
  double scale_factor = sqrt(T_0 / temperature);
  for (int n = 0; n < N; ++n) {
    vx[n] *= scale_factor;
    vy[n] *= scale_factor;
    vz[n] *= scale_factor;
  }
}

static void scale_velocity(
  int N,
  double T_0,
  std::vector<double>& mass,
  std::vector<double>& vx,
  std::vector<double>& vy,
  std::vector<double>& vz)
{
  double temperature = 0.0;
  for (int n = 0; n < N; ++n) {
    double v2 = vx[n] * vx[n] + vy[n] * vy[n] + vz[n] * vz[n];
    temperature += mass[n] * v2;
  }
  temperature /= 3.0 * K_B * N;
  double scale_factor = sqrt(T_0 / temperature);
  for (int n = 0; n < N; ++n) {
    vx[n] *= scale_factor;
    vy[n] *= scale_factor;
    vz[n] *= scale_factor;
  }
}

void initialize_position(
  int nx,
  double ax,
  std::vector<double>& box,
  std::vector<double>& x,
  std::vector<double>& y,
  std::vector<double>& z)
{
  box[0] = ax * nx;
  box[1] = ax * nx;
  box[2] = ax * nx;
  box[3] = box[0] * 0.5;
  box[4] = box[1] * 0.5;
  box[5] = box[2] * 0.5;
  double x0[4] = {0.0, 0.0, 0.5, 0.5};
  double y0[4] = {0.0, 0.5, 0.0, 0.5};
  double z0[4] = {0.0, 0.5, 0.5, 0.0};
  int n = 0;
  for (int ix = 0; ix < nx; ++ix) {
    for (int iy = 0; iy < nx; ++iy) {
      for (int iz = 0; iz < nx; ++iz) {
        for (int i = 0; i < 4; ++i) {
          x[n] = (ix + x0[i]) * ax;
          y[n] = (iy + y0[i]) * ax;
          z[n] = (iz + z0[i]) * ax;
          n++;
        }
      }
    }
  }
}

void initialize_velocity(
  int N,
  double T_0,
  std::vector<double>& mass,
  std::vector<double>& vx,
  std::vector<double>& vy,
  std::vector<double>& vz)
{
  double momentum_average[3] = {0.0, 0.0, 0.0};
  for (int n = 0; n < N; ++n) {
    vx[n] = -1.0 + (rand() * 2.0) / RAND_MAX;
    vy[n] = -1.0 + (rand() * 2.0) / RAND_MAX;
    vz[n] = -1.0 + (rand() * 2.0) / RAND_MAX;

    momentum_average[0] += mass[n] * vx[n] / N;
    momentum_average[1] += mass[n] * vy[n] / N;
    momentum_average[2] += mass[n] * vz[n] / N;
  }
  for (int n = 0; n < N; ++n) {
    vx[n] -= momentum_average[0] / mass[n];
    vy[n] -= momentum_average[1] / mass[n];
    vz[n] -= momentum_average[2] / mass[n];
  }
  scale_velocity(N, T_0, mass, vx, vy, vz);
}

static void apply_mic(double* box, double* x12, double* y12, double* z12)
{
  if (*x12 < -box[3]) {
    *x12 += box[0];
  } else if (*x12 > +box[3]) {
    *x12 -= box[0];
  }
  if (*y12 < -box[4]) {
    *y12 += box[1];
  } else if (*y12 > +box[4]) {
    *y12 -= box[1];
  }
  if (*z12 < -box[5]) {
    *z12 += box[2];
  } else if (*z12 > +box[5]) {
    *z12 -= box[2];
  }
}

void find_neighbor(
  int N,
  int MN,
  std::vector<double>& box,
  std::vector<double>& x,
  std::vector<double>& y,
  std::vector<double>& z,
  std::vector<int>& NN,
  std::vector<int>& NL)
{
  double cutoff = 11.0;
  double cutoff_square = cutoff * cutoff;

  for (int n = 0; n < N; n++) {
    NN[n] = 0;
  }

  for (int n1 = 0; n1 < N - 1; n1++) {
    for (int n2 = n1 + 1; n2 < N; n2++) {
      double x12 = x[n2] - x[n1];
      double y12 = y[n2] - y[n1];
      double z12 = z[n2] - z[n1];
      apply_mic(box.data(), &x12, &y12, &z12);
      double d_square = x12 * x12 + y12 * y12 + z12 * z12;

      if (d_square < cutoff_square) {
        NL[n1 * MN + NN[n1]++] = n2;
        NL[n2 * MN + NN[n2]++] = n1;
      }
    }
  }

  for (int n1 = 0; n1 < N - 1; n1++) {
    if (NN[n1] > MN) {
      printf("Error: MN is too small.\n");
      exit(1);
    }
  }
}

void find_force(
  int N,
  int MN,
  std::vector<double>& box,
  std::vector<double>& x,
  std::vector<double>& y,
  std::vector<double>& z,
  std::vector<int>& NN,
  std::vector<int>& NL,
  std::vector<double>& fx,
  std::vector<double>& fy,
  std::vector<double>& fz,
  std::vector<double>& pe)
{
  const double epsilon = 1.032e-2;
  const double sigma = 3.405;
  const double cutoff = 10.0;
  const double cutoff_square = cutoff * cutoff;
  const double sigma_3 = sigma * sigma * sigma;
  const double sigma_6 = sigma_3 * sigma_3;
  const double sigma_12 = sigma_6 * sigma_6;
  const double e24s6 = 24.0 * epsilon * sigma_6;
  const double e48s12 = 48.0 * epsilon * sigma_12;
  const double e4s6 = 4.0 * epsilon * sigma_6;
  const double e4s12 = 4.0 * epsilon * sigma_12;
  for (int n = 0; n < N; ++n) {
    fx[n] = fy[n] = fz[n] = pe[n] = 0.0;
  }
  for (int i = 0; i < N; ++i) {
    for (int k = 0; k < NN[i]; k++) {
      int j = NL[i * MN + k];
      if (j < i)
        continue;
      double x_ij = x[j] - x[i];
      double y_ij = y[j] - y[i];
      double z_ij = z[j] - z[i];
      apply_mic(box.data(), &x_ij, &y_ij, &z_ij);
      double r2 = x_ij * x_ij + y_ij * y_ij + z_ij * z_ij;
      if (r2 > cutoff_square) {
        continue;
      }
      double r2inv = 1.0 / r2;
      double r4inv = r2inv * r2inv;
      double r6inv = r2inv * r4inv;
      double r8inv = r4inv * r4inv;
      double r12inv = r4inv * r8inv;
      double r14inv = r6inv * r8inv;
      double f_ij = e24s6 * r8inv - e48s12 * r14inv;
      pe[i] += e4s12 * r12inv - e4s6 * r6inv;
      fx[i] += f_ij * x_ij;
      fx[j] -= f_ij * x_ij;
      fy[i] += f_ij * y_ij;
      fy[j] -= f_ij * y_ij;
      fz[i] += f_ij * z_ij;
      fz[j] -= f_ij * z_ij;
    }
  }
}

static void integrate(
  int N,
  double time_step,
  std::vector<double>& mass,
  std::vector<double>& fx,
  std::vector<double>& fy,
  std::vector<double>& fz,
  std::vector<double>& x,
  std::vector<double>& y,
  std::vector<double>& z,
  std::vector<double>& vx,
  std::vector<double>& vy,
  std::vector<double>& vz,
  std::vector<double>& ke,
  int flag)
{
  double time_step_half = time_step * 0.5;
  for (int n = 0; n < N; ++n) {
    double mass_inv = 1.0 / mass[n];
    double ax = fx[n] * mass_inv;
    double ay = fy[n] * mass_inv;
    double az = fz[n] * mass_inv;
    vx[n] += ax * time_step_half;
    vy[n] += ay * time_step_half;
    vz[n] += az * time_step_half;
    if (flag == 1) {
      x[n] += vx[n] * time_step;
      y[n] += vy[n] * time_step;
      z[n] += vz[n] * time_step;
    } else {
      double v2 = vx[n] * vx[n] + vy[n] * vy[n] + vz[n] * vz[n];
      ke[n] = mass[n] * v2 * 0.5;
    }
  }
}

int main(int argc, char** argv)
{
  int nx = 5;
  int Ne = 20000;
  int Np = 20000;

  if (argc != 3) {
    printf("Usage: %s nx Ne\n", argv[0]);
    exit(1);
  } else {
    nx = atoi(argv[1]);
    Ne = atoi(argv[2]);
    Np = Ne;
  }

  const int N = 4 * nx * nx * nx;
  int Ns = 100;
  const int MN = 200;
  double T_0 = 60.0;
  double ax = 5.385;
  double time_step = 5.0 / TIME_UNIT_CONVERSION;

  std::vector<int> NN(N);
  std::vector<int> NL(N * MN);
  std::vector<double> mass(N);
  std::vector<double> x(N);
  std::vector<double> y(N);
  std::vector<double> z(N);
  std::vector<double> vx(N);
  std::vector<double> vy(N);
  std::vector<double> vz(N);
  std::vector<double> fx(N);
  std::vector<double> fy(N);
  std::vector<double> fz(N);
  std::vector<double> pe(N);
  std::vector<double> ke(N);
  std::vector<double> box(6);

  for (int n = 0; n < N; ++n) {
    mass[n] = 40.0;
  }

  initialize_position(nx, ax, box, x, y, z);
  initialize_velocity(N, T_0, mass, vx, vy, vz);
  find_neighbor(N, MN, box, x, y, z, NN, NL);

  find_force(N, MN, box, x, y, z, NN, NL, fx, fy, fz, pe);
  for (int step = 0; step < Ne; ++step) {
    integrate(N, time_step, mass, fx, fy, fz, x, y, z, vx, vy, vz, ke, 1);
    find_force(N, MN, box, x, y, z, NN, NL, fx, fy, fz, pe);
    integrate(N, time_step, mass, fx, fy, fz, x, y, z, vx, vy, vz, ke, 2);
    scale_velocity1(N, T_0, vx, vy, vz, ke);
  }

  clock_t t_total_start = clock();

  FILE* fid = fopen("energy.txt", "w");
  for (int step = 0; step < Np; ++step) {
    integrate(N, time_step, mass, fx, fy, fz, x, y, z, vx, vy, vz, ke, 1);
    find_force(N, MN, box, x, y, z, NN, NL, fx, fy, fz, pe);
    integrate(N, time_step, mass, fx, fy, fz, x, y, z, vx, vy, vz, ke, 2);
    if (0 == step % Ns) {
      const double keTotal = std::accumulate(ke.begin(), ke.end(), 0.0);
      const double peTotal = std::accumulate(pe.begin(), pe.end(), 0.0);
      fprintf(fid, "%g %g\n", keTotal, peTotal);
    }
  }
  fclose(fid);

  clock_t t_total_stop = clock();

  float t_total = float(t_total_stop - t_total_start) / CLOCKS_PER_SEC;
  printf("Time used for production = %g s\n", t_total);

  return 0;
}
