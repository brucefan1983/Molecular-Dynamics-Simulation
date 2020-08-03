#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <numeric>
#include <vector>

const double K_B = 8.617343e-5;
const double TIME_UNIT_CONVERSION = 1.018051e+1;

void scaleVelocity(
  const int N,
  const double T0,
  const std::vector<double>& mass,
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
  double scaleFactor = sqrt(T0 / temperature);
  for (int n = 0; n < N; ++n) {
    vx[n] *= scaleFactor;
    vy[n] *= scaleFactor;
    vz[n] *= scaleFactor;
  }
}

void initializePosition(
  const int nx,
  const double ax,
  double box[6],
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

void initializeVelocity(
  const int N,
  const double T0,
  const std::vector<double>& mass,
  std::vector<double>& vx,
  std::vector<double>& vy,
  std::vector<double>& vz)
{
  double momentumAverage[3] = {0.0, 0.0, 0.0};
  for (int n = 0; n < N; ++n) {
    vx[n] = -1.0 + (rand() * 2.0) / RAND_MAX;
    vy[n] = -1.0 + (rand() * 2.0) / RAND_MAX;
    vz[n] = -1.0 + (rand() * 2.0) / RAND_MAX;

    momentumAverage[0] += mass[n] * vx[n] / N;
    momentumAverage[1] += mass[n] * vy[n] / N;
    momentumAverage[2] += mass[n] * vz[n] / N;
  }
  for (int n = 0; n < N; ++n) {
    vx[n] -= momentumAverage[0] / mass[n];
    vy[n] -= momentumAverage[1] / mass[n];
    vz[n] -= momentumAverage[2] / mass[n];
  }
  scaleVelocity(N, T0, mass, vx, vy, vz);
}

void applyMicOne(const double length, const double halfLength, double& x12)
{
  if (x12 < -halfLength)
    x12 += length;
  else if (x12 > +halfLength)
    x12 -= length;
}

void applyMic(const double box[6], double& x12, double& y12, double& z12)
{
  applyMicOne(box[0], box[3], x12);
  applyMicOne(box[1], box[4], y12);
  applyMicOne(box[2], box[5], z12);
}

void findNeighbor(
  const int N,
  const int MN,
  const double box[6],
  const std::vector<double>& x,
  const std::vector<double>& y,
  std::vector<double>& z,
  std::vector<int>& NN,
  std::vector<int>& NL)
{
  double cutoff = 11.0;
  double cutoffSquare = cutoff * cutoff;

  for (int n = 0; n < N; n++)
    NN[n] = 0;

  for (int n1 = 0; n1 < N - 1; n1++) {
    for (int n2 = n1 + 1; n2 < N; n2++) {
      double x12 = x[n2] - x[n1];
      double y12 = y[n2] - y[n1];
      double z12 = z[n2] - z[n1];
      applyMic(box, x12, y12, z12);
      double dSquare = x12 * x12 + y12 * y12 + z12 * z12;

      if (dSquare < cutoffSquare) {
        NL[n1 * MN + NN[n1]++] = n2;
        NL[n2 * MN + NN[n2]++] = n1;
      }
    }
  }

  for (int n1 = 0; n1 < N; n1++) {
    if (NN[n1] > MN) {
      printf("Error: MN is too small.\n");
      exit(1);
    }
  }
}

void find_force(
  const int N,
  const int MN,
  const double box[6],
  const std::vector<double>& x,
  const std::vector<double>& y,
  const std::vector<double>& z,
  const std::vector<int>& NN,
  const std::vector<int>& NL,
  std::vector<double>& fx,
  std::vector<double>& fy,
  std::vector<double>& fz,
  std::vector<double>& pe)
{
  const double epsilon = 1.032e-2;
  const double sigma = 3.405;
  const double cutoff = 10.0;
  const double cutoffSquare = cutoff * cutoff;
  const double sigma3 = sigma * sigma * sigma;
  const double sigma6 = sigma3 * sigma3;
  const double sigma12 = sigma6 * sigma6;
  const double e24s6 = 24.0 * epsilon * sigma6;
  const double e48s12 = 48.0 * epsilon * sigma12;
  const double e4s6 = 4.0 * epsilon * sigma6;
  const double e4s12 = 4.0 * epsilon * sigma12;
  for (int n = 0; n < N; ++n)
    fx[n] = fy[n] = fz[n] = pe[n] = 0.0;

  for (int i = 0; i < N; ++i) {
    for (int k = 0; k < NN[i]; k++) {
      int j = NL[i * MN + k];
      if (j < i)
        continue;
      double xij = x[j] - x[i];
      double yij = y[j] - y[i];
      double zij = z[j] - z[i];
      applyMic(box, xij, yij, zij);
      double r2 = xij * xij + yij * yij + zij * zij;
      if (r2 > cutoffSquare)
        continue;

      double r2inv = 1.0 / r2;
      double r4inv = r2inv * r2inv;
      double r6inv = r2inv * r4inv;
      double r8inv = r4inv * r4inv;
      double r12inv = r4inv * r8inv;
      double r14inv = r6inv * r8inv;
      double f_ij = e24s6 * r8inv - e48s12 * r14inv;
      pe[i] += e4s12 * r12inv - e4s6 * r6inv;
      fx[i] += f_ij * xij;
      fx[j] -= f_ij * xij;
      fy[i] += f_ij * yij;
      fy[j] -= f_ij * yij;
      fz[i] += f_ij * zij;
      fz[j] -= f_ij * zij;
    }
  }
}

void integrate(
  const int N,
  const double timeStep,
  const std::vector<double>& mass,
  const std::vector<double>& fx,
  const std::vector<double>& fy,
  const std::vector<double>& fz,
  std::vector<double>& x,
  std::vector<double>& y,
  std::vector<double>& z,
  std::vector<double>& vx,
  std::vector<double>& vy,
  std::vector<double>& vz,
  const int flag)
{
  double timeStepHalf = timeStep * 0.5;
  for (int n = 0; n < N; ++n) {
    double mass_inv = 1.0 / mass[n];
    double ax = fx[n] * mass_inv;
    double ay = fy[n] * mass_inv;
    double az = fz[n] * mass_inv;
    vx[n] += ax * timeStepHalf;
    vy[n] += ay * timeStepHalf;
    vz[n] += az * timeStepHalf;
    if (flag == 1) {
      x[n] += vx[n] * timeStep;
      y[n] += vy[n] * timeStep;
      z[n] += vz[n] * timeStep;
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
  const int Ns = 100;
  const int MN = 200;
  const double T0 = 60.0;
  const double ax = 5.385;
  const double timeStep = 5.0 / TIME_UNIT_CONVERSION;

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
  double box[6];

  for (int n = 0; n < N; ++n)
    mass[n] = 40.0;

  initializePosition(nx, ax, box, x, y, z);
  initializeVelocity(N, T0, mass, vx, vy, vz);
  findNeighbor(N, MN, box, x, y, z, NN, NL);

  find_force(N, MN, box, x, y, z, NN, NL, fx, fy, fz, pe);
  for (int step = 0; step < Ne; ++step) {
    integrate(N, timeStep, mass, fx, fy, fz, x, y, z, vx, vy, vz, 1);
    find_force(N, MN, box, x, y, z, NN, NL, fx, fy, fz, pe);
    integrate(N, timeStep, mass, fx, fy, fz, x, y, z, vx, vy, vz, 2);
    scaleVelocity(N, T0, mass, vx, vy, vz);
  }

  const clock_t tStart = clock();

  FILE* fid = fopen("energy.txt", "w");
  for (int step = 0; step < Np; ++step) {
    integrate(N, timeStep, mass, fx, fy, fz, x, y, z, vx, vy, vz, 1);
    find_force(N, MN, box, x, y, z, NN, NL, fx, fy, fz, pe);
    integrate(N, timeStep, mass, fx, fy, fz, x, y, z, vx, vy, vz, 2);
    if (0 == step % Ns) {
      double keTotal = 0.0;
      for (int n = 0; n < N; ++n) {
        double v2 = vx[n] * vx[n] + vy[n] * vy[n] + vz[n] * vz[n];
        keTotal += mass[n] * v2;
      }
      keTotal *= 0.5;
      const double peTotal = std::accumulate(pe.begin(), pe.end(), 0.0);
      fprintf(fid, "%g %g\n", keTotal, peTotal);
    }
  }
  fclose(fid);

  const clock_t tStop = clock();
  const float tElapsed = float(tStop - tStart) / CLOCKS_PER_SEC;
  printf("Time used for production = %g s\n", tElapsed);

  return 0;
}
