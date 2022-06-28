/*----------------------------------------------------------------------------80
    Copyright 2022 Zheyong Fan
Compile:
    using O(N^2) neighbor list algorithm:
        g++ neighbor.cpp -O3 -o neighbor
    using O(N) neighbor list algorithm:
        g++ neighbor.cpp -O3 -DUSE_ON1 -o neighbor
Run:
    command:
        neighbor numCells numSteps temperature timeStep
    example:
        neighbor 6 20000 60 5
------------------------------------------------------------------------------*/

#include <cmath>    // sqrt() function
#include <ctime>    // for timing
#include <fstream>  // file
#include <iomanip>  // std::setprecision
#include <iostream> // input/output
#include <numeric>  // summation
#include <vector>   // vector

const int Ns = 100;             // output frequency
const double K_B = 8.617343e-5; // Boltzmann's constant in natural unit
const double TIME_UNIT_CONVERSION = 1.018051e+1; // from natural unit to fs

struct Atom {
  int number;
  int numUpdates = 0;
  const int MN = 1000;
  double cutoffNeighbor = 10.0;
  double box[18];
  std::vector<int> NN, NL;
  std::vector<double> mass, x0, y0, z0, x, y, z, vx, vy, vz, fx, fy, fz, pe;
};

double findKineticEnergy(const Atom& atom)
{
  double kineticEnergy = 0.0;
  for (int n = 0; n < atom.number; ++n) {
    double v2 = atom.vx[n] * atom.vx[n] + atom.vy[n] * atom.vy[n] +
                atom.vz[n] * atom.vz[n];
    kineticEnergy += atom.mass[n] * v2;
  }
  return kineticEnergy * 0.5;
}

void scaleVelocity(const double T0, Atom& atom)
{
  const double temperature =
    findKineticEnergy(atom) * 2.0 / (3.0 * K_B * atom.number);
  double scaleFactor = sqrt(T0 / temperature);
  for (int n = 0; n < atom.number; ++n) {
    atom.vx[n] *= scaleFactor;
    atom.vy[n] *= scaleFactor;
    atom.vz[n] *= scaleFactor;
  }
}

void allocateMemory(const int numCells, Atom& atom)
{
  const int numAtomsPerCell = 4;
  atom.number = numCells * numCells * numCells * numAtomsPerCell;
  atom.NN.resize(atom.number, 0);
  atom.NL.resize(atom.number * atom.MN, 0);
  atom.mass.resize(atom.number, 40.0); // argon mass
  atom.x0.resize(atom.number, 0.0);
  atom.y0.resize(atom.number, 0.0);
  atom.z0.resize(atom.number, 0.0);
  atom.x.resize(atom.number, 0.0);
  atom.y.resize(atom.number, 0.0);
  atom.z.resize(atom.number, 0.0);
  atom.vx.resize(atom.number, 0.0);
  atom.vy.resize(atom.number, 0.0);
  atom.vz.resize(atom.number, 0.0);
  atom.fx.resize(atom.number, 0.0);
  atom.fy.resize(atom.number, 0.0);
  atom.fz.resize(atom.number, 0.0);
  atom.pe.resize(atom.number, 0.0);
}

double getDet(const double* box)
{
  return box[0] * (box[4] * box[8] - box[5] * box[7]) +
         box[1] * (box[5] * box[6] - box[3] * box[8]) +
         box[2] * (box[3] * box[7] - box[4] * box[6]);
}

void getInverseBox(double* box)
{
  box[9] = box[4] * box[8] - box[5] * box[7];
  box[10] = box[2] * box[7] - box[1] * box[8];
  box[11] = box[1] * box[5] - box[2] * box[4];
  box[12] = box[5] * box[6] - box[3] * box[8];
  box[13] = box[0] * box[8] - box[2] * box[6];
  box[14] = box[2] * box[3] - box[0] * box[5];
  box[15] = box[3] * box[7] - box[4] * box[6];
  box[16] = box[1] * box[6] - box[0] * box[7];
  box[17] = box[0] * box[4] - box[1] * box[3];
  double det = getDet(box);
  for (int n = 9; n < 18; ++n) {
    box[n] /= det;
  }
}

float getArea(const double* a, const double* b)
{
  const double s1 = a[1] * b[2] - a[2] * b[1];
  const double s2 = a[2] * b[0] - a[0] * b[2];
  const double s3 = a[0] * b[1] - a[1] * b[0];
  return sqrt(s1 * s1 + s2 * s2 + s3 * s3);
}

void getThickness(const Atom& atom, double* thickness)
{
  double volume = abs(getDet(atom.box));
  const double a[3] = {atom.box[0], atom.box[3], atom.box[6]};
  const double b[3] = {atom.box[1], atom.box[4], atom.box[7]};
  const double c[3] = {atom.box[2], atom.box[5], atom.box[8]};
  thickness[0] = volume / getArea(b, c);
  thickness[1] = volume / getArea(c, a);
  thickness[2] = volume / getArea(a, b);
}

void initializePosition(const int numCells, Atom& atom)
{
  const int numAtomsPerCell = 4;
  const double latticeConstant = 5.385;
  for (int d1 = 0; d1 < 3; ++d1) {
    for (int d2 = 0; d2 < 3; ++d2) {
      if (d1 == d2) {
        atom.box[d1 * 3 + d2] = latticeConstant * numCells;
      } else {
        atom.box[d1 * 3 + d2] = 0.0;
      }
    }
  }
  getInverseBox(atom.box);
  const double x0[numAtomsPerCell] = {0.0, 0.0, 0.5, 0.5};
  const double y0[numAtomsPerCell] = {0.0, 0.5, 0.0, 0.5};
  const double z0[numAtomsPerCell] = {0.0, 0.5, 0.5, 0.0};
  int n = 0;
  for (int ix = 0; ix < numCells; ++ix) {
    for (int iy = 0; iy < numCells; ++iy) {
      for (int iz = 0; iz < numCells; ++iz) {
        for (int i = 0; i < numAtomsPerCell; ++i) {
          atom.x[n] = (ix + x0[i]) * latticeConstant;
          atom.y[n] = (iy + y0[i]) * latticeConstant;
          atom.z[n] = (iz + z0[i]) * latticeConstant;
          ++n;
        }
      }
    }
  }
}

void initializeVelocity(const double T0, Atom& atom)
{
#ifndef DEBUG
  srand(time(NULL));
#endif
  double centerOfMassVelocity[3] = {0.0, 0.0, 0.0};
  double totalMass = 0.0;
  for (int n = 0; n < atom.number; ++n) {
    totalMass += atom.mass[n];
    atom.vx[n] = -1.0 + (rand() * 2.0) / RAND_MAX;
    atom.vy[n] = -1.0 + (rand() * 2.0) / RAND_MAX;
    atom.vz[n] = -1.0 + (rand() * 2.0) / RAND_MAX;
    centerOfMassVelocity[0] += atom.mass[n] * atom.vx[n];
    centerOfMassVelocity[1] += atom.mass[n] * atom.vy[n];
    centerOfMassVelocity[2] += atom.mass[n] * atom.vz[n];
  }
  centerOfMassVelocity[0] /= totalMass;
  centerOfMassVelocity[1] /= totalMass;
  centerOfMassVelocity[2] /= totalMass;
  for (int n = 0; n < atom.number; ++n) {
    atom.vx[n] -= centerOfMassVelocity[0];
    atom.vy[n] -= centerOfMassVelocity[1];
    atom.vz[n] -= centerOfMassVelocity[2];
  }
  scaleVelocity(T0, atom);
}

void applyMicOne(double& x12)
{
  if (x12 < -0.5)
    x12 += 1.0;
  else if (x12 > +0.5)
    x12 -= 1.0;
}

void applyMic(const double* box, double& x12, double& y12, double& z12)
{
  double sx12 = box[9] * x12 + box[10] * y12 + box[11] * z12;
  double sy12 = box[12] * x12 + box[13] * y12 + box[14] * z12;
  double sz12 = box[15] * x12 + box[16] * y12 + box[17] * z12;
  applyMicOne(sx12);
  applyMicOne(sy12);
  applyMicOne(sz12);
  x12 = box[0] * sx12 + box[1] * sy12 + box[2] * sz12;
  y12 = box[3] * sx12 + box[4] * sy12 + box[5] * sz12;
  z12 = box[6] * sx12 + box[7] * sy12 + box[8] * sz12;
}

bool checkIfNeedUpdate(const Atom& atom)
{
  bool needUpdate = false;
  for (int n = 0; n < atom.number; ++n) {
    double dx = atom.x[n] - atom.x0[n];
    double dy = atom.y[n] - atom.y0[n];
    double dz = atom.z[n] - atom.z0[n];
    if (dx * dx + dy * dy + dz * dz > 0.25) {
      needUpdate = true;
      break;
    }
  }
  return needUpdate;
}

void applyPbcOne(double& sx)
{
  if (sx < 0.0) {
    sx += 1.0;
  } else if (sx > 1.0) {
    sx -= 1.0;
  }
}

void applyPbc(Atom& atom)
{
  for (int n = 0; n < atom.number; ++n) {
    double sx = atom.box[9] * atom.x[n] + atom.box[10] * atom.y[n] +
                atom.box[11] * atom.z[n];
    double sy = atom.box[12] * atom.x[n] + atom.box[13] * atom.y[n] +
                atom.box[14] * atom.z[n];
    double sz = atom.box[15] * atom.x[n] + atom.box[16] * atom.y[n] +
                atom.box[17] * atom.z[n];
    applyPbcOne(sx);
    applyPbcOne(sy);
    applyPbcOne(sz);
    atom.x[n] = atom.box[0] * sx + atom.box[1] * sy + atom.box[2] * sz;
    atom.y[n] = atom.box[3] * sx + atom.box[4] * sy + atom.box[5] * sz;
    atom.z[n] = atom.box[6] * sx + atom.box[7] * sy + atom.box[8] * sz;
  }
}

void updateXyz0(Atom& atom)
{
  for (int n = 0; n < atom.number; ++n) {
    atom.x0[n] = atom.x[n];
    atom.y0[n] = atom.y[n];
    atom.z0[n] = atom.z[n];
  }
}

void findNeighborON2(Atom& atom)
{
  const double cutoffSquare = atom.cutoffNeighbor * atom.cutoffNeighbor;
  std::fill(atom.NN.begin(), atom.NN.end(), 0);

  for (int i = 0; i < atom.number - 1; ++i) {
    const double x1 = atom.x[i];
    const double y1 = atom.y[i];
    const double z1 = atom.z[i];
    for (int j = i + 1; j < atom.number; ++j) {
      double xij = atom.x[j] - x1;
      double yij = atom.y[j] - y1;
      double zij = atom.z[j] - z1;
      applyMic(atom.box, xij, yij, zij);
      const double distanceSquare = xij * xij + yij * yij + zij * zij;
      if (distanceSquare < cutoffSquare) {
        atom.NL[i * atom.MN + atom.NN[i]++] = j;
        if (atom.NN[i] > atom.MN) {
          std::cout << "Error: number of neighbors for atom " << i
                    << " exceeds " << atom.MN << std::endl;
          exit(1);
        }
      }
    }
  }
}

void findCell(
  const double* box,
  const double* thickness,
  const double* r,
  double cutoffInverse,
  const int* numCells,
  int* cell)
{
  double s[3];
  s[0] = box[9] * r[0] + box[10] * r[1] + box[11] * r[2];
  s[1] = box[12] * r[0] + box[13] * r[1] + box[14] * r[2];
  s[2] = box[15] * r[0] + box[16] * r[1] + box[17] * r[2];
  for (int d = 0; d < 3; ++d) {
    cell[d] = floor(s[d] * thickness[d] * cutoffInverse);
    if (cell[d] < 0)
      cell[d] += numCells[d];
    if (cell[d] >= numCells[d])
      cell[d] -= numCells[d];
  }
  cell[3] = cell[0] + numCells[0] * (cell[1] + numCells[1] * cell[2]);
}

void findNeighborON1(Atom& atom)
{
  const double cutoffInverse = 1.0 / atom.cutoffNeighbor;
  double cutoffSquare = atom.cutoffNeighbor * atom.cutoffNeighbor;
  double thickness[3];
  getThickness(atom, thickness);

  int numCells[4];

  for (int d = 0; d < 3; ++d) {
    numCells[d] = floor(thickness[d] * cutoffInverse);
  }

  numCells[3] = numCells[0] * numCells[1] * numCells[2];
  int cell[4];

  std::vector<int> cellCount(numCells[3], 0);
  std::vector<int> cellCountSum(numCells[3], 0);

  for (int n = 0; n < atom.number; ++n) {
    const double r[3] = {atom.x[n], atom.y[n], atom.z[n]};
    findCell(atom.box, thickness, r, cutoffInverse, numCells, cell);
    ++cellCount[cell[3]];
  }

  for (int i = 1; i < numCells[3]; ++i) {
    cellCountSum[i] = cellCountSum[i - 1] + cellCount[i - 1];
  }

  std::fill(cellCount.begin(), cellCount.end(), 0);

  std::vector<int> cellContents(atom.number, 0);

  for (int n = 0; n < atom.number; ++n) {
    const double r[3] = {atom.x[n], atom.y[n], atom.z[n]};
    findCell(atom.box, thickness, r, cutoffInverse, numCells, cell);
    cellContents[cellCountSum[cell[3]] + cellCount[cell[3]]] = n;
    ++cellCount[cell[3]];
  }

  std::fill(atom.NN.begin(), atom.NN.end(), 0);

  for (int n1 = 0; n1 < atom.number; ++n1) {
    const double r1[3] = {atom.x[n1], atom.y[n1], atom.z[n1]};
    findCell(atom.box, thickness, r1, cutoffInverse, numCells, cell);
    for (int k = -1; k <= 1; ++k) {
      for (int j = -1; j <= 1; ++j) {
        for (int i = -1; i <= 1; ++i) {
          int neighborCell = cell[3] + (k * numCells[1] + j) * numCells[0] + i;
          if (cell[0] + i < 0)
            neighborCell += numCells[0];
          if (cell[0] + i >= numCells[0])
            neighborCell -= numCells[0];
          if (cell[1] + j < 0)
            neighborCell += numCells[1] * numCells[0];
          if (cell[1] + j >= numCells[1])
            neighborCell -= numCells[1] * numCells[0];
          if (cell[2] + k < 0)
            neighborCell += numCells[3];
          if (cell[2] + k >= numCells[2])
            neighborCell -= numCells[3];

          for (int m = 0; m < cellCount[neighborCell]; ++m) {
            const int n2 = cellContents[cellCountSum[neighborCell] + m];
            if (n1 < n2) {
              double x12 = atom.x[n2] - r1[0];
              double y12 = atom.y[n2] - r1[1];
              double z12 = atom.z[n2] - r1[2];
              applyMic(atom.box, x12, y12, z12);
              const double d2 = x12 * x12 + y12 * y12 + z12 * z12;
              if (d2 < cutoffSquare) {
                atom.NL[n1 * atom.MN + atom.NN[n1]++] = n2;
                if (atom.NN[n1] > atom.MN) {
                  std::cout << "Error: number of neighbors for atom " << n1
                            << " exceeds " << atom.MN << std::endl;
                  exit(1);
                }
              }
            }
          }
        }
      }
    }
  }
}

void findNeighbor(Atom& atom)
{
  if (checkIfNeedUpdate(atom)) {
    atom.numUpdates++;
    applyPbc(atom);
#ifdef USE_ON1
    findNeighborON1(atom);
#else
    findNeighborON2(atom);
#endif
    updateXyz0(atom);
  }
}

void findForce(Atom& atom)
{
  const double epsilon = 1.032e-2;
  const double sigma = 3.405;
  const double cutoff = 9.0;
  const double cutoffSquare = cutoff * cutoff;
  const double sigma3 = sigma * sigma * sigma;
  const double sigma6 = sigma3 * sigma3;
  const double sigma12 = sigma6 * sigma6;
  const double e24s6 = 24.0 * epsilon * sigma6;
  const double e48s12 = 48.0 * epsilon * sigma12;
  const double e4s6 = 4.0 * epsilon * sigma6;
  const double e4s12 = 4.0 * epsilon * sigma12;
  for (int n = 0; n < atom.number; ++n)
    atom.fx[n] = atom.fy[n] = atom.fz[n] = atom.pe[n] = 0.0;

  for (int i = 0; i < atom.number; ++i) {
    const double xi = atom.x[i];
    const double yi = atom.y[i];
    const double zi = atom.z[i];
    for (int jj = 0; jj < atom.NN[i]; ++jj) {
      const int j = atom.NL[i * atom.MN + jj];
      double xij = atom.x[j] - xi;
      double yij = atom.y[j] - yi;
      double zij = atom.z[j] - zi;
      applyMic(atom.box, xij, yij, zij);
      const double r2 = xij * xij + yij * yij + zij * zij;
      if (r2 > cutoffSquare)
        continue;

      const double r2inv = 1.0 / r2;
      const double r4inv = r2inv * r2inv;
      const double r6inv = r2inv * r4inv;
      const double r8inv = r4inv * r4inv;
      const double r12inv = r4inv * r8inv;
      const double r14inv = r6inv * r8inv;
      const double f_ij = e24s6 * r8inv - e48s12 * r14inv;
      atom.pe[i] += e4s12 * r12inv - e4s6 * r6inv;
      atom.fx[i] += f_ij * xij;
      atom.fx[j] -= f_ij * xij;
      atom.fy[i] += f_ij * yij;
      atom.fy[j] -= f_ij * yij;
      atom.fz[i] += f_ij * zij;
      atom.fz[j] -= f_ij * zij;
    }
  }
}

void integrate(const bool isStepOne, const double timeStep, Atom& atom)
{
  const double timeStepHalf = timeStep * 0.5;
  for (int n = 0; n < atom.number; ++n) {
    const double mass_inv = 1.0 / atom.mass[n];
    const double ax = atom.fx[n] * mass_inv;
    const double ay = atom.fy[n] * mass_inv;
    const double az = atom.fz[n] * mass_inv;
    atom.vx[n] += ax * timeStepHalf;
    atom.vy[n] += ay * timeStepHalf;
    atom.vz[n] += az * timeStepHalf;
    if (isStepOne) {
      atom.x[n] += atom.vx[n] * timeStep;
      atom.y[n] += atom.vy[n] * timeStep;
      atom.z[n] += atom.vz[n] * timeStep;
    }
  }
}

int main(int argc, char** argv)
{
  if (argc != 5) {
    printf("usage: %s numCells numSteps temperature timeStep\n", argv[0]);
    exit(1);
  }
  const int numCells = atoi(argv[1]);
  const int numSteps = atoi(argv[2]);
  const double temperature = atof(argv[3]);
  double timeStep = atof(argv[4]);
  timeStep /= TIME_UNIT_CONVERSION; // from fs to natural unit

  Atom atom;
  allocateMemory(numCells, atom);
  initializePosition(numCells, atom);
  initializeVelocity(temperature, atom);

  const clock_t tStart = clock();
  std::ofstream ofile("energy.txt");
  ofile << std::fixed << std::setprecision(16);

  for (int step = 0; step < numSteps; ++step) {
    findNeighbor(atom);
    integrate(true, timeStep, atom);  // step 1 in the book
    findForce(atom);                  // step 2 in the book
    integrate(false, timeStep, atom); // step 3 in the book
    if (step % Ns == 0) {
      ofile << findKineticEnergy(atom) << " "
            << std::accumulate(atom.pe.begin(), atom.pe.end(), 0.0)
            << std::endl;
    }
  }
  ofile.close();
  const clock_t tStop = clock();
  const float tElapsed = float(tStop - tStart) / CLOCKS_PER_SEC;
  std::cout << atom.numUpdates << " neighbor list updates" << std::endl;
  std::cout << "Time used = " << tElapsed << " s" << std::endl;

  return 0;
}
