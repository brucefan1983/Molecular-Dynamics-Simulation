/*----------------------------------------------------------------------------80
    Copyright 2022 Zheyong Fan
Compile:
    g++ ljmd.cpp -O3 -o ljmd
Run:
    ljmd numCells numSteps temperature timeStep
    such as
        ljmd 4 20000 60 5
------------------------------------------------------------------------------*/

#include <cmath>    // sqrt() function
#include <ctime>    // for timing
#include <fstream>  // file
#include <iomanip>  // std::setprecision
#include <iostream> // input/output
#include <iterator>
#include <numeric> // std::accumulate
#include <sstream> // std::istringstream
#include <string>  // string
#include <vector>  // vector

const int Ns = 100;             // output frequency
const double K_B = 8.617343e-5; // Boltzmann's constant in natural unit
const double TIME_UNIT_CONVERSION = 1.018051e+1; // from natural unit to fs

struct Atom {
  int number;
  double box[6];
  std::vector<double> mass, x, y, z, vx, vy, vz, fx, fy, fz, pe;
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
  atom.mass.resize(atom.number, 40.0); // argon mass
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

void initializePosition(const int numCells, Atom& atom)
{
  const int numAtomsPerCell = 4;
  const double latticeConstant = 5.385;
  atom.box[0] = atom.box[1] = atom.box[2] = latticeConstant * numCells;
  atom.box[3] = atom.box[0] * 0.5;
  atom.box[4] = atom.box[1] * 0.5;
  atom.box[5] = atom.box[2] * 0.5;
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

  for (int i = 0; i < atom.number - 1; ++i) {
    for (int j = i + 1; j < atom.number; ++j) {
      double xij = atom.x[j] - atom.x[i];
      double yij = atom.y[j] - atom.y[i];
      double zij = atom.z[j] - atom.z[i];
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

std::vector<std::string> getTokens(std::ifstream& input)
{
  std::string line;
  std::getline(input, line);
  std::istringstream iss(line);
  std::vector<std::string> tokens{
    std::istream_iterator<std::string>{iss},
    std::istream_iterator<std::string>{}};
  return tokens;
}

int getInt(std::string& token)
{
  int value = 0;
  try {
    value = std::stoi(token);
  } catch (const std::exception& e) {
    std::cout << "Standard exception:" << e.what() << std::endl;
    exit(1);
  }
  return value;
}

double getDouble(std::string& token)
{
  float value = 0;
  try {
    value = std::stod(token);
  } catch (const std::exception& e) {
    std::cout << "Standard exception:" << e.what() << std::endl;
    exit(1);
  }
  return value;
}

void readRun(int& numSteps, double& timeStep, double& temperature)
{
  std::ifstream input("run.in");
  if (!input.is_open()) {
    std::cout << "Failed to open run.in." << std::endl;
    exit(1);
  }

  while (input.peek() != EOF) {
    std::vector<std::string> tokens = getTokens(input);
    if (tokens.size() > 0) {
      if (tokens[0] == "time_step") {
        timeStep = getDouble(tokens[1]);
        if (timeStep < 0) {
          std::cout << "timeStep should >= 0." << std::endl;
          exit(1);
        }
        std::cout << "timeStep = " << timeStep << " fs." << std::endl;
      } else if (tokens[0] == "run") {
        numSteps = getInt(tokens[1]);
        if (numSteps < 1) {
          std::cout << "numSteps should >= 1." << std::endl;
          exit(1);
        }
        std::cout << "numSteps = " << numSteps << std::endl;
      } else if (tokens[0] == "velocity") {
        temperature = getDouble(tokens[1]);
        if (temperature < 0) {
          std::cout << "temperature >= 0." << std::endl;
          exit(1);
        }
        std::cout << "temperature = " << temperature << " K." << std::endl;
      } else if (tokens[0][0] != '#') {
        std::cout << tokens[0] << " is not a valid keyword." << std::endl;
        exit(1);
      }
    }
  }

  input.close();
}

int main(int argc, char** argv)
{
  if (argc != 2) {
    printf("usage: %s numCells\n", argv[0]);
    exit(1);
  }
  const int numCells = atoi(argv[1]);
  int numSteps;
  double temperature;
  double timeStep;
  readRun(numSteps, timeStep, temperature);

  timeStep /= TIME_UNIT_CONVERSION; // from fs to natural unit

  Atom atom;
  allocateMemory(numCells, atom);
  initializePosition(numCells, atom);
  initializeVelocity(temperature, atom);

  const clock_t tStart = clock();
  std::ofstream ofile("energy.txt");
  ofile << std::fixed << std::setprecision(16);

  for (int step = 0; step < numSteps; ++step) {
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
  std::cout << "Time used = " << tElapsed << " s" << std::endl;

  return 0;
}
