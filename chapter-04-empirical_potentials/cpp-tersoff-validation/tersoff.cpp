/*----------------------------------------------------------------------------80
    Copyright 2022 Zheyong Fan
Compile:
    g++ tersoff.cpp -O3 -o tersoff
Run:
    path/to/tersoff
Inputs:
    none
Outputs:
    force.out
------------------------------------------------------------------------------*/

#include <cmath>    // sqrt() function
#include <ctime>    // for timing
#include <fstream>  // file
#include <iomanip>  // std::setprecision
#include <iostream> // input/output
#include <iterator>
#include <sstream> // std::istringstream
#include <string>  // string
#include <vector>  // vector

struct Atom {
  int number;
  const int MN = 4;
  const double cutoffNeighbor = 3.1;
  double box[6];
  double pe;
  std::vector<int> NN, NL;
  std::vector<double> x0, y0, z0, x, y, z, fx, fy, fz, b,
    bp;
};

void applyMicOne(
  const double length, const double halfLength, double& x12)
{
  if (x12 < -halfLength)
    x12 += length;
  else if (x12 > +halfLength)
    x12 -= length;
}

void applyMic(
  const double box[6],
  double& x12,
  double& y12,
  double& z12)
{
  applyMicOne(box[0], box[3], x12);
  applyMicOne(box[1], box[4], y12);
  applyMicOne(box[2], box[5], z12);
}

inline void
find_fr_and_frp(double d12, double& fr, double& frp)
{
  const double a = 1.8308e3;
  const double lambda = 2.4799;
  fr = a * exp(-lambda * d12);
  frp = -lambda * fr;
}

inline void
find_fa_and_fap(double d12, double& fa, double& fap)
{
  const double b = 471.18;
  const double mu = 1.7322;
  fa = b * exp(-mu * d12);
  fap = -mu * fa;
}

inline void find_fa(double d12, double& fa)
{
  const double b = 471.18;
  const double mu = 1.7322;
  fa = b * exp(-mu * d12);
}

inline void
find_fc_and_fcp(double d12, double& fc, double& fcp)
{
  const double r1 = 2.7;
  const double r2 = 3.0;
  const double pi = 3.141592653589793;
  const double pi_factor = pi / (r2 - r1);
  if (d12 < r1) {
    fc = 1.0;
    fcp = 0.0;
  } else if (d12 < r2) {
    fc = cos(pi_factor * (d12 - r1)) * 0.5 + 0.5;
    fcp = -sin(pi_factor * (d12 - r1)) * pi_factor * 0.5;
  } else {
    fc = 0.0;
    fcp = 0.0;
  }
}

inline void find_fc(double d12, double& fc)
{
  const double r1 = 2.7;
  const double r2 = 3.0;
  const double pi = 3.141592653589793;
  const double pi_factor = pi / (r2 - r1);
  if (d12 < r1) {
    fc = 1.0;
  } else if (d12 < r2) {
    fc = cos(pi_factor * (d12 - r1)) * 0.5 + 0.5;
  } else {
    fc = 0.0;
  }
}

inline void find_g_and_gp(double cos, double& g, double& gp)
{
  const double c = 1.0039e5;
  const double d = 16.217;
  const double h = -0.59825;
  const double c2 = c * c;
  const double d2 = d * d;
  const double c2overd2 = c2 / d2;
  double temp = d2 + (cos - h) * (cos - h);
  g = 1.0 + c2overd2 - c2 / temp;
  gp = 2.0 * c2 * (cos - h) / (temp * temp);
}

inline void find_g(double cos, double& g)
{
  const double c = 1.0039e5;
  const double d = 16.217;
  const double h = -0.59825;
  const double c2 = c * c;
  const double d2 = d * d;
  const double c2overd2 = c2 / d2;
  double temp = d2 + (cos - h) * (cos - h);
  g = 1.0 + c2overd2 - c2 / temp;
}

void findNeighborON2(Atom& atom)
{
  const double cutoffSquare =
    atom.cutoffNeighbor * atom.cutoffNeighbor;
  std::fill(atom.NN.begin(), atom.NN.end(), 0);

  for (int i = 0; i < atom.number; ++i) {
    const double x1 = atom.x[i];
    const double y1 = atom.y[i];
    const double z1 = atom.z[i];
    for (int j = 0; j < atom.number; ++j) {
      if (j == i)
        continue;
      double xij = atom.x[j] - x1;
      double yij = atom.y[j] - y1;
      double zij = atom.z[j] - z1;
      applyMic(atom.box, xij, yij, zij);
      const double distanceSquare =
        xij * xij + yij * yij + zij * zij;
      if (distanceSquare < cutoffSquare) {
        atom.NL[i * atom.MN + atom.NN[i]++] = j;
        if (atom.NN[i] > atom.MN) {
          std::cout
            << "Error: number of neighbors for atom " << i
            << " exceeds " << atom.MN << std::endl;
          exit(1);
        }
      }
    }
  }
}

void find_b_and_bp(Atom& atom)
{
  const double beta = 1.1000e-6;
  const double n = 0.78734;
  const double minus_half_over_n = -0.5 / n;
  for (int n1 = 0; n1 < atom.number; ++n1) {
    for (int i1 = 0; i1 < atom.NN[n1]; ++i1) {
      int n2 = atom.NL[n1 * atom.MN + i1];
      double x12, y12, z12;
      x12 = atom.x[n2] - atom.x[n1];
      y12 = atom.y[n2] - atom.y[n1];
      z12 = atom.z[n2] - atom.z[n1];
      applyMic(atom.box, x12, y12, z12);
      double d12 = sqrt(x12 * x12 + y12 * y12 + z12 * z12);

      double zeta = 0.0;
      for (int i2 = 0; i2 < atom.NN[n1]; ++i2) {
        int n3 = atom.NL[n1 * atom.MN + i2];
        if (n3 == n2) {
          continue;
        }
        double x13, y13, z13;
        x13 = atom.x[n3] - atom.x[n1];
        y13 = atom.y[n3] - atom.y[n1];
        z13 = atom.z[n3] - atom.z[n1];
        applyMic(atom.box, x13, y13, z13);

        double d13 =
          sqrt(x13 * x13 + y13 * y13 + z13 * z13);
        double cos =
          (x12 * x13 + y12 * y13 + z12 * z13) / (d12 * d13);
        double fc13, g123;
        find_fc(d13, fc13);
        find_g(cos, g123);
        zeta += fc13 * g123;
      }
      double bzn = pow(beta * zeta, n);
      double b12 = pow(1.0 + bzn, minus_half_over_n);
      atom.b[n1 * atom.MN + i1] = b12;
      atom.bp[n1 * atom.MN + i1] =
        -b12 * bzn * 0.5 / ((1.0 + bzn) * zeta);
    }
  }
}

void find_force_tersoff(Atom& atom)
{
  atom.pe = 0.0;
  for (int n = 0; n < atom.number; ++n) {
    atom.fx[n] = atom.fy[n] = atom.fz[n] = 0.0;
  }

  for (int n1 = 0; n1 < atom.number; ++n1) {
    for (int i1 = 0; i1 < atom.NN[n1]; ++i1) {
      int n2 = atom.NL[n1 * atom.MN + i1];
      double x12, y12, z12;
      x12 = atom.x[n2] - atom.x[n1];
      y12 = atom.y[n2] - atom.y[n1];
      z12 = atom.z[n2] - atom.z[n1];
      applyMic(atom.box, x12, y12, z12);

      double d12 = sqrt(x12 * x12 + y12 * y12 + z12 * z12);
      double d12inv = 1.0 / d12;
      double d12inv_square = d12inv * d12inv;

      double fc12, fcp12;
      double fa12, fap12;
      double fr12, frp12;
      find_fc_and_fcp(d12, fc12, fcp12);
      find_fa_and_fap(d12, fa12, fap12);
      find_fr_and_frp(d12, fr12, frp12);

      double b12 = atom.b[n1 * atom.MN + i1];
      double bp12 = atom.bp[n1 * atom.MN + i1];

      double f12[3] = {0.0, 0.0, 0.0};
      double factor1 = -b12 * fa12 + fr12;
      double factor2 = -b12 * fap12 + frp12;
      double factor3 =
        (fcp12 * factor1 + fc12 * factor2) / d12;
      f12[0] += x12 * factor3 * 0.5;
      f12[1] += y12 * factor3 * 0.5;
      f12[2] += z12 * factor3 * 0.5;
      
      for (int i2 = 0; i2 < atom.NN[n1]; ++i2) {
        int n3 = atom.NL[n1 * atom.MN + i2];
        if (n3 == n2) {
          continue;
        }
        double x13 = atom.x[n3] - atom.x[n1];
        double y13 = atom.y[n3] - atom.y[n1];
        double z13 = atom.z[n3] - atom.z[n1];
        applyMic(atom.box, x13, y13, z13);

        double d13 =
          sqrt(x13 * x13 + y13 * y13 + z13 * z13);
        double fc13, fa13;
        find_fc(d13, fc13);
        find_fa(d13, fa13);
        double bp13 = atom.bp[n1 * atom.MN + i2];

        double cos123 =
          (x12 * x13 + y12 * y13 + z12 * z13) / (d12 * d13);
        double g123, gp123;
        find_g_and_gp(cos123, g123, gp123);
        double cos_x =
          x13 / (d12 * d13) - x12 * cos123 / (d12 * d12);
        double cos_y =
          y13 / (d12 * d13) - y12 * cos123 / (d12 * d12);
        double cos_z =
          z13 / (d12 * d13) - z12 * cos123 / (d12 * d12);
        double factor123a = (-bp12 * fc12 * fa12 * fc13 -
                             bp13 * fc13 * fa13 * fc12) *
                            gp123;
        double factor123b =
          -bp13 * fc13 * fa13 * fcp12 * g123 * d12inv;
        f12[0] +=
          (x12 * factor123b + factor123a * cos_x) * 0.5;
        f12[1] +=
          (y12 * factor123b + factor123a * cos_y) * 0.5;
        f12[2] +=
          (z12 * factor123b + factor123a * cos_z) * 0.5;
      }

      atom.pe += factor1 * fc12 * 0.5;
      atom.fx[n1] += f12[0];
      atom.fy[n1] += f12[1];
      atom.fz[n1] += f12[2];
      atom.fx[n2] -= f12[0];
      atom.fy[n2] -= f12[1];
      atom.fz[n2] -= f12[2];
    }
  }
}

void findForce(Atom& atom)
{
  find_b_and_bp(atom);
  find_force_tersoff(atom);
}

void createXyz(Atom& atom)
{
  const int M = 8;
  double x0[M] = {0.0,  0.0,  0.5,  0.5,
                  0.25, 0.25, 0.75, 0.75};
  double y0[M] = {0.0,  0.5,  0.0,  0.5,
                  0.25, 0.75, 0.25, 0.75};
  double z0[M] = {0.0,  0.5,  0.5,  0.0,
                  0.25, 0.75, 0.75, 0.25};
  const int num_cells = 2;
  atom.number = num_cells * num_cells * num_cells * M;
  double a = 5.45;
  atom.box[0] = atom.box[1] = atom.box[2] = a * num_cells;
  atom.box[3] = atom.box[4] = atom.box[5] =
    a * num_cells * 0.5;
  atom.x0.resize(atom.number, 0.0);
  atom.y0.resize(atom.number, 0.0);
  atom.z0.resize(atom.number, 0.0);
  atom.x.resize(atom.number, 0.0);
  atom.y.resize(atom.number, 0.0);
  atom.z.resize(atom.number, 0.0);
  atom.fx.resize(atom.number, 0.0);
  atom.fy.resize(atom.number, 0.0);
  atom.fz.resize(atom.number, 0.0);
  atom.b.resize(atom.number * atom.MN, 0.0);
  atom.bp.resize(atom.number * atom.MN, 0.0);
  atom.NN.resize(atom.number, 0);
  atom.NL.resize(atom.number * atom.MN, 0);

  int n = 0;
  for (int nx = 0; nx < num_cells; ++nx) {
    for (int ny = 0; ny < num_cells; ++ny) {
      for (int nz = 0; nz < num_cells; ++nz) {
        for (int m = 0; m < M; ++m) {
          atom.x[n] = atom.x0[n] =
            (x0[m] + nx) * a +
            (double(rand()) / RAND_MAX - 0.5) * 0.5;
          atom.y[n] = atom.y0[n] =
            (y0[m] + ny) * a +
            (double(rand()) / RAND_MAX - 0.5) * 0.5;
          atom.z[n] = atom.z0[n] =
            (z0[m] + nz) * a +
            (double(rand()) / RAND_MAX - 0.5) * 0.5;
          ++n;
        }
      }
    }
  }

  std::ofstream ofile("model.xyz");
  ofile << std::fixed << std::setprecision(16);
  ofile << atom.number << std::endl;
  ofile << "Lattice=\"" << atom.box[0] << " 0 0 0 "
        << atom.box[1] << " 0 0 0 " << atom.box[2] << "\" ";
  ofile << "Properties=species:S:1:pos:R:3" << std::endl;
  for (int n = 0; n < atom.number; ++n) {
    ofile << "Si " << atom.x[n] << " " << atom.y[n] << " "
          << atom.z[n] << std::endl;
  }
  ofile.close();

  std::cout << "model.xyz is created." << std::endl;
}

int main(int argc, char** argv)
{
  Atom atom;
  createXyz(atom);
  findNeighborON2(atom);
  std::cout << "neighbor list is built." << std::endl;

  std::ofstream ofile("force.out");
  ofile << std::fixed << std::setprecision(16);

  findForce(atom);
  for (int n = 0; n < atom.number; ++n) {
    ofile << atom.fx[n] << " " << atom.fy[n] << " "
          << atom.fz[n] << std::endl;
  }
  std::cout << "analytical force is calculated."
            << std::endl;

  const double delta = 2.0e-5;
  for (int n = 0; n < atom.number; ++n) {
    atom.x[n] = atom.x0[n] + delta;
    findForce(atom);
    double pePositive = atom.pe;
    atom.x[n] = atom.x0[n] - delta;
    findForce(atom);
    double peNegative = atom.pe;
    atom.x[n] = atom.x0[n];
    ofile << (peNegative - pePositive) / (delta * 2.0)
          << " ";

    atom.y[n] = atom.y0[n] + delta;
    findForce(atom);
    pePositive = atom.pe;
    atom.y[n] = atom.y0[n] - delta;
    findForce(atom);
    peNegative = atom.pe;
    atom.y[n] = atom.y0[n];
    ofile << (peNegative - pePositive) / (delta * 2.0)
          << " ";

    atom.z[n] = atom.z0[n] + delta;
    findForce(atom);
    pePositive = atom.pe;
    atom.z[n] = atom.z0[n] - delta;
    findForce(atom);
    peNegative = atom.pe;
    atom.z[n] = atom.z0[n];
    ofile << (peNegative - pePositive) / (delta * 2.0)
          << std::endl;
  }

  std::cout << "finite-difference force is calculated."
            << std::endl;

  ofile.close();

  return 0;
}
