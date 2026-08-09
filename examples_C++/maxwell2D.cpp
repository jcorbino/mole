/**
 * 2D TMz Maxwell solver using mimetic operators and leapfrog.
 * C++ counterpart of examples_MATLAB/maxwell2D.m
 *
 * E lives at cell centres, B on the faces, so Gradient takes E to B's home and
 * Divergence brings B back -- the staggering the scheme wants, for free.
 *
 * Prints the electric field as a matrix:
 *   gnuplot> set pm3d map
 *   gnuplot> splot 'solution' matrix
 */

#include "mole.h"
#include <iostream>

using namespace std;
using namespace arma;

int main() {
  int nt = 500; // Number of time steps
  int mx = 100;
  int my = 100;

  Real west = 0, east = 1;
  Real south = 0, north = 1;

  Real dx = (east - west) / mx;
  Real dy = (north - south) / my;
  Real dt = 0.5 * std::min(dx, dy); // CFL

  // Grids
  vec xE(mx + 2), yE(my + 2);
  xE(0) = west;
  for (int i = 1; i <= mx; ++i)
    xE(i) = west + (i - 0.5) * dx;
  xE(mx + 1) = east;
  yE(0) = south;
  for (int j = 1; j <= my; ++j)
    yE(j) = south + (j - 0.5) * dy;
  yE(my + 1) = north;

  // Mimetic operators
  int k = 2;
  Gradient Gop(k, mx, my, dx, dy);
  Divergence Dop(k, mx, my, dx, dy);
  sp_mat G = dt * (sp_mat)Gop;
  sp_mat D = dt * (sp_mat)Dop;

  // Field sizes
  int NB = (mx + 1) * my + mx * (my + 1);

  // Initial condition: Gaussian pulse
  vec E((mx + 2) * (my + 2));
  for (int j = 0; j < my + 2; ++j)
    for (int i = 0; i < mx + 2; ++i)
      E(j * (mx + 2) + i) =
          exp(-100 * (pow(xE(i) - 0.5, 2) + pow(yE(j) - 0.5, 2)));

  // Magnetic field
  vec B(NB, fill::zeros);

  // Leapfrog start (half step for B)
  B = B - 0.5 * (G * E);

  for (int n = 0; n < nt; ++n) {
    B = B - G * E;
    E = E - D * B;
  }

  for (int j = 0; j < my + 2; ++j) {
    for (int i = 0; i < mx + 2; ++i)
      cout << E(j * (mx + 2) + i) << (i == mx + 1 ? "" : " ");
    cout << "\n";
  }

  return 0;
}
