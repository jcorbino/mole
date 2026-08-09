/**
 * Solves the 1D wave equation with the position Verlet algorithm.
 * C++ counterpart of examples_MATLAB/wave1D.m
 *
 * Prints "x u" per line:
 *   gnuplot> plot 'solution' with linespoints
 */

#include "mole.h"
#include <cmath>
#include <iostream>

using namespace std;
using namespace arma;

int main() {
  int k = 2;    // Order of accuracy (spatial)
  int m = 50;   // Number of cells
  Real a = 0;   // Left boundary
  Real b = 1;   // Right boundary
  Real dx = (b - a) / m;

  Laplacian L(k, m, dx);
  sp_mat Lsp = L;

  Real c = 2;      // Wave propagation speed
  Real TIME = 1;   // Simulation time
  Real dt = dx / (2 * c); // CFL

  // 1D staggered grid
  vec xgrid(m + 2);
  xgrid(0) = a;
  for (int i = 1; i <= m; ++i)
    xgrid(i) = a + (i - 0.5) * dx;
  xgrid(m + 1) = b;

  // Initial condition: position sin(pi x), velocity zero
  vec uold(m + 2), vold(m + 2, fill::zeros);
  for (int i = 0; i < m + 2; ++i)
    uold(i) = sin(datum::pi * xgrid(i));

  // "Force" function, c^2 DivGrad u
  auto F = [&](const vec &x) { return vec(c * c * (Lsp * x)); };

  // MATLAB's colon operator treats a limit that is a hair under an integer as
  // that integer, so "0 : TIME/dt" runs one step further than a plain cast
  // would. Match it, or the run stops one step early.
  Real last = TIME / dt;
  int nsteps = (int)std::floor(last);
  if (last - nsteps > 1.0 - 1e-8)
    nsteps += 1;

  for (int t = 0; t <= nsteps; ++t) {
    // Position Verlet, 2nd-order in time
    uold = uold + 0.5 * dt * vold;
    vec vnew = vold + dt * F(uold);
    vec unew = uold + 0.5 * dt * vnew;

    uold = unew;
    vold = vnew;
  }

  for (int i = 0; i < m + 2; ++i)
    cout << xgrid(i) << " " << uold(i) << "\n";

  return 0;
}
