/**
 * Solves the 1D advection equation with periodic boundary conditions, using
 * the leapfrog scheme.
 * C++ counterpart of examples_MATLAB/hyperbolic1D.m
 *
 * The periodic operators are required here, not merely convenient. Leapfrog's
 * region of absolute stability is the open segment (-i, i) of the imaginary
 * axis and nothing else, so it tolerates no real part at all: with the default
 * interpolator the boundary rows read the phantom boundary unknowns and leave
 * max Re lambda(D*I) slightly positive, which grows however small dt is.
 *
 * Prints "x u" per line:
 *   gnuplot> plot 'solution' with linespoints
 */

#include "mole.h"
#include <iostream>

using namespace std;
using namespace arma;

int main() {
  Real a = 1;    // Velocity
  Real west = 0; // Domain's limits
  Real east = 1;

  int k = 2;  // Operator's order of accuracy
  int m = 50; // Number of cells
  Real dx = (east - west) / m;

  Real t = 1; // Simulation time

  Divergence D(k, m, dx, true); // periodic
  Interpol I(m, 0.5, true);     // periodic, centered

  sp_mat DI = (sp_mat)D * (sp_mat)I;

  // Spectral radius of the spatial operator
  cx_vec eigval;
  eig_gen(eigval, mat(DI));
  Real rho = max(abs(eigval));

  // Leapfrog is stable only while every dt*eig(-a*D*I) stays on the OPEN
  // segment (-i, i), so dt < 1/(|a|*rho). The 0.9 keeps us inside it: at
  // dt = 1/(|a|*rho) exactly the two characteristic roots coalesce at +-i, and
  // that defective double root makes the solution grow with the step count.
  Real dt = 0.9 / (fabs(a) * rho);

  // Shrink dt so the last iteration lands exactly on t. Rounding the step count
  // up can only make dt smaller, so this never violates the CFL bound above.
  int n_steps = (int)ceil(t / dt);
  dt = t / n_steps;

  // 1D staggered grid
  vec grid(m + 2);
  grid(0) = west;
  for (int i = 1; i <= m; ++i)
    grid(i) = west + (i - 0.5) * dx;
  grid(m + 1) = east;

  // IC
  vec U(m + 2);
  for (int i = 0; i < m + 2; ++i)
    U(i) = sin(2 * datum::pi * grid(i));

  // Premultiply out of the time loop, since it does not change
  sp_mat A = -a * dt * 2 * DI;

  vec U2 = U + 0.5 * (A * U); // one Euler step to start
  for (int i = 1; i <= n_steps; ++i) {
    vec U3 = U + A * U2;      // leapfrog
    U = U2;
    U2 = U3;
  }

  for (int i = 0; i < m + 2; ++i)
    cout << grid(i) << " " << U2(i) << "\n";

  return 0;
}
