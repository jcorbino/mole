/**
 * Solves the 1D heat equation with Dirichlet boundary conditions.
 * C++ counterpart of examples_MATLAB/parabolic1D.m
 *
 * The Laplacian's boundary rows are empty, so the two boundary entries keep
 * whatever the initial condition put there -- that is what holds the Dirichlet
 * values fixed here, with no BC operator involved.
 *
 * Prints "x T" per line:
 *   gnuplot> plot 'solution' with linespoints
 */

#include "mole.h"
#include <cmath>
#include <iostream>

using namespace std;
using namespace arma;

int main() {
  Real alpha = 1;    // Thermal diffusivity
  Real west = 0;     // Domain's limits
  Real east = 1;
  int k = 2;         // Operator's order of accuracy
  int m = 2 * k + 1; // Minimum number of cells to attain the desired accuracy
  Real dx = (east - west) / m;

  Real t = 1;                    // Simulation time
  Real dt = dx * dx / (3 * alpha); // von Neumann stability criterion

  Laplacian L(k, m, dx);

  // Explicit scheme: one operator application per step
  sp_mat A = alpha * dt * (sp_mat)L + speye<sp_mat>(m + 2, m + 2);

  vec U(m + 2, fill::zeros);
  U(0) = 100;     // BC
  U(m + 1) = 100;

  // 1D staggered grid
  vec grid(m + 2);
  grid(0) = west;
  for (int i = 1; i <= m; ++i)
    grid(i) = west + (i - 0.5) * dx;
  grid(m + 1) = east;

  // MATLAB's colon operator treats a limit that is a hair under an integer as
  // that integer, so "0 : t/dt+1" runs one step further than a plain cast
  // would. Match it, or the run stops one step early.
  Real last = t / dt + 1;
  int nsteps = (int)std::floor(last);
  if (last - nsteps > 1.0 - 1e-8)
    nsteps += 1;

  for (int i = 0; i <= nsteps; ++i)
    U = A * U;

  for (int i = 0; i < m + 2; ++i)
    cout << grid(i) << " " << U(i) << "\n";

  return 0;
}
