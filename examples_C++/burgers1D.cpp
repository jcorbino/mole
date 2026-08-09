/**
 * Solves the 1D inviscid Burgers' equation, written in conservative form and
 * discretised with an upwind interpolator.
 * C++ counterpart of examples_MATLAB/burgers1D.m
 *
 * Interpol(m, 1) is the upwind choice, which is 1st-order by construction --
 * that truncation term is exactly its numerical diffusion. Use Interpol(m, 0)
 * for downwinding if the fluid propagates to the left.
 *
 * Prints "x u" per line:
 *   gnuplot> plot 'solution' with lines
 */

#include "mole.h"
#include <cmath>
#include <iostream>

using namespace std;
using namespace arma;

int main() {
  Real west = -15; // Domain's limits
  Real east = 15;
  int k = 2;       // Operator's order of accuracy
  int m = 300;     // Number of cells
  Real dx = (east - west) / m;

  Real t = 10;   // Simulation time
  Real dt = dx;  // CFL condition with max(|u|) <= 1

  Divergence D(k, m, dx);
  Interpol I(m, 1); // upwind

  // 1D staggered grid
  vec xgrid(m + 2);
  xgrid(0) = west;
  for (int i = 1; i <= m; ++i)
    xgrid(i) = west + (i - 0.5) * dx;
  xgrid(m + 1) = east;

  // Impose IC
  vec U = exp(-square(xgrid) / 50);

  // Premultiply out of the time loop since it does not change
  sp_mat A = -dt / 2 * ((sp_mat)D * (sp_mat)I);

  // MATLAB's colon operator treats a limit that is a hair under an integer as
  // that integer, so "0 : t/dt" runs one step further than a plain cast
  // would. Match it, or the run stops one step early.
  Real last = t / dt;
  int nsteps = (int)std::floor(last);
  if (last - nsteps > 1.0 - 1e-8)
    nsteps += 1;

  for (int i = 0; i <= nsteps; ++i)
    U = U + A * square(U);

  for (int i = 0; i < m + 2; ++i)
    cout << xgrid(i) << " " << U(i) << "\n";

  return 0;
}
