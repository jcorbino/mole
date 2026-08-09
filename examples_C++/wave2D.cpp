/**
 * Solves the 2D wave equation with the position Verlet algorithm.
 * C++ counterpart of examples_MATLAB/wave2D.m
 *
 * Uses both interpolators: Interpol maps cell centres to faces, and the
 * second-type Interpol(true, ...) maps faces back to centres, which is what
 * lets the position and velocity live on their natural staggered locations.
 *
 * Prints the membrane as a matrix:
 *   gnuplot> set pm3d map
 *   gnuplot> splot 'solution' matrix
 */

#include "mole.h"
#include <iostream>

using namespace std;
using namespace arma;

int main() {
  int k = 2;  // Order of accuracy
  int m = 50; // Number of cells along the x-axis
  int n = m;  // Number of cells along the y-axis
  Real a = 0; // West
  Real b = 1; // East
  Real c = 0; // South
  Real d = 1; // North
  Real dx = (b - a) / m;
  Real dy = (d - c) / n;

  // Mimetic operators
  Laplacian L(k, m, n, dx, dy);
  RobinBC BC(k, m, dx, n, dy, 1, 0); // Dirichlet BC
  sp_mat Lsp = (sp_mat)L + (sp_mat)BC;

  Interpol I(m, n, 0.5, 0.5);        // centres -> faces
  Interpol I2(true, m, n, 0.5, 0.5); // faces -> centres

  Real speed = 1;  // (T/p) Tension over density
  Real TIME = 1;   // Simulation time
  Real dt = dx / (2 * speed); // CFL

  // 2D staggered grid
  vec xgrid(m + 2), ygrid(n + 2);
  xgrid(0) = a;
  for (int i = 1; i <= m; ++i)
    xgrid(i) = a + (i - 0.5) * dx;
  xgrid(m + 1) = b;
  ygrid(0) = c;
  for (int j = 1; j <= n; ++j)
    ygrid(j) = c + (j - 0.5) * dy;
  ygrid(n + 1) = d;

  // Initial condition: position sin(pi x) sin(pi y), velocity zero
  vec uold((m + 2) * (n + 2));
  for (int j = 0; j < n + 2; ++j)
    for (int i = 0; i < m + 2; ++i)
      uold(j * (m + 2) + i) =
          sin(datum::pi * xgrid(i)) * sin(datum::pi * ygrid(j));

  vec vold(2 * m * n + m + n, fill::zeros);

  // Premultiply I and I2
  sp_mat Idt = dt * (sp_mat)I;
  sp_mat I2dt = 0.5 * dt * (sp_mat)I2;

  for (int t = 1; t <= (int)(TIME / dt); ++t) {
    // Position Verlet
    uold = uold + I2dt * vold;
    vec vnew = vold + Idt * (speed * speed * (Lsp * uold));
    vec unew = uold + I2dt * vnew;

    uold = unew;
    vold = vnew;
  }

  for (int j = 0; j < n + 2; ++j) {
    for (int i = 0; i < m + 2; ++i)
      cout << uold(j * (m + 2) + i) << (i == m + 1 ? "" : " ");
    cout << "\n";
  }

  return 0;
}
