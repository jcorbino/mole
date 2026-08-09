/**
 * 3D staggering example using a 3D mimetic Laplacian.
 * C++ counterpart of examples_MATLAB/elliptic3D.m
 *
 * Prints one page of the cube as a matrix:
 *   gnuplot> plot 'solution' matrix with image
 */

#include "mole.h"
#include <iostream>

using namespace std;
using namespace arma;

int main() {
  int k = 2; // Order of accuracy
  int m = 5;
  int n = 6;
  int o = 7;

  Laplacian L(k, m, n, o, 1, 1, 1);
  RobinBC BC(k, m, 1, n, 1, o, 1, 1, 0); // Dirichlet BC
  sp_mat A = (sp_mat)L + (sp_mat)BC;

  // Known value at the cube's front face, the p = 0 plane.
  // MOLE flattens with x fastest, then y, then z.
  vec rhs((m + 2) * (n + 2) * (o + 2), fill::zeros);
  for (int j = 0; j < n + 2; ++j)
    for (int i = 0; i < m + 2; ++i)
      rhs(j * (m + 2) + i) = 100;

  // Solve the system of linear equations
#ifdef EIGEN
  // Use Eigen only if SuperLU (faster) is not available
  vec sol = Utils::spsolve_eigen(A, rhs);
#else
  vec sol = spsolve(A, rhs); // Will use SuperLU
#endif

  const int p = 1; // Page to be displayed (0-based; the MATLAB script shows 2)

  for (int j = 0; j < n + 2; ++j) {
    for (int i = 0; i < m + 2; ++i)
      cout << sol(p * (m + 2) * (n + 2) + j * (m + 2) + i)
           << (i == m + 1 ? "" : " ");
    cout << "\n";
  }

  return 0;
}
