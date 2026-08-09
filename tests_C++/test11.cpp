/**
 * Tests for the 3D operators, which the rest of the suite does not reach
 *
 *   1. Correct dimensions for Divergence, Gradient and Laplacian
 *   2. Constants are annihilated
 *   3. div(grad(.)) agrees with the Laplacian constructor
 *   4. The 3D Laplacian recovers lap(x^2 + y^2 + z^2) = 6 on the interior,
 *      which is exact for every order of accuracy
 */

#include "mole.h"
#include <iostream>

static void fail() {
  cout << "\033[1;31mTest FAILED!\033[0m\n";
  exit(1);
}

int main() {
  Real tol = 1e-10;

  for (int k : {2, 4, 6}) {
    int m = 2 * k + 1;
    int n = 2 * k + 2;
    int o = 2 * k + 3;
    Real dx = 0.7, dy = 1.3, dz = 0.9;

    Divergence D(k, m, n, o, dx, dy, dz);
    Gradient G(k, m, n, o, dx, dy, dz);
    Laplacian L(k, m, n, o, dx, dy, dz);

    int nfaces = (m + 1) * n * o + m * (n + 1) * o + m * n * (o + 1);
    int ncells = (m + 2) * (n + 2) * (o + 2);

    // ---- dimensions ----
    if (D.n_rows != (u32)ncells || D.n_cols != (u32)nfaces)
      fail();
    if (G.n_rows != (u32)nfaces || G.n_cols != (u32)ncells)
      fail();
    if (L.n_rows != (u32)ncells || L.n_cols != (u32)ncells)
      fail();

    // ---- nullity ----
    if (norm((sp_mat)D * vec(nfaces, fill::ones)) > tol)
      fail();
    if (norm((sp_mat)G * vec(ncells, fill::ones)) > tol)
      fail();
    if (norm((sp_mat)L * vec(ncells, fill::ones)) > tol)
      fail();

    // ---- L is exactly div(grad(.)) ----
    sp_mat DG = (sp_mat)D * (sp_mat)G;
    if (norm(mat(DG - (sp_mat)L), "inf") > tol)
      fail();
  }

  // ---- lap(x^2 + y^2 + z^2) = 6, exact at any order ----
  {
    int k = 4;
    int m = 12, n = 12, o = 12;
    Real dx = 1.0 / m, dy = 1.0 / n, dz = 1.0 / o;

    Laplacian L(k, m, n, o, dx, dy, dz);

    // Cell-space coordinates: boundary node, m centres, boundary node
    vec xc(m + 2), yc(n + 2), zc(o + 2);
    xc(0) = 0.0;
    for (int i = 1; i <= m; ++i)
      xc(i) = (i - 0.5) * dx;
    xc(m + 1) = 1.0;
    yc(0) = 0.0;
    for (int j = 1; j <= n; ++j)
      yc(j) = (j - 0.5) * dy;
    yc(n + 1) = 1.0;
    zc(0) = 0.0;
    for (int p = 1; p <= o; ++p)
      zc(p) = (p - 0.5) * dz;
    zc(o + 1) = 1.0;

    // x running fastest, then y, then z
    vec f((m + 2) * (n + 2) * (o + 2));
    for (int p = 0; p < o + 2; ++p)
      for (int j = 0; j < n + 2; ++j)
        for (int i = 0; i < m + 2; ++i)
          f(p * (m + 2) * (n + 2) + j * (m + 2) + i) =
              xc(i) * xc(i) + yc(j) * yc(j) + zc(p) * zc(p);

    vec r = (sp_mat)L * f;

    // Interior only: the boundary rows are placeholders for a BC operator
    Real worst = 0;
    for (int p = 1; p <= o; ++p)
      for (int j = 1; j <= n; ++j)
        for (int i = 1; i <= m; ++i)
          worst = std::max(
              worst,
              std::abs(r(p * (m + 2) * (n + 2) + j * (m + 2) + i) - 6.0));

    if (worst > 1e-8)
      fail();
  }

  cout << "\033[1;32mTest PASSED!\033[0m\n";

  return 0;
}
