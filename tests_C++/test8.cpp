/**
 * Conservation test for the PERIODIC operators
 *
 * Summing a divergence over the interior cells must telescope to zero around
 * the seam, so the total is conserved for any field. Equivalently, the interior
 * column sums of the composite operators vanish:
 *
 *   sum_interior (D * G)              = 0   (diffusion conserves)
 *   sum_interior (D * diag(v) * I)    = 0   (advection conserves)
 *
 * The advective flux must itself be periodic: faces 0 and m are the same
 * physical face, so the velocity has to agree there. Sampling any periodic
 * function on the faces gives that for free; an arbitrary vector does not, and
 * then nothing is conserved -- correctly so.
 */

#include "mole.h"
#include <iostream>

static void fail() {
  cout << "\033[1;31mTest FAILED!\033[0m\n";
  exit(1);
}

int main() {
  Real tol = 1e-10;

  // ---- 1D ----
  for (int k : {2, 4, 6}) {
    int m = 30;
    Real L = 2 * datum::pi;
    Real dx = L / m;

    Divergence D(k, m, dx, true);
    Gradient G(k, m, dx, true);
    Interpol I(m, 0.5, true);

    // A periodic velocity, so v(0) == v(m) as the seam requires
    vec v(m + 1);
    for (int i = 0; i <= m; ++i)
      v(i) = 2 + cos(i * dx);

    sp_mat V(m + 1, m + 1);
    V.diag() = v;

    sp_mat diffusion = (sp_mat)D * (sp_mat)G;
    sp_mat advection = (sp_mat)D * V * (sp_mat)I;

    if (norm(sum(mat(diffusion.rows(1, m)), 0), "inf") > tol)
      fail();
    if (norm(sum(mat(advection.rows(1, m)), 0), "inf") > tol)
      fail();
  }

  // ---- 2D: the same statement over the interior cells of the grid ----
  for (int k : {2, 4}) {
    int m = 12, n = 14;
    Real Lx = 2 * datum::pi, Ly = 2 * datum::pi;
    Real dx = Lx / m, dy = Ly / n;

    Divergence D(k, m, n, dx, dy, true);
    Gradient G(k, m, n, dx, dy, true);
    Interpol I(m, n, 0.5, 0.5, true);

    int nu = (m + 1) * n;
    int nv = m * (n + 1);

    // Taylor-Green: divergence free, and periodic on the faces it lives on
    vec vel(nu + nv);
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < m + 1; ++i)
        vel(j * (m + 1) + i) = sin(i * dx) * cos((j + 0.5) * dy);
    for (int j = 0; j < n + 1; ++j)
      for (int i = 0; i < m; ++i)
        vel(nu + j * m + i) = -cos((i + 0.5) * dx) * sin(j * dy);

    sp_mat V(nu + nv, nu + nv);
    V.diag() = vel;

    sp_mat diffusion = (sp_mat)D * (sp_mat)G;
    sp_mat advection = (sp_mat)D * V * (sp_mat)I;

    // Interior cells of the (m+2) x (n+2) layout, x running fastest
    uvec interior(m * n);
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < m; ++i)
        interior(j * m + i) = (j + 1) * (m + 2) + (i + 1);

    // sp_mat::rows() only takes a contiguous range, so densify to pick out the
    // interior rows; these operators are small enough for that to be free.
    mat dif(diffusion);
    mat adv(advection);

    if (norm(sum(dif.rows(interior), 0), "inf") > tol)
      fail();
    if (norm(sum(adv.rows(interior), 0), "inf") > tol)
      fail();
  }

  cout << "\033[1;32mTest PASSED!\033[0m\n";

  return 0;
}
