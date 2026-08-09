/**
 * Structure and nullity of the PERIODIC operators (1D, 2D and 3D)
 *
 *   1. D, G and L annihilate a constant field, as in the default operators
 *   2. The two boundary-node rows of a periodic Divergence are identical --
 *      they are the same physical point
 *   3. A periodic Gradient leaves the boundary columns empty: a periodic
 *      problem carries no independent boundary unknowns
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
    int m = 2 * k + 1;
    Real dx = 0.7;

    Divergence D(k, m, dx, true);
    Gradient G(k, m, dx, true);
    Laplacian L(k, m, dx, true);

    if (norm((sp_mat)D * vec(m + 1, fill::ones)) > tol)
      fail();
    if (norm((sp_mat)G * vec(m + 2, fill::ones)) > tol)
      fail();
    if (norm((sp_mat)L * vec(m + 2, fill::ones)) > tol)
      fail();

    // The west and east boundary nodes are the same point, so the rows agree
    sp_mat Ds = D;
    if (norm(vec(mat(Ds.row(0) - Ds.row(m + 1)).t())) > tol)
      fail();

    // No independent boundary unknowns, so those columns are empty
    sp_mat Gs = G;
    if (norm(vec(mat(Gs.col(0)))) > tol || norm(vec(mat(Gs.col(m + 1)))) > tol)
      fail();
  }

  // ---- 2D ----
  for (int k : {2, 4, 6}) {
    int m = 2 * k + 1;
    int n = 2 * k + 3;
    Real dx = 0.7, dy = 1.3;

    Divergence D(k, m, n, dx, dy, true);
    Gradient G(k, m, n, dx, dy, true);
    Laplacian L(k, m, n, dx, dy, true);

    int nfaces = (m + 1) * n + m * (n + 1);
    int ncells = (m + 2) * (n + 2);

    if (D.n_rows != (u32)ncells || D.n_cols != (u32)nfaces)
      fail();
    if (G.n_rows != (u32)nfaces || G.n_cols != (u32)ncells)
      fail();

    if (norm((sp_mat)D * vec(nfaces, fill::ones)) > tol)
      fail();
    if (norm((sp_mat)G * vec(ncells, fill::ones)) > tol)
      fail();
    if (norm((sp_mat)L * vec(ncells, fill::ones)) > tol)
      fail();
  }

  // ---- 3D ----
  for (int k : {2, 4}) {
    int m = 2 * k + 1;
    int n = 2 * k + 2;
    int o = 2 * k + 3;
    Real dx = 0.7, dy = 1.3, dz = 0.9;

    Divergence D(k, m, n, o, dx, dy, dz, true);
    Gradient G(k, m, n, o, dx, dy, dz, true);

    int nfaces = (m + 1) * n * o + m * (n + 1) * o + m * n * (o + 1);
    int ncells = (m + 2) * (n + 2) * (o + 2);

    if (D.n_rows != (u32)ncells || D.n_cols != (u32)nfaces)
      fail();
    if (G.n_rows != (u32)nfaces || G.n_cols != (u32)ncells)
      fail();

    if (norm((sp_mat)D * vec(nfaces, fill::ones)) > tol)
      fail();
    if (norm((sp_mat)G * vec(ncells, fill::ones)) > tol)
      fail();
  }

  cout << "\033[1;32mTest PASSED!\033[0m\n";

  return 0;
}
