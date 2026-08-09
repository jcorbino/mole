/**
 * Tests for the Interpol operators
 *
 *   1. Correct dimensions, both directions
 *   2. Constants are reproduced exactly (the weights sum to one)
 *   3. Second-order accurate at c = 0.5, but only first-order for any other
 *      coefficient -- upwinding trades that order for its numerical diffusion
 *   4. The periodic interpolator wraps: its two boundary-node rows agree, and
 *      it never reads the boundary unknowns
 */

#include "mole.h"
#include <iostream>

static void fail() {
  cout << "\033[1;31mTest FAILED!\033[0m\n";
  exit(1);
}

// Interpolate sin(x) from cell centres to faces on [0,2pi]
static Real interp_err(int m, Real c) {
  Real L = 2 * datum::pi;
  Real dx = L / m;
  Interpol I(m, c, true);
  vec f(m + 2, fill::zeros), ex(m + 1);
  for (int i = 1; i <= m; ++i)
    f(i) = sin((i - 0.5) * dx);
  for (int i = 0; i <= m; ++i)
    ex(i) = sin(i * dx);
  return norm((sp_mat)I * f - ex, "inf");
}

int main() {
  Real tol = 1e-10;

  // ---- dimensions ----
  {
    int m = 10, n = 12, o = 14;
    Interpol I1(m, 0.5);
    Interpol I2(m, n, 0.5, 0.5);
    Interpol I3(m, n, o, 0.5, 0.5, 0.5);

    if (I1.n_rows != (u32)(m + 1) || I1.n_cols != (u32)(m + 2))
      fail();
    if (I2.n_rows != (u32)((m + 1) * n + m * (n + 1)) ||
        I2.n_cols != (u32)((m + 2) * (n + 2)))
      fail();
    if (I3.n_rows != (u32)((m + 1) * n * o + m * (n + 1) * o + m * n * (o + 1)) ||
        I3.n_cols != (u32)((m + 2) * (n + 2) * (o + 2)))
      fail();

    // Faces -> centres, the second type
    Interpol J1(true, m, 0.5);
    if (J1.n_rows != (u32)(m + 2) || J1.n_cols != (u32)(m + 1))
      fail();
  }

  // ---- constants are reproduced ----
  {
    int m = 20;
    Interpol I(m, 0.5);
    vec ones_field(m + 2, fill::ones);
    if (norm((sp_mat)I * ones_field - vec(m + 1, fill::ones), "inf") > tol)
      fail();

    Interpol Ip(m, 0.5, true);
    if (norm((sp_mat)Ip * ones_field - vec(m + 1, fill::ones), "inf") > tol)
      fail();
  }

  // ---- second order at c = 0.5 ----
  {
    Real e1 = interp_err(40, 0.5);
    Real e2 = interp_err(80, 0.5);
    if (std::log2(e1 / e2) < 1.8)
      fail();
  }

  // ---- first order otherwise: the rate must be ~1, not ~2 ----
  for (Real c : {0.0, 1.0}) {
    Real e1 = interp_err(40, c);
    Real e2 = interp_err(80, c);
    Real rate = std::log2(e1 / e2);
    if (rate < 0.9 || rate > 1.5)
      fail();
  }

  // ---- the periodic interpolator wraps around the seam ----
  {
    int m = 15;
    Interpol I(m, 0.5, true);
    sp_mat Is = I;

    // Nodes 0 and m are the same physical point
    if (norm(vec(mat(Is.row(0) - Is.row(m)).t())) > tol)
      fail();

    // and no boundary unknowns are read
    if (norm(vec(mat(Is.col(0)))) > tol || norm(vec(mat(Is.col(m + 1)))) > tol)
      fail();
  }

  cout << "\033[1;32mTest PASSED!\033[0m\n";

  return 0;
}
