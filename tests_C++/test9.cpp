/**
 * Order-of-accuracy test for Divergence and Gradient
 *
 * Each operator is applied to a smooth field sampled on the staggered grid it
 * expects, and compared against the analytic derivative. Refining the grid by
 * two must drop the error by about 2^k.
 *
 * The default operators are checked on their interior rows: the boundary rows
 * of a 1D Divergence are empty by construction, since a BC operator is meant to
 * fill them. The periodic operators are checked everywhere, having no such rows.
 *
 * k = 8 is exercised through the periodic Gradient only. The default k = 8
 * closures in gradient.cpp use rationalised decimals where the MATLAB version
 * uses exact fractions, which floors the error near 1e-10 and stops it
 * converging. That is a known, pre-existing discrepancy, not a regression.
 */

#include "mole.h"
#include <iostream>

static void fail() {
  cout << "\033[1;31mTest FAILED!\033[0m\n";
  exit(1);
}

// d/dx exp(x) on [0,1], interior cell rows only
static Real div_default_err(int k, int m) {
  Real dx = 1.0 / m;
  Divergence D(k, m, dx);
  vec f(m + 1), ex(m + 2, fill::zeros);
  for (int i = 0; i <= m; ++i)
    f(i) = exp(i * dx);
  for (int i = 1; i <= m; ++i)
    ex(i) = exp((i - 0.5) * dx);
  vec r = (sp_mat)D * f;
  return norm(r.subvec(1, m) - ex.subvec(1, m), "inf");
}

static Real grad_default_err(int k, int m) {
  Real dx = 1.0 / m;
  Gradient G(k, m, dx);
  vec f(m + 2), ex(m + 1);
  f(0) = exp(0.0);
  f(m + 1) = exp(1.0);
  for (int i = 1; i <= m; ++i)
    f(i) = exp((i - 0.5) * dx);
  for (int i = 0; i <= m; ++i)
    ex(i) = exp(i * dx);
  return norm((sp_mat)G * f - ex, "inf");
}

// d/dx sin(x) on [0,2pi], every row
static Real div_periodic_err(int k, int m) {
  Real L = 2 * datum::pi;
  Real dx = L / m;
  Divergence D(k, m, dx, true);
  vec f(m + 1), ex(m + 2);
  for (int i = 0; i <= m; ++i)
    f(i) = sin(i * dx);
  ex(0) = cos(0.0);
  ex(m + 1) = cos(L);
  for (int i = 1; i <= m; ++i)
    ex(i) = cos((i - 0.5) * dx);
  return norm((sp_mat)D * f - ex, "inf");
}

static Real grad_periodic_err(int k, int m) {
  Real L = 2 * datum::pi;
  Real dx = L / m;
  Gradient G(k, m, dx, true);
  vec f(m + 2, fill::zeros), ex(m + 1);
  for (int i = 1; i <= m; ++i)
    f(i) = sin((i - 0.5) * dx);
  for (int i = 0; i <= m; ++i)
    ex(i) = cos(i * dx);
  return norm((sp_mat)G * f - ex, "inf");
}

// Observed rate must reach 90% of the nominal order
static void check(Real coarse, Real fine, int k) {
  Real rate = std::log2(coarse / fine);
  if (rate < 0.9 * k)
    fail();
}

int main() {
  // Grid pairs chosen so the fine error stays well above roundoff
  for (int k : {2, 4, 6}) {
    int m = (k == 2 ? 40 : k == 4 ? 20 : 16);
    check(div_default_err(k, m), div_default_err(k, 2 * m), k);
    check(grad_default_err(k, m), grad_default_err(k, 2 * m), k);
  }

  for (int k : {2, 4, 6}) {
    int m = (k <= 4 ? 24 : 32);
    check(div_periodic_err(k, m), div_periodic_err(k, 2 * m), k);
    check(grad_periodic_err(k, m), grad_periodic_err(k, 2 * m), k);
  }

  // Gradient reaches k = 8; Divergence asserts k < 7, so it is not exercised here
  check(grad_periodic_err(8, 32), grad_periodic_err(8, 64), 8);

  cout << "\033[1;32mTest PASSED!\033[0m\n";

  return 0;
}
