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
 * A convergence rate is a blunt instrument for the boundary closures: on a field
 * whose truncation error is large, a closure coefficient wrong in the 11th digit
 * changes nothing measurable. So the closures are checked separately, by
 * polynomial exactness -- a k-th order operator must differentiate every
 * polynomial up to degree k exactly, with no truncation error to hide behind.
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

// Polynomial exactness. Sampling x^p on the grid an operator expects, its output
// must be the exact derivative at every point it writes -- closures included.
static Real grad_poly_err(int k, int m) {
  Real dx = 1.0 / m;
  Gradient G(k, m, dx);
  vec xc(m + 2), xf(m + 1);
  xc(0) = 0;
  for (int i = 1; i <= m; ++i)
    xc(i) = (i - 0.5) * dx;
  xc(m + 1) = 1;
  for (int i = 0; i <= m; ++i)
    xf(i) = i * dx;

  Real worst = 0;
  for (int p = 0; p <= k; ++p) {
    vec f(m + 2), ex(m + 1);
    for (int i = 0; i < m + 2; ++i)
      f(i) = pow(xc(i), p);
    for (int i = 0; i <= m; ++i)
      ex(i) = (p == 0) ? 0.0 : p * pow(xf(i), p - 1);
    worst = std::max(worst, norm((sp_mat)G * f - ex, "inf"));
  }
  return worst;
}

static Real div_poly_err(int k, int m) {
  Real dx = 1.0 / m;
  Divergence D(k, m, dx);
  vec xc(m + 2), xf(m + 1);
  for (int i = 1; i <= m; ++i)
    xc(i) = (i - 0.5) * dx;
  for (int i = 0; i <= m; ++i)
    xf(i) = i * dx;

  Real worst = 0;
  for (int p = 0; p <= k; ++p) {
    vec f(m + 1), ex(m + 2, fill::zeros);
    for (int i = 0; i <= m; ++i)
      f(i) = pow(xf(i), p);
    for (int i = 1; i <= m; ++i)
      ex(i) = (p == 0) ? 0.0 : p * pow(xc(i), p - 1);
    // The boundary rows of a default Divergence are empty by construction
    vec r = (sp_mat)D * f;
    worst = std::max(worst, norm(r.subvec(1, m) - ex.subvec(1, m), "inf"));
  }
  return worst;
}

// Curl and RobinBC build Gradients internally, so they inherit whatever the
// closures are. Checking them directly pins that down: neither can drift from
// the Gradient, but nothing else would notice if the Gradient itself changed.

// RobinBC with (a, b) = (0, 1) is pure Neumann, so its two boundary rows must
// return the OUTWARD normal derivative: -f'(west) and +f'(east).
static Real robin_poly_err(int k, int m) {
  Real dx = 1.0 / m;
  RobinBC BC(k, m, dx, 0, 1);
  vec xc(m + 2);
  xc(0) = 0;
  for (int i = 1; i <= m; ++i)
    xc(i) = (i - 0.5) * dx;
  xc(m + 1) = 1;

  Real worst = 0;
  for (int p = 0; p <= k; ++p) {
    vec f(m + 2);
    for (int i = 0; i < m + 2; ++i)
      f(i) = pow(xc(i), p);
    vec r = (sp_mat)BC * f;
    Real dwest = (p == 0) ? 0.0 : p * pow(0.0, p - 1);
    Real deast = (p == 0) ? 0.0 : p * pow(1.0, p - 1);
    worst = std::max(worst, std::fabs(r(0) - (-dwest)));
    worst = std::max(worst, std::fabs(r(m + 1) - deast));
  }
  return worst;
}

// curl of (Fx, Fy) = (y^p, x^p) is p*x^(p-1) - p*y^(p-1)
static Real curl_poly_err(int k, int m) {
  int n = m;
  Real dx = 1.0 / m, dy = 1.0 / n;
  Curl C(k, m, n, dx, dy);

  vec xn(m + 1), yn(n + 1), xcb(m + 2), ycb(n + 2);
  for (int i = 0; i <= m; ++i)
    xn(i) = i * dx;
  for (int j = 0; j <= n; ++j)
    yn(j) = j * dy;
  xcb(0) = 0;
  for (int i = 1; i <= m; ++i)
    xcb(i) = (i - 0.5) * dx;
  xcb(m + 1) = 1;
  ycb(0) = 0;
  for (int j = 1; j <= n; ++j)
    ycb(j) = (j - 0.5) * dy;
  ycb(n + 1) = 1;

  Real worst = 0;
  for (int p = 0; p <= k; ++p) {
    vec Fx((m + 1) * (n + 2)), Fy((m + 2) * (n + 1));
    for (int j = 0; j < n + 2; ++j)
      for (int i = 0; i < m + 1; ++i)
        Fx(j * (m + 1) + i) = pow(ycb(j), p);
    for (int j = 0; j < n + 1; ++j)
      for (int i = 0; i < m + 2; ++i)
        Fy(j * (m + 2) + i) = pow(xcb(i), p);

    vec r = (sp_mat)C * join_vert(Fx, Fy);
    for (int j = 0; j <= n; ++j)
      for (int i = 0; i <= m; ++i) {
        Real ex =
            (p == 0) ? 0.0 : p * pow(xn(i), p - 1) - p * pow(yn(j), p - 1);
        worst = std::max(worst, std::fabs(r(j * (m + 1) + i) - ex));
      }
  }
  return worst;
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

  // Polynomial exactness of the boundary closures, at every supported order.
  // This is what catches a closure coefficient that is merely close.
  for (int k : {2, 4, 6, 8})
    if (grad_poly_err(k, 2 * k + 4) > 1e-12)
      fail();
  for (int k : {2, 4, 6})
    if (div_poly_err(k, 2 * k + 4) > 1e-12)
      fail();
  for (int k : {2, 4, 6, 8}) {
    if (robin_poly_err(k, 2 * k + 4) > 1e-12)
      fail();
    if (curl_poly_err(k, 2 * k + 4) > 1e-12)
      fail();
  }

  cout << "\033[1;32mTest PASSED!\033[0m\n";

  return 0;
}
