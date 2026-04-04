/**
 * Tests for the 2D mimetic Curl operator:
 *   1. Correct matrix dimensions
 *   2. curl(-y, x) = 2 (linear field, should be exact)
 *   3. curl(grad(phi)) = 0 (fundamental identity)
 */

#include "mole.h"
#include <iostream>

int main() {
  Real tol = 1e-10;
  int k = 2;
  int m = 20;
  int n = 20;
  Real dx = 1.0 / m;
  Real dy = 1.0 / n;

  // Build curl and 1D gradient operators
  Curl C(k, m, n, dx, dy);
  Gradient Gx(k, m, dx);
  Gradient Gy(k, n, dy);

  // Test 1: Dimensions
  if (C.n_rows != (u32)((m + 1) * (n + 1)) ||
      C.n_cols != (u32)((m + 1) * (n + 2) + (m + 2) * (n + 1))) {
    cout << "\033[1;31mTest FAILED!\033[0m\n";
    exit(1);
  }

  // Test 2: curl(-y, x) = 2
  vec Fx((m + 1) * (n + 2));
  for (int j = 0; j < (int)(n + 2); j++) {
    Real y;
    if (j == 0)
      y = 0.0;
    else if (j == (int)(n + 1))
      y = n * dy;
    else
      y = (j - 0.5) * dy;
    for (int i = 0; i < (int)(m + 1); i++)
      Fx(j * (m + 1) + i) = -y;
  }

  vec Fy((m + 2) * (n + 1));
  for (int j = 0; j < (int)(n + 1); j++) {
    for (int i = 0; i < (int)(m + 2); i++) {
      Real x;
      if (i == 0)
        x = 0.0;
      else if (i == (int)(m + 1))
        x = m * dx;
      else
        x = (i - 0.5) * dx;
      Fy(j * (m + 2) + i) = x;
    }
  }

  vec curlF = C * join_vert(Fx, Fy);
  if (norm(curlF - 2.0 * ones(curlF.n_elem), "inf") > tol) {
    cout << "\033[1;31mTest FAILED!\033[0m\n";
    exit(1);
  }

  // Test 3: curl(grad(phi)) = 0
  sp_mat G_full =
      Utils::spjoin_cols(Utils::spkron(speye(n + 2, n + 2), (sp_mat)Gx),
                         Utils::spkron((sp_mat)Gy, speye(m + 2, m + 2)));

  arma_rng::set_seed(42);
  vec phi = randu((m + 2) * (n + 2));

  if (norm(C * (G_full * phi), "inf") > tol) {
    cout << "\033[1;31mTest FAILED!\033[0m\n";
    exit(1);
  }

  cout << "\033[1;32mTest PASSED!\033[0m\n";

  return 0;
}
