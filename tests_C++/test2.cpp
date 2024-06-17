/**
 * Nullity test of Gradient operator
 */

#include "mole.h"
#include <iostream>

void run_nullity_test(int k, Real tol) {
  int m = 2 * k + 1;
  Real dx = 1;

  Gradient G(k, m, dx);
  vec field(m + 2, fill::ones);

  vec sol = G * field;

  if (norm(sol) > tol) {
    cout << "\033[1;31mTest FAILED!\033[0m\n";
    exit(1);
  }
}

int main() {
  Real tol = 1e-12;

  for (int k : {2, 4, 6})
    run_nullity_test(k, tol);

  cout << "\033[1;32mTest PASSED!\033[0m\n";

  return 0;
}
