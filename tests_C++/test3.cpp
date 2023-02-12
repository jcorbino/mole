/**
 * Nullity test of Laplacian operator
 */

#include <iostream>
#include "mole.h"

int main() {
    int k = 2;
    int m = 2*k+1;
    Real dx = 1;
    Real tol = 1e-12;
    
    Laplacian L(k, m, dx);
    vec field(m+2, fill::ones);
    
    vec sol = L*field;
    
    norm(sol) < tol ? cout << "\033[1;32mTest PASSED!\033[0m\n" : cout << "\033[1;31mTest FAILED!\033[0m\n";
    
    return 0;
}
