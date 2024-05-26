% Nullity test of Divergence operator

addpath('../mole_MATLAB')

k = 2;
m = 2 * k + 1;
dx = 1;
tol = 1e-12;

D = div(k, m, dx);

field = ones(m + 1, 1);

sol = D * field;

if (norm(sol) < tol)
    fprintf("Test PASSED!\n");
else
    fprintf("Test FAILED!");
end