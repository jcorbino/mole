% Nullity test of Gradient operator

addpath('../mole_MATLAB')

k = 2;
m = 2 * k + 1;
dx = 1;
tol = 1e-12;

G = grad(k, m, dx);

field = ones(m + 2, 1);

sol = G * field;

if (norm(sol) < tol)
    fprintf("Test PASSED!\n");
else
    fprintf("Test FAILED!");
end