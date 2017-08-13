% Compact mimetic operators
% Factorizes higher order mimetic divergence or gradient operators
% STAR (D2*STAR*G2) operator has dimensions: m+1 by m+1
% Where STAR = RKD*LKG

addpath('../mole_MATLAB')

k = 4;
m = 2*k+1;
dx = 1;

% Factorize Divergence ----------------------------------------------------
D2 = div(2, m, dx);

DK = div(k, m, dx);

RKD = full(D2(2:end-1, :)\DK(2:end-1, :));

LKD = full(DK(2:end-1, :)/D2(2:end-1, :));

% Factorize Gradient ------------------------------------------------------
G2 = grad(2, m, dx);

GK = grad(k, m, dx);

RKG = full(G2\GK);

LKG = full(GK/G2);

STAR = RKD*LKG
