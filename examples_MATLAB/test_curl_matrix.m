% test_curl_matrix.m
%
% Two test cases are plotted side by side (computed, analytical, error):
%   1. F = (-y, x) on [-10,10]^2 — analytical curl = 2 everywhere.
%   2. F = (sin(pi*y), cos(pi*x)) on [0,1]^2 — analytical curl = -pi*sin(pi*x) - pi*cos(pi*y).

clc
close all

addpath('../mole_MATLAB')

order =   2;
west  = -10;
east  =  10;
south = -10;
north =  10;
m = 20;
n = 20;
dx = (east - west) / m;
dy = (north - south) / n;

% --- Staggered grid positions (ndgrid ordering: x varies fastest) ---

% Node positions (cell boundaries)
xnodes = linspace(west, east, m + 1);
ynodes = linspace(south, north, n + 1);

% Cell-center + boundary positions
xcb = [west, (west + dx/2) : dx : (east - dx/2), east];
ycb = [south, (south + dy/2) : dy : (north - dy/2), north];

% --- Sample vector field on face-staggered positions ---

% Fx at x-faces: x = xnodes (m+1), y = ycb (n+2)
[Xfx, Yfx] = ndgrid(xnodes, ycb);
Fx = P(Xfx, Yfx);

% Fy at y-faces: x = xcb (m+2), y = ynodes (n+1)
[Xfy, Yfy] = ndgrid(xcb, ynodes);
Fy = Q(Xfy, Yfy);

% Assemble face vector F = [Fx; Fy]
F = [Fx(:); Fy(:)];

% --- Compute mimetic curl ---
C = curl2DMatrix(order, m, dx, n, dy);
curlF = C * F;
curlF = reshape(curlF, m + 1, n + 1);

% --- Analytical solution at nodes ---
[Xn, Yn] = ndgrid(xnodes, ynodes);
exact = 2 * ones(m + 1, n + 1);

% --- Plot computed vs expected ---
figure('Position', [100 100 1200 500], 'Color', 'w');

% Computed solution
subplot(1, 3, 1);
pcolor(Xn', Yn', curlF');
shading interp;
colorbar;
title('Computed curl (curl2D)');
xlabel('x'); ylabel('y');
clim([min(curlF(:)) - 0.1, max(curlF(:)) + 0.1]);
axis equal tight;

% Analytical solution
subplot(1, 3, 2);
pcolor(Xn', Yn', exact');
shading interp;
colorbar;
title('Analytical curl (= 2)');
xlabel('x'); ylabel('y');
clim([min(curlF(:)) - 0.1, max(curlF(:)) + 0.1]);
axis equal tight;

% Error
subplot(1, 3, 3);
err = curlF - exact;
pcolor(Xn', Yn', err');
shading interp;
colorbar;
title('Error (computed - analytical)');
xlabel('x'); ylabel('y');
axis equal tight;

sgtitle(sprintf('Mimetic curl2D test: F = (-y, x), order = %d, %dx%d grid', order, m, n));

fprintf('Analytical curl = 2\n');
fprintf('Computed curl:   min = %.10f, max = %.10f\n', min(curlF(:)), max(curlF(:)));
fprintf('Max abs error:   %.2e\n', max(abs(err(:))));

% --- Also test with a non-trivial field on [0,1]x[0,1] ---
figure('Position', [100 650 1200 500], 'Color', 'w');

% Separate grid for the non-trivial test
m2 = 40; n2 = 40;
west2 = 0; east2 = 1; south2 = 0; north2 = 1;
dx2 = (east2 - west2) / m2;
dy2 = (north2 - south2) / n2;

xnodes2 = linspace(west2, east2, m2 + 1);
ynodes2 = linspace(south2, north2, n2 + 1);
xcb2 = [west2, (west2 + dx2/2) : dx2 : (east2 - dx2/2), east2];
ycb2 = [south2, (south2 + dy2/2) : dy2 : (north2 - dy2/2), north2];

[Xfx2, Yfx2] = ndgrid(xnodes2, ycb2);
[Xfy2, Yfy2] = ndgrid(xcb2, ynodes2);

% Fx = sin(pi*y), Fy = cos(pi*x)
% curl = dFy/dx - dFx/dy = -pi*sin(pi*x) - pi*cos(pi*y)
Fx2 = sin(pi * Yfx2(:));
Fy2 = cos(pi * Xfy2(:));
F2 = [Fx2; Fy2];

C2 = curl2DMatrix(order, m2, dx2, n2, dy2);
curlF2 = C2 * F2;
curlF2 = reshape(curlF2, m2 + 1, n2 + 1);

[Xn2, Yn2] = ndgrid(xnodes2, ynodes2);
exact2 = -pi*sin(pi*Xn2) - pi*cos(pi*Yn2);

subplot(1, 3, 1);
pcolor(Xn2', Yn2', curlF2');
shading interp;
colorbar;
title('Computed curl');
xlabel('x'); ylabel('y');
axis equal tight;

subplot(1, 3, 2);
pcolor(Xn2', Yn2', exact2');
shading interp;
colorbar;
title('Analytical curl');
xlabel('x'); ylabel('y');
axis equal tight;

subplot(1, 3, 3);
err2 = curlF2 - exact2;
pcolor(Xn2', Yn2', err2');
shading interp;
colorbar;
title('Error');
xlabel('x'); ylabel('y');
axis equal tight;

sgtitle(sprintf('curl(sin(\\pi y), cos(\\pi x)) on [0,1]^2, order = %d, %dx%d grid', order, m2, n2));

fprintf('\nNon-trivial field: F = (sin(pi*y), cos(pi*x)) on [0,1]^2, %dx%d\n', m2, n2);
fprintf('Max abs error:   %.2e\n', max(abs(err2(:))));

% --- Vector field components ---
function U = P(~, Y)
    U = -Y;
end

function V = Q(X, ~)
    V = X;
end
