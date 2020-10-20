clc
close all

addpath('../mole_MATLAB')

% Parameters
k = 2;
m = 50; % Number of nodes along x-axis
n = 50; % Number of nodes along y-axis

[X, Y] = genCurvGrid(n, m);
% [X, Y] = meshgrid(1:m, 1:n);

mesh(X, Y, zeros(n, m), 'Marker', '.', 'MarkerSize', 10)
view([0 90])
axis equal
set(gcf, 'Color', 'w')

C = X.^2+Y.^2; % Given scalar field (on a nodal grid)
C_ = reshape(C.', [], 1);

[J, Xe, Xn, Ye, Yn] = jacobian2D(k, X, Y);
mn = numel(J);

N = nodal2D(k, m, 1, n, 1);
Ne = N(1:m*n, :);
Nn = N(m*n+1:end, :);

J = spdiags(1./J, 0, mn, mn);
Nx = J*(spdiags(Yn, 0, mn, mn)*Ne-spdiags(Ye, 0, mn, mn)*Nn);
Ny = J*(spdiags(-Xn, 0, mn, mn)*Ne+spdiags(Xe, 0, mn, mn)*Nn);

N = [Nx; Ny];

TMP = N*C_;
Nx = TMP(1:mn);
Ny = TMP(mn+1:end);

Nx = reshape(Nx, m, n)';
Ny = reshape(Ny, m, n)';

figure
set(gcf, 'Color', 'w')
subplot(3, 1, 1)
surf(X, Y, C, 'EdgeColor', 'none');
colorbar
xlabel('x')
ylabel('y')
title('C')
axis equal
view([0 90])
set(gcf, 'Color', 'w')
subplot(3, 1, 2)
surf(X, Y, Nx, 'EdgeColor', 'none');
colorbar
xlabel('x')
ylabel('y')
title('U')
axis equal
view([0 90])
subplot(3, 1, 3)
surf(X, Y, Ny, 'EdgeColor', 'none');
colorbar
xlabel('x')
ylabel('y')
title('V')
axis equal
view([0 90])