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

[J, Xe, Xn, Ye, Yn] = jacobian2D(k, X, Y);

G = nodal2D(k, m, 1, n, 1);
Ge = G(1:m*n, :);
Gn = G(m*n+1:end, :);

C = X.^2+Y.^2;

C_ = reshape(C.', [], 1);

Ce = Ge*C_;
Cn = Gn*C_;

Gx = 1./J.*(Yn.*Ce-Ye.*Cn);
Gy = 1./J.*(-Xn.*Ce+Xe.*Cn);

Gx = reshape(Gx, m, n)';
Gy = reshape(Gy, m, n)';

figure
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
surf(X, Y, Gx, 'EdgeColor', 'none');
colorbar
xlabel('x')
ylabel('y')
title('U')
axis equal
view([0 90])
subplot(3, 1, 3)
surf(X, Y, Gy, 'EdgeColor', 'none');
colorbar
xlabel('x')
ylabel('y')
title('V')
axis equal
view([0 90])