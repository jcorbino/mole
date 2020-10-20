% Tests the 2D curvilinear gradient
clc
close all

addpath('../mole_MATLAB')

% Parameters
k = 2;
m = 10;
n = 10;

[X, Y] = genCurvGrid(n, m);
%[X, Y] = meshgrid(1:m, 1:n);

mesh(X, Y, zeros(n, m), 'Marker', '.', 'MarkerSize', 10)
view([0 90])
axis equal
set(gcf, 'Color', 'w')
hold on

[n, m] = size(X);
n = n-1;
m = m-1;

Ux = (X(1:end-1, :) + X(2:end, :))/2;
Uy = (Y(1:end-1, :) + Y(2:end, :))/2;
scatter3(Ux(:), Uy(:), zeros(n*(m+1), 1), '+', 'MarkerEdgeColor', 'k')

Vx = (X(:, 1:end-1) + X(:, 2:end))/2;
Vy = (Y(:, 1:end-1) + Y(:, 2:end))/2;
scatter3(Vx(:), Vy(:), zeros((n+1)*m, 1), '*', 'MarkerEdgeColor', 'k')

Cx = (Vx(1:end-1, :) + Vx(2:end, :))/2;
Cy = (Uy(:, 1:end-1) + Uy(:, 2:end))/2;
scatter3(Cx(:), Cy(:), zeros(n*m, 1), '.', 'MarkerEdgeColor', 'r')

% West-East sides
Cx = [Ux(:, 1) Cx];
Cy = [Uy(:, 1) Cy];
Cx = [Cx Ux(:, end)];
Cy = [Cy Uy(:, end)];

% South-North sides
Cx = [[0 Vx(1, :) 0]; Cx];
Cy = [[0 Vy(1, :) 0]; Cy];
Cx = [Cx; [0 Vx(end, :) 0]];
Cy = [Cy; [0 Vy(end, :) 0]];

% Corners
Cx(1, 1) = X(1, 1);
Cy(1, 1) = Y(1, 1);
Cx(1, end) = X(1, end);
Cy(1, end) = Y(1, end);
Cx(end, 1) = X(end, 1);
Cy(end, 1) = Y(end, 1);
Cx(end, end) = X(end, end);
Cy(end, end) = Y(end, end);

scatter3(Cx(:), Cy(:), zeros((m+2)*(n+2), 1), 'o', 'MarkerEdgeColor', 'r')
legend('Nodal points', 'u', 'v', 'Centers', 'All centers')
hold off

tic
G = grad2DCurv(k, X, Y);
toc

Cgiven = Cx.^2 + Cy.^2;

UV = G*reshape(Cgiven', [], 1);

U = UV(1:n*(m+1));
U = reshape(U, m+1, n)';
V = UV(n*(m+1)+1:end);
V = reshape(V, m, n+1)';

figure
subplot(5, 1, 4)
surf(Vx, Vy, V, 'EdgeColor', 'none');
colorbar
view([0 90])
xlabel('x')
ylabel('y')
title('Numerical V')
axis equal
subplot(5, 1, 2)
surf(Ux, Uy, U, 'EdgeColor', 'none');
colorbar
view([0 90])
xlabel('x')
ylabel('y')
title('Numerical U')
axis equal
subplot(5, 1, 3)
surf(Vx, Vy, 2*Vy, 'EdgeColor', 'none');
colorbar
view([0 90])
xlabel('x')
ylabel('y')
title('Analytical V')
axis equal
set(gcf, 'Color', 'w')
subplot(5, 1, 1)
surf(Ux, Uy, 2*Ux, 'EdgeColor', 'none');
colorbar
view([0 90])
xlabel('x')
ylabel('y')
title('Analytical U')
axis equal
set(gcf, 'Color', 'w')
subplot(5, 1, 5)
surf(Cx, Cy, Cgiven, 'EdgeColor', 'none');
colorbar
view([0 90])
xlabel('x')
ylabel('y')
title('Scalar field')
axis equal
set(gcf, 'Color', 'w')

max(max(abs(2*Ux - U)))
max(max(abs(2*Vy - V)))
