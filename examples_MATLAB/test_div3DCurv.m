% Tests the 3D curvilinear divergence
clc
close all

addpath('../mole_MATLAB')

% Parameters
k = 2;  % Order of accuracy
m = 20; % Number of nodes along x-axis
n = 20; % Number of nodes along y-axis
o = 20; % Number of nodes along z-axis

[X, Y] = genCurvGrid(n, m);
X = repmat(X, [1 1 o]);
Y = repmat(Y, [1 1 o]);
[~, ~, Z] = meshgrid(1:m, 1:n, 1:o);
% [X, Y, Z] = meshgrid(1:m, 1:n, 1:o);

% Plot the physical grid
scatter3(X(:), Y(:), Z(:), 50, 'Filled');
title('Given scalar field')
xlabel('x')
ylabel('y')
zlabel('z')
axis equal

Ux = (X(1:end-1, :, :) + X(2:end, :, :))/2;
Ux = (Ux(:, :, 1:end-1) + Ux(:, :, 2:end))/2;
Uy = (Y(1:end-1, :, :) + Y(2:end, :, :))/2;
Uy = (Uy(:, :, 1:end-1) + Uy(:, :, 2:end))/2;
Uz = (Z(1:end-1, :, :) + Z(2:end, :, :))/2;
Uz = (Uz(:, :, 1:end-1) + Uz(:, :, 2:end))/2;

Vx = (X(:, 1:end-1, :) + X(:, 2:end, :))/2;
Vx = (Vx(:, :, 1:end-1) + Vx(:, :, 2:end))/2;
Vy = (Y(:, 1:end-1, :) + Y(:, 2:end, :))/2;
Vy = (Vy(:, :, 1:end-1) + Vy(:, :, 2:end))/2;
Vz = (Z(:, 1:end-1, :) + Z(:, 2:end, :))/2;
Vz = (Vz(:, :, 1:end-1) + Vz(:, :, 2:end))/2;

Wx = (X(1:end-1, :, :) + X(2:end, :, :))/2;
Wx = (Wx(:, 1:end-1, :) + Wx(:, 2:end, :))/2;
Wy = (Y(1:end-1, :, :) + Y(2:end, :, :))/2;
Wy = (Wy(:, 1:end-1, :) + Wy(:, 2:end, :))/2;
Wz = (Z(1:end-1, :, :) + Z(2:end, :, :))/2;
Wz = (Wz(:, 1:end-1, :) + Wz(:, 2:end, :))/2;

% Interpolate U values
Ugiven = sin(X);
interpolant = scatteredInterpolant([X(:) Y(:) Z(:)], Ugiven(:));
U = interpolant(Ux, Uy, Uz);
% Interpolate V values
Vgiven = cos(Y);
interpolant = scatteredInterpolant([X(:) Y(:) Z(:)], Vgiven(:));
V = interpolant(Vx, Vy, Vz);
% Interpolate W values
Wgiven = cos(Y);
interpolant = scatteredInterpolant([X(:) Y(:) Z(:)], Wgiven(:));
W = interpolant(Wx, Wy, Wz);

U = reshape(permute(U, [2, 1, 3]), [], 1);
V = reshape(permute(V, [2, 1, 3]), [], 1);
W = reshape(permute(W, [2, 1, 3]), [], 1);

% Get 3D curvilinear mimetic divergence
D = div3DCurv(k, X, Y, Z);

spy(D-div3D(k, m-1, 1, n-1, 1, o-1, 1))
D*[U; V; W];