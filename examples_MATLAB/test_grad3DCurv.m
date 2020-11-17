clc
close all

addpath('../mole_MATLAB')

% Parameters
k = 4;  % Order of accuracy
m = 20; % Number of nodes along x-axis
n = 20; % Number of nodes along y-axis
o = 20; % Number of nodes along z-axis

[X, Y, Z] = meshgrid(1:m, 1:n, 1:o);
X = X.*X;
Y = Y.*Y;
Z = Z.*Z;

C = X.^2+Y.^2+Z.^2; % Given scalar field (on a nodal grid)

% Plot the physical grid
scatter3(X(:), Y(:), Z(:), 50, C(:), 'Filled');
title('Given scalar field')
xlabel('x')
ylabel('y')
zlabel('z')

% Staggered logical grid
[Xs, Ys, Zs] = meshgrid([1 1.5 : 1 : m-0.5 m], [1 1.5 : 1 : n-0.5 n], [1 1.5 : 1 : o-0.5 o]);
% Interpolate scalar field from nodal logical to staggered logical
Cs = interp3(C, Xs, Ys, Zs);
% Reshape the field so it can be multiplied by the operator later on
C_ = reshape(permute(Cs, [2, 1, 3]), [], 1);

% Get 3D curvilinear mimetic gradient
G = grad3DCurv(k, X, Y, Z);

% Apply the operator to the field
TMP = G*C_;
Gx = TMP(1:m*(n-1)*(o-1));
Gy = TMP(m*(n-1)*(o-1)+1:m*(n-1)*(o-1)+(m-1)*n*(o-1));
Gz = TMP(m*(n-1)*(o-1)+(m-1)*n*(o-1)+1:end);

% Reshape for visualization
Gx = permute(reshape(Gx, m, n-1, o-1), [2, 1, 3]);
Gy = permute(reshape(Gy, m-1, n, o-1), [2, 1, 3]);
Gz = permute(reshape(Gz, m-1, n-1, o), [2, 1, 3]);

xslice = 10;
yslice = 10;
zslice = 10;

figure
slice(Gx, xslice, yslice, zslice)
title('Gx')
xlabel('x')
ylabel('y')
zlabel('z')

figure
slice(Gy, xslice, yslice, zslice)
title('Gy')
xlabel('x')
ylabel('y')
zlabel('z')

figure
slice(Gz, xslice, yslice, zslice)
title('Gz')
xlabel('x')
ylabel('y')
zlabel('z')