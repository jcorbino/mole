% https://en.wikipedia.org/wiki/Transfinite_interpolation
clc 
close all

addpath('grids/Horseshoe')

% Grid resolution
m = 20;
n = 20;

% Logical grid
Xl = linspace(0, 1, m);
Yl = linspace(0, 1, n);

% Allocate space for physical grid
X = zeros(m, n);
Y = zeros(m, n);

for i = 1 : m
    u = Xl(i);
    for j = 1 : n
        v = Yl(j);
        % Transfinite interpolation
        XY = (1-v)*bottom(u)+v*top(u)+(1-u)*left(v)+u*right(v)-...
            (u*v*top(1)+u*(1-v)*bottom(1)+v*(1-u)*top(0)+(1-u)*(1-v)*bottom(0));
        X(i, j) = XY(1);
        Y(i, j) = XY(2);
    end
end

mesh(X, Y, zeros(n, m), 'Marker', '.', 'MarkerSize', 10, 'EdgeColor', 'b')
title(['Physical grid. m = ' num2str(m) ', n = ' num2str(n)])
set(gcf, 'color', 'w')
axis equal
axis off
view([0 90])