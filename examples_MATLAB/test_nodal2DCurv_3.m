clc
close all

addpath('../mole_MATLAB')

% Parameters
k = 2;  % Order of accuracy
m = 40; % Number of nodes along x-axis
n = 40; % Number of nodes along y-axis

[X, Y] = genCurvGrid(n, m);
% [X, Y] = meshgrid(1:m, 1:n);
[Xl, Yl] = meshgrid(1:m, 1:n); % Logical grids

% Plot the physical grid
mesh(X, Y, zeros(n, m), 'Marker', '.', 'MarkerSize', 10)
view([0 90])
axis equal
set(gcf, 'Color', 'w')

C = X.^2+Y.^2; % Given scalar field (on a nodal grid)
% Staggered logical grid
[Xs, Ys] = meshgrid([1 1.5 : 1 : m-0.5 m], [1 1.5 : 1 : n-0.5 n]);
% Interpolate scalar field from nodal logical to staggered logical
Cs = interp2(Xl, Yl, C, Xs, Ys);
% Reshape the field so it can be multiplied by the operator later on
C_ = reshape(Cs', [], 1);

% Get the determinant of the jacobian and the metrics
[J, Xe, Xn, Ye, Yn] = jacobian2D(k, X, Y);

% Make them surfaces so they can be interpolated
J = reshape(J, m, n)';
Xe = reshape(Xe, m, n)';
Xn = reshape(Xn, m, n)';
Ye = reshape(Ye, m, n)';
Yn = reshape(Yn, m, n)';

% Interpolate the metrics on the logical grid for positions u and v
Ju = interp2(Xl, Yl, J, (Xl(1:end-1, :)+Xl(2:end, :))/2,...
                                        (Yl(1:end-1, :)+Yl(2:end, :))/2);
Jv = interp2(Xl, Yl, J, (Xl(:, 1:end-1)+Xl(:, 2:end))/2,...
                                        (Yl(:, 1:end-1)+Yl(:, 2:end))/2);
Xeu = interp2(Xl, Yl, Xe, (Xl(1:end-1, :)+Xl(2:end, :))/2,...
                                        (Yl(1:end-1, :)+Yl(2:end, :))/2);
Xev = interp2(Xl, Yl, Xe, (Xl(:, 1:end-1)+Xl(:, 2:end))/2,...
                                        (Yl(:, 1:end-1)+Yl(:, 2:end))/2);
Xnu = interp2(Xl, Yl, Xn, (Xl(1:end-1, :)+Xl(2:end, :))/2,...
                                        (Yl(1:end-1, :)+Yl(2:end, :))/2);
Xnv = interp2(Xl, Yl, Xn, (Xl(:, 1:end-1)+Xl(:, 2:end))/2,...
                                        (Yl(:, 1:end-1)+Yl(:, 2:end))/2);
Yeu = interp2(Xl, Yl, Ye, (Xl(1:end-1, :)+Xl(2:end, :))/2,...
                                        (Yl(1:end-1, :)+Yl(2:end, :))/2);
Yev = interp2(Xl, Yl, Ye, (Xl(:, 1:end-1)+Xl(:, 2:end))/2,...
                                        (Yl(:, 1:end-1)+Yl(:, 2:end))/2);
Ynu = interp2(Xl, Yl, Yn, (Xl(1:end-1, :)+Xl(2:end, :))/2,...
                                        (Yl(1:end-1, :)+Yl(2:end, :))/2);
Ynv = interp2(Xl, Yl, Yn, (Xl(:, 1:end-1)+Xl(:, 2:end))/2,...
                                        (Yl(:, 1:end-1)+Yl(:, 2:end))/2);

% Convert metrics to diagonal matrices so they can be multiplied by the logical
% operators
Ju = spdiags(1./reshape(Ju', [], 1), 0, numel(Ju), numel(Ju));
Jv = spdiags(1./reshape(Jv', [], 1), 0, numel(Jv), numel(Jv));
Xeu = spdiags(reshape(Xeu', [], 1), 0, numel(Xeu), numel(Xeu));
Xev = spdiags(reshape(Xev', [], 1), 0, numel(Xev), numel(Xev));
Xnu = spdiags(reshape(Xnu', [], 1), 0, numel(Xnu), numel(Xnu));
Xnv = spdiags(reshape(Xnv', [], 1), 0, numel(Xnv), numel(Xnv));
Yeu = spdiags(reshape(Yeu', [], 1), 0, numel(Yeu), numel(Yeu));
Yev = spdiags(reshape(Yev', [], 1), 0, numel(Yev), numel(Yev));
Ynu = spdiags(reshape(Ynu', [], 1), 0, numel(Ynu), numel(Ynu));
Ynv = spdiags(reshape(Ynv', [], 1), 0, numel(Ynv), numel(Ynv));

% Construct 2D uniform mimetic gradient operator (d/de, d/dn)
G = grad2D(k, m-1, 1, n-1, 1);
Ge = G(1:m*(n-1), :);
Gn = G(m*(n-1)+1:end, :);

% Apply transformation
Gx = Ju*(Ynu*Ge-Yeu*I(Gn, m-1, n-1, 'Gn'));
Gy = Jv*(-Xnv*I(Ge, m-1, n-1, 'Ge')+Xev*Gn);

% Final 2D curvilinear mimetic gradient operator (d/dx, d/dy)
G = [Gx; Gy];

% Apply the operator to the field
TMP = G*C_;
Gx = TMP(1:m*(n-1));
Gy = TMP(m*(n-1)+1:end);

% Reshape for visualization
Gx = reshape(Gx, m, n-1)';
Gy = reshape(Gy, m-1, n)';

% Plot results
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
subplot(3, 1, 2)
surf((X(1:end-1, :)+X(2:end, :))/2, (Y(1:end-1, :)+Y(2:end, :))/2, Gx,...
    'EdgeColor', 'none');
colorbar
xlabel('x')
ylabel('y')
title('U')
axis equal
view([0 90])
subplot(3, 1, 3)
surf((X(:, 1:end-1)+X(:, 2:end))/2, (Y(:, 1:end-1)+Y(:, 2:end))/2, Gy,...
    'EdgeColor', 'none');
colorbar
xlabel('x')
ylabel('y')
title('V')
axis equal
view([0 90])
