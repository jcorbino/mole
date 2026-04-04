% test_curl_convergence.m
%
% Convergence study for curl2DMatrix applied to F = (sin(pi*y), cos(pi*x))
% on [0,1]^2.

clc
close all

addpath('../mole_MATLAB')

order = 4;
west = 0; east = 1;
south = 0; north = 1;

ms = [10 20 40 80 160];
errs_all = zeros(size(ms));
errs_int = zeros(size(ms));
dxs = zeros(size(ms));

for idx = 1:length(ms)
    m = ms(idx); n = m;
    dx = (east - west) / m;
    dy = (north - south) / n;
    dxs(idx) = dx;

    % Staggered grid positions
    xnodes = linspace(west, east, m + 1);
    ynodes = linspace(south, north, n + 1);
    xcb = [west, (west + dx/2) : dx : (east - dx/2), east];
    ycb = [south, (south + dy/2) : dy : (north - dy/2), north];

    % Sample field on faces
    [~, Yfx] = ndgrid(xnodes, ycb);
    [Xfy, ~] = ndgrid(xcb, ynodes);
    Fx = sin(pi * Yfx(:));
    Fy = cos(pi * Xfy(:));
    F = [Fx; Fy];

    % Compute curl
    C = curl2DMatrix(order, m, dx, n, dy);
    curlF = C * F;
    curlF = reshape(curlF, m + 1, n + 1);

    % Analytical curl at nodes
    [Xn, Yn] = ndgrid(xnodes, ynodes);
    exact = -pi*sin(pi*Xn) - pi*cos(pi*Yn);

    % Error over all nodes
    errs_all(idx) = norm(curlF(:) - exact(:), inf);

    % Error over interior nodes only (skip boundary band)
    bw = order / 2;
    interior = false(m + 1, n + 1);
    interior(1+bw:end-bw, 1+bw:end-bw) = true;
    errs_int(idx) = norm(curlF(interior(:)) - exact(interior(:)), inf);
end

% Convergence rates
rates_all = log2(errs_all(1:end-1) ./ errs_all(2:end));
rates_int = log2(errs_int(1:end-1) ./ errs_int(2:end));

% Print results
fprintf('Order = %d\n\n', order);
fprintf('%6s  %12s  %6s  %12s  %6s\n', 'm', 'err_all', 'rate', 'err_interior', 'rate');
fprintf('%s\n', repmat('-', 1, 50));
for idx = 1:length(ms)
    if idx == 1
        fprintf('%6d  %12.4e  %6s  %12.4e  %6s\n', ms(idx), errs_all(idx), '--', errs_int(idx), '--');
    else
        fprintf('%6d  %12.4e  %6.2f  %12.4e  %6.2f\n', ms(idx), errs_all(idx), rates_all(idx-1), errs_int(idx), rates_int(idx-1));
    end
end

% Plot
figure('Color', 'w');
loglog(dxs, errs_all, 'o-', 'LineWidth', 1.5, 'DisplayName', 'All nodes');
hold on;
loglog(dxs, errs_int, 's-', 'LineWidth', 1.5, 'DisplayName', 'Interior only');
loglog(dxs, dxs.^order * errs_int(1)/dxs(1)^order, 'k--', 'DisplayName', sprintf('O(dx^%d)', order));
xlabel('dx');
ylabel('Max abs error');
title(sprintf('curl2DMatrix convergence (order %d)', order));
legend('Location', 'northwest');
grid on;
