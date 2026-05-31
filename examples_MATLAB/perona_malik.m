% Solves the 2D Perona-Malik nonlinear diffusion model with MOLE operators
%
% The equation is
%
%     u_t = div(c(|grad u|) grad u),
%
% where u is the cameraman image and c is an edge-stopping diffusivity.
% Small gradients have c close to 1 and are smoothed strongly. Large
% gradients have c close to 0 and are smoothed weakly, preserving edges.
%
% This example uses the mimetic gradient and divergence from ../mole_MATLAB.
% The gradient maps cell-centered values to staggered faces, then the
% nonlinear diffusivity is evaluated on those faces and passed back through
% the mimetic divergence.

clc
close all
clear

addpath('../mole_MATLAB')

% Load the standard cameraman image. MATLAB ships cameraman.tif with the
% Image Processing Toolbox; many Octave installations include it too. The
% image is converted to double precision and normalized to [0, 1].
try
    U0 = double(imread('cameraman.tif'));
catch
    error(['Could not find cameraman.tif. Make sure MATLAB can resolve ' ...
           'the standard cameraman image, or place cameraman.tif in this directory.'])
end
if ndims(U0) == 3
    % If the input image is RGB/color, convert it to grayscale before
    % applying the scalar Perona-Malik model.
    U0 = rgb2gray(U0);
end
U0 = U0 - min(U0(:));
if max(U0(:)) > 0
    U0 = U0/max(U0(:));
end

% Spatial discretization. MOLE stores boundary centers explicitly, so the
% image resolution is used for the interior cells and one layer of boundary
% values is added by replication below.
k = 2;          % Operator order of accuracy
[n, m] = size(U0);
dx = 1;         % Pixel spacing along x
dy = 1;         % Pixel spacing along y
assert(m >= 2*k+1 && n >= 2*k+1, ...
       ['Image dimensions must be at least ' num2str(2*k+1) ...
        ' pixels in each direction for this MOLE operator order.'])

% Time discretization. Since 0 <= c <= 1, this conservative explicit step is
% based on the 2D heat-equation restriction.
tFinal = 20;
dt = 0.15*min(dx, dy)^2;
nSteps = ceil(tFinal/dt);
dt = tFinal/nSteps;

% Perona-Malik parameters.
lambda = 0.08;          % Edge threshold: smaller values preserve more edges
useExponential = true;  % true: exp(-(s/lambda)^2), false: 1/(1+(s/lambda)^2)

% MOLE operators on the staggered grid.
D = div2D(k, m, dx, n, dy);
G = grad2D(k, m, dx, n, dy);

% Coordinates used only for displaying the original image pixels. The
% replicated boundary layer added below is not shown.
xImage = 1:m;
yImage = 1:n;

% Add the boundary layer required by grad2D/div2D. Replicating the nearest
% image pixels is consistent with the no-flux boundary condition imposed on
% the fluxes inside the time loop.
U0 = [U0(1, 1) U0(1, :) U0(1, end); ...
      U0(:, 1) U0 U0(:, end); ...
      U0(end, 1) U0(end, :) U0(end, end)];

% MOLE uses x as the fastest varying index. Store the matrix transposed when
% converting to and from a column vector.
u = reshape(U0.', [], 1);

% Index where the x-face gradient block ends in the vector produced by grad2D.
nxFaces = n*(m+1);

figure
set(gcf, 'Color', 'w')

for step = 1:nSteps
    gradU = G*u;
    Gx = reshape(gradU(1:nxFaces), m+1, n).';
    Gy = reshape(gradU(nxFaces+1:end), m, n+1).';

    % Estimate |grad u| on each staggered face. The normal derivative is
    % already on that face; the transverse derivative is averaged from the
    % surrounding transverse faces.
    GyOnX = zeros(n, m+1);
    GyOnX(:, 1) = 0.5*(Gy(1:n, 1) + Gy(2:n+1, 1));
    GyOnX(:, end) = 0.5*(Gy(1:n, end) + Gy(2:n+1, end));
    GyOnX(:, 2:end-1) = 0.25*( ...
        Gy(1:n, 1:end-1) + Gy(2:n+1, 1:end-1) + ...
        Gy(1:n, 2:end) + Gy(2:n+1, 2:end));

    GxOnY = zeros(n+1, m);
    GxOnY(1, :) = 0.5*(Gx(1, 1:m) + Gx(1, 2:m+1));
    GxOnY(end, :) = 0.5*(Gx(end, 1:m) + Gx(end, 2:m+1));
    GxOnY(2:end-1, :) = 0.25*( ...
        Gx(1:end-1, 1:m) + Gx(1:end-1, 2:m+1) + ...
        Gx(2:end, 1:m) + Gx(2:end, 2:m+1));

    if useExponential
        % Perona-Malik diffusivity 1:
        % c(s) = exp(-(s/lambda)^2). This choice shuts diffusion down more
        % sharply near strong gradients, so it preserves edges aggressively.
        cX = exp(-(sqrt(Gx.^2 + GyOnX.^2)/lambda).^2);
        cY = exp(-(sqrt(GxOnY.^2 + Gy.^2)/lambda).^2);
    else
        % Perona-Malik diffusivity 2:
        % c(s) = 1/(1 + (s/lambda)^2). This choice decays more gradually,
        % so it is often a more forgiving default for denoising.
        cX = 1./(1 + (sqrt(Gx.^2 + GyOnX.^2)/lambda).^2);
        cY = 1./(1 + (sqrt(GxOnY.^2 + Gy.^2)/lambda).^2);
    end

    % Build the nonlinear flux q = c(|grad u|) grad u on staggered faces.
    % Boundary fluxes are set to zero, giving homogeneous Neumann/no-flux
    % behavior, a natural boundary condition for image denoising.
    Fx = cX.*Gx;
    Fy = cY.*Gy;
    Fx(:, 1) = 0;
    Fx(:, end) = 0;
    Fy(1, :) = 0;
    Fy(end, :) = 0;
    flux = [reshape(Fx.', [], 1); reshape(Fy.', [], 1)];

    u = u + dt*(D*flux);
    u = min(max(u, 0), 1);

    if mod(step, 25) == 0 || step == 1 || step == nSteps
        U = reshape(u, m+2, n+2).';
        imagesc(xImage, yImage, U(2:end-1, 2:end-1))
        axis image
        colormap(gca, gray)
        colorbar
        title(sprintf('Perona-Malik diffusion, t = %.4f', step*dt))
        xlabel('x')
        ylabel('y')
        drawnow
    end
end

U = reshape(u, m+2, n+2).';
gradU = G*u;
Gx = reshape(gradU(1:nxFaces), m+1, n).';
Gy = reshape(gradU(nxFaces+1:end), m, n+1).';
GyOnX = zeros(n, m+1);
GyOnX(:, 1) = 0.5*(Gy(1:n, 1) + Gy(2:n+1, 1));
GyOnX(:, end) = 0.5*(Gy(1:n, end) + Gy(2:n+1, end));
GyOnX(:, 2:end-1) = 0.25*( ...
    Gy(1:n, 1:end-1) + Gy(2:n+1, 1:end-1) + ...
    Gy(1:n, 2:end) + Gy(2:n+1, 2:end));
GxOnY = zeros(n+1, m);
GxOnY(1, :) = 0.5*(Gx(1, 1:m) + Gx(1, 2:m+1));
GxOnY(end, :) = 0.5*(Gx(end, 1:m) + Gx(end, 2:m+1));
GxOnY(2:end-1, :) = 0.25*( ...
    Gx(1:end-1, 1:m) + Gx(1:end-1, 2:m+1) + ...
    Gx(2:end, 1:m) + Gx(2:end, 2:m+1));
if useExponential
    % Same Perona-Malik diffusivity 1 used during the time integration.
    cX = exp(-(sqrt(Gx.^2 + GyOnX.^2)/lambda).^2);
    cY = exp(-(sqrt(GxOnY.^2 + Gy.^2)/lambda).^2);
else
    % Same Perona-Malik diffusivity 2 used during the time integration.
    cX = 1./(1 + (sqrt(Gx.^2 + GyOnX.^2)/lambda).^2);
    cY = 1./(1 + (sqrt(GxOnY.^2 + Gy.^2)/lambda).^2);
end

% Interpolate the staggered edge indicator to centers for visualization.
edgeCenter = zeros(n+2, m+2);
edgeCenter(2:end-1, 2:end-1) = 0.5*( ...
    0.5*(cX(:, 1:end-1) + cX(:, 2:end)) + ...
    0.5*(cY(1:end-1, :) + cY(2:end, :)));

figure
set(gcf, 'Color', 'w')

subplot(1, 3, 1)
imagesc(xImage, yImage, U0(2:end-1, 2:end-1))
axis image
colormap(gca, gray)
title('Original')
xlabel('x')
ylabel('y')
colorbar

subplot(1, 3, 2)
imagesc(xImage, yImage, U(2:end-1, 2:end-1))
axis image
colormap(gca, gray)
title('Perona-Malik')
xlabel('x')
ylabel('y')
colorbar

subplot(1, 3, 3)
imagesc(xImage, yImage, edgeCenter(2:end-1, 2:end-1))
axis image
colormap(gca, gray)
title('Edge indicator c')
xlabel('x')
ylabel('y')
colorbar
