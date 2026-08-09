% Solves the 2D advection-diffusion equation for a passive scalar carried by a
% static field of Taylor-Green vortices, using Corbino-Castillo mimetic operators
% and RK4.
%
%   ds/dt = div(D grad s) - div(v s)    on [-pi, pi]^2, periodic
%
%   v(x, y) = [10 sin(x) cos(y), -10 cos(x) sin(y)]   Taylor-Green, div v = 0
%   s(x, y, 0) = exp(-r^2/(2 sigma^2))                r from (0, 0)
%
%   Every operator carries the 'periodic' flag, which wraps the
%   interior staggered stencil around the seam, and mass is then conserved to
%   roundoff. That is the BC this problem wants anyway: the domain, the vortices
%   and the velocity field are all 2*pi-periodic, and the Taylor-Green field is
%   tangent to the boundary (v.n = 0 on all four sides).
%
%   dt is derived from the spectral radius of the assembled operator, so it
%   stays correct for any k, grid, Diff. Coeff. or velocity field.

clc
close all

addpath('../mole_MATLAB')

% Parameters
k = 2;          % Order of accuracy of the mimetic operators
m = 51;         % Number of cells along the x-axis
n = 51;         % Number of cells along the y-axis
a = -pi;        % West
b = pi;         % East
c = -pi;        % South
d = pi;         % North

D_coeff = 0.5;  % Diffusivity
TIME = 5;       % s Simulated time

dx = (b-a)/m;
dy = (d-c)/n;

% 2D staggered grid
% Cell centers
xc = [a a+dx/2 : dx : b-dx/2 b];   % m+2 long
yc = [c c+dy/2 : dy : d-dy/2 d];   % n+2 long
% Faces
xf = a : dx : b;                   % m+1 long
yf = c : dy : d;                   % n+1 long

% MOLE flattens 2D fields with x running fastest, which is exactly what ndgrid
% produces, so no transposing is needed anywhere except when plotting
[X, Y] = ndgrid(xc, yc);

% Mimetic operators
Dm = div2D(k, m, dx, n, dy, 'periodic');       % Divergence
G  = grad2D(k, m, dx, n, dy, 'periodic');      % Gradient
I  = interpol2D(m, n, 0.5, 0.5, 'periodic');   % Centers -> faces, centered
% c1 = c2 = 0.5 keeps the advection term free of numerical diffusion. Upwinding
% (c = 1 or 0) is not an option for a single interpolator here anyway, since the
% Taylor-Green velocity changes sign inside the domain.
%
% Worth knowing before raising k: MOLE's interpolators are 2nd-order,
% so the composite div*interpol is 2nd-order whatever k is.

% Velocity field, sampled analytically on the faces it lives on
% u sits on vertical faces: (m+1) x n,  v sits on horizontal faces: m x (n+1)
[Xu, Yu] = ndgrid(xf, yc(2:n+1));
[Xv, Yv] = ndgrid(xc(2:m+1), yf);
u = 10*sin(Xu).*cos(Yu);
v = -10*cos(Xv).*sin(Yv);
vel = [u(:); v(:)];

% Assemble the (constant, linear) right-hand side  ds/dt = A s
% A = div(D grad .) - div(v .). Both terms are premultiplied once, so the time
% loop is four sparse matrix-vector products per step.
A_diffusion = D_coeff*(Dm*G);       % scalar D; use Dm*spdiags(Dv,...)*G if variable
A_advection = Dm*spdiags(vel, 0, numel(vel), numel(vel))*I;
A = A_diffusion - A_advection;

% Time step
% The absolute-stability region of classical RK4 contains the whole left
% half-disk |z| <= 2.6155, Re z <= 0 (it cannot contain a full disk, since
% |R(eps)| > 1 for real eps > 0). A is dissipative, so its spectrum lies in
% Re <= 0, and dt*rho(A) <= 2.6155 then guarantees |R(dt*lambda)| <= 1 for
% every eigenvalue, whatever the mix of advection (imaginary) and diffusion
% (real).
% Any induced norm bounds the spectral radius, and the infinity-norm (largest
% absolute row sum) is exact, instant, and lands within 0.3-6% of rho(A) for
% these operators -- no knowledge of the velocity field or of the stencil
% constants needed, so refining the grid or raising k just works. Use inf and
% not 1: the column sums are ~2.5x looser here.
safety = 0.9;
dt = safety*2.6155/norm(A, inf);
steps = ceil(TIME/dt);
dt = TIME/steps;                    % trim so the run lands exactly on TIME

% Initial condition: Gaussian sitting on the center of the domain
x0 = 0;
y0 = 0;
sigma = pi/8;
s = exp(-((X-x0).^2 + (Y-y0).^2)/(2*sigma^2));
s = s(:);

% The static velocity field over the initial scalar
% MOLE's layout has x running fastest; transpose only here, where MATLAB's
% plotting routines expect the meshgrid convention.
figure('Name', 'Initial condition and velocity field')
pcolor(X', Y', reshape(s, m+2, n+2)')
shading interp
hold on
sk = max(1, round(m/16));               % subsample so arrows stay readable
[Xq, Yq] = meshgrid(xc(2:sk:m+1), yc(2:sk:n+1));
quiver(Xq, Yq, 10*sin(Xq).*cos(Yq), -10*cos(Xq).*sin(Yq), 'w')
hold off
axis equal
axis([a b c d])
title('s(x, y, 0) and the Taylor-Green velocity field')
xlabel('x')
ylabel('y')
colorbar
colormap jet
set(gcf, 'color', 'w')

% Time loop: classical RK4, from mole_MATLAB/rk4.m
% Only the operator is handed over -- rk4 needs nothing but A*s, so the same
% call works for any problem, and a handle would work too if A were nonlinear
% or never assembled:  @(t, s) D_coeff*(Dm*(G*s)) - Dm*(vel.*(I*s))
% The callback fires after every step. It keeps the memory flat: without it rk4
% returns the whole trajectory, which here would be (m+2)*(n+2) by steps.
figure('Name', 'Advection-diffusion of a passive scalar')
[~, s] = rk4(A, [0 TIME], dt, s, ...
             @(t, s, step) tgvPlotFrame(t, s, step, steps, X, Y, m, n, a, b, c, d));
