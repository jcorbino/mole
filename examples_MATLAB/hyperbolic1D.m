% Solves the 1D Advection equation with periodic boundary conditions

clc
close all

addpath('../mole_MATLAB')

a = 1; % Velocity
west = 0; % Domain's limits
east = 1;

k = 2; % Operator's order of accuracy
m = 50; % Number of cells
dx = (east-west)/m;

t = 1; % Simulation time

D = div(k, m, dx, 'periodic'); % 1D Mimetic divergence operator
I = interpol(m, 0.5, 'periodic'); % 1D 2nd order interpolator

rho = max(abs(eig(full(D*I))));   % Spectral radius of the spatial operator

% CFL condition for explicit schemes. Leapfrog is stable only while every
% dt*eig(-a*D*I) stays on the OPEN segment (-i, i) of the imaginary axis, so
% dt < 1/(|a|*rho). The 0.9 keeps us inside it: at dt = 1/(|a|*rho) exactly the
% two characteristic roots coalesce at +-i, and that defective double root
% makes the solution grow linearly with the step count.
dt = 0.9/(abs(a)*rho);

% Shrink dt so the last iteration lands exactly on t. Rounding the step count
% up can only make dt smaller, so this never violates the CFL bound above.
n_steps = ceil(t/dt);
dt = t/n_steps;

% 1D Staggered grid
grid = [west west+dx/2: dx :east-dx/2 east];

% IC
U = sin(2*pi*grid)';

% Premultiply out of the time loop (since it doesn't change)
D = -a*dt*2*D*I;
% One could also have said: D = -a*dt*2*I*D if the grid 
% was defined as: grid = west : dx : east (nodal)

U2 = U + D/2*U; % Compute one step using Euler's method

% Time integration loop
for i = 1 : n_steps
    plot(grid, U2, 'o-') % Plot approximation
    hold on
    plot(grid, sin(2*pi*(grid - a*i*dt))) % Plot exact solution
    hold off
    str = sprintf('t = %.2f', i*dt);
    title(str)
    xlabel('x')
    ylabel('u(x, t)')
    axis([west east -1.5 1.5])
    pause(0.04)
    U3 = U + D*U2; % Compute next step using Leapfrog scheme
    U = U2;
    U2 = U3;
end
