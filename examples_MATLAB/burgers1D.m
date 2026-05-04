% Solves the 1D Inviscid Burgers' equation.
% Upwind scheme is used and the equation is written in conservative form.
% Initial condition: exp(-x^2/50)

clc
close all

addpath('../mole_MATLAB')

west = -15; % Domain's limits
east = 15;

k = 2; % Operator's order of accuracy
m = 300; % Number of cells
dx = (east-west)/m;

t = 10; % Simulation time
dt = dx; % CFL condition with max(|u|) <= 1

D = div(k, m, dx); % 1D Mimetic divergence operator
I = interpol(m, 1); % 1D interpolator
% Use I = interpol(m, 0) (downwind) if the fluid propagates to the left

% 1D Staggered grid
xgrid = [west west+dx/2: dx :east-dx/2 east];

% Impose IC
U = exp(-(xgrid.^2)/50)';

% Premultiply out of the time loop since it does not change
D = -dt/2*D*I;

% Set up plot
figure
h = plot(xgrid, U, 'LineWidth', 2);
title(sprintf('t = %.2f', 0))
xlabel('x')
ylabel('u(x, t)')
grid on

% Time integration loop
for i = 0 : t/dt

    trapz(xgrid, U) % Check for conservation

    % Update plot
    set(h, 'YData', U)
    title(sprintf('t = %.2f', i*dt))
    drawnow

    % Advance solution
    U = U + D*U.^2;
end
