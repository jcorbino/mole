% 1D Inviscid Burgers' equation
% Conservative form: u_t + (u^2/2)_x = 0
% Upwind mimetic discretization with adaptive CFL timestep
% Initial condition: exp(-x^2/50)

clc
close all

addpath('../mole_MATLAB')

% Domain's limits
west = -15;
east = 15;

k = 2;              % Mimetic operator accuracy
m = 300;            % Number of cells
dx = (east - west)/m;

CFL = 0.5;          % Stable CFL factor
tfinal = 10;        % Simulation time

% Mimetic operators
D = div(k, m, dx);    % Divergence operator
I = interpol(m, 1);   % Upwind interpolation operator

% Premultiply out of the time loop
DI = D * I;

% Staggered grid
xgrid = [west west+dx/2 : dx : east-dx/2 east]';

% Initial condition
U = exp(-(xgrid.^2)/50);

% Time initialization
t = 0;
step = 0;

% Plot setup
figure
h = plot(xgrid, U, 'LineWidth', 2);
grid on
xlabel('x')
ylabel('u(x,t)')

% Time integration loop
while t < tfinal
    % Maximum wave speed (Burgers characteristic speed)
    umax = max(abs(U));

    % Adaptive CFL timestep
    dt = CFL * dx / max(umax, 1e-10);

    % Prevent overshoot of final time
    if t + dt > tfinal
        dt = tfinal - t;
    end

    % Explicit conservative update
    F = 0.5 * U.^2;
    U = U - dt * (DI * F);

    % Update time
    t = t + dt;
    step = step + 1;

    if mod(step, 10) == 0 || t >= tfinal
        % Update plot
        set(h, 'YData', U);
        title(sprintf('t = %.2f', t));
        drawnow;

        % Conservation check
        area = trapz(xgrid, U);
        fprintf('area = %.4f\n', area)
    end
end
