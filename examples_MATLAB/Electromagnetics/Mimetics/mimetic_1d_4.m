%% 1D Maxwell with mimetic operators (Sinusoidal wave, dielectric, ABC)
clear;
clc;
close all;

addpath("../../mole_MATLAB");

m = 200; % Number of cells
west = 0;
east = 1;
dx = (east - west) / m;
dt = 0.5 * dx; % CFL

% Mimetic operators
k = 2; % Order of accuracy
D = dt * div(k, m, dx); % Discrete divergence (maps H -> E)
G = dt * grad(k, m, dx); % Discrete gradient (maps E -> H)

% Fields
ex = zeros(m+2, 1); % Electric
hy = zeros(m+1, 1); % Magnetic

% Wave parameters (mirroring Sullivan's 4th example)
freq = 700e6; % Frequency in Hz
ddt = 0.01 / (2 * 3e8); % Source-phase time increment

% Number of time steps
NSTEPS = 425;

% Dielectric slab
kstart = 100; % Index where slab begins (in E-grid index)
epsilon = 4; % Relative permittivity in slab

% Relative permittivity
epsilon_r = ones(m+2, 1);
epsilon_r(kstart:end) = epsilon;

% Modify divergence operator
D = D ./ epsilon_r;

% Leapfrog start: stagger Hy by half time step
% Since ex = 0 initially, G*ex = 0 and this does nothing in this example.
% hy = hy - 0.5 * G * ex;

% Set up figure with subplots
figure('Color', 'w');
subplot(2, 1, 1);
hEx = plot(1:m, ex(2:end-1), 'b', 'LineWidth', 1.5);
title('Electric field');
ylabel('E_x');
ylim([-2, 2]);
xline(kstart, '--k')
grid on;

subplot(2, 1, 2);
hHy = plot(1:m, hy(1:end-1), 'r', 'LineWidth', 1.5);
title('Magnetic field');
xlabel('Cell index');
ylabel('H_y');
ylim([-2, 2]);
xline(kstart, '--k')
grid on;

% Time loop
for n = 1:NSTEPS
    % Store boundary values
    El = ex(2);
    Er = ex(end-1);

    % Update Ex
    ex = ex - D * hy;

    % Sinusoidal source
    ex(5) = ex(5) + sin(2*pi*freq*ddt*n); % Soft

    % Absorbing boundary conditions
    ex(1) = El;
    ex(end) = Er;

    % Update Hy
    hy = hy - G * ex;

    % Update plots (restrict to interior m points to mirror Sullivan)
    set(hEx, 'YData', ex(2:end-1));
    set(hHy, 'YData', hy(1:end-1));
    sgtitle(sprintf('1D Mimetic, sinusoidal wave, striking dielectric, ABC, step %.0f', n));
    drawnow;
end
