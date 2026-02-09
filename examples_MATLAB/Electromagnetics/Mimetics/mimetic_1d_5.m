%% 1D Maxwell with mimetic operators (Sinusoidal wave, lossy dielectric, ABC)
clear;
clc;
close all;

addpath("../../mole_MATLAB");

m = 200; % Number of cells

% Mimetic operators
k = 2; % Order of accuracy
D = div(k, m, 1); % Discrete divergence (maps H -> E)
G = grad(k, m, 1); % Discrete gradient (maps E -> H)

% Fields
ex = zeros(m+2, 1); % Electric
hy = zeros(m+1, 1); % Magnetic

% Wave parameters (mirroring Sullivan's 5th example)
freq = 700e6; % Frequency in Hz
dt = 0.01 / (2 * 3e8); % Source-phase time increment

% Dielectric lossy slab
kstart = 100; % Start of the dielectric medium
epsilon = 4; % Relative permitivity
cb = 0.5 * ones(m+2, 1); % Initialize to free space
ca = ones(m+2, 1);
sigma = 0.04; % Conductivity
epsz = 8.85419e-12; % Permitivity of free space
loss_term = dt * sigma / (2 * epsz * epsilon);
ca(kstart:end) = (1 - loss_term) / (1 + loss_term);
cb(kstart:end) = (0.5 / epsilon) / (1 + loss_term);

% Number of time steps
NSTEPS = 500;

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
    ex = ca .* ex - cb .* D * hy;

    % Sinusoidal source
    ex(5) = ex(5) + sin(2*pi*freq*dt*n); % Soft

    % Absorbing boundary conditions
    ex(1) = El;
    ex(end) = Er;

    % Update Hy
    hy = hy - 0.5 * (G * ex);

    % Update plots (restrict to interior m points to mirror Sullivan)
    set(hEx, 'YData', ex(2:end-1));
    set(hHy, 'YData', hy(1:end-1));
    sgtitle(sprintf('1D Mimetic, sinusoidal wave, striking lossy dielectric, ABC, \nstep %.0f', n));
    drawnow;
end
