%% 5th example from Sullivan's book
clear; clc; close all;

KE = 200; % Number of cells

% Fields
ex = zeros(1, KE); % Electric
hy = zeros(1, KE); % Magnetic

% Wave parameters
freq = 700e6;        % Frequency in Hz
dt   = 0.01/(2*3e8); % Source-phase time increment

% Dielectric lossy slab
kstart  = 100;          % Start of the dielectric medium
epsilon = 4;            % Relative permitivity
cb = 0.5 * ones(1, KE); % Initialize to free space
ca = ones(1, KE);
sigma = 0.04;           % Conductivity
epsz = 8.85419e-12;     % Permitivity of free space
loss_term = dt * sigma/ (2 * epsz * epsilon);
ca(kstart:KE) = (1 - loss_term) / (1 + loss_term);
cb(kstart:KE) = (0.5 / epsilon) / (1 + loss_term);

% Number of time steps
NSTEPS = 500;

% Set up figure with subplots
figure('Color','w');
subplot(2,1,1);
hEx = plot(ex, 'b', 'LineWidth', 1.5);
title('Electric field');
ylabel('E_x');
ylim([-2 2]);
xline(kstart, '--k')
grid on;

subplot(2,1,2);
hHy = plot(hy, 'r', 'LineWidth', 1.5);
title('Magnetic field');
xlabel('Cell index');
ylabel('H_y');
ylim([-2 2]);
xline(kstart, '--k')
grid on;

% Boundary memory variables
ex_low_m1  = 0;
ex_low_m2  = 0;
ex_high_m1 = 0;
ex_high_m2 = 0;

% Time loop
for n = 1:NSTEPS
    % Update Ex
    ex(2:KE) = ca(2:KE) .* ex(2:KE) + cb(2:KE) .* (hy(1:KE-1) - hy(2:KE));

    % Sinusoidal source
    ex(5) = ex(5) + sin(2*pi*freq*dt*n); % Soft

    % Absorbing boundary conditions
    % Left
    ex(1) = ex_low_m2;
    ex_low_m2 = ex_low_m1;
    ex_low_m1 = ex(2);
    % Right
    ex(KE) = ex_high_m2;
    ex_high_m2 = ex_high_m1;
    ex_high_m1 = ex(KE-1);

    % Update Hy
    hy(1:KE-1) = hy(1:KE-1) + 0.5 * (ex(1:KE-1) - ex(2:KE));

    % Update plots
    set(hEx, 'YData', ex);
    set(hHy, 'YData', hy);
    sgtitle(sprintf('1D FDTD, sinusoidal wave, striking lossy dielectric, ABC, \nstep %.0f', n));
    drawnow;
end
