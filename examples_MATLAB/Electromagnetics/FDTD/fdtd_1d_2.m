% 2nd example from Sullivan's book
clear; clc; close all;

KE = 200; % Number of cells

% Fields
ex = zeros(1, KE); % Electric
hy = zeros(1, KE); % Magnetic

% Pulse parameters
kc     = KE/2; % Center of the problem space
t0     = 40;   % Temporal center of the incident pulse
spread = 12;   % Width of the pulse

% Number of time steps
NSTEPS = 250;

% Set up figure with subplots
figure('Color','w');
subplot(2,1,1);
hEx = plot(ex, 'b', 'LineWidth', 1.5);
title('Electric field');
ylabel('E_x');
ylim([-2 2]);
grid on;

subplot(2,1,2);
hHy = plot(hy, 'r', 'LineWidth', 1.5);
title('Magnetic field');
xlabel('Cell index');
ylabel('H_y');
ylim([-2 2]);
grid on;

% Boundary memory variables
ex_low_m1  = 0;
ex_low_m2  = 0;
ex_high_m1 = 0;
ex_high_m2 = 0;

% Time loop
for n = 1:NSTEPS
    % Update Ex
    ex(2:KE) = ex(2:KE) + 0.5 .* (hy(1:KE-1) - hy(2:KE));

    % Gaussian pulse
    ex(kc) = ex(kc) + exp(-0.5 * ((t0 - n) / spread)^2); % Soft

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
    sgtitle(sprintf('1D FDTD, Gaussian pulse, free space, ABC, step %.0f', n));
    drawnow;
end
