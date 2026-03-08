% Solves ODE using backward Euler method

clc;
clear;
close all;

h = 0.1;                                    % Step-size
t = 0 : h : 5;                              % Calculates up to y(5)
y = zeros(1, length(t));
y(1) = 2;                                   % Initial condition
f = @(t, y) sin(t)^2*y;                     % f(t, y)

for i = 1 : length(t) - 1                   % Stages
    old_y = y(i);
    y(i+1) = fzero(@(y) y - h*f(t(i+1), y) - old_y, old_y); % Backward Euler
end

% Exact solution
y_exact = 2*exp(t/2 - sin(2*t)/4);

% Plot both
plot(t, y, 'ro-');
hold on;
plot(t, y_exact, 'b--');
title('IVP solved with backward Euler');
xlabel('x');
ylabel('y');
legend('Approximation', 'Exact', 'Location', 'northwest');
grid on;
