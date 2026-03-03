% Solves the initial value problem y' = y - x^2 + 1 using forward Euler's method.
% This example is provided for completeness, Euler's method should not be used in practice.

clc;
clear;
close all;

% Define function dy/dx = f(x,y)
f = @(x,y) y - x.^2 + 1;

% Initial conditions
x0 = 0;
y0 = 0.5;

% Step size
h = 0.1;

% Interval
xn = 2;

% Number of steps
N = (xn - x0)/h;

% Initialize arrays
x = zeros(1, N+1);
y = zeros(1, N+1);

% Set initial values
x(1) = x0;
y(1) = y0;

% Compute next x and y
for i = 1:N
    x(i+1) = x(i) + h;
    y(i+1) = y(i) + h * f(x(i), y(i));
end

% Exact solution
x_exact = linspace(x0, xn, 100);
y_exact = (x_exact + 1).^2 - 0.5*exp(x_exact);

% Plot results
plot(x, y, 'ro-', 'LineWidth', 1.5);
hold on;
plot(x_exact, y_exact, 'b--', 'LineWidth', 1.5);
legend('Approximation', 'Exact', 'Location', 'northwest');
xlabel('x');
ylabel('y');
xlim([x0 xn]);
title('IVP solved with forward Euler');
grid on;

% Final error
exact_final = (xn + 1)^2 - 0.5*exp(xn);
error = abs(y(end) - exact_final);
disp(['Error at x = ', num2str(xn), ' is ', num2str(error)]);
