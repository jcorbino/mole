% Solves a scalar ODE with the explicit RK4 method, written out in full
%
% NOT a mimetic example, and not the way to integrate in time with MOLE. There
% is no PDE and no mimetic operator here: this integrates y' = sin(t)^2*y on
% [0, 5] purely to show the four RK4 stages, and it hand-rolls the stepper
% instead of calling the library's own integrator, mole_MATLAB/rk4.m.
%
% Kept as a minimal reference for the scheme itself. For the intended usage see
%   van_der_pol.m          - calls mole_MATLAB/rk4.m on a system of ODEs
%   taylor_green_vortex.m  - calls it on a PDE built from mimetic operators

h = .1;                                     % Step-size
t = 0 : h : 5;                              % Calculates up to y(5)
y = zeros(1, length(t));
y(1) = 2;                                   % Initial condition
f = @(t, y) sin(t)^2*y;                     % f(t, y)

for i = 1 : length(t) - 1                   % Stages
    k1 = f(t(i),       y(i));
    k2 = f(t(i) + h/2, y(i) + h/2*k1);
    k3 = f(t(i) + h/2, y(i) + h/2*k2);
    k4 = f(t(i) + h,   y(i) + h*k3);
    
    y(i + 1) = y(i) + h/6*(k1 + 2*k2 + 2*k3 + k4);  % y(i + 1)
end

plot(t, y)
title('4th-order approximation to y(t)')
xlabel('t')
ylabel('y')
grid on
