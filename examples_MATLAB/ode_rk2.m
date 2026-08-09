% Solves a scalar ODE with the explicit RK2 (midpoint) method, written out in full
%
% NOT a mimetic example. There is no PDE and no mimetic operator here: this
% integrates y' = sin(t)^2*y on [0, 5] purely to show the two RK2 stages. MOLE
% ships an RK4 integrator (mole_MATLAB/rk4.m) but no RK2, so there is no library
% routine this could call.
%
% Kept as a minimal reference for the scheme itself. For time integration with
% mimetic operators see taylor_green_vortex.m, which uses mole_MATLAB/rk4.m.

h = .1;                                     % Step-size
t = 0 : h : 5;                              % Calculates up to y(5)
y = zeros(1, length(t));
y(1) = 2;                                   % Initial condition
f = @(t, y) sin(t)^2*y;                     % f(t, y)

for i = 1 : length(t) - 1                   % Stages
    k1 = f(t(i),       y(i));
    k2 = f(t(i) + h/2, y(i) + h/2*k1);
    
    y(i + 1) = y(i) + h*k2;                 % y(i + 1)
end

plot(t, y)
title('2nd-order approximation to y(t)')
xlabel('t')
ylabel('y')
grid on
