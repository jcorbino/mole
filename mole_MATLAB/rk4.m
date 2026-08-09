function [t, y] = rk4(func, tspan, dt, y0, callback)
% Explicit Runge-Kutta 4th-order method
%
% Returns: t (evaluation points) and y (solutions) of the specified ODE
%
% Parameters:
%                func : Function handler f(t, y) returning dy/dt, or, for the
%                       linear case dy/dt = A*y, the matrix A itself. Only the
%                       product A*y is ever needed, so the operator may be
%                       assembled however the problem requires -- or not
%                       assembled at all:
%                           f = @(t, s) D*(k*(G*s)) - D*(vel.*(I*s));
%               tspan : [t0 tf]
%                  dt : Step size
%                  y0 : Initial conditions
%            callback : (optional) handler called as callback(t, y, step) after
%                       every step, for plotting or diagnostics. When given,
%                       only the current state is kept and y returns just the
%                       final one, so the memory stays flat. A PDE needs this:
%                       storing the whole trajectory costs
%                       numel(y0)-by-number-of-steps.

    % Linear case: wrap the operator so the loop below only sees a handle
    if ~isa(func, 'function_handle')
        A = func;
        func = @(t, y) A*y;
    end

    t = tspan(1) : dt : tspan(2);

    streaming = nargin >= 5 && ~isempty(callback);

    if streaming
        y = y0(:);
    else
        y = zeros(length(y0), length(t));
        y(:, 1) = y0;
    end

    for i = 1 : length(t) - 1
        if streaming
            yi = y;
        else
            yi = y(:, i);
        end

        k1 = func(t(i),        yi);
        k2 = func(t(i) + dt/2, yi + dt/2*k1);
        k3 = func(t(i) + dt/2, yi + dt/2*k2);
        k4 = func(t(i) + dt,   yi + dt*k3);

        ynew = yi + dt/6*(k1 + 2*k2 + 2*k3 + k4);

        if streaming
            y = ynew;
            callback(t(i + 1), y, i);
        else
            y(:, i + 1) = ynew;
        end
    end
end
