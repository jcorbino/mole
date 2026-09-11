%% 2D Lock exchange test case using mimetic methods
% It uses Boussinesq's approximation to solve for nonisothermal flow: two
% fluids of slightly different density, initially separated by a vertical
% interface in a closed tank, are released and slump past each other.
% Numerical diffusion acts to limit the formation of the billows along the
% interface, so sharper billows are obtained by increasing the resolution
% (m and n). The time step is derived from the grid, so nothing else needs
% to change.
%
% Layout: x runs along the first array index, and therefore fastest in the
% flattened vectors, which is what the MOLE operators expect, so A(:) can be
% handed to them directly. Scalars (T, p) live at the cell centers plus the
% boundary nodes, (m+2)-by-(n+2); u at the vertical faces, (m+1)-by-n; v at
% the horizontal faces, m-by-(n+1).

clc
close all

%% MOLE's path so mimetic methods can be used
addpath('../mole_MATLAB')

tic  % Start the timer
%--------------------------------------------------------------------------

%% Grid [a b]x[c d] in meters
a = 0;
b = 100;
c = 0;
d = 20;
m = 1000;      % Number of cells along x-axis
n = 200;       % Number of cells along y-axis
dx = (b-a)/m;  % Step size along x-axis
dy = (d-c)/n;  % Step size along y-axis

% Cell centers
xc = a+dx*((1:m)-0.5);
yc = c+dy*((1:n)-0.5);

% Physical coordinates of the scalars (cell centers plus boundary nodes), of
% the u-velocities (vertical faces) and of the v-velocities (horizontal faces)
[X, Y] = ndgrid([a xc b], [c yc d]);
[Xu, Yu] = ndgrid(a+dx*(0:m), yc);
[Xv, Yv] = ndgrid(xc, c+dy*(0:n));

%% Parameters (for seawater @ S = 35)
alpha = 1.664e-4;                     % 1/°C Coefficient of thermal expansion
T_0 = 10;                             % °C Reference temperature
rho_0 = 1027;                         % kg/m3 Reference density
interface_width = 0.1;                % m Width of the interface
g = 9.806;                            % m/s2 Gravity acceleration
reduced_gravity = 0.01;               % m/s2 Reduced gravity
delta_rho = reduced_gravity*rho_0/g;  % kg/m3 Density jump between the fluids
mu = 0.00141;                         % kg/(m⋅s) Dynamic viscosity
time = 200;                           % s Simulation time

%% Initial condition (heavy fluid on the left, fluid at rest)
% Linear equation of state: rho = rho_0(1-alpha(T-T_0))
rho = rho_0+delta_rho/2*(1-erf((X-(a+b)/2)/interface_width));
T = (1-rho/rho_0)/alpha+T_0;
u = zeros(m+1, n);
v = zeros(m, n+1);

%% Middle points and viscosities
T_middle = (T(1)+T(end))/2;
rho_middle = rho_0*(1-alpha*(T_middle-T_0));  % Same as (rho(1)+rho(end))/2
nu = mu/rho_middle;                           % m2/s Kinematic viscosity

%% Time step
% Both explicit (forward Euler) pieces of the scheme bound dt:
%   upwind advection:   max(|u|)dt/dx + max(|v|)dt/dy <= 1
%   centered diffusion: 2 nu dt (1/dx^2 + 1/dy^2) <= 1
% The velocity scale of a lock exchange is the buoyancy velocity sqrt(g'H)
% (the fronts travel at about half of it), so it is used for both components
% together with the safety factor CFL. This makes dt scale with the grid, so
% the run stays stable when m and n are changed.
H = d-c;
U = sqrt(reduced_gravity*H);
CFL = 0.5;
dt = CFL/(U/dx+U/dy);
dt = min(dt, 0.25/(nu*(1/dx^2+1/dy^2)));
iterations = ceil(time/dt);  % Number of iterations
dt = time/iterations;        % s Time step, landing exactly on t = time
fprintf('dt = %g s, %d iterations\n', dt, iterations)

%% Mimetic operators
k = 2;                        % Spatial order of accuracy
D = div2D(k, m, dx, n, dy);   % Divergence
G = grad2D(k, m, dx, n, dy);  % Gradient
L = D*G;                      % Laplacian

% Apply Neumann BCs to Laplacian
L = L+robinBC2D(k, m, dx, n, dy, 0, 1);

% A pure-Neumann Laplacian is singular (the constants are in its null space),
% so the pressure is pinned at one cell. The right-hand side D*u_s is
% compatible, since u_s has no normal component on the walls, hence the
% solution still satisfies the divergence-free condition dropped at that cell.
pin = sub2ind([m+2 n+2], 2, 2);  % First interior cell
L(pin, :) = 0;
L(pin, pin) = 1;

% Factorize once instead of solving from scratch at every iteration (this is
% the most time-consuming part). P*(R\L)*Q = Ll*Lu, so that
% L\b = Q*(Lu\(Ll\(P*(R\b)))).
[Ll, Lu, Lp, Lq, Lr] = lu(L);
solve = @(b) Lq*(Lu\(Ll\(Lp*(Lr\b))));

% Interpolators (cell centers to faces) for upwinding
I1 = interpol2D(m, n, 1, 1);  % Takes the left/bottom cell: upwind for UV >= 0
I0 = interpol2D(m, n, 0, 0);  % Takes the right/top cell:   upwind for UV < 0

u_length = (m+1)*n;  % The first u_length entries of a face vector are u's

%% Iterate over time
for t = 1 : iterations

    %% Predictor step
    % At k = 2 the interior mimetic stencils coincide with the centered
    % differences used here, so the operators are not needed for this step.
    % Forward Euler with centered differencing for the diffusion of momentum
    % and sided differencing (based on the sign of the velocity) for its
    % advection. The normal velocity is zero on the walls and is never updated
    % there. For the tangential velocity, no-slip is imposed through ghost
    % values that mirror the first interior value, so that the velocity
    % vanishes on the wall itself, which lies halfway between them.
    up = [-u(:, 1) u -u(:, end)];    % (m+1)-by-(n+2)
    vp = [-v(1, :); v; -v(end, :)];  % (m+2)-by-(n+1)

    % u* at the interior vertical faces (i = 2:m), with v averaged onto them
    uC = up(2:m, 2:n+1);
    uW = up(1:m-1, 2:n+1);
    uE = up(3:m+1, 2:n+1);
    uS = up(2:m, 1:n);
    uN = up(2:m, 3:n+2);
    vf = 0.25*(v(1:m-1, 1:n)+v(2:m, 1:n)+v(1:m-1, 2:n+1)+v(2:m, 2:n+1));
    d2u = (uW-2*uC+uE)/dx^2+(uS-2*uC+uN)/dy^2;
    udu_dx = max(uC, 0).*(uC-uW)/dx+min(uC, 0).*(uE-uC)/dx;
    vdu_dy = max(vf, 0).*(uC-uS)/dy+min(vf, 0).*(uN-uC)/dy;
    u_s = u;
    u_s(2:m, :) = uC+dt*(nu*d2u-(udu_dx+vdu_dy));

    % v* at the interior horizontal faces (j = 2:n), with u and T averaged
    % onto them
    vC = vp(2:m+1, 2:n);
    vW = vp(1:m, 2:n);
    vE = vp(3:m+2, 2:n);
    vS = vp(2:m+1, 1:n-1);
    vN = vp(2:m+1, 3:n+1);
    uf = 0.25*(u(1:m, 1:n-1)+u(2:m+1, 1:n-1)+u(1:m, 2:n)+u(2:m+1, 2:n));
    Tf = 0.5*(T(2:m+1, 2:n)+T(2:m+1, 3:n+1));
    d2v = (vW-2*vC+vE)/dx^2+(vS-2*vC+vN)/dy^2;
    udv_dx = max(uf, 0).*(vC-vW)/dx+min(uf, 0).*(vE-vC)/dx;
    vdv_dy = max(vC, 0).*(vC-vS)/dy+min(vC, 0).*(vN-vC)/dy;
    v_s = v;
    v_s(:, 2:n) = vC+dt*(nu*d2v-(udv_dx+vdv_dy)+g*alpha*(Tf-T_middle));

    %% Solve for pressure
    % Projection: L phi = D u*, with phi = dt p/rho_middle
    UV = [u_s(:); v_s(:)];
    rhs = D*UV;
    rhs(pin) = 0;
    phi = solve(rhs);

    %% Corrector step
    UV = UV-G*phi;
    u = reshape(UV(1:u_length), m+1, n);
    v = reshape(UV(u_length+1:end), m, n+1);

    % The corrected normal velocity vanishes on the walls only up to the
    % accuracy of the linear solve, so it is zeroed exactly
    u([1 end], :) = 0;
    v(:, [1 end]) = 0;
    UV = [u(:); v(:)];

    %% Advection of heat
    % Flux form with the divergence-free velocity, taking T from the upwind
    % cell (I1 where UV >= 0, I0 where UV < 0). Forward Euler's used here, but
    % can easily be replaced by any RK-method.
    T = T(:);
    T = T-dt*(D*(max(UV, 0).*(I1*T)+min(UV, 0).*(I0*T)));
    T = reshape(T, m+2, n+2);

    %% Stability check
    cfl = dt*(max(abs(u(:)))/dx+max(abs(v(:)))/dy);
    if ~(cfl <= 1)  % Also catches NaN
        error('CFL = %g at t = %g s: the solution blew up, decrease CFL.', ...
              cfl, t*dt)
    end
    if mod(t, ceil(iterations/10)) == 0 || t == iterations
        fprintf('t = %6.f s   CFL = %.2f\n', t*dt, cfl)
    end
end

%% Plotting
in = {2:m+1, 2:n+1};  % Cell centers (the boundary nodes of T are never updated)

% Pressure (defined up to the constant fixed by the pinned cell)
p = reshape(rho_middle/dt*phi, m+2, n+2);
subplot(2, 2, 1)
surf(X(in{:}), Y(in{:}), p(in{:}))
title('p')
shading interp
axis equal
h = colorbar;
h.Label.String = 'Pa';  % N/m2 = kg/(m⋅s2)
view([0 90])
xlabel('x (m)')
ylabel('y (m)')

% Temperature
subplot(2, 2, 2)
surf(X(in{:}), Y(in{:}), T(in{:}))
title('T')
shading interp
axis equal
h = colorbar;
h.Label.String = '°C';
view([0 90])
xlabel('x (m)')
ylabel('y (m)')

% u-velocity
subplot(2, 2, 3)
surf(Xu, Yu, u)
title('u')
shading interp
axis equal
h = colorbar;
h.Label.String = 'm/s';
view([0 90])
xlabel('x (m)')
ylabel('y (m)')

% v-velocity
subplot(2, 2, 4)
surf(Xv, Yv, v)
title('v')
shading interp
axis equal
h = colorbar;
h.Label.String = 'm/s';
view([0 90])
xlabel('x (m)')
ylabel('y (m)')

if exist('sgtitle', 'file')  % MATLAB >= R2018b
    sgtitle(['t = ' num2str(time) 's'])
end
set(gcf, 'color', 'w')
colormap jet

% Get final density profile from temperatures using equation of state
rho = rho_0*(1-alpha*(T-T_0));

% Plot the density
figure
surf(X(in{:}), Y(in{:}), rho(in{:}))
title(['\rho @ t = ' num2str(time) 's'])
shading interp
axis equal
h = colorbar;
h.Label.String = 'kg/m^3';
view([0 90])
xlabel('x (m)')
ylabel('y (m)')
set(gcf, 'color', 'w')
colormap jet

%--------------------------------------------------------------------------
toc  % Stop the timer
