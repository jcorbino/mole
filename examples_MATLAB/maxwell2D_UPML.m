%% 2D TMz Maxwell Solver using Mimetic Operators, Leapfrog, and UPML
clear; clc; close all;

addpath("../mole_MATLAB");

% Parameters
nt = 200; 
mx = 100;
my = 100;

west  = 0; east  = 1;
south = 0; north = 1;

dx = (east - west)/mx;
dy = (north - south)/my;
dt = 0.5 * min(dx, dy);   % CFL condition

% Grids
xE = [west, west+dx/2:dx:east-dx/2, east];
yE = [south, south+dy/2:dy:north-dy/2, north];
[XE, YE] = meshgrid(xE, yE);

% Mimetic operators
k = 2;
G = dt * grad2D(k, mx, dx, my, dy);
D = dt * div2D(k, mx, dx, my, dy);

% Field sizes
NBx = mx*(my+1);
NBy = (mx+1)*my;
NB  = NBx + NBy;

% Initial condition for E_z
E = exp(-100*((XE-0.5).^2 + (YE-0.5).^2));
E = E(:);

% Magnetic field (B = [Bx; By])
B = zeros(NB, 1);

% PML parameters (UPML-style damping on E and B)
pml_size  = 30;   % number of cells at each boundary
sigma_max = 100;  % max conductivity
m         = 4;    % polynomial grading
eps0      = 1;    % normalized permittivity
mu0       = 1;    % normalized permeability

% 1D conductivity profiles in x and y for E-grid indices
sigma_x = zeros(mx+2,1);
sigma_y = zeros(my+2,1);

for i = 1:mx+2
    if i <= pml_size
        sigma_x(i) = sigma_max * ((pml_size - i + 1)/pml_size)^m;
    elseif i >= mx+3 - pml_size
        sigma_x(i) = sigma_max * ((i - (mx+2 - pml_size))/pml_size)^m;
    end
end

for j = 1:my+2
    if j <= pml_size
        sigma_y(j) = sigma_max * ((pml_size - j + 1)/pml_size)^m;
    elseif j >= my+3 - pml_size
        sigma_y(j) = sigma_max * ((j - (my+2 - pml_size))/pml_size)^m;
    end
end

% 2D conductivity map for E field
[SIGX, SIGY] = meshgrid(sigma_x, sigma_y);
SIGE = SIGX + SIGY;   % total electric damping

% Damping factor for E update
aE = exp(-SIGE * dt / eps0);
aE = aE(:);   % vectorize to match E

% Build B-field damping
sigmaBx = zeros(my+1, mx);
sigmaBy = zeros(my, mx+1);

for j = 1:my+1
    for i = 1:mx
        sx = sigma_x(i+1);
        sy = sigma_y(j);
        sigmaBx(j,i) = sx + sy;
    end
end

for j = 1:my
    for i = 1:mx+1
        sx = sigma_x(i);
        sy = sigma_y(j+1);
        sigmaBy(j,i) = sx + sy;
    end
end

sigmaB = [sigmaBx(:); sigmaBy(:)];
aB = exp(-sigmaB * dt / mu0);

% Leapfrog start (half step for B with damping)
B = aB .* (B - 0.5 * G * E);

% Pre-plot
E2D = reshape(E, my+2, mx+2);
figure('Color','w');
hSurf = surf(xE, yE, E2D);
shading interp;
colormap(hot);
xlabel('x'); ylabel('y'); zlabel('E_z');
lighting phong; material shiny;
zlim([-1 1]); clim([-1 1]);
camlight('headlight');
titleHandle = title('Gaussian Pulse Propagation with PML. Time step: 0');

% Time-stepping loop
for n = 1:nt
    % Update magnetic field with damping
    B = aB .* (B - G * E);

    % Update electric field with damping
    E = aE .* (E - D * B);

    % Update surface plot
    E2D = reshape(E, my+2, mx+2);
    set(hSurf, 'ZData', E2D);
    set(titleHandle, 'String', ['Gaussian Pulse Propagation with PML. Time step: ' num2str(n)]);
    drawnow;
    pause(0.01)
end
