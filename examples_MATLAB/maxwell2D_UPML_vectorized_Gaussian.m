%% 2D TMz Maxwell Solver with Mimetic Operators + Vectorized UPML
clear; clc; close all;
addpath("../mole_MATLAB");

% Parameters
nt = 200;
mx = 100;  my = 100;
west = 0; east = 1;  south = 0; north = 1;

dx = (east-west)/mx;
dy = (north-south)/my;
dt = 0.5 * min(dx,dy);

% Grids
xE = [west, west+dx/2:dx:east-dx/2, east];
yE = [south, south+dy/2:dy:north-dy/2, north];
[XE,YE] = meshgrid(xE,yE);

% Mimetic operators
k = 2;
G = dt * grad2D(k,mx,dx,my,dy);
D = dt * div2D(k,mx,dx,my,dy);

% Fields
E = exp(-400*((XE-0.5).^2 + (YE-0.5).^2));
E = E(:);
NBx = mx*(my+1);  NBy = (mx+1)*my;
B = zeros(NBx+NBy,1);

% ============================
%   Vectorized UPML parameters
%   (Uniaxial Perfectly Matched Layer)
%   The PML is implemented as a lossy anisotropic medium
%   with conductivity profiles σx(x) and σy(y).
%   These profiles smoothly ramp up inside the PML region
%   to absorb outgoing waves without reflection.
% ============================

pml = 30;          % PML thickness in grid cells
sigma_max = 100;   % Maximum conductivity at outer PML edge
m = 4;             % Polynomial grading exponent (smooth ramp)

% Grid indices for E-field nodes
ix = (1:mx+2)'; 
iy = (1:my+2)';

% Allocate 1D conductivity profiles
sigma_x = zeros(mx+2,1);
sigma_y = zeros(my+2,1);

% Logical masks for left/right and bottom/top PML regions
Lx = ix <= pml;                     Rx = ix >= mx+3-pml;
Ly = iy <= pml;                     Ry = iy >= my+3-pml;

% Polynomial ramp σx(x)
% Left PML
sigma_x(Lx) = sigma_max*((pml - ix(Lx) + 1)/pml).^m;
% Right PML
sigma_x(Rx) = sigma_max*((ix(Rx) - (mx+2-pml))/pml).^m;

% Polynomial ramp σy(y)
% Bottom PML
sigma_y(Ly) = sigma_max*((pml - iy(Ly) + 1)/pml).^m;
% Top PML
sigma_y(Ry) = sigma_max*((iy(Ry) - (my+2-pml))/pml).^m;

% ============================
%   2D damping for E-field
%   SIGE(i,j) = σx(i) + σy(j)
%   This creates a separable anisotropic loss profile.
%   aE is the multiplicative damping factor applied each step.
% ============================

SIGE = sigma_y + sigma_x.';   % Outer sum → full 2D conductivity map
aE = exp(-SIGE(:)*dt);        % Convert conductivity to damping factor

% ============================
%   2D damping for B-field
%   Bx and By live on staggered grids, so we use shifted indices.
%   Same idea: anisotropic loss profile for magnetic field.
% ============================

sigmaBx = sigma_y(1:my+1) + sigma_x(2:mx+1).';   % Bx grid
sigmaBy = sigma_y(2:my+1) + sigma_x(1:mx+1).';   % By grid

aB = exp(-[sigmaBx(:); sigmaBy(:)] * dt);        % Magnetic damping

% Leapfrog start (half-step for B)
B = aB .* (B - 0.5 * G * E);

% Plot
figure('Color','w');
hSurf = surf(xE,yE,reshape(E,my+2,mx+2));
shading interp; colormap(hot);
lighting phong; material shiny; camlight('headlight');

xlabel x; ylabel y; zlabel('E_z');
zlim([-1 1]); clim([-1 1]);
titleHandle = title('Gaussian Pulse with UPML. Step: 0');

% Time stepping
for n = 1:nt
    B = aB .* (B - G*E);   % Magnetic update with PML damping
    E = aE .* (E - D*B);   % Electric update with PML damping

    set(hSurf,'ZData',reshape(E,my+2,mx+2));
    set(titleHandle,'String',['Gaussian Pulse with UPML. Step: ' num2str(n)]);
    drawnow;
    pause(0.01);
end
