% 2D TMz Maxwell Solver with Mimetic Operators + Vectorized UPML
clear; clc; close all;
addpath("../mole_MATLAB");

% Domain
nt = 500;
mx = 100;
my = 100;
west = 0;
east = 1;
south = 0;
north = 1;

dx = (east-west)/mx;
dy = (north-south)/my;
dt = 0.5 * min(dx,dy);

% Grids
xE = [west, west+dx/2:dx:east-dx/2, east];
yE = [south, south+dy/2:dy:north-dy/2, north];
[XE,YE] = meshgrid(xE,yE);

% Mimetic operators
k = 4;
G = dt * grad2D(k,mx,dx,my,dy);
D = dt * div2D(k,mx,dx,my,dy);

% Fields
E = zeros((mx+2)*(my+2),1);
NBx = mx*(my+1);
NBy = (mx+1)*my;
B = zeros(NBx+NBy,1);

% Vectorized UPML parameters
pml = 30;
sigma_max = 100;
m = 4;

ix = (1:mx+2)';
iy = (1:my+2)';

sigma_x = zeros(mx+2,1);
sigma_y = zeros(my+2,1);

Lx = ix <= pml;
Rx = ix >= mx+3-pml;
Ly = iy <= pml;
Ry = iy >= my+3-pml;

sigma_x(Lx) = sigma_max*((pml - ix(Lx) + 1)/pml).^m;
sigma_x(Rx) = sigma_max*((ix(Rx) - (mx+2-pml))/pml).^m;

sigma_y(Ly) = sigma_max*((pml - iy(Ly) + 1)/pml).^m;
sigma_y(Ry) = sigma_max*((iy(Ry) - (my+2-pml))/pml).^m;

SIGE = sigma_y + sigma_x.';
aE = exp(-SIGE(:)*dt); % Damping profile for E

sigmaBx = sigma_y(1:my+1) + sigma_x(2:mx+1).';
sigmaBy = sigma_y(2:my+1) + sigma_x(1:mx+1).';

aB = exp(-[sigmaBx(:); sigmaBy(:)] * dt); % Damping profile for B

% Leapfrog start
B = aB .* (B - 0.5 * G * E);

% Fancy plotting
figure('Color','w');
hSurf = surf(xE,yE,reshape(E,my+2,mx+2));
shading interp; colormap(hot);
lighting phong; material shiny; camlight('headlight');

xlabel x; ylabel y; zlabel('E_z');
zlim([-1 1]); clim([-1 1]);
titleHandle = title('Sinusoidal pulse with UPML. Step: 0');

% Sinusoidal source setup
f = 10; % Frequency
ix0 = round((mx+2)/2);
iy0 = round((my+2)/2);
src_index = sub2ind([my+2, mx+2], iy0, ix0);

% Time stepping
for n = 1:nt
    % Magnetic update
    B = aB .* (B - G*E);

    % Electric update
    E = aE .* (E - D*B);

    % Inject source
    src = sin(2*pi*f*n*dt);
    E(src_index) = E(src_index) + src;

    % Update plot
    set(hSurf,'ZData',reshape(E,my+2,mx+2));
    set(titleHandle,'String',['Sinusoidal pulse with UPML. Step: ' num2str(n)]);
    drawnow;
    pause(0.01);
end
