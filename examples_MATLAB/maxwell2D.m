%% 2D TMz Maxwell Solver using Mimetic Operators and Leapfrog
clear; clc; close all;

addpath("../mole_MATLAB");

% Parameters
nt = 500;
mx = 100;
my = 100;

west = 0; east = 1;
south = 0; north = 1;

dx = (east - west)/mx;
dy = (north - south)/my;
dt = 0.5 * min(dx, dy);   % CFL

% Grids
xE = [west, west+dx/2:dx:east-dx/2, east];
yE = [south, south+dy/2:dy:north-dy/2, north];
[XE, YE] = meshgrid(xE, yE);

% Mimetic operators
k = 2;
G = dt * grad2D(k, mx, dx, my, dy);
D = dt * div2D(k, mx, dx, my, dy);

% Field sizes
NE = (mx+2)*(my+2);
NBx = mx*(my+1);
NBy = (mx+1)*my;

% Initial condition
E = exp(-100*((XE-0.5).^2 + (YE-0.5).^2)); % Gaussian pulse
E = E(:);

% Magnetic field
B = zeros(NBx + NBy, 1);

% Leapfrog start (half step for B)
B = B - 0.5 * G * E;

% Pre-plot
E2D = reshape(E, my+2, mx+2);
figure('Color','w');
hSurf = surf(xE, yE, E2D);
shading interp;
colormap(hot);
axis equal
xlabel('x'); ylabel('y'); zlabel('E_z');
lighting phong; material shiny;
zlim([-1 1]); clim([-1 1]);
hLight = camlight('headlight');
titleHandle = title('Time step: 0');

% Time-stepping loop
for n = 1:nt
    % Update fields
    B = B - G * E;
    E = E - D * B;

    % Update surface
    E2D = reshape(E, my+2, mx+2);
    set(hSurf, 'ZData', E2D);
    set(titleHandle, 'String', ['Time step: ', num2str(n)]);
    drawnow;
    pause(0.01)
end
