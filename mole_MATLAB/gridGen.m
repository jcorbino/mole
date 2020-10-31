% https://en.wikipedia.org/wiki/Transfinite_interpolation
function [X, Y] = gridGen(grid_name, m, n, plot_grid)
% Returns X and Y which are both m by n matrices that contains the physical
% coordinates
%
% Parameters:
%        grid_name : String with the name of the grid folder
%                m : Number of nodes along the horizontal axis
%                n : Number of nodes along the vertical axis
%        plot_grid : If defined -> grid will be plotted
    
    assert(nargin == 4, 'Must specify all input arguments')
    assert(m > 4 && n > 4, 'm and n must be greater than 4')
    
    addpath(['grids/' grid_name])
    
    % Logical grid
    Xl = linspace(0, 1, m);
    Yl = linspace(0, 1, n);
    
    % Allocate space for physical grid
    X = zeros(m, n);
    Y = zeros(m, n);
    
    for i = 1 : m
        u = Xl(i);
        for j = 1 : n
            v = Yl(j);
            % Transfinite interpolation
            XY = (1-v)*bottom(u)+v*top(u)+(1-u)*left(v)+u*right(v)-...
                (u*v*top(1)+u*(1-v)*bottom(1)+v*(1-u)*top(0)+(1-u)*(1-v)*bottom(0));
            X(i, j) = XY(1);
            Y(i, j) = XY(2);
        end
    end
    
    if plot_grid
        figure
        mesh(X, Y, zeros(m, n), 'Marker', '.', 'MarkerSize', 10, 'EdgeColor', 'b')
        title(['Physical grid. m = ' num2str(m) ', n = ' num2str(n)])
        set(gcf, 'color', 'w')
        axis equal
        axis off
        view([0 90])
    end
end