function N = nodal(k, m, dx)
% Returns a m+1 by m+1 one-dimensional operator that approximates the 
% first-order derivatives on a nodal grid
%
% Parameters:
%                k : Order of accuracy
%                m : Number of cells
%               dx : Step size
    
    N = grad(k, m, dx, 'nodal');
end