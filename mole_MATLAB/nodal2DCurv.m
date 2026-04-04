function [Nx, Ny] = nodal2DCurv(k, X, Y, nu)
% Parameters:
%                k : Order of accuracy
%                X : x-coordinates (physical) of meshgrid
%                Y : y-coordinates (physical) of meshgrid
%               nu : (Optional) Order of the derivative (default = 1)
%
% Note: the curvilinear chain-rule transformation is exact only for nu = 1.
% For nu > 1, the operator applies nu-th order computational-space
% differences and then multiplies by the first-order metric coefficients;
% this does NOT yield the nu-th physical partial derivatives in general.

    if nargin < 4
        nu = 1;
    end

    % Dimensions of nodal grid
    [n, m] = size(X);

    % Validate grid dimensions (k, nu compatibility is checked inside nodal)
    assert(k + 2*ceil(nu/2) - 1 <= m, ...
        'Not enough cells along x for desired order of accuracy k=%d and derivative order nu=%d (got %d cells).', ...
        k, nu, m-1);
    assert(k + 2*ceil(nu/2) - 1 <= n, ...
        'Not enough cells along y for desired order of accuracy k=%d and derivative order nu=%d (got %d cells).', ...
        k, nu, n-1);

    % Get the determinant of the jacobian and the metrics
    [J, Xe, Xn, Ye, Yn] = jacobian2D(k, X, Y);
    
    len = n*m;
    
    % Convert metrics to diagonal matrices
    J = spdiags(1./J, 0, len, len);
    Xe = spdiags(Xe, 0, len, len);
    Xn = spdiags(Xn, 0, len, len);
    Ye = spdiags(Ye, 0, len, len);
    Yn = spdiags(Yn, 0, len, len);
    
    % Construct 2D uniform nodal operator
    N = nodal2D(k, m, 1, n, 1, nu); % N is tall and skinny
    Ne = N(1:len, :);
    Nn = N(len+1:end, :);
    
    % Apply transformation
    Nx = J*(Yn*Ne-Ye*Nn);
    Ny = J*(-Xn*Ne+Xe*Nn);
    
    % [J*Yn -J*Ye; J*-Xn J*Xe]*[Ne; Nn]
end
