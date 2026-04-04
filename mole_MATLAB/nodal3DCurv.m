function [Nx, Ny, Nz] = nodal3DCurv(k, X, Y, Z, nu)
% Parameters:
%                k : Order of accuracy
%                X : x-coordinates (physical) of meshgrid
%                Y : y-coordinates (physical) of meshgrid
%                Z : z-coordinates (physical) of meshgrid
%               nu : (Optional) Order of the derivative (default = 1)
%
% Note: the curvilinear chain-rule transformation is exact only for nu = 1.
% For nu > 1, the operator applies nu-th order computational-space
% differences and then multiplies by the first-order metric coefficients;
% this does NOT yield the nu-th physical partial derivatives in general.

    if nargin < 5
        nu = 1;
    end

    % Dimensions of nodal grid
    [n, m, o] = size(X);

    % Validate grid dimensions (k, nu compatibility is checked inside nodal)
    assert(k + 2*ceil(nu/2) - 1 <= m, ...
        'Not enough cells along x for desired order of accuracy k=%d and derivative order nu=%d (got %d cells).', ...
        k, nu, m-1);
    assert(k + 2*ceil(nu/2) - 1 <= n, ...
        'Not enough cells along y for desired order of accuracy k=%d and derivative order nu=%d (got %d cells).', ...
        k, nu, n-1);
    assert(k + 2*ceil(nu/2) - 1 <= o, ...
        'Not enough cells along z for desired order of accuracy k=%d and derivative order nu=%d (got %d cells).', ...
        k, nu, o-1);

    % Get the determinant of the jacobian and the metrics
    [J, Xe, Xn, Xc, Ye, Yn, Yc, Ze, Zn, Zc] = jacobian3D(k, X, Y, Z);
    
    len = n*m*o;
    
    % Convert metrics to diagonal matrices
    J = spdiags(1./J, 0, len, len);
    A = spdiags(Yn.*Zc-Zn.*Yc, 0, len, len);
    B = spdiags(Zn.*Xc-Xn.*Zc, 0, len, len);
    C = spdiags(Xn.*Yc-Yn.*Xc, 0, len, len);
    D = spdiags(Ze.*Yc-Ye.*Zc, 0, len, len);
    E = spdiags(Xe.*Zc-Ze.*Xc, 0, len, len);
    F = spdiags(Ye.*Xc-Xe.*Yc, 0, len, len);
    G = spdiags(Ye.*Zn-Ze.*Yn, 0, len, len);
    H = spdiags(Ze.*Xn-Xe.*Zn, 0, len, len);
    I = spdiags(Xe.*Yn-Ye.*Xn, 0, len, len);
    
    % Construct 3D uniform nodal operator
    N = nodal3D(k, m, 1, n, 1, o, 1, nu); % N is tall and skinny
    Ne = N(1:len, :);
    Nn = N(len+1:2*len, :);
    Nc = N(2*len+1:end, :);
    
    % Apply transformation
    Nx = J*(A*Ne+D*Nn+G*Nc);
    Ny = J*(B*Ne+E*Nn+H*Nc);
    Nz = J*(C*Ne+F*Nn+I*Nc);
    
    %  [J*A J*D J*G;
    %   J*B J*E J*H;
    %   J*C J*F J*I]*[Ne; Nn; Nc]
end
