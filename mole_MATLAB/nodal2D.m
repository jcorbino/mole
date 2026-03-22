function N = nodal2D(k, m, dx, n, dy, nu)
% Returns a two-dimensional operator that approximates the nu-th order
% partial derivatives on a uniform nodal grid
%
% Parameters:
%                k : Order of accuracy
%                m : Number of nodes along x-axis
%               dx : Step size along x-axis
%                n : Number of nodes along y-axis
%               dy : Step size along y-axis
%               nu : (Optional) Order of the derivative (default = 1)

    if nargin < 6
        nu = 1;
    end

    % Validate grid dimensions (k, nu compatibility is checked inside nodal)
    assert(k + 1 <= m, ...
        'Not enough cells along x for desired order of accuracy k=%d and derivative order nu=%d (got %d cells).', ...
        k, nu, m-1);
    assert(k + 1 <= n, ...
        'Not enough cells along y for desired order of accuracy k=%d and derivative order nu=%d (got %d cells).', ...
        k, nu, n-1);

    Nx = nodal(k, m, dx, nu);
    Ny = nodal(k, n, dy, nu);
    
    Im = speye(m, m);
    In = speye(n, n);
    
    N = [kron(In, Nx); kron(Ny, Im)];
end