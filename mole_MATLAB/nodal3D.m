function N = nodal3D(k, m, dx, n, dy, o, dz, nu)
% Returns a three-dimensional operator that approximates the nu-th order
% partial derivatives on a uniform nodal grid
%
% Parameters:
%                k : Order of accuracy
%                m : Number of nodes along x-axis
%               dx : Step size along x-axis
%                n : Number of nodes along y-axis
%               dy : Step size along y-axis
%                o : Number of nodes along z-axis
%               dz : Step size along z-axis
%               nu : (Optional) Order of the derivative (default = 1)

    if nargin < 8
        nu = 1;
    end

    % Validate grid dimensions (k, nu compatibility is checked inside nodal)
    assert(k + 1 <= m, ...
        'Not enough cells along x for desired order of accuracy k=%d and derivative order nu=%d (got %d cells).', ...
        k, nu, m-1);
    assert(k + 1 <= n, ...
        'Not enough cells along y for desired order of accuracy k=%d and derivative order nu=%d (got %d cells).', ...
        k, nu, n-1);
    assert(k + 1 <= o, ...
        'Not enough cells along z for desired order of accuracy k=%d and derivative order nu=%d (got %d cells).', ...
        k, nu, o-1);

    Nx = nodal(k, m, dx, nu);
    Ny = nodal(k, n, dy, nu);
    Nz = nodal(k, o, dz, nu);
    
    Im = speye(m, m);
    In = speye(n, n);
    Io = speye(o, o);
    
    Sx = kron(kron(Io, In), Nx);
    Sy = kron(kron(Io, Ny), Im);
    Sz = kron(kron(Nz, In), Im);
    
    N = [Sx; Sy; Sz];
end