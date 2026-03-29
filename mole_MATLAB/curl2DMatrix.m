function C = curl2DMatrix(k, m, dx, n, dy)
% curl2DMatrix  2D Mimetic Curl Operator
%
%   C = curl2DMatrix(k, m, dx, n, dy)
%
%   Returns a sparse matrix C such that C*F computes the scalar curl of
%   a face-staggered 2D vector field F = [Fx; Fy].
%
%   Parameters:
%     k  : Order of accuracy (2, 4, 6, or 8)
%     m  : Number of cells along x-axis
%     dx : Step size along x-axis
%     n  : Number of cells along y-axis
%     dy : Step size along y-axis
%
%   Staggering:
%     Input vector F = [Fx; Fy] where:
%       Fx : x-component at x-faces, size (m+1)*(n+2)
%            (m+1 x-positions at cell boundaries,
%             n+2 y-positions at cell centers + boundary points)
%       Fy : y-component at y-faces, size (m+2)*(n+1)
%            (m+2 x-positions at cell centers + boundary points,
%             n+1 y-positions at cell boundaries)
%
%     Output: scalar curl at grid nodes (vertices), size (m+1)*(n+1)
%       curl(F) = dFy/dx - dFx/dy
%
%   Vectorization convention:
%     x varies fastest (column-major with x as first dimension).
%     E.g., Fx is stored as Fx(:) from an (m+1) x (n+2) array using ndgrid.
%
%   Key property: C * G_full = 0  (curl of gradient is identically zero)
%     where G_full = [kron(speye(n+2), Gx); kron(Gy, speye(m+2))]
%     is the full 2D gradient without boundary injection.
%
%   Construction:
%     C = [-kron(Gy, I_{m+1}),  kron(I_{n+1}, Gx)]
%     where Gx, Gy are 1D mimetic gradient operators from MOLE.
%
%   Requires: grad.m
%

    % 1D mimetic gradient operators
    Gx = grad(k, m, dx);  % (m+1) x (m+2)
    Gy = grad(k, n, dy);  % (n+1) x (n+2)

    % Identity matrices at node dimensions
    Im1 = speye(m + 1);
    In1 = speye(n + 1);

    % 2D curl: curl(F) = dFy/dx - dFx/dy
    %
    % dFx/dy: apply Gy along y-dimension of Fx (x-size m+1, y-size n+2)
    %   kron(Gy, Im1) : (n+1)(m+1) x (n+2)(m+1)
    %
    % dFy/dx: apply Gx along x-dimension of Fy (x-size m+2, y-size n+1)
    %   kron(In1, Gx) : (n+1)(m+1) x (n+1)(m+2)

    C = [-kron(Gy, Im1),  kron(In1, Gx)];
end
