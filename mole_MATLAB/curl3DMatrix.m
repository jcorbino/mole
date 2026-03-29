function C = curl3DMatrix(k, m, dx, n, dy, o, dz)
% curl3DMatrix  3D Mimetic Curl Operator
%
%   C = curl3DMatrix(k, m, dx, n, dy, o, dz)
%
%   Returns a sparse matrix C such that C*F computes the vector curl of
%   a face-staggered 3D vector field F = [Fx; Fy; Fz].
%
%   Parameters:
%     k  : Order of accuracy (2, 4, 6, or 8)
%     m  : Number of cells along x-axis
%     dx : Step size along x-axis
%     n  : Number of cells along y-axis
%     dy : Step size along y-axis
%     o  : Number of cells along z-axis
%     dz : Step size along z-axis
%
%   Staggering:
%     Input vector F = [Fx; Fy; Fz] where:
%       Fx : x-component at x-faces, size (m+1)*(n+2)*(o+2)
%       Fy : y-component at y-faces, size (m+2)*(n+1)*(o+2)
%       Fz : z-component at z-faces, size (m+2)*(n+2)*(o+1)
%
%     Output: vector curl at grid edges [Cx; Cy; Cz] where:
%       Cx = dFz/dy - dFy/dz  at x-edges, size (m+2)*(n+1)*(o+1)
%       Cy = dFx/dz - dFz/dx  at y-edges, size (m+1)*(n+2)*(o+1)
%       Cz = dFy/dx - dFx/dy  at z-edges, size (m+1)*(n+1)*(o+2)
%
%   Vectorization convention:
%     x varies fastest, then y, then z (column-major with ndgrid ordering).
%     E.g., Fx is stored as Fx(:) from an (m+1) x (n+2) x (o+2) array.
%
%   Key property: C * G_full = 0  (curl of gradient is identically zero)
%     where G_full is the full 3D gradient:
%       G_full = [kron(kron(Io2, In2), Gx);
%                 kron(kron(Io2, Gy), Im2);
%                 kron(kron(Gz, In2), Im2)]
%
%   Construction (block structure):
%     C = [  0     -dz    +dy  ]   where dx = I_o1 ⊗ I_n2 ⊗ Gx
%         [ +dz     0     -dx  ]         dy = I_o1 ⊗ Gy  ⊗ I_m2
%         [ -dy    +dx     0   ]         dz = Gz  ⊗ I_n? ⊗ I_m?
%     with appropriate identity sizes for each block.
%
%   Requires: grad.m
%

    % 1D mimetic gradient operators
    Gx = grad(k, m, dx);  % (m+1) x (m+2)
    Gy = grad(k, n, dy);  % (n+1) x (n+2)
    Gz = grad(k, o, dz);  % (o+1) x (o+2)

    % Identity matrices
    Im1 = speye(m + 1);  Im2 = speye(m + 2);
    In1 = speye(n + 1);  In2 = speye(n + 2);
    Io1 = speye(o + 1);  Io2 = speye(o + 2);

    % --- Face-based input sizes ---
    sFx = (m + 1) * (n + 2) * (o + 2);
    sFy = (m + 2) * (n + 1) * (o + 2);
    sFz = (m + 2) * (n + 2) * (o + 1);

    % --- Edge-based output sizes ---
    sCx = (m + 2) * (n + 1) * (o + 1);
    sCy = (m + 1) * (n + 2) * (o + 1);
    sCz = (m + 1) * (n + 1) * (o + 2);

    % === Row 1: Cx = dFz/dy - dFy/dz ===
    % dFz/dy: Fz has grid (m+2, n+2, o+1), apply Gy along y
    dFz_dy = kron(kron(Io1, Gy), Im2);    % (o+1)(n+1)(m+2) x (o+1)(n+2)(m+2)
    % dFy/dz: Fy has grid (m+2, n+1, o+2), apply Gz along z
    dFy_dz = kron(kron(Gz, In1), Im2);    % (o+1)(n+1)(m+2) x (o+2)(n+1)(m+2)

    % === Row 2: Cy = dFx/dz - dFz/dx ===
    % dFx/dz: Fx has grid (m+1, n+2, o+2), apply Gz along z
    dFx_dz = kron(kron(Gz, In2), Im1);    % (o+1)(n+2)(m+1) x (o+2)(n+2)(m+1)
    % dFz/dx: Fz has grid (m+2, n+2, o+1), apply Gx along x
    dFz_dx = kron(kron(Io1, In2), Gx);    % (o+1)(n+2)(m+1) x (o+1)(n+2)(m+2)

    % === Row 3: Cz = dFy/dx - dFx/dy ===
    % dFy/dx: Fy has grid (m+2, n+1, o+2), apply Gx along x
    dFy_dx = kron(kron(Io2, In1), Gx);    % (o+2)(n+1)(m+1) x (o+2)(n+1)(m+2)
    % dFx/dy: Fx has grid (m+1, n+2, o+2), apply Gy along y
    dFx_dy = kron(kron(Io2, Gy), Im1);    % (o+2)(n+1)(m+1) x (o+2)(n+2)(m+1)

    % === Zero blocks ===
    Z1 = sparse(sCx, sFx);
    Z2 = sparse(sCy, sFy);
    Z3 = sparse(sCz, sFz);

    % === Assemble the 3x3 block curl matrix ===
    %          [   Fx       Fy        Fz    ]
    C = [ Z1,    -dFy_dz,   dFz_dy;   ...  % Cx = dFz/dy - dFy/dz
          dFx_dz, Z2,      -dFz_dx;   ...  % Cy = dFx/dz - dFz/dx
         -dFx_dy,  dFy_dx,  Z3      ];     % Cz = dFy/dx - dFx/dy
end
