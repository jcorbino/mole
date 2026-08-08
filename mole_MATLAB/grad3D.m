function G = grad3D(k, m, dx, n, dy, o, dz, bc)
% Returns a three-dimensional mimetic gradient operator
%
% Parameters:
%                k : Order of accuracy
%                m : Number of cells along x-axis
%               dx : Step size along x-axis
%                n : Number of cells along y-axis
%               dy : Step size along y-axis
%                o : Number of cells along z-axis
%               dz : Step size along z-axis
%               bc : (optional) boundary-condition flag. Pass 'periodic' (or a
%                    nonzero/true value) to build the PERIODIC operator, which
%                    wraps the interior staggered stencil around the domain
%                    instead of using one-sided boundary closures. Omit (or pass
%                    'none'/0) for the standard mimetic operator.
%
%                    The flag may also be given per axis, as a three-entry cell
%                    array or vector ordered {x, y, z}, e.g.
%                    {'periodic', 'periodic', 'none'} for a duct that is
%                    periodic in x and y only.
%
%                    Along a periodic axis G ignores the boundary entries of the
%                    scalar field (the corresponding columns come out empty),
%                    since a periodic problem carries no independent boundary
%                    unknowns -- the cell values alone close the system.

    if nargin < 8
        bc = 'none';
    end

    [periodic_x, periodic_y, periodic_z] = mole_periodic_flags3D(bc);

    Gx = grad(k, m, dx, periodic_x);
    Gy = grad(k, n, dy, periodic_y);
    Gz = grad(k, o, dz, periodic_z);

    % Im', In' and Io' only select the interior cells of the transverse indices,
    % which is what the faces need in every case, so these stay untouched by bc.
    Im = sparse(m + 2, m);
    In = sparse(n + 2, n);
    Io = sparse(o + 2, o);

    Im(2:(m + 2) - 1, :) = speye(m, m);
    In(2:(n + 2) - 1, :) = speye(n, n);
    Io(2:(o + 2) - 1, :) = speye(o, o);

    Sx = kron(kron(Io', In'), Gx);
    Sy = kron(kron(Io', Gy), Im');
    Sz = kron(kron(Gz, In'), Im');

    G = [Sx; Sy; Sz];
end

function [periodic_x, periodic_y, periodic_z] = mole_periodic_flags3D(bc)
% Normalizes the bc argument into one logical flag per axis.
    if iscell(bc)
        assert(numel(bc) == 3, 'a cell-array bc needs one flag per axis');
        periodic_x = mole_isperiodic3D(bc{1});
        periodic_y = mole_isperiodic3D(bc{2});
        periodic_z = mole_isperiodic3D(bc{3});
    elseif ~ischar(bc) && numel(bc) == 3
        periodic_x = mole_isperiodic3D(bc(1));
        periodic_y = mole_isperiodic3D(bc(2));
        periodic_z = mole_isperiodic3D(bc(3));
    else
        periodic_x = mole_isperiodic3D(bc);
        periodic_y = periodic_x;
        periodic_z = periodic_x;
    end
end

function tf = mole_isperiodic3D(flag)
    if ischar(flag) || isa(flag, 'string')
        assert(any(strcmpi(flag, {'periodic', 'none'})), ...
               'bc must be ''periodic'' or ''none''');
        tf = strcmpi(flag, 'periodic');
    else
        tf = ~isempty(flag) && all(logical(flag));
    end
end
