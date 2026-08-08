function D = div3D(k, m, dx, n, dy, o, dz, bc)
% Returns a three-dimensional mimetic divergence operator
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
%                    Along a periodic axis the two boundary nodes are the same
%                    physical point, so the transverse components reaching them
%                    are taken as the order-k staggered midpoint interpolation
%                    of the cells straddling the seam rather than being left at
%                    zero. That makes every row of D -- boundary faces, edges
%                    and corners included -- an actual divergence, which is what
%                    a periodic problem needs, since there is no BC operator
%                    coming afterwards to overwrite those rows.

    if nargin < 8
        bc = 'none';
    end

    [periodic_x, periodic_y, periodic_z] = mole_periodic_flags3D(bc);

    Dx = div(k, m, dx, periodic_x);
    Dy = div(k, n, dy, periodic_y);
    Dz = div(k, o, dz, periodic_z);

    Im = mole_lift3D(m, k, periodic_x);
    In = mole_lift3D(n, k, periodic_y);
    Io = mole_lift3D(o, k, periodic_z);

    Sx = kron(kron(Io, In), Dx);
    Sy = kron(kron(Io, Dy), Im);
    Sz = kron(kron(Dz, In), Im);

    D = [Sx Sy Sz];
end

function P = mole_lift3D(m, k, is_periodic)
% Returns the (m+2) by m matrix that lifts the m interior cell values of a
% transverse index into the (m+2)-long cell-space layout of the divergence
% output: the interior lands in rows 2:m+1 and the two boundary nodes are
% rows 1 and m+2.
    P = sparse(m+2, m);
    P(2:(m+2)-1, :) = speye(m, m);

    if is_periodic
        % Both boundary nodes are the same point, sitting halfway between the
        % last cell and the first one, so each takes the order-k staggered
        % midpoint interpolation of the cells straddling the seam. A plain
        % [1/2 1/2] average would be second-order and would cap the whole
        % boundary shell at second order however large k is, while the interior
        % stayed order k.
        switch k
            case 2
                w = [1/2 1/2];
            case 4
                w = [-1/16 9/16 9/16 -1/16];
            case 6
                w = [3/256 -25/256 150/256 150/256 -25/256 3/256];
            case 8
                w = [-5/2048 49/2048 -245/2048 1225/2048 ...
                     1225/2048 -245/2048 49/2048 -5/2048];
        end

        % Cells m-p+1..m and 1..p, straddling the seam. div already asserts
        % m >= 2k+1, so the stencil never wraps onto itself.
        p = numel(w)/2;
        off = [(-p+1):0, 1:p];
        for t = 1 : numel(w)
            col = mod(m + off(t) - 1, m) + 1;
            P(1,   col) = P(1,   col) + w(t);
            P(m+2, col) = P(m+2, col) + w(t);
        end
    end
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
