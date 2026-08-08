function G = grad2D(k, m, dx, n, dy, bc)
% Returns a two-dimensional mimetic gradient operator
%
% Parameters:
%                k : Order of accuracy
%                m : Number of cells along x-axis
%               dx : Step size along x-axis
%                n : Number of cells along y-axis
%               dy : Step size along y-axis
%               bc : (optional) boundary-condition flag. Pass 'periodic' (or a
%                    nonzero/true value) to build the PERIODIC operator, which
%                    wraps the interior staggered stencil around the domain
%                    instead of using one-sided boundary closures. Omit (or pass
%                    'none'/0) for the standard mimetic operator.
%
%                    The flag may also be given per axis, as a two-entry cell
%                    array or vector ordered {x, y}, e.g. {'periodic', 'none'}
%                    for a channel that is periodic along x only.
%
%                    Along a periodic axis G ignores the boundary entries of the
%                    scalar field (the corresponding columns come out empty),
%                    since a periodic problem carries no independent boundary
%                    unknowns -- the cell values alone close the system.

    if nargin < 6
        bc = 'none';
    end

    [periodic_x, periodic_y] = mole_periodic_flags2D(bc);

    Gx = grad(k, m, dx, periodic_x);
    Gy = grad(k, n, dy, periodic_y);

    % Im' and In' only select the interior cells of the transverse index, which
    % is what the faces need in every case, so these stay untouched by bc.
    Im = sparse(m + 2, m);
    In = sparse(n + 2, n);

    Im(2:(m+2)-1, :) = speye(m, m);
    In(2:(n+2)-1, :) = speye(n, n);

    Sx = kron(In', Gx);
    Sy = kron(Gy, Im');

    G = [Sx; Sy];
end

function [periodic_x, periodic_y] = mole_periodic_flags2D(bc)
% Normalizes the bc argument into one logical flag per axis.
    if iscell(bc)
        assert(numel(bc) == 2, 'a cell-array bc needs one flag per axis');
        periodic_x = mole_isperiodic2D(bc{1});
        periodic_y = mole_isperiodic2D(bc{2});
    elseif ~ischar(bc) && numel(bc) == 2
        periodic_x = mole_isperiodic2D(bc(1));
        periodic_y = mole_isperiodic2D(bc(2));
    else
        periodic_x = mole_isperiodic2D(bc);
        periodic_y = periodic_x;
    end
end

function tf = mole_isperiodic2D(flag)
    if ischar(flag) || isa(flag, 'string')
        assert(any(strcmpi(flag, {'periodic', 'none'})), ...
               'bc must be ''periodic'' or ''none''');
        tf = strcmpi(flag, 'periodic');
    else
        tf = ~isempty(flag) && all(logical(flag));
    end
end
