function I = interpol2D(m, n, c1, c2, bc)
% Returns a two-dimensional interpolator of 2nd-order
%                m : Number of cells along x-axis
%                n : Number of cells along y-axis
%               c1 : Left interpolation coeff.
%               c2 : Bottom interpolation coeff.
%               bc : (optional) boundary-condition flag. Pass 'periodic' (or a
%                    nonzero/true value) to build the PERIODIC interpolator, in
%                    which the boundary faces take the wrapped average of the
%                    last and the first cell instead of copying the boundary
%                    entries of the scalar field, which a periodic problem does
%                    not carry. Omit (or pass 'none'/0) for the standard one.
%
%                    The flag may also be given per axis, as a two-entry cell
%                    array or vector ordered {x, y}, e.g. {'periodic', 'none'}
%                    for a channel that is periodic along x only.
%
%                    Use this together with div2D(..., 'periodic'): the default
%                    interpolator reads the boundary unknowns instead of
%                    wrapping, which leaves the eigenvalues of D*I off the
%                    imaginary axis and slowly feeds energy into the solution.

    if nargin < 5
        bc = 'none';
    end

    [periodic_x, periodic_y] = mole_periodic_flags2D(bc);

    Ix = interpol(m, c1, periodic_x);
    Iy = interpol(n, c2, periodic_y);

    % Im' and In' only select the interior cells of the transverse index, which
    % is what the faces need in every case, so these stay untouched by bc.
    Im = sparse(m + 2, m);
    In = sparse(n + 2, n);

    Im(2:(m+2)-1, :) = speye(m, m);
    In(2:(n+2)-1, :) = speye(n, n);

    Sx = kron(In', Ix);
    Sy = kron(Iy, Im');

    I = [Sx; Sy];
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
