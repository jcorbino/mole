function I = interpol3D(m, n, o, c1, c2, c3, bc)
% Returns a three-dimensional interpolator of 2nd-order
%                m : Number of cells along x-axis
%                n : Number of cells along y-axis
%                o : Number of cells along z-axis
%               c1 : Left interpolation coeff.
%               c2 : Bottom interpolation coeff.
%               c3 : Front interpolation coeff.
%               bc : (optional) boundary-condition flag. Pass 'periodic' (or a
%                    nonzero/true value) to build the PERIODIC interpolator, in
%                    which the boundary faces take the wrapped average of the
%                    last and the first cell instead of copying the boundary
%                    entries of the scalar field, which a periodic problem does
%                    not carry. Omit (or pass 'none'/0) for the standard one.
%
%                    The flag may also be given per axis, as a three-entry cell
%                    array or vector ordered {x, y, z}, e.g.
%                    {'periodic', 'periodic', 'none'} for a duct that is
%                    periodic in x and y only.
%
%                    Use this together with div3D(..., 'periodic'): the default
%                    interpolator reads the boundary unknowns instead of
%                    wrapping, which leaves the eigenvalues of D*I off the
%                    imaginary axis and slowly feeds energy into the solution.

    if nargin < 7
        bc = 'none';
    end

    [periodic_x, periodic_y, periodic_z] = mole_periodic_flags3D(bc);

    Ix = interpol(m, c1, periodic_x);
    Iy = interpol(n, c2, periodic_y);
    Iz = interpol(o, c3, periodic_z);

    % Im', In' and Io' only select the interior cells of the transverse indices,
    % which is what the faces need in every case, so these stay untouched by bc.
    Im = sparse(m + 2, m);
    In = sparse(n + 2, n);
    Io = sparse(o + 2, o);

    Im(2:(m + 2) - 1, :) = speye(m, m);
    In(2:(n + 2) - 1, :) = speye(n, n);
    Io(2:(o + 2) - 1, :) = speye(o, o);

    Sx = kron(kron(Io', In'), Ix);
    Sy = kron(kron(Io', Iy), Im');
    Sz = kron(kron(Iz, In'), Im');

    I = [Sx; Sy; Sz];
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
