function I = interpol(m, c, bc)
% Returns a m+1 by m+2 one-dimensional interpolator of 2nd-order
%
% Parameters:
%               m : Number of cells
%               c : Left interpolation coeff.
%              bc : (optional) boundary-condition flag. Pass 'periodic'
%                   (or a nonzero/true value) to build the PERIODIC
%                   interpolator, in which the two boundary nodes -- which
%                   are the same physical point -- carry the same wrapped
%                   average of the last and first cells that every interior
%                   node carries. Omit (or pass 'none'/0) for the standard
%                   interpolator, which copies the boundary values verbatim.
%
%                   Use this together with div(k, m, dx, 'periodic'). The
%                   default interpolator reads the two boundary unknowns
%                   instead of wrapping, which leaves eigenvalues of D*I off
%                   the imaginary axis; non-dissipative time integrators such
%                   as leapfrog are then weakly unstable no matter how small
%                   the time step is.

    if nargin < 3
        bc = 'none';
    end

    % Assertions:
    assert(m >= 4, 'm >= 4');
    assert(c >= 0 && c <= 1, '0 <= c <= 1');

    if ischar(bc) || isa(bc, 'string')
        assert(any(strcmpi(bc, {'periodic', 'none'})), ...
               'bc must be ''periodic'' or ''none''');
        is_periodic = strcmpi(bc, 'periodic');
    else
        is_periodic = ~isempty(bc) && all(logical(bc));
    end

    % Dimensions of I:
    n_rows = m+1;
    n_cols = m+2;

    I = sparse(n_rows, n_cols);

    % Average between two continuous cells
    avg = [c 1-c];

    if is_periodic
        % Nodes 1 and m+1 coincide, halfway between the last cell (column
        % m+1) and the first one (column 2), so both rows get the interior
        % stencil wrapped around the seam. Columns 1 and m+2 are left empty:
        % a periodic problem carries no independent boundary unknowns.
        I(1,   [m+1 2]) = avg;
        I(end, [m+1 2]) = avg;
    else
        I(1, 1) = 1;
        I(end, end) = 1;
    end

    j = 2;
    for i = 2 : n_rows - 1
        I(i, j:j+2-1) = avg;
        j = j + 1;
    end
end
