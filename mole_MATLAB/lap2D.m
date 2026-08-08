function L = lap2D(k, m, dx, n, dy, bc)
% Returns a two-dimensional mimetic laplacian operator
%
% Parameters:
%                k : Order of accuracy
%                m : Number of cells along x-axis
%               dx : Step size along x-axis
%                n : Number of cells along y-axis
%               dy : Step size along y-axis
%               bc : (optional) boundary-condition flag, forwarded to div2D
%                    and grad2D. Pass 'periodic' (or a nonzero/true value) to
%                    build the PERIODIC Laplacian, which wraps the interior
%                    staggered stencils around the domain instead of using
%                    one-sided boundary closures. Omit (or pass 'none'/0) for
%                    the standard operator, whose boundary rows are placeholders
%                    meant to be overwritten by a BC operator such as
%                    robinBC2D. A periodic problem has no such operator, so
%                    every row of the periodic Laplacian is a real Laplacian.
%                    The flag may also be given per axis, as a two-entry cell
%                    array or vector ordered {x, y}.

    if nargin < 6
        bc = 'none';
    end

    D = div2D(k, m, dx, n, dy, bc);
    G = grad2D(k, m, dx, n, dy, bc);
    
    L = D*G;
end
