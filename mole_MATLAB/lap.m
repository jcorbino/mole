function L = lap(k, m, dx, bc)
% Returns a m+2 by m+2 one-dimensional mimetic laplacian operator
%
% Parameters:
%                k : Order of accuracy
%                m : Number of cells
%               dx : Step size
%               bc : (optional) boundary-condition flag, forwarded to div
%                    and grad. Pass 'periodic' (or a nonzero/true value) to
%                    build the PERIODIC Laplacian, which wraps the interior
%                    staggered stencils around the domain instead of using
%                    one-sided boundary closures. Omit (or pass 'none'/0) for
%                    the standard operator, whose boundary rows are placeholders
%                    meant to be overwritten by a BC operator such as
%                    robinBC. A periodic problem has no such operator, so
%                    every row of the periodic Laplacian is a real Laplacian.

    if nargin < 4
        bc = 'none';
    end

    D = div(k, m, dx, bc);
    G = grad(k, m, dx, bc);
    
    L = D*G;
end
