function L = lap3D(k, m, dx, n, dy, o, dz, bc)
% Returns a three-dimensional mimetic laplacian operator
%
% Parameters:
%                k : Order of accuracy
%                m : Number of cells along x-axis
%               dx : Step size along x-axis
%                n : Number of cells along y-axis
%               dy : Step size along y-axis
%                o : Number of cells along z-axis
%               dz : Step size along z-axis
%               bc : (optional) boundary-condition flag, forwarded to div3D
%                    and grad3D. Pass 'periodic' (or a nonzero/true value) to
%                    build the PERIODIC Laplacian, which wraps the interior
%                    staggered stencils around the domain instead of using
%                    one-sided boundary closures. Omit (or pass 'none'/0) for
%                    the standard operator, whose boundary rows are placeholders
%                    meant to be overwritten by a BC operator such as
%                    robinBC3D. A periodic problem has no such operator, so
%                    every row of the periodic Laplacian is a real Laplacian.
%                    The flag may also be given per axis, as a three-entry cell
%                    array or vector ordered {x, y, z}.

    if nargin < 8
        bc = 'none';
    end

    D = div3D(k, m, dx, n, dy, o, dz, bc);
    G = grad3D(k, m, dx, n, dy, o, dz, bc);
    
    L = D*G;
end
