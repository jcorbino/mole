function D = div3DNonUniform(k, xticks, yticks, zticks)
% Returns a three-dimensional non-uniform mimetic divergence operator
%
% Parameters:
%                k : Order of accuracy
%                xticks : Edges' ticks (x-axis)
%                yticks : Edges' ticks (y-axis)
%                zticks : Edges' ticks (z-axis)

    Im = sparse(m + 2, m);
    Im(2:(m + 2) - 1, :) = speye(m, m);
    
    Dx = divNonUniform(k, xticks);
    
    In = sparse(n + 2, n);
    In(2:(n + 2) - 1, :) = speye(n, n);
    
    Dy = divNonUniform(k, yticks);
    
    Io = sparse(o + 2, o);
    Io(2:(o + 2) - 1, :) = speye(o, o);
    
    Dz = divNonUniform(k, zticks);
    
    Sx = kron(kron(Io, In), Dx);
    Sy = kron(kron(Io, Dy), Im);
    Sz = kron(kron(Dz, In), Im);
    
    D = [Sx Sy Sz];
end
