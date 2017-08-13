function G = grad3DNonUniform(k, xticks, yticks, zticks)
% Returns a three-dimensional non-uniform mimetic gradient operator
%
% Parameters:
%                k : Order of accuracy
%                xticks : Centers' ticks (x-axis)
%                yticks : Centers' ticks (y-axis)
%                zticks : Centers' ticks (z-axis)
%                         (including the boundaries!)

    Im = sparse(m + 2, m);
    Im(2:(m + 2) - 1, :) = speye(m, m);
    
    Gx = gradNonUniform(k, xticks);
    
    In = sparse(n + 2, n);
    In(2:(n + 2) - 1, :) = speye(n, n);
    
    Gy = gradNonUniform(k, yticks);
    
    Io = sparse(o + 2, o);
    Io(2:(o + 2) - 1, :) = speye(o, o);
    
    Gz = gradNonUniform(k, zticks);
    
    Sx = kron(kron(Io', In'), Gx);
    Sy = kron(kron(Io', Gy), Im');
    Sz = kron(kron(Gz, In'), Im');
    
    G = [Sx; Sy; Sz];
end
