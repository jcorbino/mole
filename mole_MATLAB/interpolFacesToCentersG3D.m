function I = interpolFacesToCentersG3D(k, m, n, o)
% 3D interpolation from faces to centers
% centers logical coordinates [1,1.5:m-0.5,m]x[1,1.5:n-0.5,n]x[1,1.5:o-0.5,o]
% m, n, o, are the number of cells in the logical x-, y-, z- axes

    cells = (o+2)*(n+2)*(m+2);

    Ix = interpolFacesToCentersG1D(k, m);
    Iy = interpolFacesToCentersG1D(k, n);
    Iz = interpolFacesToCentersG1D(k, o);

    Im = sparse(m + 2, m);
    In = sparse(n + 2, n);
    Io = sparse(o + 2, o);

    Im(2:(m+2)-1, :) = speye(m, m);
    In(2:(n+2)-1, :) = speye(n, n);
    Io(2:(o+2)-1, :) = speye(o, o);

    Sx = kron(kron(Io, In), Ix);
    Sy = kron(kron(Io, Iy), Im);
    Sz = kron(kron(Iz, In), Im);

    I = sparse(3*cells, o*n*(m+1)+o*(n+1)*m+(o+1)*n*m);
    
    I(1:cells, 1:o*n*(m+1)) = Sx; 
    I(cells+1:2*cells, o*n*(m+1)+1:o*n*(m+1)+o*(n+1)*m) = Sy;  
    I(2*cells+1:end, o*n*(m+1)+o*(n+1)*m+1:end) = Sz;
end