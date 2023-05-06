function I = interpolCentersToNodes3D(k, m, n, o)
% interpolation operator from staggered to nodes
% m, n, o, are the number of cells in the logical x-, y-, z- axes
% nodal logical coordinates are [1:1:m]x[1:1:n]x[1:1:o]
% centers logical coordinates [1,1.5:m-0.5,m]x[1,1.5:n-0.5,n]x[1,1.5:o-0.5,o]

    I1 = interpolCentersToFacesD1D(k, m);
    I2 = interpolCentersToFacesD1D(k, n);
    I3 = interpolCentersToFacesD1D(k, o);

    I = kron(I3, kron(I2, I1));
end