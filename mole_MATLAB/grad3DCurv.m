function G = grad3DCurv(k, X, Y, Z)
    % Get the determinant of the jacobian and the metrics
    [J, Xe, Xn, Xc, Ye, Yn, Yc, Ze, Zn, Zc] = jacobian3D(k, X, Y, Z);
    
    % Dimensions of nodal grid
    [n, m, o] = size(X);
    
    % Make them volumes so they can be interpolated
    J = reshape(J, m, n, o)';
    Xe = reshape(Xe, m, n, o)';
    Xn = reshape(Xn, m, n, o)';
    Xc = reshape(Xc, m, n, o)';
    Ye = reshape(Ye, m, n, o)';
    Yn = reshape(Yn, m, n, o)';
    Yc = reshape(Yc, m, n, o)';
    Ze = reshape(Ze, m, n, o)';
    Zn = reshape(Zn, m, n, o)';
    Zc = reshape(Zc, m, n, o)';
    
    % Logical grids
    [Xl, Yl, Zl] = meshgrid(1:m, 1:n, 1:o);
    
    % Interpolate the metrics on the logical grid for positions u, v and w
    Ju = interp3(Xl, Yl, Zl, J, (Xl(1:end-1, :)+Xl(2:end, :))/2,...
                                            (Yl(1:end-1, :)+Yl(2:end, :))/2);
    Jv = interp3(Xl, Yl, Zl, J, (Xl(:, 1:end-1)+Xl(:, 2:end))/2,...
                                            (Yl(:, 1:end-1)+Yl(:, 2:end))/2);
    Xev = interp3(Xl, Yl, Zl, Xe, (Xl(:, 1:end-1)+Xl(:, 2:end))/2,...
                                            (Yl(:, 1:end-1)+Yl(:, 2:end))/2);
    Xnv = interp3(Xl, Yl, Zl, Xn, (Xl(:, 1:end-1)+Xl(:, 2:end))/2,...
                                            (Yl(:, 1:end-1)+Yl(:, 2:end))/2);
    Yeu = interp3(Xl, Yl, Zl, Ye, (Xl(1:end-1, :)+Xl(2:end, :))/2,...
                                            (Yl(1:end-1, :)+Yl(2:end, :))/2);
    Ynu = interp3(Xl, Yl, Zl, Yn, (Xl(1:end-1, :)+Xl(2:end, :))/2,...
                                            (Yl(1:end-1, :)+Yl(2:end, :))/2);
    
    % Convert metrics to diagonal matrices so they can be multiplied by the 
    % logical operators
    Ju = spdiags(1./reshape(Ju', [], 1), 0, numel(Ju), numel(Ju));
    Jv = spdiags(1./reshape(Jv', [], 1), 0, numel(Jv), numel(Jv));
    Xev = spdiags(reshape(Xev', [], 1), 0, numel(Xev), numel(Xev));
    Xnv = spdiags(reshape(Xnv', [], 1), 0, numel(Xnv), numel(Xnv));
    Yeu = spdiags(reshape(Yeu', [], 1), 0, numel(Yeu), numel(Yeu));
    Ynu = spdiags(reshape(Ynu', [], 1), 0, numel(Ynu), numel(Ynu));
    
    % Construct 3D uniform mimetic gradient operator (d/de, d/dn, d/dc)
    G = grad3D(k, m-1, 1, n-1, 1, o-1, 1);
    Ge = G(1:m*(n-1)*o, :);
    Gn = G(m*(n-1)*o+1:m*(n-1)*o+(m-1)*n*o, :);
    Gc = G(m*(n-1)*o+(m-1)*n*o+1:end, :);
    
    % Apply transformation
    Gx = Ju*(Ynu*Ge-Yeu*GI2(Gn, m-1, n-1, 'Gn'));
    Gy = Jv*(-Xnv*GI2(Ge, m-1, n-1, 'Ge')+Xev*Gn);
    Gz = Gy;
    
    % Final 3D curvilinear mimetic gradient operator (d/dx, d/dy, d/dz)
    G = [Gx; Gy; Gz];
end