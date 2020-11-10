function G = grad3DCurv(k, X, Y, Z)
    % Get the determinant of the jacobian and the metrics
    [J, Xe, Xn, Xc, Ye, Yn, Yc, Ze, Zn, Zc] = jacobian3D(k, X, Y, Z);
    
    % Dimensions of nodal grid
    [n, m, o] = size(X);
    
    % Make them volumes so they can be interpolated
    J = permute(reshape(J, m, n, o), [2, 1, 3]);
    
    % Logical grids
    [Xl, Yl, Zl] = meshgrid(1:m, 1:n, 1:o);
    
    % Interpolate the metrics on the logical grid for positions u, v and w
    Xl_ = (Xl(1:end-1, :, :)+Xl(2:end, :, :))/2;
    Xl_ = (Xl_(:, :, 1:end-1)+Xl_(:, :, 2:end))/2;
    Yl_ = (Yl(1:end-1, :, :)+Yl(2:end, :, :))/2;
    Yl_ = (Yl_(:, :, 1:end-1)+Yl_(:, :, 2:end))/2;
    Zl_ = (Zl(1:end-1, :, :)+Zl(2:end, :, :))/2;
    Zl_ = (Zl_(:, :, 1:end-1)+Zl_(:, :, 2:end))/2;
    Ju = interp3(Xl, Yl, Zl, J, Xl_, Yl_, Zl_);
    
    Xl_ = (Xl(:, 1:end-1, :)+Xl(:, 2:end, :))/2;
    Xl_ = (Xl_(:, :, 1:end-1)+Xl_(:, :, 2:end))/2;
    Yl_ = (Yl(:, 1:end-1, :)+Yl(:, 2:end, :))/2;
    Yl_ = (Yl_(:, :, 1:end-1)+Yl_(:, :, 2:end))/2;
    Zl_ = (Zl(:, 1:end-1, :)+Zl(:, 2:end, :))/2;
    Zl_ = (Zl_(:, :, 1:end-1)+Zl_(:, :, 2:end))/2;
    Jv = interp3(Xl, Yl, Zl, J, Xl_, Yl_, Zl_);
    
    Xl_ = (Xl(1:end-1, :, :)+Xl(2:end, :, :))/2;
    Xl_ = (Xl_(:, 1:end-1, :)+Xl_(:, 2:end, :))/2;
    Yl_ = (Yl(1:end-1, :, :)+Yl(2:end, :, :))/2;
    Yl_ = (Yl_(:, 1:end-1, :)+Yl_(:, 2:end, :))/2;
    Zl_ = (Zl(1:end-1, :, :)+Zl(2:end, :, :))/2;
    Zl_ = (Zl_(:, 1:end-1, :)+Zl_(:, 2:end, :))/2;
    Jw = interp3(Xl, Yl, Zl, J, Xl_, Yl_, Zl_);
    
    % Convert metrics to diagonal matrices so they can be multiplied by the 
    % logical operators
    Ju = spdiags(1./reshape(permute(Ju, [2, 1, 3]), [], 1), 0, numel(Ju), numel(Ju));
    Jv = spdiags(1./reshape(permute(Jv, [2, 1, 3]), [], 1), 0, numel(Jv), numel(Jv));
    Jw = spdiags(1./reshape(permute(Jw, [2, 1, 3]), [], 1), 0, numel(Jw), numel(Jw));
    
    % Construct 3D uniform mimetic gradient operator (d/de, d/dn, d/dc)
    G = grad3D(k, m-1, 1, n-1, 1, o-1, 1);
    Ge = G(1:m*(n-1)*o, :);
    Gn = G(m*(n-1)*o+1:m*(n-1)*o+(m-1)*n*o, :);
    Gc = G(m*(n-1)*o+(m-1)*n*o+1:end, :);
    
    % Apply transformation
    %Gx = ;
    %Gy = ;
    %Gz = ;
    
    % Final 3D curvilinear mimetic gradient operator (d/dx, d/dy, d/dz)
    G = [Ge; Gn; Gc];
end