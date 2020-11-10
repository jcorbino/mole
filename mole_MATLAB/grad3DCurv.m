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
    Grad = grad3D(k, m-1, 1, n-1, 1, o-1, 1);
    Ge = Grad(1:m*(n-1)*(o-1), :);
    Gn = Grad(m*(n-1)*(o-1)+1:m*(n-1)*(o-1)+(m-1)*n*(o-1), :);
    Gc = Grad(m*(n-1)*(o-1)+(m-1)*n*(o-1)+1:end, :);
    
    % Apply transformation
    A = Yn.*Zc-Zn.*Yc;
    A = spdiags(reshape(permute(A, [2, 1, 3]), [], 1), 0, numel(A), numel(A));
    B = Zn.*Xc-Xn.*Zc;
    B = spdiags(reshape(permute(B, [2, 1, 3]), [], 1), 0, numel(B), numel(B));
    C = Xn.*Yc-Yn.*Xc;
    C = spdiags(reshape(permute(C, [2, 1, 3]), [], 1), 0, numel(C), numel(C));
    D = Ze.*Yc-Ye.*Zc;
    D = spdiags(reshape(permute(D, [2, 1, 3]), [], 1), 0, numel(D), numel(D));
    E = Xe.*Zc-Ze.*Xc;
    E = spdiags(reshape(permute(E, [2, 1, 3]), [], 1), 0, numel(E), numel(E));
    F = Ye.*Xc-Xe.*Yc;
    F = spdiags(reshape(permute(F, [2, 1, 3]), [], 1), 0, numel(F), numel(F));
    G = Ye.*Zn-Ze.*Yn;
    G = spdiags(reshape(permute(G, [2, 1, 3]), [], 1), 0, numel(G), numel(G));
    H = Ze.*Xn-Xe.*Zn;
    H = spdiags(reshape(permute(H, [2, 1, 3]), [], 1), 0, numel(H), numel(H));
    I = Xe.*Yn-Ye.*Xn;
    I = spdiags(reshape(permute(I, [2, 1, 3]), [], 1), 0, numel(I), numel(I));
    
    Gx = Ju*Ge;
    Gy = Jv*Gn;
    Gz = Jw*Gc;
    
    % Final 3D curvilinear mimetic gradient operator (d/dx, d/dy, d/dz)
    G = [Gx; Gy; Gz];
end