function G = grad2DCurv(k, X, Y)
% Returns a two-dimensional curvilinear mimetic gradient operator
%
% Parameters:
%                k : Order of accuracy
%                X : x-coordinates of meshgrid
%                Y : y-coordinates of meshgrid

    [n, m] = size(X);
    n = n-1;
    m = m-1;
    
    Ux = [(X(1:end-1, 1) + X(2:end, 1))/2 (X(1:end-1, end) + X(2:end, end))/2];
    Vy = [(Y(1, 1:end-1) + Y(1, 2:end))/2; (Y(end, 1:end-1) + Y(end, 2:end))/2];
    Vx = (X(:, 1:end-1) + X(:, 2:end))/2;
    Uy = (Y(1:end-1, :) + Y(2:end, :))/2;
    Cx = (Vx(1:end-1, :) + Vx(2:end, :))/2;
    Cy = (Uy(:, 1:end-1) + Uy(:, 2:end))/2;
    
    Cx = [Ux(:, 1) Cx];
    Cy = [Uy(:, 1) Cy];
    Cx = [Cx Ux(:, end)];
    Cy = [Cy Uy(:, end)];
    
    Cx = [[X(1, 1) Vx(1, :) X(1, end)]; Cx];
    Cy = [[Y(1, 1) Vy(1, :) Y(1, end)]; Cy];
    Cx = [Cx; [X(end, 1) Vx(end, :) X(end, end)]];
    Cy = [Cy; [Y(end, 1) Vy(end, :) Y(end, end)]];
    
    Cx_ = reshape(Cx', [], 1);
    Cy_ = reshape(Cy', [], 1);
    
    G_ = grad2D(k, m, 1, n, 1);
    Ge = G_(1:n*m+n, :);
    Gn = G_(n*m+n+1:end, :);
    
    Gex = spdiags((Ge*Cx_).^-1, 0, size(Ge, 1), size(Ge, 1));
    Sx = Gex*Ge;
    
    Gny = spdiags((Gn*Cy_).^-1, 0, size(Gn, 1), size(Gn, 1));
    Sy = Gny*Gn;
    
    G = [Sx; Sy];
end
