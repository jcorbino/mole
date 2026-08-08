function D = div(k, m, dx, bc)
% Returns a m+2 by m+1 one-dimensional mimetic divergence operator
%
% Parameters:
%                k : Order of accuracy
%                m : Number of cells
%               dx : Step size
%               bc : (optional) boundary-condition flag. Pass 'periodic'
%                    (or a nonzero/true value) to build the PERIODIC operator,
%                    in which every cell row uses the interior (staggered)
%                    stencil wrapped around the domain, and the two boundary
%                    node rows use the order-k nodal centered difference wrapped
%                    around the seam. Omit (or pass 'none'/0) for the standard
%                    mimetic operator with one-sided boundary closures.

    if nargin < 4
        bc = 'none';
    end

    % Assertions:
    assert(k >= 2, 'k >= 2');
    assert(k <= 8, 'k <= 8');
    assert(mod(k, 2) == 0, 'k % 2 = 0');
    assert(m >= 2*k+1, ['m >= ' num2str(2*k+1) ' for k = ' num2str(k)]);

    if ischar(bc) || isa(bc, 'string')
        assert(any(strcmpi(bc, {'periodic', 'none'})), ...
               'bc must be ''periodic'' or ''none''');
        is_periodic = strcmpi(bc, 'periodic');
    else
        is_periodic = ~isempty(bc) && all(logical(bc));
    end

    D = sparse(m+2, m+1);

    if is_periodic
        % --- Periodic operator -------------------------------------------
        % S : interior (staggered) mimetic stencil, offS its column offsets
        % C : nodal centered first-difference stencil for the boundary nodes
        [S, offS, C, offC] = mole_periodic_stencils(k);

        % Interior cell rows: staggered stencil, wrapped across the seam.
        for i = 2:m+1
            for t = 1:numel(S)
                c = mole_wrapcol(i + offS(t), m);
                D(i, c) = D(i, c) + S(t);
            end
        end

        % Boundary node rows: west node = row 1 (node 1),
        %                     east node = row m+2 (node m+1).
        for t = 1:numel(C)
            cw = mole_wrapcol(1     + offC(t), m);
            ce = mole_wrapcol((m+1) + offC(t), m);
            D(1,   cw) = D(1,   cw) + C(t);
            D(m+2, ce) = D(m+2, ce) + C(t);
        end

    else
        % --- Standard mimetic operator (one-sided boundary closures) -----
        switch k
            case 2
                for i = 2:m+1
                   D(i, i-1:i) = [-1 1];
                end

            case 4
                A = [-11/12 17/24 3/8 -5/24 1/24];
                D(2, 1:5) = A;
                D(m+1, m-3:end) = -fliplr(A);
                for i = 3:m
                   D(i, i-2:i+1) = [1/24 -9/8 9/8 -1/24];
                end

            case 6
                A = [-1627/1920  211/640  59/48  -235/192 91/128 -443/1920 31/960; ...
                        31/960  -687/640 129/128   19/192 -3/32    21/640  -3/640];
                D(2:3, 1:7) = A;
                D(m:m+1, m-5:end) = -rot90(A,2);
                for i = 4:m-1
                    D(i, i-3:i+2) = [-3/640 25/384 -75/64 75/64 -25/384 3/640];
                end

            case 8
                A = [-1423/1792     -491/7168   7753/3072 -18509/5120  3535/1024 -2279/1024  953/1024 -1637/7168  2689/107520; ...
                      2689/107520 -36527/35840  4259/5120   6497/15360 -475/1024  1541/5120 -639/5120  1087/35840  -59/17920; ...
                       -59/17920    1175/21504 -1165/1024   1135/1024    25/3072  -251/5120   25/1024   -45/7168     5/7168];
                D(2:4, 1:9) = A;
                D(m-1:m+1, m-7:end) = -rot90(A,2);
                for i = 5:m-2
                    D(i, i-4:i+3) = [5/7168 -49/5120 245/3072 -1225/1024 1225/1024 -245/3072 49/5120 -5/7168];
                end
        end
    end

    D = (1/dx).*D;
end

function [S, offS, C, offC] = mole_periodic_stencils(k)
% Interior staggered mimetic stencil S (with column offsets offS relative to
% the row/cell index) and the order-k nodal centered first-difference stencil C
% (with offsets offC relative to the node index) used on the two boundary nodes.
    switch k
        case 2
            S = [-1 1];                                       offS = [-1 0];
            C = [-1/2 0 1/2];                                 offC = [-1 0 1];
        case 4
            S = [1/24 -9/8 9/8 -1/24];                        offS = [-2 -1 0 1];
            C = [1/12 -2/3 0 2/3 -1/12];                      offC = [-2 -1 0 1 2];
        case 6
            S = [-3/640 25/384 -75/64 75/64 -25/384 3/640];   offS = [-3 -2 -1 0 1 2];
            C = [-1/60 3/20 -3/4 0 3/4 -3/20 1/60];           offC = [-3 -2 -1 0 1 2 3];
        case 8
            S = [5/7168 -49/5120 245/3072 -1225/1024 1225/1024 -245/3072 49/5120 -5/7168];
            offS = [-4 -3 -2 -1 0 1 2 3];
            C = [1/280 -4/105 1/5 -4/5 0 4/5 -1/5 4/105 -1/280];
            offC = [-4 -3 -2 -1 0 1 2 3 4];
    end
end

function c = mole_wrapcol(p, m)
% Wrap a column index into 1..m+1 with period m (node m+1 aliases node 1),
% leaving in-range indices untouched so interior rows keep their native columns.
    while p < 1,     p = p + m; end
    while p > (m+1), p = p - m; end
    c = p;
end
