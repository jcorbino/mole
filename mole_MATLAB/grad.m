function G = grad(k, m, dx, bc)
% Returns a m+1 by m+2 one-dimensional mimetic gradient operator
%
% Parameters:
%                k : Order of accuracy
%                m : Number of cells
%               dx : Step size
%               bc : (optional) boundary-condition flag. Pass 'periodic'
%                    (or a nonzero/true value) to build the PERIODIC operator,
%                    in which every face row uses the interior (staggered)
%                    stencil wrapped around the domain, replacing the one-sided
%                    boundary closures. Omit (or pass 'none'/0) for the standard
%                    mimetic operator.
%
%                    In the periodic operator the first and last face rows come
%                    out identical -- x = west and x = east are the same point --
%                    and columns 1 and m+2 are left empty, since a periodic
%                    problem carries no independent boundary unknowns.

    if nargin < 4
        bc = 'none';
    end

    % Assertions:
    assert(k >= 2, 'k >= 2');
    assert(k <= 8, 'k <= 8');
    assert(mod(k, 2) == 0, 'k % 2 = 0');
    assert(m >= 2*k, ['m >= ' num2str(2*k) ' for k = ' num2str(k)]);

    if ischar(bc) || isa(bc, 'string')
        assert(any(strcmpi(bc, {'periodic', 'none'})), ...
               'bc must be ''periodic'' or ''none''');
        is_periodic = strcmpi(bc, 'periodic');
    else
        is_periodic = ~isempty(bc) && all(logical(bc));
    end

    G = sparse(m+1, m+2);

    if is_periodic
        % --- Periodic operator -------------------------------------------
        % S : interior (staggered) mimetic stencil, offS its cell offsets
        % relative to the face index. Face i lies between cells i-1 and i, so
        % offset 0 is cell i. Every face row uses S, wrapped across the seam.
        [S, offS] = mole_periodic_grad_stencil(k);

        for i = 1:m+1
            for t = 1:numel(S)
                % Column j holds cell j-1, hence the +1 on the cell index.
                c = mole_wrapcell(i + offS(t), m) + 1;
                G(i, c) = G(i, c) + S(t);
            end
        end

    else
        % --- Standard mimetic operator (one-sided boundary closures) -----
        switch k
            case 2
                A = [-8/3 3 -1/3];
                G(1, 1:3) = A;
                G(end, end-2:end) = -fliplr(A);
                for i = 2:m
                   G(i, i:i+1) = [-1 1];
                end

            case 4
                A = [-352/105  35/8  -35/24 21/40 -5/56; ...
                       16/105 -31/24  29/24 -3/40  1/168];
                G(1:2, 1:5) = A;
                G(m:m+1, m-2:end) = -rot90(A,2);
                for i = 3:m-1
                   G(i, i-1:i+2) = [1/24 -9/8 9/8 -1/24];
                end

            case 6
                A = [-13016/3465  693/128  -385/128 693/320 -495/448  385/1152 -63/1408; ...
                        496/3465 -811/640   449/384 -29/960  -11/448   13/1152 -37/21120; ...
                         -8/385   179/1920 -153/128 381/320 -101/1344   1/128   -3/7040];
                G(1:3, 1:7) = A;
                G(m-1:m+1, m-4:end) = -rot90(A,2);
                for i = 4:m-2
                    G(i, i-2:i+3) = [-3/640 25/384 -75/64 75/64 -25/384 3/640];
                end

            case 8
                A = [-182144/45045     6435/1024    -5005/1024  27027/5120  -32175/7168  25025/9216  -12285/11264  3465/13312   -143/5120; ...
                       86048/675675 -131093/107520  49087/46080 10973/76800  -4597/21504  4019/27648 -10331/168960 2983/199680 -2621/1612800; ...
                       -3776/225225    8707/107520 -17947/15360 29319/25600   -533/21504  -263/9216     903/56320  -283/66560    257/537600; ...
                          32/9009      -543/35840     265/3072  -1233/1024    8625/7168   -775/9216     639/56320  -15/13312       1/21504];
                G(1:4, 1:9) = A;
                G(m-2:m+1, m-6:end) = -rot90(A,2);
                for i = 5:m-3
                    G(i, i-3:i+4) = [5/7168 -49/5120 245/3072 -1225/1024 1225/1024 -245/3072 49/5120 -5/7168];
                end
        end
    end

    G = (1/dx).*G;
end

function [S, offS] = mole_periodic_grad_stencil(k)
% Interior staggered mimetic stencil S, with cell offsets offS relative to the
% face (row) index. These are the same coefficients the standard operator uses
% on its interior rows, just expressed as offsets so they can be wrapped.
    switch k
        case 2
            S = [-1 1];                                       offS = [-1 0];
        case 4
            S = [1/24 -9/8 9/8 -1/24];                        offS = [-2 -1 0 1];
        case 6
            S = [-3/640 25/384 -75/64 75/64 -25/384 3/640];   offS = [-3 -2 -1 0 1 2];
        case 8
            S = [5/7168 -49/5120 245/3072 -1225/1024 1225/1024 -245/3072 49/5120 -5/7168];
            offS = [-4 -3 -2 -1 0 1 2 3];
    end
end

function c = mole_wrapcell(p, m)
% Wrap a cell index into 1..m with period m.
    while p < 1, p = p + m; end
    while p > m, p = p - m; end
    c = p;
end
