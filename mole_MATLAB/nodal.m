function N = nodal(k, m, dx, nu)
% Returns a m+1 by m+1 one-dimensional operator that approximates the
% nu-th order derivative on a uniform nodal grid
%
% Parameters:
%                k : Order of accuracy
%                m : Number of nodes
%               dx : Step size
%               nu : (Optional) Order of the derivative (default = 1)

    if nargin < 4
        nu = 1;
    end

    % Validate inputs
    assert(mod(k, 2) == 0, ...
        'k must be a positive even integer (got k=%d).', k);
    assert(nu >= 1 && mod(nu, 1) == 0, ...
        'nu must be a positive integer (got nu=%d).', nu);
    assert(k >= nu, ...
        'Derivative order nu=%d cannot exceed order of accuracy k=%d.', nu, k);
    assert(k + 1 <= m, ...
        'Not enough cells for desired order of accuracy k=%d and derivative order nu=%d (got %d cells).', ...
        k, nu, m-1);

    m = m-1;
    n_rows = m+1;
    n_cols = n_rows;

    N = sparse(n_rows, n_cols);

    neighbors = zeros(1, k+1); % Bandwidth = k+1
    neighbors(1) = -k/2;
    len = k+1;
    j = 1;

    for i = 2 : len
        neighbors(i) = neighbors(i-1)+1;
    end

    % Create a k by k Vandermonde matrix based on the neighbors:
    A = vander(neighbors)';

    % nu-th order derivative
    b = zeros(len, 1);
    b(len-nu) = 1;

    % Solve the linear system to get the coefficients
    coeffs = A\b;

    for i = k/2+1 : n_rows-k/2
        N(i, j:j+len-1) = coeffs;
        j = j+1;
    end

    p = k/2;
    q = k+1;
    A = sparse(p, q);
    for i = 1 : p % For each row of A
        neighbors = zeros(1, q); % k+1 points are used for the boundaries
        neighbors(1) = 1-i; % Shifting the stencil to the right
        neighbors(2) = neighbors(1)+1;

        for j = 3 : q
            neighbors(j) = neighbors(j-1)+1;
        end

        V = vander(neighbors)';
        b = zeros(q, 1);
        b(q-nu) = 1;
        coeffs = V\b;
        A(i, 1:q) = coeffs;
    end

    % Insert A into N (upper-left corner of N)
    N(1:p, 1:q) = A;

    % Permutation matrices
    Pp = fliplr(speye(p));
    Pq = fliplr(speye(q));
    % Construct boundary operator for lower-right corner of N.
    % Sign depends on derivative order: odd → antisymmetric, even → symmetric
    A = (-1)^nu * Pp*A*Pq;

    % Insert A' into N
    N(n_rows-p+1:n_rows, n_cols-q+1:n_cols) = A;

    % Scale N: nu-th derivative requires factorial(nu)/dx^nu.
    % Effective order of accuracy for the nu-th derivative (centered symmetric stencil):
    %   nu odd  -> order = (k+1) - nu          (leading error term at power k+1 is odd; does not cancel)
    %   nu even -> order = (k+2) - nu          (leading error term at power k+1 is odd; cancels by symmetry)
    N = factorial(nu)/dx^nu * N;
end
