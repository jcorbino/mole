function I = I(M, m, n, type)
    if strcmp(type, 'Gn')
        I = speye(n);
        I1 = speye(m+1, m);
        I1(end, end) = 1;
        I = kron(I, I1);
        I = [I sparse(size(I, 1), m)];
        I = I*M;
    else
        I = speye(n+1, n);
        I(end, end) = 1;
        I1 = speye(m, m+1);
        I = kron(I, I1);
        I = I*M;
    end
end