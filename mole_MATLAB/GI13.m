function I = GI13(M, m, n, o, type)
    if strcmp(type, 'Gn')
        I = speye(n*o);
        I1 = speye(m+1, m);
        I1(end, end) = 1;
        I = kron(I, I1);
        I = [I sparse(size(I, 1), m*o)];
        I = I*M;
    elseif strcmp(type, 'Ge')
        I = speye(n+1, n);
        I(end, end) = 1;
        I1 = speye(m, m+1);
        I = kron(I, I1);
        I = kron(speye(o), I);
        I = I*M;
    elseif strcmp(type, 'Gc')
        I = speye(n*o);
        I1 = speye(m+1, m);
        I1(end, end) = 1;
        I = kron(I, I1);
        I = [I sparse(size(I, 1), m*n)];
        I = I*M;
    elseif strcmp(type, 'Gcy')
        I = speye(m*o);
        I1 = speye(n+1, n);
        I1(end, end) = 1;
        I = kron(I, I1);
        I = [I sparse(size(I, 1), m*n)];
        I = I*M;
    elseif strcmp(type, 'Gee')
        I = speye(o+1, o);
        I(end, end) = 1;
        I1 = speye(m, m+1);
        I = kron(I, I1);
        I = kron(speye(n), I);
        I = I*M;
    elseif strcmp(type, 'Gnn')
        I = speye(m*n);
        I1 = speye(o+1, o);
        I1(end, end) = 1;
        I = kron(I, I1);
        I = [I sparse(size(I, 1), m*o)];
        I = I*M;
    end
end