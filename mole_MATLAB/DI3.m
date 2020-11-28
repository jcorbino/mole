function I = DI3(m, n, o, type)
    if strcmp(type, 'Dn')
        e = ones(m, 1);
        bdry = spdiags([-0.25*e -0.25*e 0.25*e 0.25*e], [0 1 m+1 m+2], m, (m+1)*n);
        block = spdiags([[0.25*e; 0; 0] [0.25*e; 0.25; 0]], [0 1], m+2, m+1);
        pattern = spdiags([-ones(n-2, 1) ones(n-2, 1)], [0 2], n-2, n);
        middle = kron(pattern, block);
        I = [spalloc(m+3, (m+1)*n, 0); bdry; spalloc(2, (m+1)*n, 0); middle; circshift(bdry, (m+1)*(n-2), 2); spalloc(m+3, (m+1)*n, 0)];
        I = kron(speye(o), I);
        I = [spalloc((m+2)*(n+2), size(I, 2), 0); I; spalloc((m+2)*(n+2), size(I, 2), 0)];
    elseif strcmp(type, 'De')
        e = ones(m-2, 1);
        block = spdiags([[-0.25*e; -0.25; 0] [0; 0.25*e; 0.25]], [-1 1], m+2, m);
        block(1) = -0.25;
        block(m, m) = 0.25;
        pattern = spdiags([ones(n, 1) ones(n, 1)], [0 1], n, n+1);
        middle = kron(pattern, block);
        I = [spalloc(m+3, (n+1)*m, 0); middle; spalloc(m+1, (n+1)*m, 0)];
        I = kron(speye(o), I);
        I = [spalloc((m+2)*(n+2), size(I, 2), 0); I; spalloc((m+2)*(n+2), size(I, 2), 0)];
    elseif strcmp(type, 'Dc')
    end
end