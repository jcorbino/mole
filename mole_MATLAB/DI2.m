function I = DI2(m, n, type)
    if strcmp(type, 'Dn')
        padm3 = spalloc(m+3, (m+1)*n, 0);
        pad1 = spalloc(1, (m+1)*n, 0);
        e = ones(m, 1);
        top = [spdiags([-0.25*e -0.25*e], 0:1, m, m+1) spdiags([0.25*e 0.25*e], 0:1, m, (m+1)*(n-1))];
        
        I = [padm3; top; pad1];
        
        middle = [spdiags([-0.25*e -0.25*e], 0:1, m, 2*m+2) spdiags([0.25*e 0.25*e], 0:1, m, (m+1)*n - 2*m-2)];
        
        middle = kron(ones(n-2, 1), [pad1; middle; pad1]);
        
        j = m+4;
        for i = 1 : n-3
            % This line is too slow
            middle(j:j+m-1, :) = circshift(middle(j:j+m-1, :), i*(m+1), 2);
            j = j+m+2;
        end
        
        I = [I; middle; pad1; circshift(top, (m+1)*(n-2), 2); padm3];
    else
        I = spalloc((m+2)*(n+2), (n+1)*m, 0);
    end
end