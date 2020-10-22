function I = I(M, m, n, type)
    if strcmp(type, 'Gn')
        I = spalloc((m+1)*n, size(M, 2), nnz(M));
        offset = m-1;
        j = 1;
        for i = 0 : n-1
            I(j:j+offset, :) = M(j-i:j-i+offset, :);
            I(j+offset+1, :) = I(j+offset, :);
            j = j+offset+2;
        end
    else
        I = spalloc((n+1)*m, size(M, 2), nnz(M));
        offset = m-1;
        j = 1;
        for i = 0 : n-1
            I(j:j+offset, :) = M(j+i:j+i+offset, :);
            j = j+offset+1;
        end
        I(j:j+offset, :) = I(j-m:j-1, :);
    end
end