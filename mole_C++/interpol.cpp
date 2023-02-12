#include "interpol.h"

// 1-D Constructor
Interpol::Interpol(u32 m, Real c) : sp_mat(m+1, m+2)
{
    assert(m >= 4);
    assert(c >= 0 && c <= 1);

    at(0, 0) = 1;
    at(m, m+1) = 1;

    for (u32 i = 1; i < m; i++) {
        at(i, i)   = c;
        at(i, i+1) = 1-c;
    }
}

// 2-D Constructor
Interpol::Interpol(u32 m, u32 n, Real c1, Real c2)
{
    Interpol Ix(m, c1);
    Interpol Iy(n, c2);

    sp_mat Im = speye(m+2, m+2);
    sp_mat In = speye(n+2, n+2);

    Im.shed_row(0);
    Im.shed_row(m);
    In.shed_row(0);
    In.shed_row(n);

    sp_mat I1 = Utils::spkron(In, Ix);
    sp_mat I2 = Utils::spkron(Iy, Im);

    // Dimensions = 2*m*n+m+n, (m+2)*(n+2)
    if (m != n)
        *this = Utils::spjoin_cols(I1, I2);
    else {
        sp_mat A1(2, 1);
        sp_mat A2(2, 1);
        A1(0, 0) = A2(1, 0) = 1.0;
        *this = Utils::spkron(A1, I1) + Utils::spkron(A2, I2);
    }
}

// 3-D Constructor
Interpol::Interpol(u32 m, u32 n, u32 o, Real c1, Real c2, Real c3)
{
    Interpol Ix(m, c1);
    Interpol Iy(n, c2);
    Interpol Iz(o, c3);

    sp_mat Im = speye(m+2, m+2);
    sp_mat In = speye(n+2, n+2);
    sp_mat Io = speye(o+2, o+2);

    Im.shed_row(0);
    Im.shed_row(m);
    In.shed_row(0);
    In.shed_row(n);
    Io.shed_row(0);
    Io.shed_row(o);

    sp_mat I1 = Utils::spkron(Utils::spkron(Io, In), Ix);
    sp_mat I2 = Utils::spkron(Utils::spkron(Io, Iy), Im);
    sp_mat I3 = Utils::spkron(Utils::spkron(Iz, In), Im);

    // Dimensions = 3*m*n*o+m*n+m*o+n*o, (m+2)*(n+2)*(o+2)
    if ((m != n) || (n != o))
        *this = Utils::spjoin_cols(Utils::spjoin_cols(I1, I2), I3);
    else {
        sp_mat A1(3, 1);
        sp_mat A2(3, 1);
        sp_mat A3(3, 1);
        A1(0, 0) = A2(1, 0) = A3(2, 0) = 1.0;
        *this = Utils::spkron(A1, I1) + Utils::spkron(A2, I2) + Utils::spkron(A3, I3);
    }
}
