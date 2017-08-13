#ifndef INTERPOL_H
#define INTERPOL_H

#include <cassert>
#include "utils.h"

class Interpol : public sp_mat
{

public:
    using sp_mat::operator=;

    // 1-D Constructor
    Interpol(u32 m, double c);
    // 2-D Constructor
    Interpol(u32 m, u32 n, double c1, double c2);
    // 3-D Constructor
    Interpol(u32 m, u32 n, u32 o, double c1, double c2, double c3);
};

#endif // INTERPOL_H
