#ifndef OPERATORS_H
#define OPERATORS_H

#include "laplacian.h"
#include "interpol.h"

inline sp_mat operator*(const Divergence &div, const Gradient &grad)
{
    return (sp_mat)div*(sp_mat)grad;
}

inline vec operator*(const Divergence &div, const vec &v)
{
    return (sp_mat)div*v;
}

inline vec operator*(const Gradient &grad, const vec &v)
{
    return (sp_mat)grad*v;
}

inline vec operator*(const Laplacian &lap, const vec &v)
{
    return (sp_mat)lap*v;
}

inline vec operator*(const Interpol &I, const vec &v)
{
    return (sp_mat)I*v;
}

#endif // OPERATORS_H
