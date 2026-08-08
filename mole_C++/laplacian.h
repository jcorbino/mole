/**
 * @file laplacian.h
 * @author Johnny Corbino (johnnycorbino@gmail.com)
 * @brief Header file for the Laplacian operator
 * @version 0.1
 * @date 2024-05-25
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef LAPLACIAN_H
#define LAPLACIAN_H

#include "divergence.h"
#include "gradient.h"

/**
 * @brief Mimetic Laplacian operator
 *
 */
class Laplacian : public sp_mat {

public:
  using sp_mat::operator=;

  /**
   * @brief 1-D Constructor
   *
   * @param k Order of accuracy
   * @param m Number of cells
   * @param dx Spacing between cells
   * @param periodic Pass true to build the PERIODIC Laplacian, forwarded to
   *        Divergence and Gradient: the interior staggered stencils wrap around
   *        the domain instead of using one-sided boundary closures. The default
   *        operator's boundary rows are placeholders meant to be overwritten by
   *        a BC operator such as RobinBC; a periodic problem has no such
   *        operator, so every row of the periodic Laplacian is a real
   *        Laplacian. Defaults to false.
   */
  Laplacian(u16 k, u32 m, Real dx, bool periodic = false);
  
  /**
   * @brief 2-D Constructor
   *
   * @param k Order of accuracy
   * @param m Number of cells along x-axis
   * @param n Number of cells along y-axis
   * @param dx Spacing between cells along x-axis
   * @param dy Spacing between cells along y-axis
   * @param periodic Pass true to build the PERIODIC Laplacian; see the 1-D
   *        constructor. Defaults to false.
   */
  Laplacian(u16 k, u32 m, u32 n, Real dx, Real dy, bool periodic = false);
  
  /**
   * @brief 3-D Constructor
   *
   * @param k Order of accuracy
   * @param m Number of cells along x-axis
   * @param n Number of cells along y-axis
   * @param o Number of cells along z-axis
   * @param dx Spacing between cells along x-axis
   * @param dy Spacing between cells along y-axis
   * @param dz Spacing between cells along z-axis
   * @param periodic Pass true to build the PERIODIC Laplacian; see the 1-D
   *        constructor. Defaults to false.
   */
  Laplacian(u16 k, u32 m, u32 n, u32 o, Real dx, Real dy, Real dz,
            bool periodic = false);
};

#endif // LAPLACIAN_H
