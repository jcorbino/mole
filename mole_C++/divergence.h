/**
 * @file divergence.h
 * @author Johnny Corbino (johnnycorbino@gmail.com)
 * @brief Header file for the Divergence operator
 * @version 0.1
 * @date 2024-05-25
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef DIVERGENCE_H
#define DIVERGENCE_H

#include "utils.h"
#include <cassert>

/**
 * @brief Mimetic Divergence operator
 *
 */
class Divergence : public sp_mat {

public:
  using sp_mat::operator=;

  /**
   * @brief 1-D Constructor
   *
   * @param k Order of accuracy
   * @param m Number of cells
   * @param dx Spacing between cells
   * @param periodic Pass true to build the PERIODIC operator, in which every
   *        cell row uses the interior (staggered) stencil wrapped around the
   *        domain, and the two boundary node rows use the order-k nodal
   *        centered difference wrapped around the seam. Defaults to false,
   *        giving the standard operator with one-sided boundary closures.
   *        Pair it with Interpol(m, c, true).
   */
  Divergence(u16 k, u32 m, Real dx, bool periodic = false);

  /**
   * @brief 2-D Constructor
   *
   * @param k Order of accuracy
   * @param m Number of cells along x-axis
   * @param n Number of cells along y-axis
   * @param dx Spacing between cells along x-axis
   * @param dy Spacing between cells along y-axis
   */
  Divergence(u16 k, u32 m, u32 n, Real dx, Real dy);

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
   */
  Divergence(u16 k, u32 m, u32 n, u32 o, Real dx, Real dy, Real dz);

  /**
   * @brief Q Vector for the weighted inner product
   * 
   * @return vec Weights depend on the order of accuracy
   */
  vec getQ();

private:
  vec Q;
};

#endif // DIVERGENCE_H
