/**
 * @file gradient.h
 * @author Johnny Corbino (johnnycorbino@gmail.com)
 * @brief Header file for the Gradient operator
 * @version 0.1
 * @date 2024-05-25
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef GRADIENT_H
#define GRADIENT_H

#include "utils.h"
#include <cassert>

/**
 * @brief Mimetic Gradient operator
 *
 */
class Gradient : public sp_mat {

public:
  using sp_mat::operator=;

  /**
   * @brief 1-D Constructor
   *
   * @param k Order of accuracy
   * @param m Number of cells
   * @param dx Spacing between cells
   * @param periodic Pass true to build the PERIODIC operator, in which every
   *        face row uses the interior (staggered) stencil wrapped around the
   *        domain, replacing the one-sided boundary closures. Defaults to
   *        false, giving the standard operator. In the periodic operator the
   *        first and last face rows come out identical -- west and east are
   *        the same physical point -- and columns 0 and m+1 are left empty,
   *        since a periodic problem carries no independent boundary unknowns.
   */
  Gradient(u16 k, u32 m, Real dx, bool periodic = false);
  
  /**
   * @brief 2-D Constructor
   *
   * @param k Order of accuracy
   * @param m Number of cells along x-axis
   * @param n Number of cells along y-axis
   * @param dx Spacing between cells along x-axis
   * @param dy Spacing between cells along y-axis
   * @param periodic Pass true to build the PERIODIC operator, wrapping the
   *        interior staggered stencil around the domain on every axis instead
   *        of using one-sided boundary closures. G then ignores the boundary
   *        entries of the scalar field -- those columns come out empty --
   *        since a periodic problem carries no independent boundary unknowns.
   *        Defaults to false. Pair it with Divergence(k, m, n, dx, dy, true).
   */
  Gradient(u16 k, u32 m, u32 n, Real dx, Real dy, bool periodic = false);
  
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
   * @param periodic Pass true to build the PERIODIC operator; see the 2-D
   *        constructor. Defaults to false.
   */
  Gradient(u16 k, u32 m, u32 n, u32 o, Real dx, Real dy, Real dz,
           bool periodic = false);
  
  /**
   * @brief P Vector for the weighted inner product
   * 
   * @return vec Weights depend on the order of accuracy
   */
  vec getP();

private:
  vec P;
};

#endif // GRADIENT_H
