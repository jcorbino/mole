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
   */
  Laplacian(u16 k, u32 m, Real dx);
  
  /**
   * @brief 2-D Constructor
   *
   * @param k Order of accuracy
   * @param m Number of cells along x-axis
   * @param n Number of cells along y-axis
   * @param dx Spacing between cells along x-axis
   * @param dy Spacing between cells along y-axis
   */
  Laplacian(u16 k, u32 m, u32 n, Real dx, Real dy);
  
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
  Laplacian(u16 k, u32 m, u32 n, u32 o, Real dx, Real dy, Real dz);
};

#endif // LAPLACIAN_H
