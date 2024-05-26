/**
 * @file robinbc.h
 * @author Johnny Corbino (johnnycorbino@gmail.com)
 * @brief Header file for the Boundary operator
 * @version 0.1
 * @date 2024-05-25
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef ROBINBC_H
#define ROBINBC_H

#include "gradient.h"

/**
 * @brief Mimetic Robin BC operator
 *
 */
class RobinBC : public sp_mat {

public:
  using sp_mat::operator=;

  /**
   * @brief 1-D Constructor
   *
   * @param k Order of accuracy
   * @param m Number of cells
   * @param dx Spacing between cells
   * @param a Dirichlet coefficient
   * @param b Neumann coefficient
   */
  RobinBC(u16 k, u32 m, Real dx, Real a, Real b);
  
  /**
   * @brief 2-D Constructor
   *
   * @param k Order of accuracy
   * @param m Number of cells along x-axis
   * @param dx Spacing between cells along x-axis
   * @param n Number of cells along y-axis
   * @param dy Spacing between cells along y-axis
   * @param a Dirichlet coefficient
   * @param b Neumann coefficient
   */
  RobinBC(u16 k, u32 m, Real dx, u32 n, Real dy, Real a, Real b);
  
  /**
   * @brief 3-D Constructor
   *
   * @param k Order of accuracy
   * @param m Number of cells along x-axis
   * @param dx Spacing between cells along x-axis
   * @param n Number of cells along y-axis
   * @param dy Spacing between cells along y-axis
   * @param o Number of cells along z-axis
   * @param dz Spacing between cells along z-axis
   * @param a Dirichlet coefficient
   * @param b Neumann coefficient
   */
  RobinBC(u16 k, u32 m, Real dx, u32 n, Real dy, u32 o, Real dz, Real a,
          Real b);
};

#endif // ROBINBC_H
