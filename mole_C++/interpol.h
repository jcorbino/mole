/**
 * @file interpol.h
 * @author Johnny Corbino (johnnycorbino@gmail.com)
 * @brief Header file for the Interpolation operator
 * @version 0.1
 * @date 2024-05-25
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef INTERPOL_H
#define INTERPOL_H

#include "utils.h"
#include <cassert>

/**
 * @brief Interpolation operator
 *
 */
class Interpol : public sp_mat {

public:
  using sp_mat::operator=;

  /**
   * @brief 1-D Constructor
   *
   * @param m Number of cells
   * @param c Interpolation coefficient
   */
  Interpol(u32 m, Real c);
  
  /**
   * @brief 2-D Constructor
   *
   * @param m Number of cells along x-axis
   * @param n Number of cells along y-axis
   * @param c1 Interpolation coefficient along x-axis
   * @param c2 Interpolation coefficient along y-axis
   */
  Interpol(u32 m, u32 n, Real c1, Real c2);
  
  /**
   * @brief 3-D Constructor
   *
   * @param m Number of cells along x-axis
   * @param n Number of cells along y-axis
   * @param o Number of cells along z-axis
   * @param c1 Interpolation coefficient along x-axis
   * @param c2 Interpolation coefficient along y-axis
   * @param c3 Interpolation coefficient along z-axis
   */
  Interpol(u32 m, u32 n, u32 o, Real c1, Real c2, Real c3);
  
  // 1-D Constructor for second type
  Interpol(bool type, u32 m, Real c);
  // 2-D Constructor for second type
  Interpol(bool type, u32 m, u32 n, Real c1, Real c2);
  // 3-D Constructor for second type
  Interpol(bool type, u32 m, u32 n, u32 o, Real c1, Real c2, Real c3);
};

#endif // INTERPOL_H
