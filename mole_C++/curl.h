/**
 * @file curl.h
 * @brief Header file for the Curl operator
 *
 */

#ifndef CURL_H
#define CURL_H

#include "gradient.h"
#include <cassert>

/**
 * @brief Mimetic Curl operator
 *
 * Constructs the curl matrix from Kronecker products of 1D mimetic gradients.
 *
 * 2D: Returns a matrix C such that C*F computes the scalar curl of a
 * face-staggered vector field F = [Fx; Fy].
 *   - Input:  (m+1)(n+2) + (m+2)(n+1) face DOFs
 *   - Output: (m+1)(n+1) nodal values
 *
 * 3D: Returns a matrix C such that C*F computes the vector curl of a
 * face-staggered vector field F = [Fx; Fy; Fz].
 *   - Input:  (m+1)(n+2)(o+2) + (m+2)(n+1)(o+2) + (m+2)(n+2)(o+1) face DOFs
 *   - Output: (m+2)(n+1)(o+1) + (m+1)(n+2)(o+1) + (m+1)(n+1)(o+2) edge DOFs
 *
 * Key property: C * G_full = 0 (curl of gradient is identically zero)
 */
class Curl : public sp_mat {

public:
  using sp_mat::operator=;

  /**
   * @brief 2-D Constructor
   *
   * @param k Order of accuracy
   * @param m Number of cells along x-axis
   * @param n Number of cells along y-axis
   * @param dx Spacing between cells along x-axis
   * @param dy Spacing between cells along y-axis
   */
  Curl(u16 k, u32 m, u32 n, Real dx, Real dy);

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
  Curl(u16 k, u32 m, u32 n, u32 o, Real dx, Real dy, Real dz);
};

#endif // CURL_H
