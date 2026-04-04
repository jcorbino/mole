/**
 * @file curl.cpp
 * @brief Implementation of the Curl operator
 *
 */

#include "curl.h"

// 2D Curl
Curl::Curl(u16 k, u32 m, u32 n, Real dx, Real dy) {
  assert(k >= 2 && k % 2 == 0);
  assert(m >= 2 * k);
  assert(n >= 2 * k);
  assert(dx > 0 && dy > 0);

  // 1D mimetic gradient operators
  Gradient Gx(k, m, dx); // (m+1) x (m+2)
  Gradient Gy(k, n, dy); // (n+1) x (n+2)

  // Identity matrices at node dimensions
  sp_mat Im1 = speye(m + 1, m + 1);
  sp_mat In1 = speye(n + 1, n + 1);

  // 2D curl: C = [-kron(Gy, Im1),  kron(In1, Gx)]
  sp_mat dFx_dy = Utils::spkron((sp_mat)Gy, Im1);
  sp_mat dFy_dx = Utils::spkron(In1, (sp_mat)Gx);

  *this = Utils::spjoin_rows(-dFx_dy, dFy_dx);
}

// 3D Curl
Curl::Curl(u16 k, u32 m, u32 n, u32 o, Real dx, Real dy, Real dz) {
  assert(k >= 2 && k % 2 == 0);
  assert(m >= 2 * k);
  assert(n >= 2 * k);
  assert(o >= 2 * k);
  assert(dx > 0 && dy > 0 && dz > 0);

  // 1D mimetic gradient operators
  Gradient Gx(k, m, dx); // (m+1) x (m+2)
  Gradient Gy(k, n, dy); // (n+1) x (n+2)
  Gradient Gz(k, o, dz); // (o+1) x (o+2)

  // Identity matrices
  sp_mat Im1 = speye(m + 1, m + 1);
  sp_mat Im2 = speye(m + 2, m + 2);
  sp_mat In1 = speye(n + 1, n + 1);
  sp_mat In2 = speye(n + 2, n + 2);
  sp_mat Io1 = speye(o + 1, o + 1);
  sp_mat Io2 = speye(o + 2, o + 2);

  // Face-based input sizes
  u32 sFx = (m + 1) * (n + 2) * (o + 2);
  u32 sFy = (m + 2) * (n + 1) * (o + 2);
  u32 sFz = (m + 2) * (n + 2) * (o + 1);

  // Edge-based output sizes
  u32 sCx = (m + 2) * (n + 1) * (o + 1);
  u32 sCy = (m + 1) * (n + 2) * (o + 1);
  u32 sCz = (m + 1) * (n + 1) * (o + 2);

  // Row 1: Cx = dFz/dy - dFy/dz
  sp_mat dFz_dy = Utils::spkron(Utils::spkron(Io1, (sp_mat)Gy), Im2);
  sp_mat dFy_dz = Utils::spkron(Utils::spkron((sp_mat)Gz, In1), Im2);

  // Row 2: Cy = dFx/dz - dFz/dx
  sp_mat dFx_dz = Utils::spkron(Utils::spkron((sp_mat)Gz, In2), Im1);
  sp_mat dFz_dx = Utils::spkron(Utils::spkron(Io1, In2), (sp_mat)Gx);

  // Row 3: Cz = dFy/dx - dFx/dy
  sp_mat dFy_dx = Utils::spkron(Utils::spkron(Io2, In1), (sp_mat)Gx);
  sp_mat dFx_dy = Utils::spkron(Utils::spkron(Io2, (sp_mat)Gy), Im1);

  // Zero blocks
  sp_mat Z1(sCx, sFx);
  sp_mat Z2(sCy, sFy);
  sp_mat Z3(sCz, sFz);

  // Assemble 3x3 block curl matrix
  //        [   Fx       Fy        Fz    ]
  // Row 1: [  Z1     -dFy_dz    dFz_dy  ]  Cx = dFz/dy - dFy/dz
  // Row 2: [ dFx_dz    Z2      -dFz_dx  ]  Cy = dFx/dz - dFz/dx
  // Row 3: [-dFx_dy   dFy_dx     Z3     ]  Cz = dFy/dx - dFx/dy
  sp_mat row1 = Utils::spjoin_rows(Utils::spjoin_rows(Z1, -dFy_dz), dFz_dy);
  sp_mat row2 = Utils::spjoin_rows(Utils::spjoin_rows(dFx_dz, Z2), -dFz_dx);
  sp_mat row3 = Utils::spjoin_rows(Utils::spjoin_rows(-dFx_dy, dFy_dx), Z3);

  *this = Utils::spjoin_cols(Utils::spjoin_cols(row1, row2), row3);
}
