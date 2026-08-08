#include "laplacian.h"

// 1-D Constructor
Laplacian::Laplacian(u16 k, u32 m, Real dx, bool periodic) {
  Divergence div(k, m, dx, periodic);
  Gradient grad(k, m, dx, periodic);

  // Dimensions = m+2, m+2
  *this = (sp_mat)div * (sp_mat)grad;
}

// 2-D Constructor
Laplacian::Laplacian(u16 k, u32 m, u32 n, Real dx, Real dy, bool periodic) {
  Divergence div(k, m, n, dx, dy, periodic);
  Gradient grad(k, m, n, dx, dy, periodic);

  // Dimensions = (m+2)*(n+2), (m+2)*(n+2)
  *this = (sp_mat)div * (sp_mat)grad;
}

// 3-D Constructor
Laplacian::Laplacian(u16 k, u32 m, u32 n, u32 o, Real dx, Real dy, Real dz,
                     bool periodic) {
  Divergence div(k, m, n, o, dx, dy, dz, periodic);
  Gradient grad(k, m, n, o, dx, dy, dz, periodic);

  // Dimensions = (m+2)*(n+2)*(o+2), (m+2)*(n+2)*(o+2)
  *this = (sp_mat)div * (sp_mat)grad;
}
