#include "divergence.h"
#include <vector>

namespace {

// Interior (staggered) mimetic stencil S with its face offsets offS, relative
// to the cell row index, plus the order-k nodal centered first difference C
// (offsets offC) used on the two boundary node rows of the periodic operator.
void periodicStencils(u16 k, std::vector<Real> &S, std::vector<int> &offS,
                      std::vector<Real> &C, std::vector<int> &offC) {
  switch (k) {
  case 2:
    S = {-1.0, 1.0};
    offS = {-1, 0};
    C = {-1.0 / 2.0, 0.0, 1.0 / 2.0};
    offC = {-1, 0, 1};
    break;
  case 4:
    S = {1.0 / 24.0, -9.0 / 8.0, 9.0 / 8.0, -1.0 / 24.0};
    offS = {-2, -1, 0, 1};
    C = {1.0 / 12.0, -2.0 / 3.0, 0.0, 2.0 / 3.0, -1.0 / 12.0};
    offC = {-2, -1, 0, 1, 2};
    break;
  case 6:
    S = {-3.0 / 640.0, 25.0 / 384.0, -75.0 / 64.0,
         75.0 / 64.0,  -25.0 / 384.0, 3.0 / 640.0};
    offS = {-3, -2, -1, 0, 1, 2};
    C = {-1.0 / 60.0, 3.0 / 20.0, -3.0 / 4.0, 0.0,
         3.0 / 4.0,   -3.0 / 20.0, 1.0 / 60.0};
    offC = {-3, -2, -1, 0, 1, 2, 3};
    break;
  }
}

// Wrap a face index into [0, m] with period m (face m aliases face 0), leaving
// in-range indices untouched so interior rows keep their native columns.
u32 wrapFace(int p, u32 m) {
  const int mm = static_cast<int>(m);
  while (p < 0)
    p += mm;
  while (p > mm)
    p -= mm;
  return static_cast<u32>(p);
}

// The (m+2) x m matrix that lifts the m interior cell values of a transverse
// index into the (m+2)-long cell-space layout of the 2-D/3-D divergence output:
// the interior lands in rows 1..m, and rows 0 and m+1 are the two boundary
// nodes. Those two rows are empty in the default operator -- a BC operator
// overwrites them later -- but a periodic problem has no such operator, so they
// carry the order-k staggered midpoint interpolation of the cells straddling
// the seam. A plain [1/2 1/2] average would cap the whole boundary shell at
// second order however large k is, while the interior stayed order k.
sp_mat transverseLift(u32 m, u16 k, bool periodic) {
  sp_mat P = speye(m + 2, m + 2);
  P.shed_col(0);
  P.shed_col(m);

  if (!periodic)
    return P;

  std::vector<Real> w;
  switch (k) {
  case 2:
    w = {1.0 / 2.0, 1.0 / 2.0};
    break;
  case 4:
    w = {-1.0 / 16.0, 9.0 / 16.0, 9.0 / 16.0, -1.0 / 16.0};
    break;
  case 6:
    w = {3.0 / 256.0,   -25.0 / 256.0, 150.0 / 256.0,
         150.0 / 256.0, -25.0 / 256.0, 3.0 / 256.0};
    break;
  }

  // Cells m-p+1..m then 1..p in 1-based terms, straddling the seam. The 1-D
  // constructor already asserts m > 2k, so the stencil never wraps onto itself.
  const int mm = static_cast<int>(m);
  const int p = static_cast<int>(w.size()) / 2;
  for (int t = 0; t < 2 * p; t++) {
    const int off = t - p + 1;
    int col = (mm + off - 1) % mm;
    if (col < 0)
      col += mm;
    P(0, static_cast<u32>(col)) += w[static_cast<size_t>(t)];
    P(m + 1, static_cast<u32>(col)) += w[static_cast<size_t>(t)];
  }

  return P;
}

} // namespace

// 1-D Constructor
Divergence::Divergence(u16 k, u32 m, Real dx, bool periodic)
    : sp_mat(m + 2, m + 1) {
  assert(!(k % 2));
  assert(k > 1 && k < 7);
  assert(m > 2 * k);

  if (periodic) {
    std::vector<Real> S, C;
    std::vector<int> offS, offC;
    periodicStencils(k, S, offS, C, offC);

    // Interior cell rows: staggered stencil, wrapped across the seam.
    for (u32 i = 1; i < m + 1; i++)
      for (size_t t = 0; t < S.size(); t++)
        at(i, wrapFace(static_cast<int>(i) + offS[t], m)) += S[t];

    // Boundary node rows: west node = row 0, east node = row m+1. They come
    // out identical, as they must: both are the same physical point.
    for (size_t t = 0; t < C.size(); t++) {
      at(0, wrapFace(offC[t], m)) += C[t];
      at(m + 1, wrapFace(static_cast<int>(m) + offC[t], m)) += C[t];
    }

    // No one-sided closures, so the mimetic weights are uniform.
    Q = ones<vec>(2 * k + 1);

    // Scaling
    *this /= dx;
    return;
  }

  switch (k) {
  case 2:
    for (u32 i = 1; i < m + 1; i++) {
      at(i, i - 1) = -1.0;
      at(i, i) = 1.0;
    }
    // Weights
    Q = {1.0, 1.0, 1.0, 1.0, 1.0};
    break;
  case 4:
    // A
    at(1, 0) = -11.0 / 12.0;
    at(1, 1) = 17.0 / 24.0;
    at(1, 2) = 3.0 / 8.0;
    at(1, 3) = -5.0 / 24.0;
    at(1, 4) = 1.0 / 24.0;
    // A'
    at(m, m) = 11.0 / 12.0;
    at(m, m - 1) = -17.0 / 24.0;
    at(m, m - 2) = -3.0 / 8.0;
    at(m, m - 3) = 5.0 / 24.0;
    at(m, m - 4) = -1.0 / 24.0;
    // Middle
    for (u32 i = 2; i < m; i++) {
      at(i, i - 2) = 1.0 / 24.0;
      at(i, i - 1) = -9.0 / 8.0;
      at(i, i) = 9.0 / 8.0;
      at(i, i + 1) = -1.0 / 24.0;
    }
    // Weights
    Q = {2186.0 / 1943.0, 2125.0 / 2828.0, 1441.0 / 1240.0,
         648.0 / 673.0,   349.0 / 350.0,   648.0 / 673.0,
         1441.0 / 1240.0, 2125.0 / 2828.0, 2186.0 / 1943.0};
    break;
  case 6:
    // A
    at(1, 0) = -1627.0 / 1920.0;
    at(1, 1) = 211.0 / 640.0;
    at(1, 2) = 59.0 / 48.0;
    at(1, 3) = -235.0 / 192.0;
    at(1, 4) = 91.0 / 128.0;
    at(1, 5) = -443.0 / 1920.0;
    at(1, 6) = 31.0 / 960.0;
    at(2, 0) = 31.0 / 960.0;
    at(2, 1) = -687.0 / 640.0;
    at(2, 2) = 129.0 / 128.0;
    at(2, 3) = 19.0 / 192.0;
    at(2, 4) = -3.0 / 32.0;
    at(2, 5) = 21.0 / 640.0;
    at(2, 6) = -3.0 / 640.0;
    // A'
    at(m, m) = 1627.0 / 1920.0;
    at(m, m - 1) = -211.0 / 640.0;
    at(m, m - 2) = -59.0 / 48.0;
    at(m, m - 3) = 235.0 / 192.0;
    at(m, m - 4) = -91.0 / 128.0;
    at(m, m - 5) = 443.0 / 1920.0;
    at(m, m - 6) = -31.0 / 960.0;
    at(m - 1, m) = -31.0 / 960.0;
    at(m - 1, m - 1) = 687.0 / 640.0;
    at(m - 1, m - 2) = -129.0 / 128.0;
    at(m - 1, m - 3) = -19.0 / 192.0;
    at(m - 1, m - 4) = 3.0 / 32.0;
    at(m - 1, m - 5) = -21.0 / 640.0;
    at(m - 1, m - 6) = 3.0 / 640.0;
    // Middle
    for (u32 i = 3; i < m - 1; i++) {
      at(i, i - 3) = -3.0 / 640.0;
      at(i, i - 2) = 25.0 / 384.0;
      at(i, i - 1) = -75.0 / 64.0;
      at(i, i) = 75.0 / 64.0;
      at(i, i + 1) = -25.0 / 384.0;
      at(i, i + 2) = 3.0 / 640.0;
    }
    // Weights
    Q = {2383.0 / 2005.0, 929.0 / 2002.0,  887.0 / 531.0,   3124.0 / 5901.0,
         1706.0 / 1457.0, 457.0 / 467.0,   1057.0 / 1061.0, 457.0 / 467.0,
         1706.0 / 1457.0, 3124.0 / 5901.0, 887.0 / 531.0,   929.0 / 2002.0,
         2383.0 / 2005.0};
    break;
  }

  // Scaling
  *this /= dx;
}

// 2-D Constructor
Divergence::Divergence(u16 k, u32 m, u32 n, Real dx, Real dy, bool periodic) {
  Divergence Dx(k, m, dx, periodic);
  Divergence Dy(k, n, dy, periodic);

  sp_mat Im = transverseLift(m, k, periodic);
  sp_mat In = transverseLift(n, k, periodic);

  sp_mat D1 = Utils::spkron(In, Dx);
  sp_mat D2 = Utils::spkron(Dy, Im);

  // Dimensions = (m+2)*(n+2), 2*m*n+m+n
  if (m != n)
    *this = Utils::spjoin_rows(D1, D2);
  else {
    sp_mat A1(1, 2);
    sp_mat A2(1, 2);
    A1(0, 0) = A2(0, 1) = 1.0;
    *this = Utils::spkron(A1, D1) + Utils::spkron(A2, D2);
  }
}

// 3-D Constructor
Divergence::Divergence(u16 k, u32 m, u32 n, u32 o, Real dx, Real dy, Real dz,
                       bool periodic) {
  Divergence Dx(k, m, dx, periodic);
  Divergence Dy(k, n, dy, periodic);
  Divergence Dz(k, o, dz, periodic);

  sp_mat Im = transverseLift(m, k, periodic);
  sp_mat In = transverseLift(n, k, periodic);
  sp_mat Io = transverseLift(o, k, periodic);

  sp_mat D1 = Utils::spkron(Utils::spkron(Io, In), Dx);
  sp_mat D2 = Utils::spkron(Utils::spkron(Io, Dy), Im);
  sp_mat D3 = Utils::spkron(Utils::spkron(Dz, In), Im);

  // Dimensions = (m+2)*(n+2)*(o+2), 3*m*n*o+m*n+m*o+n*o
  if ((m != n) || (n != o))
    *this = Utils::spjoin_rows(Utils::spjoin_rows(D1, D2), D3);
  else {
    sp_mat A1(1, 3);
    sp_mat A2(1, 3);
    sp_mat A3(1, 3);
    A1(0, 0) = A2(0, 1) = A3(0, 2) = 1.0;
    *this =
        Utils::spkron(A1, D1) + Utils::spkron(A2, D2) + Utils::spkron(A3, D3);
  }
}

// Returns weights
vec Divergence::getQ() { return Q; }
