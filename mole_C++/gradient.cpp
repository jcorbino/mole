#include "gradient.h"
#include <vector>

namespace {

// Interior (staggered) mimetic stencil S, with cell offsets offS relative to
// the face (row) index. Face i lies between cells i-1 and i, so offset 0 is
// cell i. These are the same coefficients the standard operator uses on its
// interior rows, just expressed as offsets so they can be wrapped.
void periodicGradStencil(u16 k, std::vector<Real> &S, std::vector<int> &offS) {
  switch (k) {
  case 2:
    S = {-1.0, 1.0};
    offS = {-1, 0};
    break;
  case 4:
    S = {1.0 / 24.0, -9.0 / 8.0, 9.0 / 8.0, -1.0 / 24.0};
    offS = {-2, -1, 0, 1};
    break;
  case 6:
    S = {-3.0 / 640.0, 25.0 / 384.0, -75.0 / 64.0,
         75.0 / 64.0,  -25.0 / 384.0, 3.0 / 640.0};
    offS = {-3, -2, -1, 0, 1, 2};
    break;
  case 8:
    S = {5.0 / 7168.0,    -49.0 / 5120.0,  245.0 / 3072.0, -1225.0 / 1024.0,
         1225.0 / 1024.0, -245.0 / 3072.0, 49.0 / 5120.0,  -5.0 / 7168.0};
    offS = {-4, -3, -2, -1, 0, 1, 2, 3};
    break;
  }
}

// Wrap a cell index into [0, m-1] with period m.
u32 wrapCell(int p, u32 m) {
  const int mm = static_cast<int>(m);
  while (p < 0)
    p += mm;
  while (p >= mm)
    p -= mm;
  return static_cast<u32>(p);
}

} // namespace

// 1-D Constructor
Gradient::Gradient(u16 k, u32 m, Real dx, bool periodic)
    : sp_mat(m + 1, m + 2) {
  assert(!(k % 2));
  assert(k > 1 && k < 9);
  assert(m >= 2 * k);

  if (periodic) {
    std::vector<Real> S;
    std::vector<int> offS;
    periodicGradStencil(k, S, offS);

    // Every face row uses the interior stencil, wrapped across the seam. The
    // first and last rows come out identical (same physical point), and
    // columns 0 and m+1 are never touched: a periodic problem carries no
    // independent boundary unknowns.
    for (u32 i = 0; i < m + 1; i++)
      for (size_t t = 0; t < S.size(); t++)
        at(i, wrapCell(static_cast<int>(i) + offS[t], m) + 1) += S[t];

    // No one-sided closures, so the mimetic weights are uniform.
    P = ones<vec>(2 * k + 1);

    // Scaling
    *this /= dx;
    return;
  }

  switch (k) {
  case 2:
    // A
    at(0, 0) = -8.0 / 3.0;
    at(0, 1) = 3.0;
    at(0, 2) = -1.0 / 3.0;
    // A'
    at(m, m + 1) = 8.0 / 3.0;
    at(m, m) = -3.0;
    at(m, m - 1) = 1.0 / 3.0;
    // Middle
    for (u32 i = 1; i < m; i++) {
      at(i, i) = -1.0;
      at(i, i + 1) = 1.0;
    }
    // Weights
    P = {3.0 / 8.0, 9.0 / 8.0, 1.0, 9.0 / 8.0, 3.0 / 8.0};
    break;
  case 4:
    // A
    at(0, 0) = -352.0 / 105.0;
    at(0, 1) = 35.0 / 8.0;
    at(0, 2) = -35.0 / 24.0;
    at(0, 3) = 21.0 / 40.0;
    at(0, 4) = -5.0 / 56.0;
    at(1, 0) = 16.0 / 105.0;
    at(1, 1) = -31.0 / 24.0;
    at(1, 2) = 29.0 / 24.0;
    at(1, 3) = -3.0 / 40.0;
    at(1, 4) = 1.0 / 168.0;
    // A'
    at(m, m + 1) = 352.0 / 105.0;
    at(m, m) = -35.0 / 8.0;
    at(m, m - 1) = 35.0 / 24.0;
    at(m, m - 2) = -21.0 / 40.0;
    at(m, m - 3) = 5.0 / 56.0;
    at(m - 1, m + 1) = -16.0 / 105.0;
    at(m - 1, m) = 31.0 / 24.0;
    at(m - 1, m - 1) = -29.0 / 24.0;
    at(m - 1, m - 2) = 3.0 / 40.0;
    at(m - 1, m - 3) = -1.0 / 168.0;
    // Middle
    for (u32 i = 2; i < m - 1; i++) {
      at(i, i - 1) = 1.0 / 24.0;
      at(i, i) = -9.0 / 8.0;
      at(i, i + 1) = 9.0 / 8.0;
      at(i, i + 2) = -1.0 / 24.0;
    }
    // Weights
    P = {1606.0 / 4535.0, 941.0 / 766.0, 1384.0 / 1541.0,
         1371.0 / 1346.0, 701.0 / 700.0, 1371.0 / 1346.0,
         1384.0 / 1541.0, 941.0 / 766.0, 1606.0 / 4535.0};
    break;
  case 6:
    // A
    at(0, 0) = -13016.0 / 3465.0;
    at(0, 1) = 693.0 / 128.0;
    at(0, 2) = -385.0 / 128.0;
    at(0, 3) = 693.0 / 320.0;
    at(0, 4) = -495.0 / 448.0;
    at(0, 5) = 385.0 / 1152.0;
    at(0, 6) = -63.0 / 1408.0;
    at(1, 0) = 496.0 / 3465.0;
    at(1, 1) = -811.0 / 640.0;
    at(1, 2) = 449.0 / 384.0;
    at(1, 3) = -29.0 / 960.0;
    at(1, 4) = -11.0 / 448.0;
    at(1, 5) = 13.0 / 1152.0;
    at(1, 6) = -37.0 / 21120.0;
    at(2, 0) = -8.0 / 385.0;
    at(2, 1) = 179.0 / 1920.0;
    at(2, 2) = -153.0 / 128.0;
    at(2, 3) = 381.0 / 320.0;
    at(2, 4) = -101.0 / 1344.0;
    at(2, 5) = 1.0 / 128.0;
    at(2, 6) = -3.0 / 7040.0;
    // A'
    at(m, m + 1) = 13016.0 / 3465.0;
    at(m, m) = -693.0 / 128.0;
    at(m, m - 1) = 385.0 / 128.0;
    at(m, m - 2) = -693.0 / 320.0;
    at(m, m - 3) = 495.0 / 448.0;
    at(m, m - 4) = -385.0 / 1152.0;
    at(m, m - 5) = 63.0 / 1408.0;
    at(m - 1, m + 1) = -496.0 / 3465.0;
    at(m - 1, m) = 811.0 / 640.0;
    at(m - 1, m - 1) = -449.0 / 384.0;
    at(m - 1, m - 2) = 29.0 / 960.0;
    at(m - 1, m - 3) = 11.0 / 448.0;
    at(m - 1, m - 4) = -13.0 / 1152.0;
    at(m - 1, m - 5) = 37.0 / 21120.0;
    at(m - 2, m + 1) = 8.0 / 385.0;
    at(m - 2, m) = -179.0 / 1920.0;
    at(m - 2, m - 1) = 153.0 / 128.0;
    at(m - 2, m - 2) = -381.0 / 320.0;
    at(m - 2, m - 3) = 101.0 / 1344.0;
    at(m - 2, m - 4) = -1.0 / 128.0;
    at(m - 2, m - 5) = 3.0 / 7040.0;
    // Middle
    for (u32 i = 3; i < m - 2; i++) {
      at(i, i - 2) = -3.0 / 640.0;
      at(i, i - 1) = 25.0 / 384.0;
      at(i, i) = -75.0 / 64.0;
      at(i, i + 1) = 75.0 / 64.0;
      at(i, i + 2) = -25.0 / 384.0;
      at(i, i + 3) = 3.0 / 640.0;
    }
    // Weights
    P = {420249.0 / 1331069.0,  2590978.0 / 1863105.0, 882762.0 / 1402249.0,
         1677712.0 / 1359311.0, 239985.0 / 261097.0,   664189.0 / 657734.0,
         756049.0 / 754729.0,   664189.0 / 657734.0,   239985.0 / 261097.0,
         1677712.0 / 1359311.0, 882762.0 / 1402249.0,  2590978.0 / 1863105.0,
         420249.0 / 1331069.0};
    break;
  case 8:
    // A
    at(0, 0) = -182144.0 / 45045.0;
    at(0, 1) = 6435.0 / 1024.0;
    at(0, 2) = -5005.0 / 1024.0;
    at(0, 3) = 27027.0 / 5120.0;
    at(0, 4) = -32175.0 / 7168.0;
    at(0, 5) = 25025.0 / 9216.0;
    at(0, 6) = -12285.0 / 11264.0;
    at(0, 7) = 3465.0 / 13312.0;
    at(0, 8) = -143.0 / 5120.0;
    at(1, 0) = 86048.0 / 675675.0;
    at(1, 1) = -131093.0 / 107520.0;
    at(1, 2) = 49087.0 / 46080.0;
    at(1, 3) = 10973.0 / 76800.0;
    at(1, 4) = -4597.0 / 21504.0;
    at(1, 5) = 4019.0 / 27648.0;
    at(1, 6) = -10331.0 / 168960.0;
    at(1, 7) = 2983.0 / 199680.0;
    at(1, 8) = -2621.0 / 1612800.0;
    at(2, 0) = -3776.0 / 225225.0;
    at(2, 1) = 8707.0 / 107520.0;
    at(2, 2) = -17947.0 / 15360.0;
    at(2, 3) = 29319.0 / 25600.0;
    at(2, 4) = -533.0 / 21504.0;
    at(2, 5) = -263.0 / 9216.0;
    at(2, 6) = 903.0 / 56320.0;
    at(2, 7) = -283.0 / 66560.0;
    at(2, 8) = 257.0 / 537600.0;
    at(3, 0) = 32.0 / 9009.0;
    at(3, 1) = -543.0 / 35840.0;
    at(3, 2) = 265.0 / 3072.0;
    at(3, 3) = -1233.0 / 1024.0;
    at(3, 4) = 8625.0 / 7168.0;
    at(3, 5) = -775.0 / 9216.0;
    at(3, 6) = 639.0 / 56320.0;
    at(3, 7) = -15.0 / 13312.0;
    at(3, 8) = 1.0 / 21504.0;
    // A'
    at(m, m + 1) = 182144.0 / 45045.0;
    at(m, m) = -6435.0 / 1024.0;
    at(m, m - 1) = 5005.0 / 1024.0;
    at(m, m - 2) = -27027.0 / 5120.0;
    at(m, m - 3) = 32175.0 / 7168.0;
    at(m, m - 4) = -25025.0 / 9216.0;
    at(m, m - 5) = 12285.0 / 11264.0;
    at(m, m - 6) = -3465.0 / 13312.0;
    at(m, m - 7) = 143.0 / 5120.0;
    at(m - 1, m + 1) = -86048.0 / 675675.0;
    at(m - 1, m) = 131093.0 / 107520.0;
    at(m - 1, m - 1) = -49087.0 / 46080.0;
    at(m - 1, m - 2) = -10973.0 / 76800.0;
    at(m - 1, m - 3) = 4597.0 / 21504.0;
    at(m - 1, m - 4) = -4019.0 / 27648.0;
    at(m - 1, m - 5) = 10331.0 / 168960.0;
    at(m - 1, m - 6) = -2983.0 / 199680.0;
    at(m - 1, m - 7) = 2621.0 / 1612800.0;
    at(m - 2, m + 1) = 3776.0 / 225225.0;
    at(m - 2, m) = -8707.0 / 107520.0;
    at(m - 2, m - 1) = 17947.0 / 15360.0;
    at(m - 2, m - 2) = -29319.0 / 25600.0;
    at(m - 2, m - 3) = 533.0 / 21504.0;
    at(m - 2, m - 4) = 263.0 / 9216.0;
    at(m - 2, m - 5) = -903.0 / 56320.0;
    at(m - 2, m - 6) = 283.0 / 66560.0;
    at(m - 2, m - 7) = -257.0 / 537600.0;
    at(m - 3, m + 1) = -32.0 / 9009.0;
    at(m - 3, m) = 543.0 / 35840.0;
    at(m - 3, m - 1) = -265.0 / 3072.0;
    at(m - 3, m - 2) = 1233.0 / 1024.0;
    at(m - 3, m - 3) = -8625.0 / 7168.0;
    at(m - 3, m - 4) = 775.0 / 9216.0;
    at(m - 3, m - 5) = -639.0 / 56320.0;
    at(m - 3, m - 6) = 15.0 / 13312.0;
    at(m - 3, m - 7) = -1.0 / 21504.0;
    // Middle
    for (u32 i = 4; i < m - 3; i++) {
      at(i, i - 3) = 5.0 / 7168.0;
      at(i, i - 2) = -49.0 / 5120.0;
      at(i, i - 1) = 245.0 / 3072.0;
      at(i, i) = -1225.0 / 1024.0;
      at(i, i + 1) = 1225.0 / 1024.0;
      at(i, i + 2) = -245.0 / 3072.0;
      at(i, i + 3) = 49.0 / 5120.0;
      at(i, i + 4) = -5.0 / 7168.0;
    }
    // Weights
    P = {267425.0 / 904736.0,   2307435.0 / 1517812.0, 847667.0 / 3066027.0,
         4050911.0 / 2301238.0, 498943.0 / 1084999.0,  211042.0 / 170117.0,
         2065895.0 / 2191686.0, 1262499.0 / 1258052.0, 1314891.0 / 1312727.0,
         1262499.0 / 1258052.0, 2065895.0 / 2191686.0, 211042.0 / 170117.0,
         498943.0 / 1084999.0,  4050911.0 / 2301238.0, 847667.0 / 3066027.0,
         2307435.0 / 1517812.0, 267425.0 / 904736.0};
    break;
  }

  // Scaling
  *this /= dx;
}

// 2-D Constructor
Gradient::Gradient(u16 k, u32 m, u32 n, Real dx, Real dy, bool periodic) {
  Gradient Gx(k, m, dx, periodic);
  Gradient Gy(k, n, dy, periodic);

  sp_mat Im = speye(m + 2, m + 2);
  sp_mat In = speye(n + 2, n + 2);

  Im.shed_row(0);
  Im.shed_row(m);
  In.shed_row(0);
  In.shed_row(n);

  sp_mat G1 = Utils::spkron(In, Gx);
  sp_mat G2 = Utils::spkron(Gy, Im);

  // Dimensions = 2*m*n+m+n, (m+2)*(n+2)
  if (m != n)
    *this = Utils::spjoin_cols(G1, G2);
  else {
    sp_mat A1(2, 1);
    sp_mat A2(2, 1);
    A1(0, 0) = A2(1, 0) = 1.0;
    *this = Utils::spkron(A1, G1) + Utils::spkron(A2, G2);
  }
}

// 3-D Constructor
Gradient::Gradient(u16 k, u32 m, u32 n, u32 o, Real dx, Real dy, Real dz,
                   bool periodic) {
  Gradient Gx(k, m, dx, periodic);
  Gradient Gy(k, n, dy, periodic);
  Gradient Gz(k, o, dz, periodic);

  sp_mat Im = speye(m + 2, m + 2);
  sp_mat In = speye(n + 2, n + 2);
  sp_mat Io = speye(o + 2, o + 2);

  Im.shed_row(0);
  Im.shed_row(m);
  In.shed_row(0);
  In.shed_row(n);
  Io.shed_row(0);
  Io.shed_row(o);

  sp_mat G1 = Utils::spkron(Utils::spkron(Io, In), Gx);
  sp_mat G2 = Utils::spkron(Utils::spkron(Io, Gy), Im);
  sp_mat G3 = Utils::spkron(Utils::spkron(Gz, In), Im);

  // Dimensions = 3*m*n*o+m*n+m*o+n*o, (m+2)*(n+2)*(o+2)
  if ((m != n) || (n != o))
    *this = Utils::spjoin_cols(Utils::spjoin_cols(G1, G2), G3);
  else {
    sp_mat A1(3, 1);
    sp_mat A2(3, 1);
    sp_mat A3(3, 1);
    A1(0, 0) = A2(1, 0) = A3(2, 0) = 1.0;
    *this =
        Utils::spkron(A1, G1) + Utils::spkron(A2, G2) + Utils::spkron(A3, G3);
  }
}

// Returns weights
vec Gradient::getP() { return P; }
