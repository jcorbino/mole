// clang-format off
/**
 * Two-level Adaptive Mesh Refinement (AMR) for 2D Poisson using MOLE
 *
 * Demonstrates how to combine Corbino-Castillo mimetic operators with
 * patch-based structured AMR (Berger-Oliger approach):
 *
 *   1. Solve -∇²u = f on a coarse uniform grid using MOLE
 *   2. Estimate the error via gradient magnitude
 *   3. Tag cells exceeding a threshold; compute their bounding box
 *   4. Create a fine patch (4:1 refinement) over the tagged region
 *   5. Set fine-patch Dirichlet BCs from bilinear interpolation of coarse soln
 *   6. Solve the fine patch independently with MOLE operators
 *   7. Output coarse + fine solutions for visualization
 *
 * Problem:  -∇²u = f  on [0,1]²,  u = 0 on ∂Ω
 *           f = localized Gaussian source at (0.75, 0.75)
 *
 * The sharp source creates multi-scale features: a smooth background that the
 * coarse grid resolves well, plus a narrow peak that demands local refinement.
 *
 * Visualize with gnuplot:
 *   set view equal xyz
 *   set pm3d border lc rgb "black" lw 0.5
 *   splot 'amr_coarse.dat' w pm3d title 'Coarse', 'amr_fine.dat' w pm3d title 'Fine patch'
 */
// clang-format on

#include "mole.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

using namespace arma;
using namespace std;

// Localized Gaussian source — drives refinement near (0.75, 0.75)
static double source_term(double x, double y) {
  const double x0 = 0.75, y0 = 0.75, sigma = 0.03;
  return 200.0 * exp(-((x - x0) * (x - x0) + (y - y0) * (y - y0)) /
                     (2 * sigma * sigma));
}

// Build staggered-grid coordinates: [boundary, cell-centers, boundary]
static vec build_coords(int m, double lo, double hi) {
  double dx = (hi - lo) / m;
  vec x(m + 2);
  x(0) = lo;
  for (int i = 1; i <= m; ++i)
    x(i) = lo + (i - 0.5) * dx;
  x(m + 1) = hi;
  return x;
}

// Bilinear interpolation of a grid function U(xc, yc) at point (x, y)
static double bilinear_interp(const mat &U, const vec &xc, const vec &yc,
                              double x, double y) {
  int nx = static_cast<int>(xc.n_elem);
  int ny = static_cast<int>(yc.n_elem);

  x = max(xc(0), min(x, xc(nx - 1)));
  y = max(yc(0), min(y, yc(ny - 1)));

  int i0 = 0, j0 = 0;
  for (int i = 0; i < nx - 1; ++i)
    if (xc(i) <= x && x <= xc(i + 1)) {
      i0 = i;
      break;
    }
  for (int j = 0; j < ny - 1; ++j)
    if (yc(j) <= y && y <= yc(j + 1)) {
      j0 = j;
      break;
    }

  double dx = xc(i0 + 1) - xc(i0);
  double dy = yc(j0 + 1) - yc(j0);
  double t = (dx > 0) ? (x - xc(i0)) / dx : 0.0;
  double s = (dy > 0) ? (y - yc(j0)) / dy : 0.0;

  return (1 - t) * (1 - s) * U(i0, j0) + t * (1 - s) * U(i0 + 1, j0) +
         (1 - t) * s * U(i0, j0 + 1) + t * s * U(i0 + 1, j0 + 1);
}

// Write solution on a structured grid to a gnuplot-compatible file
static void write_solution(const string &filename, const mat &U, const vec &x,
                           const vec &y) {
  ofstream fout(filename);
  for (int j = 0; j < static_cast<int>(y.n_elem); ++j) {
    for (int i = 0; i < static_cast<int>(x.n_elem); ++i)
      fout << x(i) << " " << y(j) << " " << U(i, j) << "\n";
    fout << "\n";
  }
  fout.close();
}

int main() {
  cout << "============================================\n"
       << " 2D Poisson with AMR using MOLE (mimetic)\n"
       << " -nabla^2 u = f,  u = 0 on boundary\n"
       << " f = sharp Gaussian at (0.75, 0.75)\n"
       << "============================================\n\n";

  int k = 2; // Mimetic operator order of accuracy

  // ==========================================================
  // LEVEL 0 (Coarse): uniform grid over [0,1]^2
  // ==========================================================
  int mc = 24, nc = 24;
  double dx_c = 1.0 / mc, dy_c = 1.0 / nc;
  vec xc = build_coords(mc, 0.0, 1.0);
  vec yc = build_coords(nc, 0.0, 1.0);

  Laplacian L_c(k, mc, nc, dx_c, dy_c);
  RobinBC BC_c(k, mc, dx_c, nc, dy_c, 1, 0); // Dirichlet
  L_c = L_c + BC_c;

  // The mimetic operator discretises nabla^2, and the problem is
  // -nabla^2 u = f, so the right-hand side is -f. Without the sign the solve
  // returns -u, whose maximum is the zero Dirichlet boundary rather than the
  // peak of the solution.
  mat rhs_c(mc + 2, nc + 2, fill::zeros);
  for (int i = 1; i <= mc; ++i)
    for (int j = 1; j <= nc; ++j)
      rhs_c(i, j) = -source_term(xc(i), yc(j));

#ifdef EIGEN
  vec sol_c = Utils::spsolve_eigen(L_c, vectorise(rhs_c));
#else
  vec sol_c = spsolve(L_c, vectorise(rhs_c));
#endif
  mat U_c = reshape(sol_c, mc + 2, nc + 2);

  cout << "[Level 0] Coarse grid: " << mc << " x " << nc << " cells ("
       << (mc + 2) * (nc + 2) << " unknowns)\n"
       << "          Peak value : " << U_c.max() << "\n";

  // ==========================================================
  // ERROR ESTIMATION: gradient magnitude at cell centers
  // ==========================================================
  mat grad_err(mc + 2, nc + 2, fill::zeros);
  for (int i = 1; i <= mc; ++i) {
    for (int j = 1; j <= nc; ++j) {
      double ux =
          (U_c(min(i + 1, mc + 1), j) - U_c(max(i - 1, 0), j)) / (2 * dx_c);
      double uy =
          (U_c(i, min(j + 1, nc + 1)) - U_c(i, max(j - 1, 0))) / (2 * dy_c);
      grad_err(i, j) = sqrt(ux * ux + uy * uy);
    }
  }
  double err_max = grad_err.max();
  double threshold = 0.25 * err_max;

  // ==========================================================
  // TAG CELLS: bounding box of cells exceeding the threshold
  // ==========================================================
  int i_lo = mc, i_hi = 1, j_lo = nc, j_hi = 1;
  int tagged = 0;
  for (int i = 1; i <= mc; ++i)
    for (int j = 1; j <= nc; ++j)
      if (grad_err(i, j) > threshold) {
        i_lo = min(i_lo, i);
        i_hi = max(i_hi, i);
        j_lo = min(j_lo, j);
        j_hi = max(j_hi, j);
        ++tagged;
      }

  if (tagged == 0) {
    cout << "\nNo cells flagged — skipping refinement.\n";
    write_solution("amr_coarse.dat", U_c, xc, yc);
    return 0;
  }

  // Buffer zone: 2 coarse cells around the tagged region
  const int buf = 2;
  i_lo = max(1, i_lo - buf);
  i_hi = min(mc, i_hi + buf);
  j_lo = max(1, j_lo - buf);
  j_hi = min(nc, j_hi + buf);

  double xlo = xc(i_lo) - 0.5 * dx_c;
  double xhi = xc(i_hi) + 0.5 * dx_c;
  double ylo = yc(j_lo) - 0.5 * dy_c;
  double yhi = yc(j_hi) + 0.5 * dy_c;

  cout << "\n[Tagging] " << tagged
       << " cells flagged (threshold = " << threshold << ")\n"
       << "          Patch region: [" << xlo << ", " << xhi << "] x [" << ylo
       << ", " << yhi << "]\n";

  // ==========================================================
  // LEVEL 1 (Fine): 4:1 refined patch over the tagged region
  // ==========================================================
  const int ref = 4;
  int mf = (i_hi - i_lo + 1) * ref;
  int nf = (j_hi - j_lo + 1) * ref;
  double dx_f = (xhi - xlo) / mf;
  double dy_f = (yhi - ylo) / nf;
  vec xf = build_coords(mf, xlo, xhi);
  vec yf = build_coords(nf, ylo, yhi);

  Laplacian L_f(k, mf, nf, dx_f, dy_f);
  RobinBC BC_f(k, mf, dx_f, nf, dy_f, 1, 0); // Dirichlet
  L_f = L_f + BC_f;

  // RHS: source in interior, interpolated coarse values on boundary
  mat rhs_f(mf + 2, nf + 2, fill::zeros);
  for (int i = 1; i <= mf; ++i)
    for (int j = 1; j <= nf; ++j)
      rhs_f(i, j) = -source_term(xf(i), yf(j));

  // Dirichlet BCs from coarse solution (coarse-to-fine coupling)
  for (int j = 0; j <= nf + 1; ++j) {
    rhs_f(0, j) = bilinear_interp(U_c, xc, yc, xf(0), yf(j));
    rhs_f(mf + 1, j) = bilinear_interp(U_c, xc, yc, xf(mf + 1), yf(j));
  }
  for (int i = 1; i <= mf; ++i) {
    rhs_f(i, 0) = bilinear_interp(U_c, xc, yc, xf(i), yf(0));
    rhs_f(i, nf + 1) = bilinear_interp(U_c, xc, yc, xf(i), yf(nf + 1));
  }

#ifdef EIGEN
  vec sol_f = Utils::spsolve_eigen(L_f, vectorise(rhs_f));
#else
  vec sol_f = spsolve(L_f, vectorise(rhs_f));
#endif
  mat U_f = reshape(sol_f, mf + 2, nf + 2);

  cout << "\n[Level 1] Fine patch: " << mf << " x " << nf << " cells ("
       << (mf + 2) * (nf + 2) << " unknowns)\n"
       << "          Peak value : " << U_f.max() << "\n";

  // ==========================================================
  // SUMMARY
  // ==========================================================
  int amr_total = (mc + 2) * (nc + 2) + (mf + 2) * (nf + 2);
  int uniform_eq = (mc * ref + 2) * (nc * ref + 2);
  cout << "\n[Summary]\n"
       << "  AMR unknowns:     " << amr_total << "\n"
       << "  Uniform equiv:    " << uniform_eq << " (" << mc * ref << "x"
       << nc * ref << " cells)\n"
       << "  Efficiency ratio: " << 100.0 * amr_total / uniform_eq
       << "% of uniform cost\n";

  // ==========================================================
  // OUTPUT
  // ==========================================================
  write_solution("amr_coarse.dat", U_c, xc, yc);
  write_solution("amr_fine.dat", U_f, xf, yf);

  return 0;
}
