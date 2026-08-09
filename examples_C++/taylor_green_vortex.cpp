// clang-format off
/**
 * Solves the 2D advection-diffusion equation for a passive scalar carried by a
 * static field of Taylor-Green vortices, using Corbino-Castillo mimetic
 * operators and RK4. C++ counterpart of examples_MATLAB/taylor_green_vortex.m.
 *
 *   ds/dt = div(D grad s) - div(v s)    on [-pi, pi]^2, periodic
 *
 *   v(x, y) = [10 sin(x) cos(y), -10 cos(x) sin(y)]   Taylor-Green, div v = 0
 *   s(x, y, 0) = exp(-r^2/(2 sigma^2))                r from (0, 0)
 *
 *   Every operator carries the periodic flag, which wraps the interior
 *   staggered stencil around the seam, and mass is then conserved to roundoff.
 *   That is the BC this problem wants anyway: the domain, the vortices and the
 *   velocity field are all 2*pi-periodic, and the Taylor-Green field is
 *   tangent to the boundary (v.n = 0 on all four sides).
 *
 *   dt is derived from the spectral radius of the assembled operator, so it
 *   stays correct for any k, grid, Diff. Coeff. or velocity field.
 *
 * Output is gnuplot blocks separated by two blank lines, so each is one gnuplot
 * "index": the scalar field first, index 0 being the initial condition and the
 * rest every fifth step, then the Taylor-Green velocity field as the last index.
 *
 *   $ ./taylor_green_vortex > solution
 *
 *   gnuplot> stats 'solution' nooutput
 *   gnuplot> vel = STATS_blocks - 1              # velocity is the last index
 *   gnuplot> set size ratio -1
 *   gnuplot> set xrange [-pi:pi]; set yrange [-pi:pi]
 *   gnuplot> set palette rgbformulae 22,13,-31   # jet-like
 *   gnuplot> set cbrange [0:*]                   # rescaled every frame
 *   gnuplot> set pm3d map
 *
 *   the initial condition with the velocity field over it:
 *   gnuplot> splot 'solution' index 0 notitle, 'solution' index vel \
 *            using 1:2:(0):3:4:(0) with vectors lc rgb 'white' notitle
 *
 *   and the animation:
 *   gnuplot> do for [i=0:vel-1] { splot 'solution' index i notitle; pause 0.05 }
 */
// clang-format on

#include "mole.h"
#include <iostream>

using namespace std;
using namespace arma;

int main() {

  // Parameters
  int k = 2;           // Order of accuracy of the mimetic operators
  int m = 51;          // Number of cells along the x-axis
  int n = 51;          // Number of cells along the y-axis
  Real a = -datum::pi; // West
  Real b = datum::pi;  // East
  Real c = -datum::pi; // South
  Real d = datum::pi;  // North

  Real D_coeff = 0.5; // Diffusivity
  Real TIME = 5;      // s Simulated time

  Real dx = (b - a) / m;
  Real dy = (d - c) / n;

  // 2D staggered grid
  // Cell centers
  vec xc(m + 2); // m+2 long
  xc(0) = a;
  for (int i = 1; i <= m; ++i)
    xc(i) = a + (i - 0.5) * dx;
  xc(m + 1) = b;

  vec yc(n + 2); // n+2 long
  yc(0) = c;
  for (int j = 1; j <= n; ++j)
    yc(j) = c + (j - 0.5) * dy;
  yc(n + 1) = d;

  // Faces
  vec xf = linspace(a, b, m + 1); // m+1 long
  vec yf = linspace(c, d, n + 1); // n+1 long

  // MOLE flattens 2D fields with x running fastest, so a field indexed (i, j)
  // lives at j*(number of x positions) + i throughout.

  // Mimetic operators
  Divergence Dm(k, m, n, dx, dy, true); // Divergence
  Gradient G(k, m, n, dx, dy, true);    // Gradient
  Interpol I(m, n, 0.5, 0.5, true);     // Centers -> faces, centered
  // c1 = c2 = 0.5 keeps the advection term free of numerical diffusion.
  // Upwinding (c = 1 or 0) is not an option for a single interpolator here
  // anyway, since the Taylor-Green velocity changes sign inside the domain.
  //
  // Worth knowing before raising k: MOLE's interpolators are 2nd-order, so the
  // composite div*interpol is 2nd-order whatever k is.

  // Armadillo's expression templates do not accept the derived operator types,
  // so bind them to the sparse base before multiplying.
  sp_mat Dsp = Dm, Gsp = G, Isp = I;

  // Velocity field, sampled analytically on the faces it lives on
  // u sits on vertical faces: (m+1) x n,  v sits on horizontal faces: m x (n+1)
  int nu = (m + 1) * n;
  int nv = m * (n + 1);
  vec vel(nu + nv);
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < m + 1; ++i)
      vel(j * (m + 1) + i) = 10 * sin(xf(i)) * cos(yc(j + 1));
  for (int j = 0; j < n + 1; ++j)
    for (int i = 0; i < m; ++i)
      vel(nu + j * m + i) = -10 * cos(xc(i + 1)) * sin(yf(j));

  // Assemble the (constant, linear) right-hand side  ds/dt = A s
  // A = div(D grad .) - div(v .). Both terms are premultiplied once, so the
  // time loop is four sparse matrix-vector products per step.
  umat loc(2, nu + nv); // spdiags(vel, 0, numel(vel), numel(vel))
  for (int i = 0; i < nu + nv; ++i)
    loc(0, i) = loc(1, i) = i;
  sp_mat Vdiag(loc, vel, nu + nv, nu + nv);

  sp_mat A_diffusion = D_coeff * (Dsp * Gsp);
  sp_mat A_advection = Dsp * Vdiag * Isp;
  sp_mat A = A_diffusion - A_advection;

  // Time step
  // The absolute-stability region of classical RK4 contains the whole left
  // half-disk |z| <= 2.6155, Re z <= 0 (it cannot contain a full disk, since
  // |R(eps)| > 1 for real eps > 0). A is dissipative, so its spectrum lies in
  // Re <= 0, and dt*rho(A) <= 2.6155 then guarantees |R(dt*lambda)| <= 1 for
  // every eigenvalue, whatever the mix of advection (imaginary) and diffusion
  // (real).
  // Any induced norm bounds the spectral radius, and the infinity-norm (largest
  // absolute row sum) is exact, instant, and lands within 0.3-6% of rho(A) for
  // these operators -- no knowledge of the velocity field or of the stencil
  // constants needed, so refining the grid or raising k just works. Use inf and
  // not 1: the column sums are ~2.5x looser here.
  Real safety = 0.9;
  Real dt = safety * 2.6155 / norm(A, "inf");
  int steps = (int)ceil(TIME / dt);
  dt = TIME / steps; // trim so the run lands exactly on TIME

  // Initial condition: Gaussian sitting on the center of the domain
  Real x0 = 0;
  Real y0 = 0;
  Real sigma = datum::pi / 8;
  vec s((m + 2) * (n + 2));
  for (int j = 0; j < n + 2; ++j)
    for (int i = 0; i < m + 2; ++i)
      s(j * (m + 2) + i) =
          exp(-(pow(xc(i) - x0, 2) + pow(yc(j) - y0, 2)) / (2 * sigma * sigma));

  // One gnuplot block of the scalar field: "x y s", rows separated by a blank
  // line and the block closed by a second one, which gnuplot reads as an index.
  auto dump = [&](const vec &f) {
    for (int j = 0; j < n + 2; ++j) {
      for (int i = 0; i < m + 2; ++i)
        cout << xc(i) << " " << yc(j) << " " << f(j * (m + 2) + i) << "\n";
      cout << "\n";
    }
    cout << "\n";
  };

  dump(s); // index 0: the initial condition

  // Time loop: classical RK4, from mole_C++/rk4.h
  // Only the operator is handed over -- rk4 needs nothing but A*s, so the same
  // call works for any problem, and a function would work too if A were
  // nonlinear or never assembled:
  //   [&](Real, const vec &y) { return vec(Dsp*(D_coeff*(Gsp*y)) -
  //                                        Dsp*(vel % (Isp*y))); }
  // The callback fires after every step, so only the current state is held.
  s = rk4(A, {0.0, TIME}, dt, s, [&](Real, const vec &y, u32 step) {
    // Dump every few steps
    if (step % 5 == 0 || step == (u32)steps)
      dump(y);
  });

  // Last index: the static velocity field, subsampled so the arrows stay
  // readable. gnuplot draws vectors in data units and never auto-scales them
  // the way MATLAB's quiver does, so the components are scaled here to make the
  // longest arrow about one subsample spacing.
  int sk = max(1, (int)round(m / 16.0));
  Real arrow = sk * dx / 10.0;
  for (int j = 1; j <= n; j += sk) {
    for (int i = 1; i <= m; i += sk)
      cout << xc(i) << " " << yc(j) << " "
           << arrow * 10 * sin(xc(i)) * cos(yc(j)) << " "
           << arrow * -10 * cos(xc(i)) * sin(yc(j)) << "\n";
    cout << "\n";
  }

  return 0;
}
