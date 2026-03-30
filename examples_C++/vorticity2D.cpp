/**
 * This example uses MOLE to compute the mimetic curl of a 2D vector field
 *
 * Given F = (sin(pi*y), cos(pi*x)) on [0,1]^2, the analytical curl is:
 *   curl(F) = dFy/dx - dFx/dy = -pi*sin(pi*x) - pi*cos(pi*y)
 *
 * Three blocks are printed separated by two blank lines (gnuplot index):
 *   0 — computed curl
 *   1 — analytical curl
 *   2 — pointwise error
 *
 * Plot with:
 *   gnuplot> set multiplot layout 1,3 title "Curl of F = (sin({/Symbol p}y), cos({/Symbol p}x))"
 *   gnuplot> set pm3d map
 *   gnuplot> splot 'solution' index 0 title "Computed"
 *   gnuplot> splot 'solution' index 1 title "Analytical"
 *   gnuplot> splot 'solution' index 2 title "Error"
 *   gnuplot> unset multiplot
 */

#include "mole.h"
#include <iostream>

using namespace std;
using namespace arma;

int main() {
  int k = 4;  // Order of accuracy
  int m = 50; // Number of cells in x
  int n = 50; // Number of cells in y

  Real west = 0.0, east = 1.0;
  Real south = 0.0, north = 1.0;
  Real dx = (east - west) / m;
  Real dy = (north - south) / n;

  // Build mimetic curl operator
  Curl C(k, m, n, dx, dy);

  // Staggered grid positions

  // Node positions (cell boundaries)
  vec xnodes = linspace(west, east, m + 1);
  vec ynodes = linspace(south, north, n + 1);

  // Cell-center positions + boundaries
  vec xcb(m + 2);
  xcb(0) = west;
  for (int i = 1; i <= m; ++i)
    xcb(i) = west + (i - 0.5) * dx;
  xcb(m + 1) = east;

  vec ycb(n + 2);
  ycb(0) = south;
  for (int j = 1; j <= n; ++j)
    ycb(j) = south + (j - 0.5) * dy;
  ycb(n + 1) = north;

  // Sample F = (sin(pi*y), cos(pi*x)) on face-staggered grid

  // Fx at x-faces: (m+1) x-positions (xnodes) by (n+2) y-positions (ycb)
  vec Fx((m + 1) * (n + 2));
  for (int j = 0; j < n + 2; ++j)
    for (int i = 0; i < m + 1; ++i)
      Fx(j * (m + 1) + i) = sin(datum::pi * ycb(j));

  // Fy at y-faces: (m+2) x-positions (xcb) by (n+1) y-positions (ynodes)
  vec Fy((m + 2) * (n + 1));
  for (int j = 0; j < n + 1; ++j)
    for (int i = 0; i < m + 2; ++i)
      Fy(j * (m + 2) + i) = cos(datum::pi * xcb(i));

  // Compute mimetic curl
  vec curlF = C * join_vert(Fx, Fy);

  // Analytical curl at node positions
  vec exact((m + 1) * (n + 1));
  for (int j = 0; j < n + 1; ++j)
    for (int i = 0; i < m + 1; ++i)
      exact(j * (m + 1) + i) = -datum::pi * sin(datum::pi * xnodes(i)) -
                               datum::pi * cos(datum::pi * ynodes(j));

  // Output for gnuplot (x y z format, blocks separated by blank lines)

  // Block 0: Computed
  for (int j = 0; j < n + 1; ++j) {
    for (int i = 0; i < m + 1; ++i)
      cout << xnodes(i) << " " << ynodes(j) << " " << curlF(j * (m + 1) + i)
           << "\n";
    cout << "\n"; // blank line between rows for gnuplot
  }
  cout << "\n"; // double blank line = new index

  // Block 1: Analytical
  for (int j = 0; j < n + 1; ++j) {
    for (int i = 0; i < m + 1; ++i)
      cout << xnodes(i) << " " << ynodes(j) << " " << exact(j * (m + 1) + i)
           << "\n";
    cout << "\n";
  }
  cout << "\n";

  // Block 2: Error
  vec err = curlF - exact;
  for (int j = 0; j < n + 1; ++j) {
    for (int i = 0; i < m + 1; ++i)
      cout << xnodes(i) << " " << ynodes(j) << " " << err(j * (m + 1) + i)
           << "\n";
    cout << "\n";
  }

  // Print max error to stderr so it doesn't interfere with gnuplot
  cerr << "Max absolute error: " << norm(err, "inf") << endl;

  return 0;
}
