#include "rk4.h"
#include <cmath>

namespace {

// Number of steps spanned by tspan at this dt, matching the MATLAB colon
// operator. dt usually comes from (tf - t0)/steps, whose rounding can leave the
// ratio a hair under an integer, so absorb that rather than lose a step.
u32 stepCount(const vec &tspan, Real dt) {
  const Real ratio = (tspan(1) - tspan(0)) / dt;

  if (!(ratio > 0))
    return 0;

  u32 n = static_cast<u32>(std::floor(ratio));
  if (ratio - static_cast<Real>(n) > 1.0 - 1e-8)
    n += 1;
  return n;
}

// One classical RK4 stage sweep
vec step(const RhsFunction &f, Real t, Real dt, const vec &y) {
  const vec k1 = f(t, y);
  const vec k2 = f(t + dt / 2, y + dt / 2 * k1);
  const vec k3 = f(t + dt / 2, y + dt / 2 * k2);
  const vec k4 = f(t + dt, y + dt * k3);
  return y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
}

} // namespace

mat rk4(const RhsFunction &f, const vec &tspan, Real dt, const vec &y0) {
  const u32 nsteps = stepCount(tspan, dt);

  mat y(y0.n_elem, nsteps + 1);
  y.col(0) = y0;

  for (u32 i = 0; i < nsteps; ++i)
    y.col(i + 1) = step(f, tspan(0) + i * dt, dt, y.col(i));

  return y;
}

vec rk4(const RhsFunction &f, const vec &tspan, Real dt, const vec &y0,
        const StepCallback &callback) {
  const u32 nsteps = stepCount(tspan, dt);

  vec y = y0;

  for (u32 i = 0; i < nsteps; ++i) {
    y = step(f, tspan(0) + i * dt, dt, y);
    if (callback)
      callback(tspan(0) + (i + 1) * dt, y, i + 1);
  }

  return y;
}

mat rk4(const sp_mat &A, const vec &tspan, Real dt, const vec &y0) {
  return rk4(RhsFunction([&A](Real, const vec &y) { return vec(A * y); }),
             tspan, dt, y0);
}

vec rk4(const sp_mat &A, const vec &tspan, Real dt, const vec &y0,
        const StepCallback &callback) {
  return rk4(RhsFunction([&A](Real, const vec &y) { return vec(A * y); }),
             tspan, dt, y0, callback);
}
