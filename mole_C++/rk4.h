/**
 * @file rk4.h
 * @brief Explicit Runge-Kutta 4th-order time integrator
 *
 * The integrator itself is not mimetic -- RK4 is RK4. All it ever needs is a
 * way to evaluate dy/dt from y, so the mimetic content lives entirely in what
 * the caller hands over: either an assembled operator A, used only as A*y, or a
 * function that applies the operators in sequence and never assembles anything.
 */

#ifndef RK4_H
#define RK4_H

#include "utils.h"
#include <functional>

/**
 * @brief Right-hand side of dy/dt = f(t, y)
 */
typedef std::function<vec(Real, const vec &)> RhsFunction;

/**
 * @brief Called after every step, as callback(t, y, step) with step 1-based
 */
typedef std::function<void(Real, const vec &, u32)> StepCallback;

/**
 * @brief Integrates dy/dt = f(t, y), keeping the whole trajectory
 *
 * @param f Right-hand side
 * @param tspan {t0, tf}
 * @param dt Step size
 * @param y0 Initial condition
 * @return mat One column per evaluation point, numel(y0) by number of steps + 1
 */
mat rk4(const RhsFunction &f, const vec &tspan, Real dt, const vec &y0);

/**
 * @brief Integrates dy/dt = f(t, y), keeping only the current state
 *
 * The callback runs after every step, for plotting or diagnostics. Only the
 * current state is held, so the memory stays flat -- which is what a PDE needs,
 * the full trajectory being numel(y0) by number of steps.
 *
 * @return vec The final state
 */
vec rk4(const RhsFunction &f, const vec &tspan, Real dt, const vec &y0,
        const StepCallback &callback);

/**
 * @brief Linear case dy/dt = A*y, keeping the whole trajectory
 *
 * Only the product A*y is ever formed, so A may be assembled however the
 * problem requires.
 */
mat rk4(const sp_mat &A, const vec &tspan, Real dt, const vec &y0);

/**
 * @brief Linear case dy/dt = A*y, keeping only the current state
 */
vec rk4(const sp_mat &A, const vec &tspan, Real dt, const vec &y0,
        const StepCallback &callback);

#endif // RK4_H
