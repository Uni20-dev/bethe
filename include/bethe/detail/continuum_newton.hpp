// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/newton_backtracking.hpp>
#include <bethe/solver.hpp>
#include <uni20/linalg/ops/linear_solve.hpp>

namespace bethe::detail
{
/// c/(c*c+x*x), without squaring large/small dimensional parameters.
template <uni20::Real Real> Real rational_scattering_kernel(Real x, Real c)
{
  if (std::abs(x) > c)
  {
    Real const ratio = c / x;
    return (ratio / x) / (Real{1} + ratio * ratio);
  }
  Real const ratio = x / c;
  return (Real{1} / c) / (Real{1} + ratio * ratio);
}

/// Shared physical-domain damped Newton iteration for real-root Bethe systems.
/// The caller owns validation, seed, residual scaling, observables and the final
/// residual check after conversion back to physical units. Returns accepted steps.
template <uni20::Real Real, typename System>
std::size_t continuum_newton(System const& system, std::vector<Real>& q, SolverOptions<Real> const& options)
{
  auto const n = q.size();
  uni20::DenseMatrix<Real> jacobian(n, n), step(n, 1);
  std::size_t iterations = 0;
  for (;;)
  {
    auto const evaluation = system.evaluate(q, &jacobian);
    if (evaluation.norm <= options.residual_tolerance || iterations == options.max_iterations) break;
    for (std::size_t j = 0; j < n; ++j)
      step[j, 0] = -evaluation.residual[j];
    uni20::linalg::solve_inplace(jacobian, step);
    bool const accepted = backtrack_newton(
        q, [&](std::size_t j) { return step[j, 0]; },
        [&](auto const& trial, Real damping) {
          return system.physical(trial) &&
                 newton_decreases(system.evaluate(trial).norm, evaluation.norm, damping, options.residual_tolerance);
        });
    if (!accepted) break;
    ++iterations;
  }
  return iterations;
}
} // namespace bethe::detail
