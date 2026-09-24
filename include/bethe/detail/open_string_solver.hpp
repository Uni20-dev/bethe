// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/xxz_quantum_group.hpp>

namespace bethe::xxz::quantum_group::detail
{
template <uni20::Real Real> struct StringPhase
{
    Real value, angle, width;
};

/// Theta(beta;w), with derivatives, and its complementary phase near w=0.
template <uni20::Real Real> StringPhase<Real> string_phase(Real beta, Real w, bool complement = false)
{
  Real const s = std::sin(beta / Real{2}), c = std::cos(beta / Real{2}), t = std::tanh(w);
  Real const den = s * s + t * t * c * c;
  Real const da = t / den, dw = -Real{2} * s * c * (Real{1} - t * t) / den;
  if (complement) return {Real{2} * std::atan2(t * c, s), -da, -dw};
  return {Real{2} * std::atan2(s, t * c), da, dw};
}

template <uni20::Real Real> struct StringIteration
{
    std::vector<Real> x;
    std::size_t iterations = 0;
    bool converged = false;
    SolveStatus status = SolveStatus::iteration_limit;
};

/// Shared logarithmic-string Newton driver. The system owns the equations,
/// ideal initializer, admissible domain and optional coordinate scaling.
/// Budgets count both stages; acceptance and diagnostics use original units.
template <uni20::Real Real, typename System>
StringIteration<Real> solve_log_string(System const& system, SolverOptions<Real> const& options)
{
  if (!uni20::isfinite(options.residual_tolerance) || options.residual_tolerance <= Real{0})
    throw std::invalid_argument("residual tolerance must be finite and positive");
  check_matrix_size<Real>(system.order);
  StringIteration<Real> result{.x = system.seed()};
  auto& x = result.x;
  if (!system.physical(x)) throw std::overflow_error("string seed is not resolvable at the selected precision");
  std::vector<Real> jac, step(system.order);
  bool ideal = true;
  for (;;)
  {
    auto evaluation = system.evaluate(x, &jac, ideal);
    if (ideal && evaluation.norm <= std::sqrt(uni20::numeric_limits<Real>::epsilon()))
    {
      ideal = false;
      evaluation = system.evaluate(x, &jac);
    }
    if (!ideal && evaluation.norm <= options.residual_tolerance)
    {
      result.converged = true;
      result.status = SolveStatus::converged;
      break;
    }
    if (result.iterations == options.max_iterations) break;
    for (std::size_t i = 0; i < system.order; ++i)
      step[i] = -evaluation.residual[i];
    auto const scales = system.coordinate_scales(x);
    if (!scales.empty())
      for (std::size_t i = 0; i < system.order; ++i)
      {
        Real row_scale{};
        for (std::size_t j = 0; j < system.order; ++j)
        {
          jac[i * system.order + j] *= scales[j];
          row_scale = std::max(row_scale, std::abs(jac[i * system.order + j]));
        }
        if (row_scale > Real{0})
        {
          for (std::size_t j = 0; j < system.order; ++j)
            jac[i * system.order + j] /= row_scale;
          step[i] /= row_scale;
        }
      }
    if (!bethe::detail::newton_step(jac, step))
    {
      result.status = SolveStatus::singular_jacobian;
      break;
    }
    if (!scales.empty())
      for (std::size_t j = 0; j < system.order; ++j)
        step[j] *= scales[j];
    auto trial = x;
    Real damping{1};
    bool accepted = false;
    for (int backtrack = 0; backtrack <= uni20::numeric_limits<Real>::digits; ++backtrack)
    {
      for (std::size_t i = 0; i < system.order; ++i)
        trial[i] = x[i] + damping * step[i];
      if constexpr (requires { system.normalize(trial); }) system.normalize(trial);
      if (system.physical(trial) &&
          system.evaluate(trial, nullptr, ideal).norm < (Real{1} - damping / Real{10000}) * evaluation.norm)
      {
        accepted = true;
        break;
      }
      damping /= Real{2};
    }
    if (!accepted)
    {
      result.status = SolveStatus::stalled;
      break;
    }
    x = std::move(trial);
    ++result.iterations;
  }
  return result;
}
} // namespace bethe::xxz::quantum_group::detail
