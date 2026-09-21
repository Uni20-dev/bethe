// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2015-2023, 2026 Ian McCulloch
// Successor to MPToolkit's misc/heisenberg-energy.cpp. See CITATIONS.md.

#pragma once

#include <uni20/core/types.hpp>
#if UNI20_HAS_FLOAT128
// Supplies the configured provider's scalar math overloads and numeric limits.
#include <mplapack_binary128.h>
#endif
#include <uni20/core/math.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace bethe::heisenberg
{

/// Controls the real-root, zero-field ground-state fixed-point iteration.
template <uni20::Real Real> struct SolverOptions
{
  Real residual_tolerance = Real{32} * uni20::numeric_limits<Real>::epsilon();
  std::size_t max_iterations = 10000;
};

/// The last iterate, including on iteration-limit exhaustion.
/// Energy includes the ferromagnetic reference N/4 for J=1.
template <uni20::Real Real> struct GroundState
{
  std::vector<Real> rapidities;
  Real energy = Real{0};
  /// Infinity norm of the logarithmic Bethe equations divided by N (radians).
  Real residual_norm = Real{0};
  /// Number of simultaneous updates; evaluating the initial guess is not an update.
  std::size_t iterations = 0;
  bool converged = false;
};

namespace detail
{
// Keep roundoff in the phase and energy sums from growing with the root count.
template <typename Real> class CompensatedSum
{
public:
  void add(Real value)
  {
    Real const corrected = value - correction_;
    Real const next = sum_ + corrected;
    correction_ = (next - sum_) - corrected;
    sum_ = next;
  }

  Real value() const { return sum_; }

private:
  Real sum_ = Real{0};
  Real correction_ = Real{0};
};
} // namespace detail

/// Solve the even-N periodic spin-1/2 XXX antiferromagnet at zero field.
///
/// H = sum_i S_i . S_(i+1), with S_N = S_0 and J=1. For N=2 this
/// periodic sum contains two copies of the bond, so the ground energy is -3/2.
/// Uses N/2 real roots and I_i = i - (N/2 - 1)/2, starting from zero.
/// A sweep costs O(N^2) and storage is O(N). No energy-accuracy bound is implied
/// by the residual tolerance. Zero max_iterations evaluates only the initial guess.
///
/// Throws invalid_argument for unsupported sizes or a nonpositive/nonfinite
/// tolerance, and runtime_error for nonfinite arithmetic. Exhausting the
/// iteration budget returns converged=false, with diagnostics for the last roots.
template <uni20::Real Real = double>
[[nodiscard]] GroundState<Real> ground_state(std::size_t sites, SolverOptions<Real> const& options = {})
{
  using std::abs;
  using std::atan;
  using std::tan;

  if (sites < 2 || sites % 2 != 0)
    throw std::invalid_argument("Heisenberg ground state requires an even number of sites >= 2");
  if (!uni20::isfinite(options.residual_tolerance) || options.residual_tolerance <= Real{0})
    throw std::invalid_argument("residual tolerance must be finite and positive");

  std::size_t const roots = sites / 2;
  Real const n = static_cast<Real>(sites);
  Real const pi = Real{4} * atan(Real{1});
  GroundState<Real> result;
  result.rapidities.resize(roots, Real{0});
  std::vector<Real> angles(roots);

  for (;;)
  {
    // Evaluate both the residual and the next fixed-point angles at the same
    // current roots. Never report a residual belonging to the previous iterate.
    result.residual_norm = Real{0};
    for (std::size_t i = 0; i < roots; ++i)
    {
      detail::CompensatedSum<Real> phase;
      for (std::size_t j = 0; j < roots; ++j)
      {
        if (i != j)
          phase.add(atan((result.rapidities[i] - result.rapidities[j]) / Real{2}));
      }
      Real const quantum_number = static_cast<Real>(i) - (static_cast<Real>(roots) - Real{1}) / Real{2};
      angles[i] = (pi * quantum_number + phase.value()) / n;
      Real const residual = Real{2} * abs(atan(result.rapidities[i]) - angles[i]);
      if (!uni20::isfinite(residual))
        throw std::runtime_error("nonfinite Heisenberg equation residual");
      result.residual_norm = std::max(result.residual_norm, residual);
    }

    result.converged = result.residual_norm <= options.residual_tolerance;
    if (result.converged || result.iterations == options.max_iterations)
      break;

    // All angles were computed before any roots change: this is a Jacobi update.
    for (std::size_t i = 0; i < roots; ++i)
    {
      result.rapidities[i] = tan(angles[i]);
      if (!uni20::isfinite(result.rapidities[i]))
        throw std::runtime_error("nonfinite Heisenberg rapidity");
    }
    ++result.iterations;
  }

  detail::CompensatedSum<Real> magnons;
  for (Real const z : result.rapidities)
    magnons.add(-Real{2} / (Real{1} + z * z));
  result.energy = n / Real{4} + magnons.value();
  if (!uni20::isfinite(result.energy))
    throw std::runtime_error("nonfinite Heisenberg energy");
  return result;
}

} // namespace bethe::heisenberg
