// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/xxz.hpp>
#include <optional>

namespace bethe::xxz::open
{
enum class GroundSolveStatus
{
  converged,
  iteration_limit,
  stalled
};

/// Residuals from different coordinate systems are not interchangeable.
enum class GroundResidualConvention
{
  logarithmic_phase,   // max|F|/(2*N), 0 <= Delta <= 1
  massive_regularized, // bulk phase and optional regularized boundary
  negative_rank_scaled // rank-subtracted equations divided by N*s
};

/// The last root of a massive even, zero-magnetization chain.
/// With r=(Delta-1)/(Delta+1), y=1/z_B^2=r^2*expm1(-log_distance).
/// Negative y means a purely imaginary z_B; y=0 is a regular crossing.
/// log_distance retains exponentially small deviations from y=-r^2.
template <uni20::Real Real> struct BoundaryRoot
{
    Real log_distance = Real{0};
    Real inverse_square = Real{0};
    uni20::half_int quantum_number;
};

/// Free-end XXZ sector minimum, Delta>-1, no boundary fields or momentum.
/// Distinct from RealState: no implicit conversion can silently discard a
/// boundary root. Bulk arrays exclude the optional distinguished root/label.
template <uni20::Real Real> struct GroundState
{
    std::vector<Real> rapidities;
    /// For -1<Delta<0, lambda with z=s*tanh(lambda), s=sqrt((1+Delta)/(1-Delta)).
    /// Empty for nonnegative Delta. Retain these for the scaled equations;
    /// reconstructing lambda from rounded z can lose accuracy near an endpoint.
    std::vector<Real> log_rapidities;
    QuantumNumbers quantum_numbers;
    std::optional<BoundaryRoot<Real>> boundary_root;
    Real delta = Real{0};
    /// Continuation stage reached; equals delta for the gapless/XXX solver.
    /// Energy and residual always refer to the requested delta.
    Real root_delta = Real{0};
    Real energy = Real{0};
    /// Interpretation is given by residual_convention. See docs/xxz-negative.md
    /// for the rank-subtracted negative-Delta equations, or
    /// docs/xxz-open-massive.md for the regularized massive boundary component.
    /// Neither an energy-error bound nor a bound on exponentially small gaps.
    Real residual_norm = Real{0};
    GroundResidualConvention residual_convention = GroundResidualConvention::logarithmic_phase;
    std::size_t iterations = 0;
    bool converged = false;
    GroundSolveStatus status = GroundSolveStatus::iteration_limit;
    uni20::half_int sz;
    bool spin_reversed = false;
};
} // namespace bethe::xxz::open
