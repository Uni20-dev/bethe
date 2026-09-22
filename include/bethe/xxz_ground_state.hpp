// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/xxz_common.hpp>

namespace bethe::xxz
{
enum class GroundSolveStatus
{
  converged,
  iteration_limit,
  stalled
};

enum class GroundResidualConvention
{
  logarithmic_phase,   // max|F|/N, Delta>=0
  negative_rank_scaled // rank-subtracted equations divided by N*s
};

/// Periodic ground-state result, deliberately distinct from RealState, which
/// is the explicit-label/excitation type. No implicit conversion can discard
/// representation-specific coordinates or silently change residual units.
/// Currently supports Delta>=0, plus -1<Delta<0 on EVEN rings only.
template <uni20::Real Real> struct GroundState
{
    std::vector<Real> rapidities; // Scaled z, as in the nonnegative solver.
    /// Negative Delta: lambda with z=s*tanh(lambda), s=sqrt((1+Delta)/(1-Delta)).
    /// Empty for Delta>=0; do not reconstruct these from rounded z values.
    std::vector<Real> log_rapidities;
    QuantumNumbers quantum_numbers;
    Real delta = Real{0};
    Real energy = Real{0};
    /// Units depend on residual_convention; not an energy-error estimate.
    Real residual_norm = Real{0};
    GroundResidualConvention residual_convention = GroundResidualConvention::logarithmic_phase;
    std::size_t iterations = 0;
    bool converged = false;
    GroundSolveStatus status = GroundSolveStatus::iteration_limit;
    uni20::half_int sz;
    std::size_t momentum_index = 0;
    Real momentum = Real{0};
    bool spin_reversed = false;
};
} // namespace bethe::xxz
