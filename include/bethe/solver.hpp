// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#pragma once

#include <cstddef>
#include <uni20/core/math.hpp>

namespace bethe
{
/// Controls a residual-based iteration. A zero budget evaluates the initial
/// guess only. Tolerance is in normalized equation units, not energy units.
template <uni20::Real Real> struct SolverOptions
{
    Real residual_tolerance = Real{32} * uni20::numeric_limits<Real>::epsilon();
    std::size_t max_iterations = 10000;
};

namespace detail
{
// Compensate phase and energy sums at the selected scalar precision.
template <typename Real> class CompensatedSum {
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
} // namespace bethe
