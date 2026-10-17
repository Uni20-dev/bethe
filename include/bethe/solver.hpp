// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#pragma once

#include <cstddef>
#include <stdexcept>
#include <uni20/core/math.hpp>
#include <uni20/core/runtime_precision.hpp>

namespace bethe
{
/// Controls a residual-based iteration. A zero budget evaluates the initial
/// guess only. Tolerance is in normalized equation units, not energy units.
template <uni20::Real Real> struct SolverOptions
{
    Real residual_tolerance = Real{32} * uni20::numeric_limits<Real>::epsilon();
    std::size_t max_iterations = 10000;
};

#if UNI20_HAS_MPREAL
/// Runtime-precision solvers require a finite precision per call. No MPFR global
/// or thread-local default is read or changed. Set custom tolerances at this precision.
template <> struct SolverOptions<uni20::mpreal>
{
    explicit SolverOptions(uni20::Precision precision)
        : residual_tolerance(default_tolerance(precision)), precision_(precision)
    {}
    uni20::mpreal residual_tolerance;
    std::size_t max_iterations = 10000;
    uni20::Precision precision() const { return precision_; }

  private:
    static uni20::mpreal default_tolerance(uni20::Precision p)
    {
      if (p.is_exact()) throw std::invalid_argument("a solver requires finite working precision");
      return 32 * uni20::numeric_limits<uni20::mpreal>::epsilon(p);
    }
    uni20::Precision precision_;
};
#endif

namespace detail
{
/// A finite-precision zero, without evaluating an unset runtime scalar.
template <uni20::Real Real> Real precision_zero(Real const& exemplar)
{
#if UNI20_HAS_MPREAL
  if constexpr (uni20::has_runtime_precision_v<Real>)
    if (!exemplar.initialized() || exemplar.is_exact())
      throw std::invalid_argument("a numerical calculation requires a finite-precision scalar");
#endif
  return uni20::scalar_like(exemplar, 0);
}

template <uni20::Real Real> void require_same_precision(Real const& value, Real const& exemplar)
{
  (void)value;
  (void)exemplar;
#if UNI20_HAS_MPREAL
  if constexpr (uni20::has_runtime_precision_v<Real>)
    if (!value.initialized() || value.is_exact() || value.precision() != exemplar.precision())
      throw std::invalid_argument("mixed, exact or unset solver values require explicit precision conversion");
#endif
}

template <uni20::Real Real> Real solver_zero(SolverOptions<Real> const& options)
{
  (void)options;
#if UNI20_HAS_MPREAL
  if constexpr (uni20::has_runtime_precision_v<Real>)
  {
    Real zero(0, options.precision());
    require_same_precision(options.residual_tolerance, zero);
    return zero;
  }
  else
#endif
    return Real{0};
}

template <uni20::Real Real> Real pi_like(Real const& exemplar)
{
  (void)exemplar;
#if UNI20_HAS_MPREAL
  if constexpr (uni20::has_runtime_precision_v<Real>)
    return uni20::pi<Real>.at(exemplar.precision());
  else
#endif
    return uni20::pi<Real>;
}

// Compensate phase and energy sums at the selected scalar precision.
template <typename Real> class CompensatedSum {
  public:
    CompensatedSum() = default;
    explicit CompensatedSum(Real const& zero) : sum_(zero), correction_(zero) {}
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
