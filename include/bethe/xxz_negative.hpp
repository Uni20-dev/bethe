// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/detail/newton_backtracking.hpp>
#include <bethe/xxz_common.hpp>
#include <uni20/linalg/ops/linear_solve.hpp>

namespace bethe::xxz::detail
{
enum class NegativeSolveStatus
{
  converged,
  iteration_limit,
  stalled
};

/// Internal ground-root engine for bipartite geometries only: free ends of
/// either parity, or even periodic rings. Odd periodic sector ground states
/// can have complex roots and MUST NOT be routed through this system.
template <uni20::Real Real> struct NegativeGroundRoots
{
    std::vector<Real> rapidities;     // z=s*tanh(lambda)
    std::vector<Real> log_rapidities; // lambda, used by the regularized equations
    QuantumNumbers quantum_numbers;
    Real energy = Real{0}, residual_norm = Real{0};
    std::size_t iterations = 0;
    bool converged = false;
    NegativeSolveStatus status = NegativeSolveStatus::iteration_limit;
};

template <uni20::Real Real> class NegativeGroundSystem {
  public:
    NegativeGroundSystem(std::size_t sites, std::size_t roots, Real delta, bool open)
        : sites(sites), order(roots), delta(delta), open(open), plus(Real{1} + delta), minus(Real{1} - delta),
          a(-delta), pi(bethe::detail::pi_like(delta))
    {
      checked_sites(sites);
      if (!uni20::isfinite(delta) || delta <= -Real{1} || delta >= Real{0})
        throw std::invalid_argument("negative XXZ root engine requires -1 < Delta < 0");
      if (roots > sites / 2) throw std::invalid_argument("negative XXZ roots require M <= N/2");
      if (!open && sites % 2)
        throw std::invalid_argument("odd periodic negative XXZ requires a complex-root and sector-selection solver");
      s = uni20::math::sqrt(plus / minus);
      b = uni20::math::sqrt(plus) * uni20::math::sqrt(minus);
    }

    std::vector<Real> seed() const
    {
      std::vector<Real> x(order);
      for (std::size_t i = 0; i < order; ++i)
      {
        Real const angle = open ? pi * Real(i + 1) / (Real{2} * Real(sites + 1))
                                : pi * (Real(i) - Real(order - 1) / Real{2}) / Real(sites);
        x[i] = uni20::math::atanh(uni20::math::tan(angle));
      }
      return x;
    }

    bool physical(std::span<Real const> x) const
    {
      if (x.size() != order) return false;
      for (std::size_t i = 0; i < order; ++i)
      {
        if (!uni20::isfinite(x[i]) || (open && x[i] <= Real{0}) || (i && x[i] <= x[i - 1])) return false;
        // A rounded endpoint cannot encode a finite hyperbolic root.
        if (uni20::math::abs(uni20::math::tanh(x[i])) >= Real{1}) return false;
      }
      return true;
    }

    // Complementary half-scattering phase divided by s. Exact rank-sized
    // multiples of pi were canceled analytically before numerical evaluation.
    Real phase(Real x) const
    {
      return uni20::math::copysign(uni20::math::atan2(b, a * uni20::math::abs(uni20::math::tanh(x))) / s, x);
    }
    Real derivative(Real x) const
    {
      Real const v = uni20::math::tanh(x), c = Real{1} / uni20::math::cosh(x);
      return -minus * a * c * c / (b * b + a * a * v * v);
    }

    struct Evaluation
    {
        std::vector<Real> residual;
        Real norm = Real{0};
    };

    Evaluation evaluate(std::span<Real const> x, uni20::DenseMatrix<Real>* jacobian = nullptr) const
    {
      Real const zero = bethe::detail::precision_zero(delta);
      Evaluation out{.residual = std::vector<Real>(order, zero), .norm = zero};
      for (std::size_t i = 0; i < order; ++i)
      {
        Real const v = uni20::math::tanh(x[i]), c = Real{1} / uni20::math::cosh(x[i]);
        Real const multiplier = open ? Real{2} : Real{1};
        bethe::detail::CompensatedSum<Real> f(zero), diagonal(zero);
        f.add(multiplier * Real(sites) * uni20::math::atan(s * v) / s);
        diagonal.add(multiplier * Real(sites) * c * c / (Real{1} + s * s * v * v));
        if (open)
        {
          f.add(-Real{2} * uni20::math::atan(s / v) / s);
          diagonal.add(Real{2} * c * c / (v * v + s * s));
        }
        for (std::size_t j = 0; j < order; ++j)
          if (i != j)
          {
            f.add(-phase(x[i] - x[j]));
            Real const direct = derivative(x[i] - x[j]);
            diagonal.add(-direct);
            Real reflected = Real{0};
            if (open)
            {
              f.add(-phase(x[i] + x[j]));
              reflected = derivative(x[i] + x[j]);
              diagonal.add(-reflected);
            }
            if (jacobian) (*jacobian)[i, j] = (direct - reflected) / Real(sites);
          }
        if (jacobian) (*jacobian)[i, i] = diagonal.value() / Real(sites);
        out.residual[i] = f.value() / Real(sites);
        if (!uni20::isfinite(out.residual[i])) throw std::runtime_error("nonfinite negative XXZ residual");
        out.norm = std::max(out.norm, uni20::math::abs(out.residual[i]));
      }
      return out;
    }

    std::size_t sites, order;
    Real delta;
    bool open;
    Real plus, minus, s, b, a, pi;
};

/// Newton solves rank-subtracted equations divided by s=sqrt((1+Delta)/(1-Delta)).
/// That scaling prevents false convergence as Delta approaches -1. The output
/// residual uses these scaled equations, not the existing nonnegative solver's
/// normalization. Every accepted Newton step consumes the shared budget.
template <uni20::Real Real>
NegativeGroundRoots<Real> negative_ground_roots(std::size_t sites, Real delta, uni20::half_int sz, bool open,
                                                SolverOptions<Real> const& options = {})
{
  Real const zero = solver_zero(delta, options);
  auto const m = sector_roots(sites, sz);
  NegativeGroundSystem<Real> const system(sites, m, delta, open);
  auto const elements = std::size_t(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(Real);
  if (m && m > elements / m) throw std::length_error("negative XXZ Newton matrix is too large");
  NegativeGroundRoots<Real> state;
  auto x = system.seed();
  uni20::DenseMatrix<Real> jacobian(m, m, options.precision()), step(m, 1, options.precision());
  for (;;)
  {
    if (!system.physical(x)) throw std::runtime_error("negative XXZ roots cannot be represented at this precision");
    auto const evaluation = system.evaluate(x, &jacobian);
    state.residual_norm = evaluation.norm;
    if (evaluation.norm <= options.residual_tolerance)
    {
      state.converged = true;
      state.status = NegativeSolveStatus::converged;
      break;
    }
    if (state.iterations == options.max_iterations) break;
    for (std::size_t i = 0; i < m; ++i)
      step[i, 0] = -evaluation.residual[i];
    if (!uni20::linalg::solve_inplace_with_info(jacobian, step, options.precision()).succeeded())
    {
      state.status = NegativeSolveStatus::stalled;
      break;
    }
    bool const accepted = bethe::detail::backtrack_newton(
        x, options.precision(), [&](std::size_t j) { return step[j, 0]; },
        [&](auto const& trial, Real damping) {
          return system.physical(trial) && bethe::detail::newton_decreases(system.evaluate(trial).norm, evaluation.norm,
                                                                           damping, options.residual_tolerance);
        });
    if (!accepted)
    {
      state.status = NegativeSolveStatus::stalled;
      break;
    }
    ++state.iterations;
  }
  state.log_rapidities = x;
  bethe::detail::CompensatedSum<Real> energy(zero);
  energy.add(Real(open ? sites - 1 : sites) * delta / Real{4});
  for (std::size_t i = 0; i < m; ++i)
  {
    Real const v = uni20::math::tanh(x[i]), c = Real{1} / uni20::math::cosh(x[i]);
    Real const z = system.s * v;
    state.rapidities.push_back(z);
    energy.add(-system.plus * c * c / (Real{1} + z * z));
    state.quantum_numbers.push_back(open ? uni20::half_int(std::int64_t(i + 1))
                                         : uni20::from_twice(2 * std::int64_t(i) - std::int64_t(m - 1)));
  }
  state.energy = energy.value();
  if (!uni20::isfinite(state.energy)) throw std::runtime_error("nonfinite negative XXZ energy");
  return state;
}
} // namespace bethe::xxz::detail
