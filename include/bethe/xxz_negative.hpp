// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/xxz.hpp>
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
          a(-delta), pi(Real{4} * std::atan(Real{1}))
    {
      checked_sites(sites);
      if (!uni20::isfinite(delta) || delta <= -Real{1} || delta >= Real{0})
        throw std::invalid_argument("negative XXZ root engine requires -1 < Delta < 0");
      if (roots > sites / 2) throw std::invalid_argument("negative XXZ roots require M <= N/2");
      if (!open && sites % 2)
        throw std::invalid_argument("odd periodic negative XXZ requires a complex-root and sector-selection solver");
      s = std::sqrt(plus / minus);
      b = std::sqrt(plus) * std::sqrt(minus);
    }

    std::vector<Real> seed() const
    {
      std::vector<Real> x(order);
      for (std::size_t i = 0; i < order; ++i)
      {
        Real const angle = open ? pi * Real(i + 1) / (Real{2} * Real(sites + 1))
                                : pi * (Real(i) - Real(order - 1) / Real{2}) / Real(sites);
        x[i] = std::atanh(std::tan(angle));
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
        if (std::abs(std::tanh(x[i])) >= Real{1}) return false;
      }
      return true;
    }

    // Complementary half-scattering phase divided by s. Exact rank-sized
    // multiples of pi were canceled analytically before numerical evaluation.
    Real phase(Real x) const { return std::copysign(std::atan2(b, a * std::abs(std::tanh(x))) / s, x); }
    Real derivative(Real x) const
    {
      Real const v = std::tanh(x), c = Real{1} / std::cosh(x);
      return -minus * a * c * c / (b * b + a * a * v * v);
    }

    struct Evaluation
    {
        std::vector<Real> residual;
        Real norm = Real{0};
    };

    Evaluation evaluate(std::span<Real const> x, uni20::DenseMatrix<Real>* jacobian = nullptr) const
    {
      Evaluation out{.residual = std::vector<Real>(order)};
      for (std::size_t i = 0; i < order; ++i)
      {
        Real const v = std::tanh(x[i]), c = Real{1} / std::cosh(x[i]);
        Real const multiplier = open ? Real{2} : Real{1};
        bethe::detail::CompensatedSum<Real> f, diagonal;
        f.add(multiplier * Real(sites) * std::atan(s * v) / s);
        diagonal.add(multiplier * Real(sites) * c * c / (Real{1} + s * s * v * v));
        if (open)
        {
          f.add(-Real{2} * std::atan(s / v) / s);
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
        out.norm = std::max(out.norm, std::abs(out.residual[i]));
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
  auto const m = sector_roots(sites, sz);
  NegativeGroundSystem<Real> const system(sites, m, delta, open);
  if (!uni20::isfinite(options.residual_tolerance) || options.residual_tolerance <= Real{0})
    throw std::invalid_argument("residual tolerance must be finite and positive");
  auto const elements = std::size_t(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(Real);
  if (m && m > elements / m) throw std::length_error("negative XXZ Newton matrix is too large");
  NegativeGroundRoots<Real> state;
  auto x = system.seed();
  uni20::DenseMatrix<Real> jacobian(m, m), step(m, 1);
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
    uni20::linalg::solve_inplace(jacobian, step);
    auto trial = x;
    bool accepted = false;
    Real damping = Real{1};
    for (int k = 0; k <= uni20::numeric_limits<Real>::digits; ++k)
    {
      for (std::size_t i = 0; i < m; ++i)
        trial[i] = x[i] + damping * step[i, 0];
      if (system.physical(trial))
      {
        Real const norm = system.evaluate(trial).norm;
        if (norm <= options.residual_tolerance || norm < (Real{1} - damping / Real{10000}) * evaluation.norm)
        {
          accepted = true;
          break;
        }
      }
      damping /= Real{2};
    }
    if (!accepted)
    {
      state.status = NegativeSolveStatus::stalled;
      break;
    }
    x = std::move(trial);
    ++state.iterations;
  }
  state.log_rapidities = x;
  bethe::detail::CompensatedSum<Real> energy;
  energy.add(Real(open ? sites - 1 : sites) * delta / Real{4});
  for (std::size_t i = 0; i < m; ++i)
  {
    Real const v = std::tanh(x[i]), c = Real{1} / std::cosh(x[i]);
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
