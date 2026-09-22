// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <algorithm>
#include <array>
#include <bethe/solver.hpp>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <uni20/common/half_int.hpp>
#include <uni20/linalg/ops/linear_solve.hpp>
#include <vector>

namespace bethe::su3
{
using QuantumNumbers = std::array<std::vector<uni20::half_int>, 2>;
template <uni20::Real Real> using Rapidities = std::array<std::vector<Real>, 2>;
template <uni20::Real Real> using SolverOptions = bethe::SolverOptions<Real>;
enum class SolveStatus
{
  converged,
  iteration_limit,
  stalled
};

/// Fundamental SU(3) permutation chain, H=sum_j P_(j,j+1), J=1, periodic.
/// The currently supported state is the balanced singlet for L divisible by 3.
template <uni20::Real Real> struct State
{
    std::size_t sites = 0;
    std::array<std::size_t, 3> populations{};
    QuantumNumbers quantum_numbers;
    /// Conventional lambda and mu, NOT the XXX front end's z=2*lambda.
    Rapidities<Real> rapidities;
    Real energy = Real{0}, momentum = Real{0};
    std::size_t momentum_index = 0;
    std::array<Real, 2> level_residuals{};
    Real residual_norm = Real{0};
    std::size_t iterations = 0;
    bool converged = false;
    SolveStatus status = SolveStatus::iteration_limit;
};

namespace detail
{
inline void check_sites(std::size_t sites)
{
  if (sites < 3 || sites % 3 || sites > std::size_t(std::numeric_limits<std::int64_t>::max() / 4))
    throw std::invalid_argument("SU(3) balanced ground state requires L>=3 divisible by 3, with L<=INT64_MAX/4");
}
inline QuantumNumbers centered_numbers(std::array<std::size_t, 2> counts)
{
  QuantumNumbers labels;
  for (std::size_t a = 0; a < 2; ++a)
  {
    labels[a].resize(counts[a]);
    for (std::size_t j = 0; j < counts[a]; ++j)
      labels[a][j] = uni20::from_twice(2 * std::int64_t(j) - std::int64_t(counts[a] - 1));
  }
  return labels;
}

/// The two-level permutation-chain equations, with the second sea optionally
/// empty for the SU(2) reduction test. This is not a generic Hubbard system.
template <uni20::Real Real> class GroundSystem {
  public:
    GroundSystem(std::size_t sites, std::array<std::size_t, 2> counts)
        : sites(sites), counts(counts), positive{counts[0] / 2, counts[1] / 2}, order(positive[0] + positive[1]),
          labels(centered_numbers(counts)), pi(Real{4} * std::atan(Real{1}))
    {}

    std::vector<Real> seed() const
    {
      // Invert the filled-sea thermodynamic density at I/L. This is only a
      // starting guess; all finite-size equations are subsequently solved.
      using std::atanh;
      using std::tan;
      std::vector<Real> x(order);
      Real const colors = counts[1] ? Real{3} : Real{2};
      for (std::size_t a = 0; a < 2; ++a)
        for (std::size_t j = 0; j < positive[a]; ++j)
        {
          Real const i = Real(labels[a][first(a) + j].twice()) / Real{2};
          x[offset(a) + j] =
              colors / pi * atanh(tan(pi * i / Real(sites)) * tan(pi * Real(a + 1) / (Real{2} * colors)));
        }
      return x;
    }
    bool physical(std::span<Real const> x) const
    {
      if (x.size() != order) return false;
      for (std::size_t a = 0; a < 2; ++a)
        for (std::size_t j = 0; j < positive[a]; ++j)
        {
          Real const value = x[offset(a) + j];
          if (!uni20::isfinite(value) || value <= Real{0} || value > uni20::numeric_limits<Real>::max() / Real{4} ||
              (j && value <= x[offset(a) + j - 1]))
            return false;
        }
      return true;
    }
    Rapidities<Real> expand(std::span<Real const> x) const
    {
      Rapidities<Real> roots;
      for (std::size_t a = 0; a < 2; ++a)
      {
        roots[a].reserve(counts[a]);
        for (std::size_t j = positive[a]; j-- > 0;)
          roots[a].push_back(-x[offset(a) + j]);
        if (counts[a] % 2) roots[a].push_back(Real{0});
        for (std::size_t j = 0; j < positive[a]; ++j)
          roots[a].push_back(x[offset(a) + j]);
      }
      return roots;
    }
    static Real phase(Real x, Real width) { return Real{2} * std::atan(Real{2} * x / width); }
    static Real derivative(Real x, Real width)
    {
      using std::abs;
      Real const a = width / Real{2};
      if (abs(x) > Real{1})
      {
        Real const r = a / x;
        return (Real{2} * r / x) / (Real{1} + r * r);
      }
      return Real{2} * a / (a * a + x * x);
    }
    struct Evaluation
    {
        std::vector<Real> residual;
        std::array<Real, 2> levels{};
        Real norm() const { return std::max(levels[0], levels[1]); }
    };
    Evaluation evaluate(std::span<Real const> x, uni20::DenseMatrix<Real>* jacobian = nullptr) const
    {
      using std::abs;
      auto const roots = expand(x);
      Evaluation result{.residual = std::vector<Real>(order)};
      if (jacobian)
        for (std::size_t i = 0; i < order; ++i)
          for (std::size_t j = 0; j < order; ++j)
            (*jacobian)[i, j] = Real{0};
      for (std::size_t a = 0; a < 2; ++a)
        for (std::size_t j = 0; j < positive[a]; ++j)
        {
          auto const row = offset(a) + j, full = first(a) + j;
          Real const value = x[row];
          bethe::detail::CompensatedSum<Real> sum, diagonal;
          sum.add(-pi * Real(labels[a][full].twice()) / Real(sites));
          if (a == 0)
          {
            sum.add(phase(value, Real{1}));
            diagonal.add(derivative(value, Real{1}));
          }
          for (std::size_t b = 0; b < 2; ++b)
          {
            Real const sign = a == b ? -Real{1} : Real{1}, width = a == b ? Real{2} : Real{1};
            for (std::size_t k = 0; k < counts[b]; ++k)
            {
              if (a == b && k == full) continue;
              Real const difference = value - roots[b][k];
              sum.add(sign * phase(difference, width) / Real(sites));
              if (jacobian)
              {
                Real const d = sign * derivative(difference, width) / Real(sites);
                diagonal.add(d);
                if (k < positive[b])
                  (*jacobian)[row, offset(b) + positive[b] - 1 - k] += d;
                else if (k >= first(b))
                  (*jacobian)[row, offset(b) + k - first(b)] -= d;
              }
            }
          }
          if (jacobian) (*jacobian)[row, row] += diagonal.value();
          result.residual[row] = sum.value();
          if (!uni20::isfinite(result.residual[row])) throw std::runtime_error("nonfinite SU(3) equation residual");
          result.levels[a] = std::max(result.levels[a], abs(result.residual[row]));
        }
      return result;
    }
    std::size_t offset(std::size_t a) const { return a ? positive[0] : 0; }
    std::size_t first(std::size_t a) const { return positive[a] + counts[a] % 2; }
    std::size_t const sites;
    std::array<std::size_t, 2> const counts, positive;
    std::size_t const order;
    QuantumNumbers const labels;
    Real const pi;
};

template <uni20::Real Real> void check_options(std::size_t sites, SolverOptions<Real> const& options)
{
  check_sites(sites);
  if (!uni20::isfinite(options.residual_tolerance) || options.residual_tolerance <= Real{0})
    throw std::invalid_argument("residual tolerance must be finite and positive");
  auto const order = sites / 3 + (sites / 3) / 2;
  auto const elements = std::size_t(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(Real);
  if (order > elements / order) throw std::length_error("SU(3) Newton matrix is too large");
}

template <uni20::Real Real> State<Real> solve(GroundSystem<Real> const& system, SolverOptions<Real> const& options)
{
  auto x = system.seed();
  if (!system.physical(x)) throw std::runtime_error("SU(3) seed cannot be represented at this precision");
  State<Real> result;
  result.sites = system.sites;
  result.populations = {system.sites - system.counts[0], system.counts[0] - system.counts[1], system.counts[1]};
  result.quantum_numbers = system.labels;
  uni20::DenseMatrix<Real> jacobian(system.order, system.order), step(system.order, 1);
  for (;;)
  {
    auto const evaluation = system.evaluate(x, &jacobian);
    result.level_residuals = evaluation.levels;
    result.residual_norm = evaluation.norm();
    if (result.residual_norm <= options.residual_tolerance)
    {
      result.converged = true;
      result.status = SolveStatus::converged;
      break;
    }
    if (result.iterations == options.max_iterations) break;
    for (std::size_t j = 0; j < system.order; ++j)
      step[j, 0] = -evaluation.residual[j];
    uni20::linalg::solve_inplace(jacobian, step);
    bool accepted = false;
    Real damping = Real{1};
    auto trial = x;
    for (int backtrack = 0; backtrack <= uni20::numeric_limits<Real>::digits; ++backtrack)
    {
      for (std::size_t j = 0; j < system.order; ++j)
        trial[j] = x[j] + damping * step[j, 0];
      if (system.physical(trial))
      {
        Real const norm = system.evaluate(trial).norm();
        if (norm <= options.residual_tolerance || norm < (Real{1} - damping / Real{10000}) * result.residual_norm)
        {
          accepted = true;
          break;
        }
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
  result.rapidities = system.expand(x);
  bethe::detail::CompensatedSum<Real> energy;
  energy.add(Real(system.sites));
  for (Real lambda : result.rapidities[0])
    energy.add(-GroundSystem<Real>::derivative(lambda, Real{1}));
  result.energy = energy.value();
  if (!uni20::isfinite(result.energy)) throw std::overflow_error("SU(3) energy overflow");
  // Reflection pairs contribute 2*pi; an unpaired lambda=0 contributes pi.
  // M1=2L/3 is even in the supported SU(3) singlet, hence P=0 exactly.
  result.momentum_index = system.counts[0] % 2 ? system.sites / 2 : 0;
  result.momentum = system.counts[0] % 2 ? system.pi : Real{0};
  return result;
}
} // namespace detail

inline QuantumNumbers ground_quantum_numbers(std::size_t sites)
{
  detail::check_sites(sites);
  return detail::centered_numbers({2 * (sites / 3), sites / 3});
}

/// Balanced singlet at L>=3 divisible by 3. Arbitrary color sectors, excited
/// states, and complex roots require additional state-selection machinery.
template <uni20::Real Real> State<Real> ground_state(std::size_t sites, SolverOptions<Real> const& options = {})
{
  detail::check_options(sites, options);
  return detail::solve(detail::GroundSystem<Real>(sites, {2 * (sites / 3), sites / 3}), options);
}
} // namespace bethe::su3
