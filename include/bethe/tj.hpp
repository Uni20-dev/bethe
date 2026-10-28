// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <array>
#include <bethe/detail/newton.hpp>
#include <bethe/heisenberg.hpp>

namespace bethe::tj
{
using QuantumNumbers = std::array<std::vector<uni20::half_int>, 2>;
template <uni20::Real Real> using Rapidities = std::array<std::vector<Real>, 2>;
template <uni20::Real Real> using SolverOptions = bethe::SolverOptions<Real>;
enum class SolveStatus
{
  converged,
  iteration_limit,
  stalled,
  ill_conditioned
};
enum class Branch
{
  sutherland,
  polarized_free,
  no_holes_xxx
};

/// Periodic projected fermions, t=1, J=2, with no chemical-potential shift.
/// Mixed doped sectors require odd populations of both spins. Vacuum/full
/// polarization and the no-hole XXX reduction support all physical counts.
template <uni20::Real Real> struct State
{
    std::size_t sites = 0, num_up = 0, num_down = 0, holes = 0;
    Branch branch = Branch::sutherland;
    QuantumNumbers quantum_numbers;
    Rapidities<Real> rapidities;          // Conventional lambda, mu; NOT XXX z=2*lambda.
    std::vector<std::int64_t> free_modes; // k=2*pi*j/L; only polarized_free.
    Real energy{}, momentum{}, residual_norm{};
    std::array<Real, 2> level_residuals{};
    std::size_t momentum_index = 0, iterations = 0;
    bool converged = false, spin_reversed = false;
    SolveStatus status = SolveStatus::iteration_limit;
};

namespace detail
{
inline void check_counts(std::size_t n, std::size_t up, std::size_t down)
{
  if (n < 3 || n > std::size_t(std::numeric_limits<std::int64_t>::max() / 4))
    throw std::invalid_argument("periodic t-J requires 3 <= L <= INT64_MAX/4");
  if (up > n || down > n - up) throw std::invalid_argument("t-J populations must satisfy N_up+N_down <= L");
}

inline QuantumNumbers centered_numbers(std::array<std::size_t, 2> counts)
{
  QuantumNumbers labels;
  for (std::size_t a = 0; a < 2; ++a)
    for (std::size_t j = 0; j < counts[a]; ++j)
      labels[a].push_back(uni20::from_twice(2 * std::int64_t(j) - std::int64_t(counts[a] - 1)));
  return labels;
}

using bethe::detail::newton_step;

/// Sutherland (BFF) equations, Essler-Korepin (3.73), conventional rapidities.
/// M1=N_h+min(N_up,N_down), M2=N_h; the second level has NO self-scattering.
template <uni20::Real Real> class GroundSystem {
  public:
    GroundSystem(std::size_t n, std::size_t up, std::size_t down)
      requires(!uni20::has_runtime_precision_v<Real>)
        : GroundSystem(n, up, down, Real{0})
    {}
    GroundSystem(std::size_t n, std::size_t up, std::size_t down, Real const& exemplar)
        : sites(n), counts{n - std::max(up, down), n - up - down}, positive{counts[0] / 2, counts[1] / 2},
          order(positive[0] + positive[1]), labels(centered_numbers(counts)),
          zero(bethe::detail::precision_zero(exemplar)), pi(Real{4} * uni20::math::atan(uni20::scalar_like(zero, 1)))
    {}

    std::size_t offset(std::size_t a) const { return a ? positive[0] : 0; }
    std::size_t first(std::size_t a) const { return positive[a] + counts[a] % 2; }
    static Real phase(Real x, Real width) { return Real{2} * uni20::math::atan(Real{2} * x / width); }
    static Real derivative(Real x, Real width)
    {
      Real const a = width / Real{2};
      if (uni20::math::abs(x) > Real{1})
      {
        Real const r = a / x;
        return (Real{2} * r / x) / (Real{1} + r * r);
      }
      return Real{2} * a / (a * a + x * x);
    }
    std::vector<Real> seed() const
    {
      std::vector<Real> x(order, zero);
      for (std::size_t a = 0; a < 2; ++a)
        for (std::size_t j = 0; j < positive[a]; ++j)
          x[offset(a) + j] =
              uni20::math::tan(pi * Real(labels[a][first(a) + j].twice()) / (Real{2} * Real(a ? counts[0] : sites))) /
              Real{2};
      return x;
    }
    bool physical(std::span<Real const> x) const
    {
      if (x.size() != order) return false;
      Real const limit = [&] {
        if constexpr (uni20::has_runtime_precision_v<Real>)
          return uni20::math::next_down(uni20::scalar_like(zero, std::numeric_limits<double>::infinity())) / Real{8};
        else
          return uni20::numeric_limits<Real>::max() / Real{8};
      }();
      for (std::size_t a = 0; a < 2; ++a)
        for (std::size_t j = 0; j < positive[a]; ++j)
        {
          Real const v = x[offset(a) + j];
          if (!uni20::isfinite(v) || v <= Real{0} || v > limit || (j && v <= x[offset(a) + j - 1])) return false;
        }
      return true;
    }
    Rapidities<Real> expand(std::span<Real const> x) const
    {
      Rapidities<Real> roots;
      for (std::size_t a = 0; a < 2; ++a)
      {
        for (std::size_t j = positive[a]; j-- > 0;)
          roots[a].push_back(-x[offset(a) + j]);
        if (counts[a] % 2) roots[a].push_back(zero);
        for (std::size_t j = 0; j < positive[a]; ++j)
          roots[a].push_back(x[offset(a) + j]);
      }
      return roots;
    }
    struct Evaluation
    {
        std::vector<Real> residual;
        std::array<Real, 2> levels{};
        Real norm() const { return std::max(levels[0], levels[1]); }
    };
    Evaluation evaluate(std::span<Real const> x, std::vector<Real>* jac = nullptr) const
    {
      auto const roots = expand(x);
      Evaluation out{.residual = std::vector<Real>(order, zero), .levels = {zero, zero}};
      if (jac) jac->assign(order * order, zero);
      for (std::size_t a = 0; a < 2; ++a)
        for (std::size_t j = 0; j < positive[a]; ++j)
        {
          auto const row = offset(a) + j, full = first(a) + j;
          Real const value = x[row];
          bethe::detail::CompensatedSum<Real> sum(zero), diagonal(zero);
          sum.add(-pi * Real(labels[a][full].twice()) / Real(sites));
          if (a == 0)
          {
            sum.add(phase(value, Real{1}));
            diagonal.add(derivative(value, Real{1}));
          }
          for (std::size_t b = 0; b < 2; ++b)
          {
            if (a == 1 && b == 1) continue;
            Real const sign = a == b ? -Real{1} : Real{1}, width = a == b ? Real{2} : Real{1};
            for (std::size_t k = 0; k < counts[b]; ++k)
            {
              if (a == b && k == full) continue;
              Real const diff = value - roots[b][k];
              sum.add(sign * phase(diff, width) / Real(sites));
              if (jac)
              {
                Real const d = sign * derivative(diff, width) / Real(sites);
                diagonal.add(d);
                if (k < positive[b])
                  (*jac)[row * order + offset(b) + positive[b] - 1 - k] += d;
                else if (k >= first(b))
                  (*jac)[row * order + offset(b) + k - first(b)] -= d;
              }
            }
          }
          if (jac) (*jac)[row * order + row] += diagonal.value();
          out.residual[row] = sum.value();
          if (!uni20::isfinite(out.residual[row])) throw std::runtime_error("nonfinite t-J equation residual");
          out.levels[a] = std::max(out.levels[a], uni20::math::abs(out.residual[row]));
        }
      return out;
    }
    std::size_t const sites;
    std::array<std::size_t, 2> const counts, positive;
    std::size_t const order;
    QuantumNumbers const labels;
    Real const zero, pi;
};
} // namespace detail

template <uni20::Real Real = double>
[[nodiscard]] State<Real> ground_state(std::size_t sites, std::size_t num_up, std::size_t num_down,
                                       SolverOptions<Real> const& options = {})
{
  detail::check_counts(sites, num_up, num_down);
  Real const zero = bethe::detail::solver_zero(options);
  auto const precision = options.precision();
  if (!uni20::isfinite(options.residual_tolerance) || options.residual_tolerance <= Real{0})
    throw std::invalid_argument("t-J residual tolerance must be finite and positive");
  State<Real> out;
  out.energy = out.momentum = out.residual_norm = zero;
  out.level_residuals = {zero, zero};
  out.sites = sites;
  out.num_up = num_up;
  out.num_down = num_down;
  out.holes = sites - num_up - num_down;
  out.spin_reversed = num_down > num_up;
  Real const pi = Real{4} * uni20::math::atan(uni20::scalar_like(zero, 1));
  if (num_up == 0 || num_down == 0)
  {
    out.branch = Branch::polarized_free;
    auto const particles = num_up + num_down;
    for (std::size_t j = 0; j < particles; ++j)
      out.free_modes.push_back(std::int64_t(j) - std::int64_t((particles - 1) / 2));
    if (particles && particles < sites)
      out.energy = -Real{2} * uni20::math::sin(pi * (uni20::scalar_like(zero, particles) / Real(sites))) /
                   uni20::math::sin(pi / Real(sites)) * (particles % 2 ? Real{1} : uni20::math::cos(pi / Real(sites)));
    out.momentum_index = particles % 2 ? 0 : particles / 2;
    out.converged = true;
    out.status = SolveStatus::converged;
  }
  else if (!out.holes)
  {
    out.branch = Branch::no_holes_xxx;
    auto const sz = uni20::from_twice(std::int64_t(num_up) - std::int64_t(num_down));
    auto spin = heisenberg::sector_ground_state<Real>(sites, sz, options);
    out.quantum_numbers[0] = std::move(spin.quantum_numbers);
    out.rapidities[0] = std::move(spin.rapidities);
    for (auto& v : out.rapidities[0])
      v /= Real{2};
    // Algebraically 2*E_XXX-L/2, but do not cancel two extensive constants:
    // a nearly polarized t-J energy can remain O(1) on a very long chain.
    bethe::detail::CompensatedSum<Real> energy(zero);
    for (Real v : out.rapidities[0])
      energy.add(-detail::GroundSystem<Real>::derivative(v, Real{1}));
    out.energy = energy.value();
    out.iterations = spin.iterations;
    out.residual_norm = spin.residual_norm;
    out.level_residuals[0] = spin.residual_norm;
    out.converged = spin.converged;
    out.status = spin.converged ? SolveStatus::converged : SolveStatus::iteration_limit;
    // Fermionic translation of the filled reference contributes (-1)^(L-1).
    out.momentum_index = (spin.momentum_index + (sites % 2 ? 0 : sites / 2)) % sites;
  }
  else
  {
    if (num_up % 2 == 0 || num_down % 2 == 0)
      throw std::invalid_argument("doped mixed-spin t-J ground states currently require odd N_up AND odd N_down");
    auto const m = sites - std::max(num_up, num_down), order = m / 2 + out.holes / 2;
    auto const elements = std::size_t(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(Real);
    if (order && order > elements / order) throw std::length_error("t-J Newton matrix is too large");
    detail::GroundSystem<Real> const system(sites, num_up, num_down, zero);
    out.quantum_numbers = system.labels;
    auto x = system.seed();
    if (!system.physical(x)) throw std::runtime_error("t-J seed is not representable at this precision");
    std::vector<Real> jac;
    for (;;)
    {
      auto const evaluation = system.evaluate(x, &jac);
      out.level_residuals = evaluation.levels;
      out.residual_norm = evaluation.norm();
      if (out.residual_norm <= options.residual_tolerance)
      {
        out.converged = true;
        out.status = SolveStatus::converged;
        break;
      }
      if (out.iterations == options.max_iterations) break;
      auto step = evaluation.residual;
      for (auto& v : step)
        v = -v;
      if (!detail::newton_step(std::move(jac), step, precision))
      {
        out.status = SolveStatus::ill_conditioned;
        break;
      }
      bool accepted = false;
      Real damping = uni20::scalar_like(zero, 1);
      auto trial = x;
      for (int backtrack = 0; backtrack <= uni20::numeric_limits<Real>::digits(precision); ++backtrack)
      {
        for (std::size_t j = 0; j < x.size(); ++j)
          trial[j] = x[j] + damping * step[j];
        if (system.physical(trial))
        {
          Real const norm = system.evaluate(trial).norm();
          if (norm <= options.residual_tolerance || norm < (Real{1} - damping / Real{10000}) * out.residual_norm)
          {
            accepted = true;
            break;
          }
        }
        damping /= Real{2};
      }
      if (!accepted)
      {
        out.status = SolveStatus::stalled;
        break;
      }
      x = std::move(trial);
      ++out.iterations;
    }
    out.rapidities = system.expand(x);
    bethe::detail::CompensatedSum<Real> energy(zero);
    energy.add(Real{2} * Real(out.holes));
    for (Real v : out.rapidities[0])
      energy.add(-detail::GroundSystem<Real>::derivative(v, Real{1}));
    out.energy = energy.value();
    // The filled reference phase plus M1*pi gives zero for odd N_up,N_down.
    out.momentum_index = 0;
  }
  out.momentum = Real{2} * pi * (Real(out.momentum_index) / Real(sites));
  return out;
}
} // namespace bethe::tj
