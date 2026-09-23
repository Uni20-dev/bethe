// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/detail/newton.hpp>
#include <bethe/xxz_common.hpp>

namespace bethe::xxz::quantum_group
{
enum class SolveStatus
{
  converged,
  iteration_limit,
  singular_jacobian,
  stalled
};

/// Even-chain ground state of the quantum-group-invariant open XXZ chain.
/// H=sum(sx*sx+sy*sy+Delta*sz*sz)+sqrt(Delta^2-1)/2*(sz_1-sz_N).
/// Spin-half operators, Delta>1. NOT the zero-boundary-field open XXZ model.
template <uni20::Real Real> struct GroundState
{
    std::size_t sites = 0;
    Real delta{};
    /// x=Theta(alpha;eta/2)/2, where Delta=cosh(eta), 0<x<pi/2.
    std::vector<Real> angles;
    /// Albertini's real alpha coordinates, 0<alpha<pi, not the zero-field z roots.
    std::vector<Real> rapidities;
    QuantumNumbers quantum_numbers;
    Real energy{}, residual_norm{};
    std::size_t iterations = 0;
    bool converged = false;
    SolveStatus status = SolveStatus::iteration_limit;
};

namespace detail
{
template <uni20::Real Real> class GroundSystem {
  public:
    GroundSystem(std::size_t sites, Real delta)
        : sites(sites), order(sites / 2), pi(Real{4} * std::atan(Real{1})),
          t(std::sqrt(((delta - Real{1}) / delta) / (Real{1} + Real{1} / delta))),
          tau(std::sqrt(((delta - Real{1}) / delta) * (Real{1} + Real{1} / delta)))
    {}

    std::vector<Real> seed() const
    {
      // Bare driving phases: a compact, ordered sea before scattering is included.
      // The Ising-limit sea is too spread out near Delta=1 on longer chains.
      std::vector<Real> x(order);
      for (std::size_t j = 0; j < order; ++j)
        x[j] = pi * Real(j + 1) / (Real{2} * Real(sites));
      return x;
    }

    bool physical(std::span<Real const> x) const
    {
      if (x.size() != order) return false;
      for (std::size_t j = 0; j < order; ++j)
        if (!uni20::isfinite(x[j]) || x[j] <= Real{0} || x[j] >= pi / Real{2} || (j && x[j] <= x[j - 1])) return false;
      return true;
    }

    Real alpha(Real x) const { return Real{2} * std::atan2(t * std::sin(x), std::cos(x)); }

    struct Evaluation
    {
        std::vector<Real> residual;
        Real norm{};
    };

    Evaluation evaluate(std::span<Real const> x, std::vector<Real>* jacobian = nullptr) const
    {
      if (jacobian) jacobian->resize(order * order);
      std::vector<Real> a(order), derivative(order);
      for (std::size_t j = 0; j < order; ++j)
      {
        a[j] = alpha(x[j]);
        Real const c = std::cos(x[j]), s = std::sin(x[j]);
        derivative[j] = t / (c * c + t * t * s * s); // (d alpha/d x)/2
      }
      Evaluation out{.residual = std::vector<Real>(order)};
      for (std::size_t i = 0; i < order; ++i)
      {
        bethe::detail::CompensatedSum<Real> sum, diagonal;
        sum.add(Real{2} * x[i] - pi * Real(i + 1) / Real(sites));
        diagonal.add(Real{2});
        for (std::size_t j = 0; j < order; ++j)
          if (i != j)
          {
            Real const minus = (a[i] - a[j]) / Real{2}, plus = (a[i] + a[j]) / Real{2};
            Real const sm = std::sin(minus), cm = std::cos(minus), sp = std::sin(plus), cp = std::cos(plus);
            // atan2 retains the reflected-scattering branch when alpha_i+alpha_j>pi.
            sum.add(-std::atan2(sm, tau * cm) / Real(sites));
            sum.add(-std::atan2(sp, tau * cp) / Real(sites));
            if (jacobian)
            {
              Real const km = tau / (sm * sm + tau * tau * cm * cm);
              Real const kp = tau / (sp * sp + tau * tau * cp * cp);
              diagonal.add(-(km + kp) * derivative[i] / Real(sites));
              (*jacobian)[i * order + j] = (km - kp) * derivative[j] / Real(sites);
            }
          }
        if (jacobian) (*jacobian)[i * order + i] = diagonal.value();
        out.residual[i] = sum.value();
        if (!uni20::isfinite(out.residual[i])) throw std::runtime_error("nonfinite quantum-group XXZ residual");
        out.norm = std::max(out.norm, std::abs(out.residual[i]));
      }
      return out;
    }

    std::size_t sites, order;
    Real pi, t, tau;
};
} // namespace detail

/// Filled sea I=1,...,N/2, even N>=2 and finite Delta>1.
/// Damped analytic-Jacobian Newton solve: O(N^2) storage, O(N^3) per update.
/// Residual is max|2N*Theta_1-sum(Theta_2^-+Theta_2^+)-2pi*I|/(2N).
/// A zero budget evaluates the seed; failed solves retain consistent roots/energy/residual.
template <uni20::Real Real = double>
[[nodiscard]] GroundState<Real> ground_state(std::size_t sites, Real delta, SolverOptions<Real> const& options = {})
{
  xxz::detail::checked_sites(sites);
  if (sites % 2) throw std::invalid_argument("quantum-group XXZ ground state currently requires even sites");
  if (!uni20::isfinite(delta) || delta <= Real{1})
    throw std::invalid_argument("quantum-group XXZ ground state requires finite Delta > 1");
  if (!uni20::isfinite(options.residual_tolerance) || options.residual_tolerance <= Real{0})
    throw std::invalid_argument("residual tolerance must be finite and positive");
  auto const m = sites / 2;
  auto const elements = std::size_t(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(Real);
  if (m > elements / m) throw std::length_error("quantum-group XXZ Newton matrix is too large");
  detail::GroundSystem<Real> const system(sites, delta);
  GroundState<Real> state;
  state.sites = sites;
  state.delta = delta;
  auto x = system.seed();
  std::vector<Real> jacobian, step(m);
  for (;;)
  {
    auto const evaluation = system.evaluate(x, &jacobian);
    state.residual_norm = evaluation.norm;
    if (evaluation.norm <= options.residual_tolerance)
    {
      state.converged = true;
      state.status = SolveStatus::converged;
      break;
    }
    if (state.iterations == options.max_iterations) break;
    for (std::size_t i = 0; i < m; ++i)
      step[i] = -evaluation.residual[i];
    if (!bethe::detail::newton_step(jacobian, step))
    {
      state.status = SolveStatus::singular_jacobian;
      break;
    }
    auto trial = x;
    Real damping = Real{1};
    bool accepted = false;
    for (int backtrack = 0; backtrack <= uni20::numeric_limits<Real>::digits; ++backtrack)
    {
      for (std::size_t i = 0; i < m; ++i)
        trial[i] = x[i] + damping * step[i];
      if (system.physical(trial))
      {
        auto const norm = system.evaluate(trial).norm;
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
      state.status = SolveStatus::stalled;
      break;
    }
    x = std::move(trial);
    ++state.iterations;
  }
  state.angles = std::move(x);
  bethe::detail::CompensatedSum<Real> energy;
  energy.add((Real(sites - 1) / Real{4} - Real(m)) * delta);
  for (std::size_t i = 0; i < m; ++i)
  {
    state.rapidities.push_back(system.alpha(state.angles[i]));
    state.quantum_numbers.emplace_back(std::int64_t(i + 1));
    energy.add(-std::cos(Real{2} * state.angles[i]));
  }
  state.energy = energy.value();
  if (!uni20::isfinite(state.energy)) throw std::overflow_error("quantum-group XXZ energy overflow");
  return state;
}
} // namespace bethe::xxz::quantum_group
