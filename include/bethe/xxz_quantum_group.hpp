// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/detail/newton.hpp>
#include <bethe/xxz_common.hpp>

namespace bethe::xxz::quantum_group
{
using xxz::QuantumNumbers;
enum class SolveStatus
{
  converged,
  iteration_limit,
  singular_jacobian,
  stalled
};

/// Positive finite-real-root state of the quantum-group-invariant open XXZ chain.
/// H=sum(sx*sx+sy*sy+Delta*sz*sz)+sqrt(Delta^2-1)/2*(sz_1-sz_N).
/// Spin-half operators, Delta>1. NOT the zero-boundary-field open XXZ model.
template <uni20::Real Real> struct RealState
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
    /// E-(N-1)*Delta/4, evaluated directly rather than subtracting extensive energies.
    Real energy_shift{};
};
template <uni20::Real Real> using GroundState = RealState<Real>;

namespace detail
{
inline std::size_t real_sector_roots(std::size_t sites, std::size_t through_lines)
{
  xxz::detail::checked_sites(sites);
  if (through_lines > sites || (sites - through_lines) % 2)
    throw std::invalid_argument("TL through-lines must lie in [0,N] with the same parity as N");
  return (sites - through_lines) / 2;
}

inline std::size_t sector_roots(std::size_t sites, std::size_t through_lines)
{
  xxz::detail::checked_sites(sites);
  if (sites % 2) throw std::invalid_argument("quantum-group XXZ solver currently requires even sites");
  if (through_lines > sites || through_lines % 2)
    throw std::invalid_argument("TL through-lines must be even and lie in [0,N]");
  return (sites - through_lines) / 2;
}

template <uni20::Real Real> void check_matrix_size(std::size_t m)
{
  auto const elements = std::size_t(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(Real);
  if (m && m > elements / m) throw std::length_error("quantum-group XXZ Newton matrix is too large");
}

inline QuantumNumbers consecutive(std::size_t m)
{
  QuantumNumbers numbers;
  numbers.reserve(m);
  for (std::size_t i = 0; i < m; ++i)
    numbers.emplace_back(std::int64_t(i + 1));
  return numbers;
}

template <uni20::Real Real> class GroundSystem {
  public:
    GroundSystem(std::size_t sites, Real delta, std::span<uni20::half_int const> numbers)
        : sites(sites), order(numbers.size()), numbers(numbers.begin(), numbers.end()),
          pi(Real{4} * std::atan(Real{1})), t(std::sqrt(((delta - Real{1}) / delta) / (Real{1} + Real{1} / delta))),
          tau(std::sqrt(((delta - Real{1}) / delta) * (Real{1} + Real{1} / delta)))
    {}
    GroundSystem(std::size_t sites, Real delta) : GroundSystem(sites, delta, consecutive(sites / 2)) {}

    std::vector<Real> seed() const
    {
      // Bare driving phases: a compact, ordered sea before scattering is included.
      // The Ising-limit sea is too spread out near Delta=1 on longer chains.
      std::vector<Real> x(order);
      for (std::size_t j = 0; j < order; ++j)
        x[j] = pi * Real(numbers[j].twice() / 2) / (Real{2} * Real(sites));
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
        sum.add(Real{2} * x[i] - pi * Real(numbers[i].twice() / 2) / Real(sites));
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
    QuantumNumbers numbers;
    Real pi, t, tau;
};
} // namespace detail

/// Odd/even N>=2, finite Delta>1, M<=floor(N/2), integer labels 1<=I<=N-M.
/// This is the positive finite-real family, NOT the complete TL module spectrum.
/// Its regular Bethe states belong to the module with ell=N-2*M through-lines.
/// Damped analytic-Jacobian Newton solve: O(M^2) storage, O(M^3) per update.
/// Residual is max|2N*Theta_1-sum(Theta_2^-+Theta_2^+)-2pi*I|/(2N).
/// A zero budget evaluates the seed; failed solves retain consistent roots/energy/residual.
template <uni20::Real Real = double>
[[nodiscard]] RealState<Real> solve_real(std::size_t sites, Real delta, std::span<uni20::half_int const> numbers,
                                         SolverOptions<Real> const& options = {})
{
  xxz::detail::checked_sites(sites);
  if (!uni20::isfinite(delta) || delta <= Real{1})
    throw std::invalid_argument("quantum-group XXZ solver requires finite Delta > 1");
  if (!uni20::isfinite(options.residual_tolerance) || options.residual_tolerance <= Real{0})
    throw std::invalid_argument("residual tolerance must be finite and positive");
  auto const m = numbers.size();
  if (m > sites / 2) throw std::invalid_argument("quantum-group XXZ real states require M <= N/2");
  for (std::size_t i = 0; i < m; ++i)
    if (numbers[i].twice() % 2 || numbers[i].twice() <= 0 || numbers[i].twice() > 2 * std::int64_t(sites - m) ||
        (i && numbers[i] <= numbers[i - 1]))
      throw std::invalid_argument("quantum-group XXZ labels must be increasing integers in [1,N-M]");
  detail::check_matrix_size<Real>(m);
  detail::GroundSystem<Real> const system(sites, delta, numbers);
  RealState<Real> state;
  state.sites = sites;
  state.delta = delta;
  state.quantum_numbers.assign(numbers.begin(), numbers.end());
  auto x = system.seed();
  if (!system.physical(x)) throw std::overflow_error("real-root seed is not resolvable at the selected precision");
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
  bethe::detail::CompensatedSum<Real> energy, shift;
  energy.add((Real(sites - 1) / Real{4} - Real(m)) * delta);
  for (std::size_t i = 0; i < m; ++i)
  {
    state.rapidities.push_back(system.alpha(state.angles[i]));
    energy.add(-std::cos(Real{2} * state.angles[i]));
    Real const c = std::cos(state.angles[i]);
    shift.add(-(delta - Real{1}));
    shift.add(-Real{2} * c * c);
  }
  state.energy = energy.value();
  state.energy_shift = shift.value();
  if (!uni20::isfinite(state.energy) || !uni20::isfinite(state.energy_shift))
    throw std::overflow_error("quantum-group XXZ energy overflow");
  return state;
}

/// Lowest state in the even-chain TL module ell=through_lines: I=1,...,(N-ell)/2.
/// ell/2 is auxiliary quantum-group spin, NOT the physical spin-1 total spin.
template <uni20::Real Real = double>
[[nodiscard]] RealState<Real> sector_ground_state(std::size_t sites, Real delta, std::size_t through_lines,
                                                  SolverOptions<Real> const& options = {})
{
  auto const m = detail::sector_roots(sites, through_lines);
  detail::check_matrix_size<Real>(m);
  return solve_real<Real>(sites, delta, detail::consecutive(m), options);
}

template <uni20::Real Real = double>
[[nodiscard]] GroundState<Real> ground_state(std::size_t sites, Real delta, SolverOptions<Real> const& options = {})
{
  return sector_ground_state(sites, delta, 0, options);
}
} // namespace bethe::xxz::quantum_group
