// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/eigenvalue_continuation.hpp>
#include <bethe/solver.hpp>
#include <span>
#include <stdexcept>

namespace bethe::richardson
{
template <uni20::Real Real> struct SolverOptions : bethe::SolverOptions<Real>
{
    using bethe::SolverOptions<Real>::SolverOptions;
    std::size_t max_stages = 10000; // Attempted continuation stages, including rejected ones.
};
enum class SolveStatus
{
  converged,
  iteration_limit,
  stage_limit,
  stalled,
  ill_conditioned
};

/// H=sum epsilon_i*(n_up+n_down)-g*sum_ij b_i^dagger*b_j, g>=0.
/// Distinct ascending single-particle levels, each a time-reversed doublet.
template <uni20::Real Real> struct State
{
    std::vector<Real> levels;
    std::vector<std::size_t> blocked, active;
    std::vector<Real> eigenvalue_variables; // y_i=g*sum_alpha 1/(2*epsilon_i-E_alpha); NOT occupations.
    std::size_t pairs = 0, iterations = 0, stages = 0, rejected_stages = 0;
    Real coupling{}, reached_coupling{}, energy{}, residual_norm{}, target_residual_norm{};
    Real particle_number_error{};
    bool converged = false;
    SolveStatus status = SolveStatus::iteration_limit;
};

namespace detail
{
// Enforce sum y=M by eliminating the lowest occupied variable y_0, preserving
// tiny empty-level variables at weak coupling. All N quadratic equations drive
// a rectangular QR Newton correction in N-1 variables. Repeated levels require
// additional derivative equations, not this system.
template <uni20::Real Real> class System {
  public:
    System(std::vector<Real> e, std::size_t m, uni20::Precision precision)
        : levels(std::move(e)), pairs(m), n(levels.size()), zero(uni20::scalar_at_precision<Real>(0, precision))
    {}
    System(std::vector<Real> e, std::size_t m)
      requires(!uni20::has_runtime_precision_v<Real>)
        : System(std::move(e), m, uni20::native_precision<Real>())
    {}
    std::vector<Real> expand(std::span<Real const> x) const
    {
      std::vector<Real> y(n, zero);
      bethe::detail::CompensatedSum<Real> last(zero);
      last.add(Real(pairs));
      for (Real v : x)
        last.add(-v);
      y[0] = last.value();
      std::copy(x.begin(), x.end(), y.begin() + 1);
      return y;
    }
    struct Evaluation
    {
        std::vector<Real> residual;
        Real norm{};
    };
    Evaluation evaluate(std::span<Real const> x, Real g, std::vector<Real>* jac = nullptr) const
    {
      auto const y = expand(x);
      Evaluation out{.residual = std::vector<Real>(n, zero), .norm = zero};
      if (jac) jac->assign(n * (n - 1), zero);
      for (std::size_t i = 0; i < n; ++i)
      {
        bethe::detail::CompensatedSum<Real> sum(zero), scale(zero), diagonal(zero);
        sum.add(y[i] * (y[i] - Real{1}));
        scale.add(Real{1} + uni20::math::abs(y[i] * (y[i] - Real{1})));
        diagonal.add(Real{2} * y[i] - Real{1});
        Real last = zero;
        for (std::size_t j = 0; j < n; ++j)
          if (i != j)
          {
            Real const a = g / (levels[i] - levels[j]), term = a * (y[i] - y[j]);
            sum.add(-term);
            scale.add(uni20::math::abs(a) * (uni20::math::abs(y[i]) + uni20::math::abs(y[j])));
            if (jac)
            {
              diagonal.add(-a);
              if (j > 0)
                (*jac)[i * (n - 1) + j - 1] = a;
              else
                last = a;
            }
          }
        Real const normalization = scale.value();
        if (!uni20::isfinite(sum.value()) || !uni20::isfinite(normalization))
          throw std::runtime_error("Richardson equations are not representable at this precision");
        out.norm = std::max(out.norm, uni20::math::abs(sum.value()) / normalization);
        out.residual[i] = sum.value() / normalization;
        if (jac)
        {
          if (i)
            (*jac)[i * (n - 1) + i - 1] = diagonal.value();
          else
            last = diagonal.value();
          for (std::size_t j = 0; j < n - 1; ++j)
            (*jac)[i * (n - 1) + j] = ((*jac)[i * (n - 1) + j] - last) / normalization;
        }
      }
      return out;
    }
    // Differentiate the UNscaled equations at a converged stage. The
    // row-scaled Jacobian must receive the same scaling on this RHS.
    bool tangent(std::span<Real const> x, Real g, std::vector<Real>& slope) const
    {
      auto const y = expand(x);
      std::vector<Real> jac;
      (void)evaluate(x, g, &jac);
      slope.assign(n, zero);
      for (std::size_t i = 0; i < n; ++i)
      {
        bethe::detail::CompensatedSum<Real> rhs(zero), scale(zero);
        scale.add(Real{1} + uni20::math::abs(y[i] * (y[i] - Real{1})));
        for (std::size_t j = 0; j < n; ++j)
          if (i != j)
          {
            Real const term = (y[i] - y[j]) / (levels[i] - levels[j]);
            rhs.add(term);
            scale.add(uni20::math::abs(g / (levels[i] - levels[j])) *
                      (uni20::math::abs(y[i]) + uni20::math::abs(y[j])));
          }
        slope[i] = rhs.value() / scale.value();
      }
      return bethe::detail::least_squares_step(jac, slope, n - 1, uni20::precision_of(zero));
    }
    Real energy(std::span<Real const> x, Real g) const
    {
      auto const y = expand(x);
      bethe::detail::CompensatedSum<Real> sum(zero);
      for (std::size_t i = 0; i < n; ++i)
        sum.add(levels[i] * y[i]);
      sum.add(-g * Real(pairs) * Real(n - pairs + 1));
      return sum.value();
    }
    std::vector<Real> const levels;
    std::size_t const pairs, n;
    Real const zero;
};
} // namespace detail

template <uni20::Real Real = double>
[[nodiscard]] State<Real> ground_state(std::span<Real const> levels, std::size_t pairs, Real coupling,
                                       std::span<std::size_t const> blocked = {},
                                       SolverOptions<Real> const& options = {})
{
  auto const zero = bethe::detail::solver_zero(options);
  bethe::detail::require_same_precision(coupling, zero);
  for (auto const& level : levels)
    bethe::detail::require_same_precision(level, zero);
  if (!uni20::isfinite(coupling) || coupling < Real{0})
    throw std::invalid_argument("Richardson requires finite g>=0 (attractive pairing)");
  if (!uni20::isfinite(options.residual_tolerance) || options.residual_tolerance <= Real{0})
    throw std::invalid_argument("Richardson residual tolerance must be finite and positive");
  if (levels.empty()) throw std::invalid_argument("Richardson requires at least one single-particle level");
  for (std::size_t i = 0; i < levels.size(); ++i)
    if (!uni20::isfinite(levels[i]) || (i && levels[i] <= levels[i - 1]))
      throw std::invalid_argument("Richardson single-particle levels must be finite, distinct and strictly increasing");
  std::vector<bool> is_blocked(levels.size(), false);
  for (auto i : blocked)
  {
    if (i >= levels.size() || is_blocked[i])
      throw std::invalid_argument("blocked indices must be distinct and within the level list");
    is_blocked[i] = true;
  }
  auto const n = levels.size() - blocked.size();
  if (pairs > n) throw std::invalid_argument("pair count exceeds the number of unblocked levels");
  State<Real> out;
  out.reached_coupling = out.energy = out.residual_norm = out.target_residual_norm = out.particle_number_error = zero;
  out.levels.assign(levels.begin(), levels.end());
  out.pairs = pairs;
  out.coupling = coupling;
  bethe::detail::CompensatedSum<Real> constant(zero);
  for (std::size_t i = 0; i < levels.size(); ++i)
    if (is_blocked[i])
    {
      out.blocked.push_back(i);
      constant.add(levels[i]);
    }
    else
      out.active.push_back(i);
  if (coupling == Real{0} || !pairs || pairs == n)
  {
    out.eigenvalue_variables.assign(n, zero);
    for (std::size_t j = 0; j < pairs; ++j)
    {
      constant.add(levels[out.active[j]]);
      constant.add(levels[out.active[j]]);
      out.eigenvalue_variables[j] = uni20::scalar_like(zero, 1);
    }
    constant.add(-coupling * Real(pairs));
    out.energy = constant.value();
    out.reached_coupling = coupling;
    out.converged = true;
    out.status = SolveStatus::converged;
    if (!uni20::isfinite(out.energy))
      throw std::runtime_error("Richardson energy is not representable at this precision");
    return out;
  }
  auto const elements = std::size_t(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(Real);
  if (n > elements / (n - 1)) throw std::length_error("Richardson Newton matrix is too large");
  Real const origin = levels[out.active.front()];
  Real const spread = levels[out.active.back()] - origin, scale = std::max(spread, coupling);
  if (!uni20::isfinite(scale))
    throw std::invalid_argument("Richardson level span is not representable at this precision");
  std::vector<Real> e(n, zero);
  Real gap = uni20::scalar_like(zero, 2);
  for (std::size_t i = 0; i < n; ++i)
  {
    e[i] = Real{2} * ((levels[out.active[i]] - origin) / scale);
    if (i)
    {
      if (!(e[i] > e[i - 1]))
        throw std::invalid_argument("Richardson level spacing is not representable at this precision");
      gap = std::min(gap, e[i] - e[i - 1]);
    }
  }
  Real const target = coupling / scale;
  if (!(target > Real{0}))
    throw std::invalid_argument("Richardson coupling scale is not representable at this precision");
  detail::System<Real> const system(std::move(e), pairs, options.precision());
  std::vector<Real> x(n - 1, zero);
  for (std::size_t i = 0; i < pairs - 1; ++i)
    x[i] = uni20::scalar_like(zero, 1);
  Real reached = zero, step = std::min(target, gap / Real{16});
  Real const trust = uni20::scalar_like(zero, 1) / Real{16};
  Real const eps = uni20::numeric_limits<Real>::epsilon(options.precision());
  auto predict = [&](std::vector<Real>& slope) -> std::optional<bethe::detail::ContinuationPrediction<Real>> {
    if (!system.tangent(x, reached, slope)) return std::nullopt;
    Real max_slope = zero;
    bethe::detail::CompensatedSum<Real> last(zero);
    for (Real v : slope)
    {
      max_slope = std::max(max_slope, uni20::math::abs(v));
      last.add(v);
    }
    max_slope = std::max(max_slope, uni20::math::abs(last.value()));
    return bethe::detail::ContinuationPrediction<Real>{trust, max_slope};
  };
  auto correct = [&](std::vector<Real>& candidate, Real next, Real step, std::vector<Real> const& slope,
                     Real trust_radius) {
    bool accepted = bethe::detail::correct_eigenvalue_stage(
        candidate, trust_radius, options, out.iterations,
        [&](auto const& trial, std::vector<Real>* jac) { return system.evaluate(trial, next, jac); },
        [&](auto const& trial) { return system.expand(trial); });
    if (accepted)
    {
      // Necessary physical bounds: kinetic minimum plus the maximal pair
      // interaction; variational Fermi-sea and uniform-pair (Dicke) states.
      bethe::detail::CompensatedSum<Real> free(zero), sum(zero);
      for (std::size_t i = 0; i < n; ++i)
      {
        sum.add(system.levels[i]);
        if (i < pairs) free.add(system.levels[i]);
      }
      Real const interaction = next * Real(pairs) * Real(n - pairs + 1);
      Real const lower = free.value() - interaction;
      Real const upper = std::min(free.value() - next * Real(pairs), Real(pairs) / Real(n) * sum.value() - interaction);
      Real const energy = system.energy(candidate, next);
      Real const allowance = Real{512} * Real(n) * (eps + options.residual_tolerance) *
                             (Real{1} + uni20::math::abs(lower) + uni20::math::abs(upper));
      accepted = uni20::isfinite(energy) && energy >= lower - allowance && energy <= upper + allowance;
      // The ground energy is concave in g. Its tangent is an upper bound,
      // while the norm of S+S- bounds how far it can fall in one stage.
      bethe::detail::CompensatedSum<Real> derivative(zero);
      for (std::size_t j = 1; j < n; ++j)
        derivative.add(system.levels[j] * slope[j - 1]);
      derivative.add(-Real(pairs) * Real(n - pairs + 1));
      Real const old_energy = system.energy(x, reached);
      accepted = accepted && energy <= old_energy + step * derivative.value() + allowance &&
                 energy >= old_energy - step * Real(pairs) * Real(n - pairs + 1) - allowance;
    }
    return accepted;
  };
  bethe::detail::continue_eigenvalues(x, reached, step, target, options, out, predict, correct);
  out.reached_coupling = reached == target ? coupling : reached * scale;
  out.eigenvalue_variables = system.expand(x);
  out.residual_norm = system.evaluate(x, reached).norm;
  out.target_residual_norm = system.evaluate(x, target).norm;
  bethe::detail::CompensatedSum<Real> number(zero);
  for (Real v : out.eigenvalue_variables)
    number.add(v);
  out.particle_number_error = uni20::math::abs(number.value() - Real(pairs));
  constant.add(Real{2} * (origin * Real(pairs)));
  constant.add(scale * system.energy(x, reached));
  out.energy = constant.value();
  if (!uni20::isfinite(out.energy))
    throw std::runtime_error("Richardson energy is not representable at this precision");
  if (reached == target)
  {
    out.converged = true;
    out.status = SolveStatus::converged;
  }
  else if (out.iterations == options.max_iterations)
    out.status = SolveStatus::iteration_limit;
  else if (out.stages == options.max_stages)
    out.status = SolveStatus::stage_limit;
  return out;
}
} // namespace bethe::richardson
