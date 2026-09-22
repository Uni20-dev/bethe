// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/newton.hpp>
#include <bethe/solver.hpp>
#include <span>
#include <stdexcept>

namespace bethe::richardson
{
template <uni20::Real Real> struct SolverOptions : bethe::SolverOptions<Real>
{
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
    System(std::vector<Real> e, std::size_t m) : levels(std::move(e)), pairs(m), n(levels.size()) {}
    std::vector<Real> expand(std::span<Real const> x) const
    {
      std::vector<Real> y(n);
      bethe::detail::CompensatedSum<Real> last;
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
      Evaluation out{.residual = std::vector<Real>(n)};
      if (jac) jac->assign(n * (n - 1), Real{0});
      for (std::size_t i = 0; i < n; ++i)
      {
        bethe::detail::CompensatedSum<Real> sum, scale, diagonal;
        sum.add(y[i] * (y[i] - Real{1}));
        scale.add(Real{1} + std::abs(y[i] * (y[i] - Real{1})));
        diagonal.add(Real{2} * y[i] - Real{1});
        Real last{};
        for (std::size_t j = 0; j < n; ++j)
          if (i != j)
          {
            Real const a = g / (levels[i] - levels[j]), term = a * (y[i] - y[j]);
            sum.add(-term);
            scale.add(std::abs(a) * (std::abs(y[i]) + std::abs(y[j])));
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
        out.norm = std::max(out.norm, std::abs(sum.value()) / normalization);
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
      slope.assign(n, Real{0});
      for (std::size_t i = 0; i < n; ++i)
      {
        bethe::detail::CompensatedSum<Real> rhs, scale;
        scale.add(Real{1} + std::abs(y[i] * (y[i] - Real{1})));
        for (std::size_t j = 0; j < n; ++j)
          if (i != j)
          {
            Real const term = (y[i] - y[j]) / (levels[i] - levels[j]);
            rhs.add(term);
            scale.add(std::abs(g / (levels[i] - levels[j])) * (std::abs(y[i]) + std::abs(y[j])));
          }
        slope[i] = rhs.value() / scale.value();
      }
      return bethe::detail::least_squares_step(std::move(jac), slope, n - 1);
    }
    Real energy(std::span<Real const> x, Real g) const
    {
      auto const y = expand(x);
      bethe::detail::CompensatedSum<Real> sum;
      for (std::size_t i = 0; i < n; ++i)
        sum.add(levels[i] * y[i]);
      sum.add(-g * Real(pairs) * Real(n - pairs + 1));
      return sum.value();
    }
    std::vector<Real> const levels;
    std::size_t const pairs, n;
};
} // namespace detail

template <uni20::Real Real = double>
[[nodiscard]] State<Real> ground_state(std::span<Real const> levels, std::size_t pairs, Real coupling,
                                       std::span<std::size_t const> blocked = {},
                                       SolverOptions<Real> const& options = {})
{
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
  out.levels.assign(levels.begin(), levels.end());
  out.pairs = pairs;
  out.coupling = coupling;
  bethe::detail::CompensatedSum<Real> constant;
  for (std::size_t i = 0; i < levels.size(); ++i)
    if (is_blocked[i])
    {
      out.blocked.push_back(i);
      constant.add(levels[i]);
    }
    else
      out.active.push_back(i);
  if (!coupling || !pairs || pairs == n)
  {
    out.eigenvalue_variables.assign(n, Real{0});
    for (std::size_t j = 0; j < pairs; ++j)
    {
      constant.add(levels[out.active[j]]);
      constant.add(levels[out.active[j]]);
      out.eigenvalue_variables[j] = Real{1};
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
  std::vector<Real> e(n);
  Real gap = Real{2};
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
  detail::System<Real> const system(std::move(e), pairs);
  std::vector<Real> x(n - 1, Real{0});
  for (std::size_t i = 0; i < pairs - 1; ++i)
    x[i] = Real{1};
  Real reached{}, step = std::min(target, gap / Real{16});
  Real const trust = Real{1} / Real{16}, eps = uni20::numeric_limits<Real>::epsilon();
  while (reached < target && out.iterations < options.max_iterations && out.stages < options.max_stages)
  {
    std::vector<Real> slope;
    if (!system.tangent(x, reached, slope))
    {
      out.status = SolveStatus::ill_conditioned;
      break;
    }
    Real max_slope{};
    bethe::detail::CompensatedSum<Real> last;
    for (Real v : slope)
    {
      max_slope = std::max(max_slope, std::abs(v));
      last.add(v);
    }
    max_slope = std::max(max_slope, std::abs(last.value()));
    if (max_slope > Real{0}) step = std::min(step, trust / max_slope);
    step = std::min(step, target - reached);
    Real const next = step == target - reached ? target : reached + step;
    if (!(next > reached))
    {
      out.status = SolveStatus::stalled;
      break;
    }
    ++out.stages;
    auto candidate = x;
    for (std::size_t i = 0; i < n - 1; ++i)
      candidate[i] += step * slope[i];
    auto const predictor = system.expand(candidate);
    bool accepted = false;
    for (unsigned update = 0; update < 12; ++update)
    {
      std::vector<Real> jac;
      auto const evaluation = system.evaluate(candidate, next, &jac);
      if (evaluation.norm <= options.residual_tolerance)
      {
        accepted = true;
        break;
      }
      if (out.iterations == options.max_iterations) break;
      ++out.iterations; // Attempted corrections, including failed/rejected stages.
      auto correction = evaluation.residual;
      for (auto& v : correction)
        v = -v;
      if (!bethe::detail::least_squares_step(std::move(jac), correction, n - 1)) break;
      bool improved = false;
      Real damping{1};
      for (int backtrack = 0; backtrack <= uni20::numeric_limits<Real>::digits; ++backtrack)
      {
        auto trial = candidate;
        for (std::size_t i = 0; i < n - 1; ++i)
          trial[i] += damping * correction[i];
        auto const y = system.expand(trial);
        bool near = true;
        for (std::size_t i = 0; i < n; ++i)
          if (!uni20::isfinite(y[i]) || std::abs(y[i] - predictor[i]) > trust) near = false;
        if (near)
        {
          Real const norm = system.evaluate(trial, next).norm;
          if (norm <= options.residual_tolerance || norm < (Real{1} - damping / Real{10000}) * evaluation.norm)
          {
            candidate = std::move(trial);
            improved = true;
            break;
          }
        }
        damping /= Real{2};
      }
      if (!improved) break;
    }
    if (accepted)
    {
      // Necessary physical bounds: kinetic minimum plus the maximal pair
      // interaction; variational Fermi-sea and uniform-pair (Dicke) states.
      bethe::detail::CompensatedSum<Real> free, sum;
      for (std::size_t i = 0; i < n; ++i)
      {
        sum.add(system.levels[i]);
        if (i < pairs) free.add(system.levels[i]);
      }
      Real const interaction = next * Real(pairs) * Real(n - pairs + 1);
      Real const lower = free.value() - interaction;
      Real const upper = std::min(free.value() - next * Real(pairs), Real(pairs) / Real(n) * sum.value() - interaction);
      Real const energy = system.energy(candidate, next);
      Real const allowance =
          Real{512} * Real(n) * (eps + options.residual_tolerance) * (Real{1} + std::abs(lower) + std::abs(upper));
      accepted = uni20::isfinite(energy) && energy >= lower - allowance && energy <= upper + allowance;
      // The ground energy is concave in g. Its tangent is an upper bound,
      // while the norm of S+S- bounds how far it can fall in one stage.
      bethe::detail::CompensatedSum<Real> derivative;
      for (std::size_t j = 1; j < n; ++j)
        derivative.add(system.levels[j] * slope[j - 1]);
      derivative.add(-Real(pairs) * Real(n - pairs + 1));
      Real const old_energy = system.energy(x, reached);
      accepted = accepted && energy <= old_energy + step * derivative.value() + allowance &&
                 energy >= old_energy - step * Real(pairs) * Real(n - pairs + 1) - allowance;
    }
    if (accepted)
    {
      x = std::move(candidate);
      reached = next;
      step *= Real{3} / Real{2};
    }
    else
    {
      ++out.rejected_stages;
      step /= Real{2};
    }
  }
  out.reached_coupling = reached == target ? coupling : reached * scale;
  out.eigenvalue_variables = system.expand(x);
  out.residual_norm = system.evaluate(x, reached).norm;
  out.target_residual_norm = system.evaluate(x, target).norm;
  bethe::detail::CompensatedSum<Real> number;
  for (Real v : out.eigenvalue_variables)
    number.add(v);
  out.particle_number_error = std::abs(number.value() - Real(pairs));
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
