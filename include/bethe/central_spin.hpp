// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/eigenvalue_continuation.hpp>
#include <bethe/solver.hpp>
#include <optional>
#include <span>
#include <stdexcept>
#include <uni20/common/half_int.hpp>

namespace bethe::central_spin
{
template <uni20::Real Real> struct SolverOptions : bethe::SolverOptions<Real>
{
    std::size_t max_stages = 10000;
};
enum class SolveStatus
{
  converged,
  iteration_limit,
  stage_limit,
  stalled,
  ill_conditioned
};

/// H=B*S0^z+sum_j A_j*S0.Sj; every spin is 1/2. No bath fields.
template <uni20::Real Real> struct State
{
    std::vector<Real> couplings;
    // Compactified eigenvalue variables: central first, then ORIGINAL bath order.
    // Empty for the analytic polarized/isolated-spin cases. Not occupations.
    std::vector<Real> eigenvalue_variables;
    uni20::half_int sz{};
    std::size_t up_spins = 0, iterations = 0, stages = 0, rejected_stages = 0;
    Real field{}, continuation_parameter{}, residual_norm{}, number_error{};
    // No finite-field energy exists at the infinite-field seed (zero budget).
    std::optional<Real> reached_field, energy;
    bool spin_reversed = false, converged = false;
    SolveStatus status = SolveStatus::iteration_limit;
};

namespace detail
{
// q0=0, qj=Astar/Aj. t=Astar/(2*|B|+Astar), p=1-t.
// v_i^2-p*v_i-t*sum_j(v_i-v_j)/(q_i-q_j)=0; sum v=M*p.
// Supply p separately to retain tiny nonzero fields when t rounds to one.
// Eliminate the occupied bath variable v1, NOT the tiny central variable v0.
template <uni20::Real Real> class System {
  public:
    System(std::vector<Real> poles, std::size_t count) : q(std::move(poles)), m(count), n(q.size()) {}
    std::vector<Real> expand(std::span<Real const> x, Real p) const
    {
      std::vector<Real> v(n);
      bethe::detail::CompensatedSum<Real> sum;
      sum.add(Real(m) * p);
      for (std::size_t j = 0; j < n - 1; ++j)
      {
        v[j ? j + 1 : 0] = x[j];
        sum.add(-x[j]);
      }
      v[1] = sum.value();
      return v;
    }
    struct Evaluation
    {
        std::vector<Real> residual;
        Real norm{};
    };
    Evaluation evaluate(std::span<Real const> x, Real t, Real p, std::vector<Real>* jac = nullptr,
                        std::vector<Real>* rhs = nullptr) const
    {
      auto const v = expand(x, p);
      Evaluation out{.residual = std::vector<Real>(n)};
      if (jac) jac->assign(n * (n - 1), Real{0});
      if (rhs) rhs->assign(n, Real{0});
      for (std::size_t i = 0; i < n; ++i)
      {
        bethe::detail::CompensatedSum<Real> sum, scale, diagonal, derivative;
        sum.add(v[i] * (v[i] - p));
        scale.add(Real{1} + std::abs(v[i] * (v[i] - p)));
        diagonal.add(Real{2} * v[i] - p);
        derivative.add(-v[i]);
        Real eliminated{};
        for (std::size_t j = 0; j < n; ++j)
          if (i != j)
          {
            Real const a = Real{1} / (q[i] - q[j]), b = t * a;
            sum.add(-b * (v[i] - v[j]));
            scale.add(std::abs(b) * (std::abs(v[i]) + std::abs(v[j])));
            diagonal.add(-b);
            derivative.add(a * (v[i] - v[j]));
            if (j == 1)
              eliminated = b;
            else if (jac)
              (*jac)[i * (n - 1) + (j ? j - 1 : 0)] = b;
          }
        if (i == 1)
          eliminated = diagonal.value();
        else if (jac)
          (*jac)[i * (n - 1) + (i ? i - 1 : 0)] = diagonal.value();
        Real const normalization = scale.value();
        if (!uni20::isfinite(sum.value()) || !uni20::isfinite(normalization))
          throw std::runtime_error("central-spin equations are not representable at this precision");
        out.residual[i] = sum.value() / normalization;
        out.norm = std::max(out.norm, std::abs(out.residual[i]));
        if (jac)
          for (std::size_t j = 0; j < n - 1; ++j)
            (*jac)[i * (n - 1) + j] = ((*jac)[i * (n - 1) + j] - eliminated) / normalization;
        if (rhs) (*rhs)[i] = (derivative.value() + Real(m) * eliminated) / normalization;
      }
      return out;
    }
    bool tangent(std::span<Real const> x, Real t, Real p, std::vector<Real>& slope) const
    {
      std::vector<Real> jac;
      (void)evaluate(x, t, p, &jac, &slope);
      return bethe::detail::least_squares_step(std::move(jac), slope, n - 1);
    }
    std::vector<Real> const q;
    std::size_t const m, n;
};
} // namespace detail

template <uni20::Real Real = double>
[[nodiscard]] State<Real> sector_ground_state(std::span<Real const> couplings, Real field, uni20::half_int sz,
                                              SolverOptions<Real> const& options = {})
{
  if (!uni20::isfinite(field)) throw std::invalid_argument("central-spin field must be finite");
  if (!uni20::isfinite(options.residual_tolerance) || options.residual_tolerance <= Real{0})
    throw std::invalid_argument("central-spin residual tolerance must be finite and positive");
  if (couplings.size() >= std::size_t(std::numeric_limits<std::int64_t>::max() / 2))
    throw std::length_error("too many central-spin bath sites");
  auto const n = couplings.size() + 1;
  auto const z = sz.twice(), length = std::int64_t(n);
  if (z < -length || z > length || (z + length) % 2)
    throw std::invalid_argument("Sz must be a physical total spin projection for central plus bath spins");
  State<Real> out;
  out.couplings.assign(couplings.begin(), couplings.end());
  out.field = field;
  out.sz = sz;
  out.up_spins = std::size_t((z + length) / 2);
  out.spin_reversed = field < Real{0};
  auto const m = out.spin_reversed ? n - out.up_spins : out.up_spins;
  std::vector<std::size_t> order(couplings.size());
  Real scale{};
  bethe::detail::CompensatedSum<Real> offset;
  for (std::size_t j = 0; j < couplings.size(); ++j)
  {
    Real const a = couplings[j];
    if (!uni20::isfinite(a) || a == Real{0})
      throw std::invalid_argument("bath couplings must be finite, distinct and nonzero");
    scale = std::max(scale, std::abs(a));
    offset.add(a / Real{4});
    order[j] = j;
  }
  std::sort(order.begin(), order.end(), [&](auto i, auto j) { return couplings[i] > couplings[j]; });
  for (std::size_t j = 1; j < order.size(); ++j)
    if (couplings[order[j]] == couplings[order[j - 1]])
      throw std::invalid_argument("repeated bath couplings require grouped-spin equations and are not implemented");
  if (!m || m == n)
  {
    offset.add((out.up_spins ? field : -field) / Real{2});
    out.energy = offset.value();
    out.reached_field = field;
    out.converged = true;
    out.status = SolveStatus::converged;
    if (!uni20::isfinite(*out.energy)) throw std::runtime_error("central-spin energy is not representable");
    return out;
  }
  auto const elements = std::size_t(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(Real);
  if (n > elements / (n - 1)) throw std::length_error("central-spin Newton matrix is too large");
  Real const h = std::abs(field) / scale;
  if (!uni20::isfinite(h) || (field != Real{0} && h == Real{0}))
    throw std::invalid_argument("field/coupling ratio is not representable at this precision");
  Real const target = Real{0.5} / (h + Real{0.5}), target_p = h / (h + Real{0.5});
  if (!(target > Real{0})) throw std::invalid_argument("continuation parameter is not representable");
  std::vector<Real> q(n), a(n - 1);
  Real gap{1};
  bethe::detail::CompensatedSum<Real> sum_a, abs_a, selected;
  for (std::size_t j = 1; j < n; ++j)
  {
    a[j - 1] = couplings[order[j - 1]] / scale;
    q[j] = scale / couplings[order[j - 1]];
    if (!uni20::isfinite(q[j])) throw std::invalid_argument("bath coupling ratio is not representable");
    sum_a.add(a[j - 1]);
    abs_a.add(std::abs(a[j - 1]));
    if (j <= m) selected.add(a[j - 1]);
    for (std::size_t k = 0; k < j; ++k)
    {
      Real const distance = std::abs(q[j] - q[k]);
      if (!(distance > Real{0}) || !uni20::isfinite(distance) || !uni20::isfinite(Real{1} / distance))
        throw std::invalid_argument("distinct Gaudin poles are not resolvable at this precision");
      gap = std::min(gap, distance);
    }
  }
  detail::System<Real> const system(std::move(q), m);
  std::vector<Real> x(n - 1, Real{0});
  for (std::size_t j = 1; j < m; ++j)
    x[j] = Real{1};
  Real reached{}, p{1}, step = std::min(target, gap / Real{16});
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  Real const allowance = Real{512} * Real(n) * (eps + options.residual_tolerance) * (Real{1} + abs_a.value());
  auto reduced_energy = [&](std::span<Real const> v, Real field_scaled) {
    return field_scaled * v[0] + v[0] / Real{2} + sum_a.value() / Real{4};
  };
  auto predict = [&](std::vector<Real>& slope) -> std::optional<bethe::detail::ContinuationPrediction<Real>> {
    if (!system.tangent(x, reached, p, slope)) return std::nullopt;
    Real max_slope{}, max_value{};
    bethe::detail::CompensatedSum<Real> eliminated;
    eliminated.add(-Real(m));
    for (Real v : slope)
    {
      max_slope = std::max(max_slope, std::abs(v));
      eliminated.add(-v);
    }
    max_slope = std::max(max_slope, std::abs(eliminated.value()));
    for (Real v : system.expand(x, p))
      max_value = std::max(max_value, std::abs(v));
    return bethe::detail::ContinuationPrediction<Real>{(Real{1} + max_value) / Real{16}, max_slope};
  };
  auto correct = [&](std::vector<Real>& candidate, Real next, Real /* step */, std::vector<Real> const& slope,
                     Real trust_radius) {
    Real const next_p = next == target ? target_p : Real{1} - next;
    bool accepted = bethe::detail::correct_eigenvalue_stage(
        candidate, trust_radius, options, out.iterations,
        [&](auto const& trial, std::vector<Real>* jac) { return system.evaluate(trial, next, next_p, jac); },
        [&](auto const& trial) { return system.expand(trial, next_p); });
    if (accepted)
    {
      // Work with E+|B|/2 to avoid losing the interaction energy at high field.
      Real const next_h = next == target ? h : next_p / (Real{2} * next);
      Real const energy = reduced_energy(candidate, next_h);
      Real const lower = -Real{0.75} * abs_a.value();
      Real upper = sum_a.value() / Real{4} - selected.value() / Real{2};
      bethe::detail::CompensatedSum<Real> up;
      for (std::size_t j = 0; j < m - 1; ++j)
        up.add(a[a.size() - 1 - j]);
      upper = std::min(upper, next_h - sum_a.value() / Real{4} + up.value() / Real{2});
      accepted = uni20::isfinite(energy) && energy >= lower - allowance && energy <= upper + allowance;
      if (reached > Real{0})
      {
        Real const old_h = p / (Real{2} * reached), dh = next_h - old_h;
        Real const old_energy = reduced_energy(x, old_h), derivative = x[0] - reached * slope[0];
        accepted = accepted && energy >= old_energy + dh - allowance && energy <= old_energy + allowance &&
                   energy <= old_energy + dh * derivative + allowance;
      }
    }
    return accepted;
  };
  bethe::detail::continue_eigenvalues(x, reached, step, target, options, out, predict, correct,
                                      [&](Real next) { p = next == target ? target_p : Real{1} - next; });
  out.continuation_parameter = reached;
  auto const v = system.expand(x, p);
  out.eigenvalue_variables.resize(n);
  out.eigenvalue_variables[0] = v[0];
  for (std::size_t j = 1; j < n; ++j)
    out.eigenvalue_variables[order[j - 1] + 1] = v[j];
  out.residual_norm = system.evaluate(x, reached, p).norm;
  bethe::detail::CompensatedSum<Real> number;
  for (Real value : v)
    number.add(value);
  out.number_error = std::abs(number.value() - Real(m) * p);
  if (reached > Real{0})
  {
    Real const reached_h = reached == target ? h : p / (Real{2} * reached);
    Real const actual_field = reached == target ? std::abs(field) : scale * reached_h;
    out.reached_field = out.spin_reversed ? -actual_field : actual_field;
    out.energy = scale * reduced_energy(x, reached_h) - actual_field / Real{2};
    if (!uni20::isfinite(*out.energy) || !uni20::isfinite(*out.reached_field))
      throw std::runtime_error("central-spin energy or reached field is not representable");
  }
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
} // namespace bethe::central_spin
