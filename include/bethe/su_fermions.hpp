// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/newton.hpp>
#include <bethe/solver.hpp>
#include <cstdint>
#include <optional>
#include <span>
#include <stdexcept>
#include <uni20/common/half_int.hpp>

namespace bethe::su_fermions
{
using QuantumNumbers = std::vector<std::vector<uni20::half_int>>;
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

/// Equal-mass SU(kappa) fermions on a ring; H=-sum d_j^2+2c sum delta, c>=0.
/// The interacting centered branch requires every occupied population odd.
template <uni20::Real Real> struct State
{
    std::vector<std::size_t> populations, component_order;
    std::size_t particles = 0, iterations = 0, stages = 0;
    Real length{1}, interaction{}, residual_norm{}, target_residual_norm{}, momentum{};
    std::int64_t momentum_index = 0;
    QuantumNumbers quantum_numbers;            // level 0 is charge, remaining levels are spin.
    std::vector<std::vector<Real>> rapidities; // Physical k,lambda; only accepted stages.
    std::vector<Real> level_residuals;
    std::vector<std::vector<std::int64_t>> free_modes; // Original component order; free path only.
    std::optional<Real> energy, reached_interaction;
    bool free = false, converged = false;
    SolveStatus status = SolveStatus::iteration_limit;
};

namespace detail
{
struct Sector
{
    std::size_t n = 0, order = 0;
    std::vector<std::size_t> components, counts;
};
inline Sector sector(std::span<std::size_t const> populations)
{
  if (populations.empty()) throw std::invalid_argument("at least one component population is required");
  Sector out;
  auto const limit = std::size_t(std::numeric_limits<std::int64_t>::max() / 4);
  for (std::size_t a = 0; a < populations.size(); ++a)
  {
    if (populations[a] > limit - out.n)
      throw std::invalid_argument("SU fermion particle count exceeds the label range");
    out.n += populations[a];
    if (populations[a]) out.components.push_back(a);
  }
  std::stable_sort(out.components.begin(), out.components.end(),
                   [&](auto i, auto j) { return populations[i] > populations[j]; });
  std::size_t remaining = out.n;
  for (auto a : out.components)
  {
    out.counts.push_back(remaining);
    if (remaining / 2 > std::numeric_limits<std::size_t>::max() - out.order)
      throw std::length_error("SU fermion nested root count is too large");
    out.order += remaining / 2;
    remaining -= populations[a];
  }
  return out;
}
inline void check_branch(std::span<std::size_t const> populations, Sector const& s)
{
  if (s.components.size() < 2)
    throw std::invalid_argument("interacting nested labels require at least two occupied components");
  for (auto a : s.components)
    if (populations[a] % 2 == 0)
      throw std::invalid_argument("interacting SU fermion ground states require every occupied population odd; other "
                                  "periodic shell branches are not implemented");
}
inline QuantumNumbers centered_numbers(std::span<std::size_t const> counts)
{
  QuantumNumbers out(counts.size());
  for (std::size_t a = 0; a < counts.size(); ++a)
    for (std::size_t j = 0; j < counts[a]; ++j)
      out[a].push_back(uni20::from_twice(2 * std::int64_t(j) - std::int64_t(counts[a] - 1)));
  return out;
}
template <uni20::Real Real> class GroundSystem {
  public:
    explicit GroundSystem(std::vector<std::size_t> nested_counts)
        : counts(std::move(nested_counts)), labels(centered_numbers(counts)), pi(Real{4} * std::atan(Real{1}))
    {
      offsets.push_back(0);
      for (auto count : counts)
        offsets.push_back(offsets.back() + count / 2);
    }
    std::size_t order() const { return offsets.back(); }
    static Real scale(Real g) { return std::max(Real{1}, g); }
    std::vector<Real> seed() const
    {
      std::vector<Real> x(order());
      Real const colors = Real(counts.size());
      for (std::size_t a = 0; a < counts.size(); ++a)
        for (std::size_t j = 0; j < counts[a] / 2; ++j)
        {
          Real const label = Real(labels[a][counts[a] - counts[a] / 2 + j].twice()) / Real{2};
          // Strong-coupling charge sea and inverse balanced SU(kappa) spin
          // density. For imbalanced populations this remains only a seed.
          x[offsets[a] + j] =
              a ? colors / pi *
                      std::atanh(std::tan(pi * label / Real(counts[0])) * std::tan(pi * Real(a) / (Real{2} * colors)))
                : Real{2} * pi * label;
        }
      return x;
    }
    std::vector<std::vector<Real>> expand(std::span<Real const> x, Real g) const
    {
      std::vector<std::vector<Real>> roots(counts.size());
      for (std::size_t a = 0; a < counts.size(); ++a)
      {
        Real const factor = a ? scale(g) : Real{1};
        for (std::size_t j = counts[a] / 2; j-- > 0;)
          roots[a].push_back(-factor * x[offsets[a] + j]);
        if (counts[a] % 2) roots[a].push_back(Real{0});
        for (std::size_t j = 0; j < counts[a] / 2; ++j)
          roots[a].push_back(factor * x[offsets[a] + j]);
      }
      return roots;
    }
    bool physical(std::span<Real const> x, Real g) const
    {
      if (x.size() != order()) return false;
      for (std::size_t a = 0; a < counts.size(); ++a)
        for (std::size_t j = 0; j < counts[a] / 2; ++j)
        {
          Real const value = x[offsets[a] + j];
          if (!uni20::isfinite(value) || value <= Real{0} ||
              value > uni20::numeric_limits<Real>::max() / Real{8} / (a ? scale(g) : Real{1}) ||
              (j && value <= x[offsets[a] + j - 1]))
            return false;
        }
      return true;
    }
    static Real kernel(Real d, Real w)
    {
      if (std::abs(d) > w)
      {
        Real const r = w / d;
        return (r / d) / (Real{1} + r * r);
      }
      Real const r = d / w;
      return (Real{1} / w) / (Real{1} + r * r);
    }
    struct Evaluation
    {
        std::vector<Real> residual, levels;
        Real norm{};
    };
    Evaluation evaluate(std::span<Real const> x, Real g, std::vector<Real>* jac = nullptr) const
    {
      auto const roots = expand(x, g);
      Evaluation out{.residual = std::vector<Real>(order()), .levels = std::vector<Real>(counts.size())};
      if (jac) jac->assign(order() * order(), Real{0});
      for (std::size_t a = 0; a < counts.size(); ++a)
        for (std::size_t j = 0; j < counts[a] / 2; ++j)
        {
          auto const row = offsets[a] + j, full = counts[a] - counts[a] / 2 + j;
          Real const value = roots[a][full];
          bethe::detail::CompensatedSum<Real> sum, diagonal;
          std::int64_t rank = -labels[a][full].twice();
          if (g >= Real{1}) sum.add(pi * Real(rank));
          if (a == 0)
          {
            sum.add(value);
            diagonal.add(Real{1});
          }
          for (std::size_t b = a ? a - 1 : 1; b < std::min(counts.size(), a + 2); ++b)
          {
            Real const sign = a == b ? Real{-1} : Real{1}, w = a == b ? g : g / Real{2};
            for (std::size_t k = 0; k < counts[b]; ++k)
            {
              if (a == b && k == full) continue;
              Real const d = value - roots[b][k];
              Real phase;
              if (g < Real{1} && std::abs(d) > w)
              {
                rank += (d > Real{0} ? 1 : -1) * (sign > Real{0} ? 1 : -1);
                phase = -Real{2} * std::atan(w / d);
              }
              else
                phase = Real{2} * std::atan(d / w);
              sum.add(sign * phase);
              if (jac)
              {
                Real const derivative = sign * Real{2} * kernel(d, w);
                diagonal.add(derivative * (a ? scale(g) : Real{1}));
                auto const positives = counts[b] / 2, first = counts[b] - positives;
                if (k < positives)
                  (*jac)[row * order() + offsets[b] + positives - 1 - k] += derivative * (b ? scale(g) : Real{1});
                else if (k >= first)
                  (*jac)[row * order() + offsets[b] + k - first] -= derivative * (b ? scale(g) : Real{1});
              }
            }
          }
          if (g < Real{1}) sum.add(pi * Real(rank));
          // All weak-coupling central clusters, including nested spin seas,
          // need relative sqrt(g) resolution, not just the charge roots.
          Real const normalization = g < Real{1} ? std::max(std::sqrt(g), std::abs(value)) : Real(counts[0]);
          if (!uni20::isfinite(sum.value())) throw std::runtime_error("nonfinite SU fermion equation residual");
          out.residual[row] = sum.value() / normalization;
          out.levels[a] = std::max(out.levels[a], std::abs(out.residual[row]));
          if (jac)
          {
            (*jac)[row * order() + row] += diagonal.value();
            for (std::size_t col = 0; col < order(); ++col)
              (*jac)[row * order() + col] /= normalization;
          }
        }
      for (Real v : out.levels)
        out.norm = std::max(out.norm, v);
      return out;
    }
    void rescale(std::vector<Real>& x, Real from, Real to) const
    {
      for (std::size_t j = offsets[1]; j < order(); ++j)
        x[j] *= scale(from) / scale(to);
    }
    std::vector<std::size_t> const counts;
    std::vector<std::size_t> offsets;
    QuantumNumbers const labels;
    Real const pi;
};
template <uni20::Real Real> Real physical_root(Real value, Real length)
{
  Real const r = value / length;
  if (!uni20::isfinite(r) || (value != Real{0} && r == Real{0}))
    throw std::runtime_error("SU fermion physical root is not representable; rescale the units");
  return r;
}
template <uni20::Real Real> Real energy(std::span<Real const> momenta)
{
  bethe::detail::CompensatedSum<Real> sum;
  for (Real k : momenta)
  {
    if (k != Real{0} && k * k == Real{0})
      throw std::runtime_error("SU fermion kinetic energy underflow; rescale the units");
    sum.add(k * k);
  }
  if (!uni20::isfinite(sum.value())) throw std::runtime_error("SU fermion energy overflow; rescale the units");
  return sum.value();
}
} // namespace detail

inline QuantumNumbers ground_quantum_numbers(std::span<std::size_t const> populations)
{
  auto const s = detail::sector(populations);
  detail::check_branch(populations, s);
  return detail::centered_numbers(s.counts);
}

template <uni20::Real Real = double>
[[nodiscard]] State<Real> ground_state(std::span<std::size_t const> populations, Real length, Real c,
                                       SolverOptions<Real> const& options = {})
{
  auto const sector = detail::sector(populations);
  if (!uni20::isfinite(length) || length <= Real{0} || !uni20::isfinite(c) || c < Real{0})
    throw std::invalid_argument("repulsive SU fermions require finite length>0 and c>=0");
  if (!uni20::isfinite(options.residual_tolerance) || options.residual_tolerance <= Real{0})
    throw std::invalid_argument("SU fermion residual tolerance must be finite and positive");
  State<Real> out;
  out.populations.assign(populations.begin(), populations.end());
  out.component_order = sector.components;
  out.particles = sector.n;
  out.length = length;
  out.interaction = c;
  Real const pi = Real{4} * std::atan(Real{1});
  if (c == Real{0} || sector.components.size() < 2)
  {
    out.free_modes.resize(populations.size());
    out.rapidities.resize(1);
    for (std::size_t a = 0; a < populations.size(); ++a)
    {
      auto const count = populations[a];
      for (std::size_t j = 0; j < count; ++j)
      {
        auto const mode = std::int64_t(j) - std::int64_t((count - 1) / 2);
        out.free_modes[a].push_back(mode);
        out.rapidities[0].push_back(detail::physical_root(Real{2} * pi * Real(mode), length));
      }
      if (count % 2 == 0) out.momentum_index += std::int64_t(count / 2);
    }
    std::sort(out.rapidities[0].begin(), out.rapidities[0].end());
    out.energy = detail::energy<Real>(out.rapidities[0]);
    out.momentum = detail::physical_root(Real{2} * pi * Real(out.momentum_index), length);
    out.reached_interaction = c;
    out.free = out.converged = true;
    out.status = SolveStatus::converged;
    return out;
  }
  detail::check_branch(populations, sector);
  Real const target = c * length;
  if (!uni20::isfinite(target) || target / Real{2} <= Real{0})
    throw std::invalid_argument("c*length/2 must be finite and nonzero at the selected precision");
  auto const elements = std::size_t(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(Real);
  if (sector.order > elements / sector.order) throw std::length_error("SU fermion Newton matrix is too large");
  detail::GroundSystem<Real> const system(sector.counts);
  out.quantum_numbers = system.labels;
  auto x = system.seed();
  Real g = std::max(target, Real{64} * Real(sector.n)), accepted_g{};
  std::vector<Real> accepted_x;
  bool finished = false;
  while (!finished && out.iterations < options.max_iterations && out.stages < options.max_stages)
  {
    ++out.stages;
    bool accepted = false;
    for (;;)
    {
      if (!system.physical(x, g)) throw std::runtime_error("SU fermion roots are not representable at this precision");
      std::vector<Real> jac;
      auto const evaluation = system.evaluate(x, g, &jac);
      if (evaluation.norm <= options.residual_tolerance)
      {
        accepted = true;
        break;
      }
      if (out.iterations == options.max_iterations) break;
      ++out.iterations;
      auto correction = evaluation.residual;
      for (auto& v : correction)
        v = -v;
      if (!bethe::detail::newton_step(std::move(jac), correction))
      {
        out.status = SolveStatus::ill_conditioned;
        break;
      }
      auto trial = x;
      Real damping{1};
      bool improved = false;
      for (int backtrack = 0; backtrack <= uni20::numeric_limits<Real>::digits; ++backtrack)
      {
        for (std::size_t j = 0; j < x.size(); ++j)
          trial[j] = x[j] + damping * correction[j];
        if (system.physical(trial, g))
        {
          Real const norm = system.evaluate(trial, g).norm;
          if (norm <= options.residual_tolerance || norm < (Real{1} - damping / Real{10000}) * evaluation.norm)
          {
            improved = true;
            break;
          }
        }
        damping /= Real{2};
      }
      if (!improved)
      {
        out.status = SolveStatus::stalled;
        break;
      }
      x = std::move(trial);
    }
    if (!accepted) break;
    accepted_x = x;
    accepted_g = g;
    if (g == target)
    {
      out.converged = true;
      out.status = SolveStatus::converged;
      finished = true;
    }
    else
    {
      Real const next = std::max(target, g / Real{2});
      system.rescale(x, g, next);
      g = next;
    }
  }
  if (!accepted_x.empty())
  {
    out.rapidities = system.expand(accepted_x, accepted_g);
    for (auto& level : out.rapidities)
      for (auto& r : level)
        r = detail::physical_root(r, length);
    out.energy = detail::energy<Real>(out.rapidities[0]);
    out.reached_interaction = accepted_g == target ? c : detail::physical_root(accepted_g, length);
    auto const evaluation = system.evaluate(accepted_x, accepted_g);
    out.residual_norm = evaluation.norm;
    out.level_residuals = evaluation.levels;
    system.rescale(accepted_x, accepted_g, target);
    out.target_residual_norm = system.evaluate(accepted_x, target).norm;
  }
  if (!out.converged)
  {
    if (out.iterations == options.max_iterations)
      out.status = SolveStatus::iteration_limit;
    else if (out.stages == options.max_stages)
      out.status = SolveStatus::stage_limit;
  }
  return out;
}
} // namespace bethe::su_fermions
