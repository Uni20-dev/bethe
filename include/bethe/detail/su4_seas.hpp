// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <array>
#include <bethe/detail/newton.hpp>
#include <bethe/solver.hpp>
#include <cstdint>
#include <functional>
#include <span>
#include <stdexcept>
#include <uni20/common/half_int.hpp>

namespace bethe::ladder::detail
{
using Shape = std::array<std::size_t, 4>;
using Labels = std::vector<std::vector<uni20::half_int>>;
enum class BranchStatus
{
  converged,
  iteration_limit,
  stalled,
  ill_conditioned
};
template <uni20::Real Real> struct Branch
{
    Shape shape{};
    Labels labels;
    std::vector<std::vector<Real>> rapidities;
    Real energy{}, residual{};
    std::size_t momentum_index = 0, iterations = 0;
    bool converged = false;
    BranchStatus status = BranchStatus::iteration_limit;
};

// Full (not reflection-reduced) real nested seas. Mixed row parities require
// half-step displaced seas and can have nonzero lattice momentum.
template <uni20::Real Real> class System {
  public:
    System(Shape const& shape, std::array<int, 3> const& twice_shift)
      requires(!uni20::has_runtime_precision_v<Real>)
        : System(shape, twice_shift, Real{0})
    {}
    System(Shape const& shape, std::array<int, 3> const& twice_shift, Real const& exemplar)
        : zero(bethe::detail::precision_zero(exemplar)), pi(Real{4} * uni20::math::atan(uni20::scalar_like(zero, 1)))
    {
      for (auto p : shape)
        sites += p;
      while (colors < 4 && shape[colors])
        ++colors;
      std::size_t remaining = sites;
      offsets.push_back(0);
      for (std::size_t a = 0; a + 1 < colors; ++a)
      {
        remaining -= shape[a];
        counts.push_back(remaining);
        offsets.push_back(offsets.back() + remaining);
        labels.emplace_back();
        for (std::size_t j = 0; j < remaining; ++j)
          labels.back().push_back(
              uni20::from_twice(2 * std::int64_t(j) - std::int64_t(remaining - 1) + twice_shift[a]));
      }
    }
    std::size_t order() const { return offsets.back(); }
    std::vector<Real> seed() const
    {
      std::vector<Real> x(order(), zero);
      for (std::size_t a = 0; a < counts.size(); ++a)
        for (std::size_t j = 0; j < counts[a]; ++j)
        {
          Real const label = Real(labels[a][j].twice()) / Real{2};
          Real const argument = uni20::math::tan(pi * label / Real(sites)) *
                                uni20::math::tan(pi * Real(a + 1) / (Real{2} * Real(colors)));
          // The balanced thermodynamic density supplies only a starting
          // guess, including for shifted and imbalanced finite-size seas.
          Real const bound = uni20::scalar_like(zero, 97) / Real{100};
          x[offsets[a] + j] = Real(colors) / pi * uni20::math::atanh(std::clamp(argument, -bound, bound));
        }
      return x;
    }
    bool physical(std::span<Real const> x) const
    {
      if (x.size() != order()) return false;
      Real const limit = [&] {
        if constexpr (uni20::has_runtime_precision_v<Real>)
          return uni20::math::next_down(uni20::scalar_like(zero, std::numeric_limits<double>::infinity())) / Real{8};
        else
          return uni20::numeric_limits<Real>::max() / Real{8};
      }();
      for (std::size_t a = 0; a < counts.size(); ++a)
        for (std::size_t j = 0; j < counts[a]; ++j)
        {
          Real const v = x[offsets[a] + j];
          if (!uni20::isfinite(v) || uni20::math::abs(v) > limit || (j && v <= x[offsets[a] + j - 1])) return false;
        }
      return true;
    }
    static Real phase(Real d, Real w) { return Real{2} * uni20::math::atan(d / w); }
    static Real derivative(Real d, Real w)
    {
      if (uni20::math::abs(d) > w)
      {
        Real const r = w / d;
        return (Real{2} * r / d) / (Real{1} + r * r);
      }
      Real const r = d / w;
      return (Real{2} / w) / (Real{1} + r * r);
    }
    struct Evaluation
    {
        std::vector<Real> residual;
        Real norm{};
    };
    Evaluation evaluate(std::span<Real const> x, std::vector<Real>* jac = nullptr) const
    {
      Evaluation out{.residual = std::vector<Real>(order(), zero), .norm = zero};
      Real const half = uni20::scalar_like(zero, 1) / Real{2};
      if (jac) jac->assign(order() * order(), zero);
      for (std::size_t a = 0; a < counts.size(); ++a)
        for (std::size_t j = 0; j < counts[a]; ++j)
        {
          auto const row = offsets[a] + j;
          Real const value = x[row];
          bethe::detail::CompensatedSum<Real> sum(zero), diagonal(zero);
          sum.add(-pi * Real(labels[a][j].twice()) / Real(sites));
          if (a == 0)
          {
            sum.add(phase(value, half));
            diagonal.add(derivative(value, half));
          }
          for (std::size_t b = a ? a - 1 : 0; b < std::min(counts.size(), a + 2); ++b)
          {
            Real const sign = a == b ? Real{-1} : Real{1}, w = a == b ? uni20::scalar_like(zero, 1) : half;
            for (std::size_t k = 0; k < counts[b]; ++k)
            {
              auto const col = offsets[b] + k;
              if (col == row) continue;
              Real const diff = value - x[col];
              sum.add(sign * phase(diff, w) / Real(sites));
              if (jac)
              {
                Real const d = sign * derivative(diff, w) / Real(sites);
                diagonal.add(d);
                (*jac)[row * order() + col] -= d;
              }
            }
          }
          if (!uni20::isfinite(sum.value())) throw std::runtime_error("nonfinite SU(4) sea residual");
          out.residual[row] = sum.value();
          out.norm = std::max(out.norm, uni20::math::abs(sum.value()));
          if (jac) (*jac)[row * order() + row] += diagonal.value();
        }
      return out;
    }
    std::size_t sites = 0, colors = 0;
    std::vector<std::size_t> counts, offsets;
    Labels labels;
    Real const zero, pi;
};

template <uni20::Real Real>
Branch<Real> solve(Shape const& shape, std::array<int, 3> const& shift, Real tolerance, std::size_t budget)
{
  Real const zero = bethe::detail::precision_zero(tolerance);
  auto const precision = uni20::precision_of(tolerance);
  System<Real> const system(shape, shift, zero);
  Branch<Real> out;
  out.energy = out.residual = zero;
  out.shape = shape;
  out.labels = system.labels;
  auto x = system.seed();
  for (;;)
  {
    if (!system.physical(x)) throw std::runtime_error("SU(4) sea roots are not representable at this precision");
    std::vector<Real> jac;
    auto const evaluation = system.evaluate(x, &jac);
    out.residual = evaluation.norm;
    if (evaluation.norm <= tolerance)
    {
      out.converged = true;
      out.status = BranchStatus::converged;
      break;
    }
    if (out.iterations == budget) break;
    ++out.iterations;
    auto correction = evaluation.residual;
    for (auto& v : correction)
      v = -v;
    if (!bethe::detail::newton_step(std::move(jac), correction, precision))
    {
      out.status = BranchStatus::ill_conditioned;
      break;
    }
    auto trial = x;
    bool improved = false;
    Real damping = uni20::scalar_like(zero, 1);
    for (int backtrack = 0; backtrack <= uni20::numeric_limits<Real>::digits(precision); ++backtrack)
    {
      for (std::size_t j = 0; j < x.size(); ++j)
        trial[j] = x[j] + damping * correction[j];
      if (system.physical(trial))
      {
        Real const norm = system.evaluate(trial).norm;
        if (norm <= tolerance || norm < (Real{1} - damping / Real{10000}) * evaluation.norm)
        {
          improved = true;
          break;
        }
      }
      damping /= Real{2};
    }
    if (!improved)
    {
      out.status = BranchStatus::stalled;
      break;
    }
    x = std::move(trial);
  }
  if (!out.converged) return out; // Never expose an unconverged branch as a candidate eigenstate.
  for (std::size_t a = 0; a < system.counts.size(); ++a)
    out.rapidities.emplace_back(x.begin() + system.offsets[a], x.begin() + system.offsets[a + 1]);
  bethe::detail::CompensatedSum<Real> energy(zero);
  energy.add(Real(system.sites));
  Real const half = uni20::scalar_like(zero, 1) / Real{2};
  if (!out.rapidities.empty())
    for (Real v : out.rapidities[0])
      energy.add(-System<Real>::derivative(v, half));
  out.energy = energy.value();
  // P=pi*M1-(2*pi/L)*sum(all labels). Wider arithmetic avoids overflow
  // before the exact modulo-L reduction. Infinite-root descendants add no P.
  __int128 twice_sum = 0;
  for (auto const& sea : out.labels)
    for (auto label : sea)
      twice_sum += label.twice();
  __int128 const numerator = __int128(system.sites) * (system.counts.empty() ? 0 : system.counts[0]) - twice_sum;
  if (numerator % 2) throw std::logic_error("invalid SU(4) sea momentum parity");
  auto index = (numerator / 2) % system.sites;
  if (index < 0) index += system.sites;
  out.momentum_index = std::size_t(index);
  if (!uni20::isfinite(out.energy)) throw std::runtime_error("nonfinite SU(4) sea energy");
  return out;
}

inline bool dominates(Shape const& highest, Shape const& population)
{
  auto target = population;
  std::sort(target.begin(), target.end(), std::greater<>{});
  std::size_t lhs = 0, rhs = 0;
  for (std::size_t a = 0; a < 4; ++a)
  {
    lhs += highest[a];
    rhs += target[a];
    if (lhs < rhs) return false;
  }
  return lhs == rhs;
}
// Independent shifted-sea branches modulo simultaneous reflection. At each
// level, centering is admissible iff the two adjacent root counts sum to even.
inline std::vector<std::array<int, 3>> shifts(Shape const& shape)
{
  std::array<std::size_t, 5> m{};
  for (std::size_t a = 0; a < 4; ++a)
    for (std::size_t b = a; b < 4; ++b)
      m[a] += shape[b];
  std::vector<std::array<int, 3>> out(1);
  bool first = true;
  for (std::size_t a = 0; a < 3 && m[a + 1]; ++a)
    if ((m[a] + m[a + 2]) % 2)
    {
      auto const count = out.size();
      if (!first)
        for (std::size_t j = 0; j < count; ++j)
        {
          auto copy = out[j];
          copy[a] = 1;
          out.push_back(copy);
        }
      for (std::size_t j = 0; j < count; ++j)
        out[j][a] = -1;
      first = false;
    }
  return out;
}
} // namespace bethe::ladder::detail
