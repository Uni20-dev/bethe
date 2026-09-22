// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <algorithm>
#include <bethe/solver.hpp>
#include <cstdint>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <vector>

namespace bethe::sutherland
{
inline constexpr std::size_t max_particles = 1000000;
namespace detail
{
template <uni20::Real Real> void validate(std::size_t n, Real length, Real lambda)
{
  if (n == 0 || n > max_particles) throw std::invalid_argument("Sutherland requires 1<=particles<=1000000");
  if (!uni20::isfinite(length) || length <= Real{0}) throw std::invalid_argument("length must be finite and positive");
  if (!uni20::isfinite(lambda) || lambda < Real{0})
    throw std::invalid_argument("lambda must be finite and nonnegative");
}
inline std::optional<std::size_t> window_count(std::size_t n, std::size_t window, std::size_t cap)
{
  if (cap == 0) return {};
  std::size_t value = 1, total = n + 2 * window, k = std::min(n, 2 * window);
  for (std::size_t j = 1; j <= k; ++j)
  {
    auto a = total - k + j, b = j;
    auto const g = std::gcd(a, b);
    a /= g;
    b /= g;
    value /= b; // Binomial recurrence: the remaining denominator divides value.
    if (value > cap / a) return {};
    value *= a;
  }
  return value;
}
} // namespace detail

/// Scalar periodic bosons, collision branch psi ~ |x_i-x_j|^lambda.
/// Labels are ascending integers, repeated labels allowed, including boosts.
template <uni20::Real Real> struct State
{
    std::size_t particles = 0;
    Real length{}, lambda{}, energy{}, ground_energy{}, gap{}, momentum{};
    std::int64_t momentum_index = 0; // P=2*pi*index/L, NOT modulo anything.
    std::vector<std::int64_t> labels;
    std::vector<Real> pseudomomenta;
};

/// H=-sum d_i^2+2*lambda*(lambda-1)*(pi/L)^2 sum_{i<j} csc^2(pi*(x_i-x_j)/L).
/// The branch parameter, not just the potential coefficient, fixes the domain.
template <uni20::Real Real>
State<Real> evaluate(std::size_t n, Real length, Real lambda, std::vector<std::int64_t> labels)
{
  detail::validate(n, length, lambda);
  if (labels.size() != n || !std::is_sorted(labels.begin(), labels.end()))
    throw std::invalid_argument("supply exactly N nondecreasing integer labels");
  // N<=10^6: a wide accumulator safely sums any int64 labels, even when
  // intermediate prefixes overflow int64 but the final momentum does not.
  __int128_t index = 0;
  for (auto label : labels)
    index += label;
  if (index < std::numeric_limits<std::int64_t>::min() || index > std::numeric_limits<std::int64_t>::max())
    throw std::overflow_error("total momentum index exceeds int64 range");
  State<Real> out;
  out.particles = n;
  out.length = length;
  out.lambda = lambda;
  out.momentum_index = static_cast<std::int64_t>(index);
  out.labels = std::move(labels);
  Real const pi = Real{4} * std::atan(Real{1}), unit = Real{2} * (pi / length);
  bethe::detail::CompensatedSum<Real> free_weight, interaction_weight;
  out.pseudomomenta.resize(n);
  for (std::size_t j = 0; j < n; ++j)
  {
    Real const label = Real(out.labels[j]);
    out.pseudomomenta[j] = unit * (label + lambda * (Real(j) - Real(n - 1) / Real{2}));
    // Positive expansion avoids cancellation of small gaps against E0, or of
    // large common boosts in the interaction term. Unsigned subtraction is
    // the exact nonnegative int64 label difference, including across zero.
    // Sum the integer coefficients before introducing q and lambda. For modest
    // labels these are exactly representable, preserving reflection degeneracies
    // instead of accumulating differently rounded momentum-square terms.
    free_weight.add(label * label);
    if (j && lambda != Real{0})
    {
      auto const difference = static_cast<std::uint64_t>(out.labels[j]) - static_cast<std::uint64_t>(out.labels[j - 1]);
      if (difference) interaction_weight.add(Real(j * (n - j)) * Real(difference));
    }
  }
  Real const edge = (pi / length) * lambda;
  out.ground_energy = n == 1 ? Real{0} : edge * edge * Real(n) * Real(n - 1) * Real(n + 1) / Real{3};
  out.gap = unit * (unit * free_weight.value());
  if (interaction_weight.value() != Real{0}) out.gap += (unit * lambda) * (unit * interaction_weight.value());
  out.energy = out.ground_energy + out.gap;
  out.momentum = unit * Real(out.momentum_index);
  if (!uni20::isfinite(out.energy) || !uni20::isfinite(out.momentum) || !uni20::isfinite(out.gap) ||
      !std::all_of(out.pseudomomenta.begin(), out.pseudomomenta.end(), [](Real p) { return uni20::isfinite(p); }))
    throw std::overflow_error("Sutherland energy or momentum is outside the selected scalar range");
  if ((n > 1 && lambda > Real{0} && out.ground_energy == Real{0}) ||
      (out.gap == Real{0} && std::any_of(out.labels.begin(), out.labels.end(), [](auto label) { return label != 0; })))
    throw std::underflow_error("positive Sutherland energy or gap underflowed in the selected scalar type");
  return out;
}

template <uni20::Real Real> State<Real> ground_state(std::size_t n, Real length, Real lambda)
{
  detail::validate(n, length, lambda);
  return evaluate(n, length, lambda, std::vector<std::int64_t>(n, 0));
}

template <uni20::Real Real> struct Spectrum
{
    std::vector<State<Real>> states;
    std::size_t total_states = 0;
    bool complete = false; // Complete enumeration of the window, NOT the infinite spectrum.
};

/// Lowest count states in [-window,window]^N (nullopt means all in that window).
/// No global low-energy completeness outside the specified window is claimed.
template <uni20::Real Real>
Spectrum<Real> spectrum(std::size_t n, Real length, Real lambda, std::size_t window,
                        std::optional<std::size_t> count = {}, std::size_t max_states = 100000)
{
  detail::validate(n, length, lambda);
  if (window > 1000000) throw std::invalid_argument("label window must be <=1000000");
  if (count && *count == 0) throw std::invalid_argument("level count must be positive");
  Spectrum<Real> out;
  auto const total = detail::window_count(n, window, max_states);
  if (!total) return out;
  out.total_states = *total;
  out.states.reserve(*total);
  auto const w = static_cast<std::int64_t>(window);
  std::vector<std::int64_t> labels(n, -w);
  for (;;)
  {
    out.states.push_back(evaluate(n, length, lambda, labels));
    auto i = n;
    while (i && labels[i - 1] == w)
      --i;
    if (i == 0) break;
    std::fill(labels.begin() + (i - 1), labels.end(), labels[i - 1] + 1);
  }
  std::sort(out.states.begin(), out.states.end(), [](auto const& a, auto const& b) {
    if (a.gap != b.gap) return a.gap < b.gap;
    if (a.momentum_index != b.momentum_index) return a.momentum_index < b.momentum_index;
    return a.labels < b.labels;
  });
  if (count && *count < out.states.size()) out.states.resize(*count);
  out.complete = true;
  return out;
}
} // namespace bethe::sutherland
