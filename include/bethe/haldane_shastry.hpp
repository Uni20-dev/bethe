// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <algorithm>
#include <bethe/solver.hpp>
#include <bethe/spin_multiplets.hpp>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <uni20/common/half_int.hpp>
#include <vector>

namespace bethe::haldane_shastry
{
// This bound keeps all exact energy numerators within int64 and bounds the
// allocation for a single packed motif. Spectrum enumeration has its own cap.
inline constexpr std::size_t max_sites = 1000000;
namespace detail
{
inline void validate_sites(std::size_t n)
{
  if (n < 2 || n > max_sites) throw std::invalid_argument("Haldane-Shastry requires 2<=sites<=1000000");
}
inline void validate_motif(std::size_t n, std::vector<std::size_t> const& motif)
{
  validate_sites(n);
  std::size_t previous = 0;
  for (auto m : motif)
  {
    if (m == 0 || m >= n || (previous && (m <= previous || m - previous < 2)))
      throw std::invalid_argument("motif positions must increase in [1,N-1] with gaps >=2");
    previous = m;
  }
}
inline std::int64_t ground_numerator(std::size_t n)
{
  auto const m = static_cast<std::int64_t>(n);
  return -m * (m * m + (n % 2 ? -1 : 5));
}
// Number of motifs, with saturation represented by nullopt, not overflow.
inline std::optional<std::size_t> motif_count(std::size_t n, std::size_t cap)
{
  std::size_t previous = 1, current = 1;
  for (std::size_t k = 2; k <= n; ++k)
  {
    if (current > cap || previous > cap - current) return std::nullopt;
    auto const next = previous + current;
    previous = current;
    current = next;
  }
  return current;
}
} // namespace detail

/// One Yangian multiplet, not necessarily a single SU(2) spin multiplet.
template <uni20::Real Real> struct Level
{
    std::size_t sites = 0;
    std::vector<std::size_t> motif;
    std::int64_t energy_numerator = 0; // E = (pi/N)^2 * numerator/24.
    std::size_t momentum_index = 0, spinons = 0;
    uni20::half_int maximum_spin{0};
    std::optional<std::uint64_t> degeneracy; // Exact, or absent on count overflow.
    Real energy{}, gap{}, momentum{};        // Gap above the global AF ground energy.
};

/// Evaluate H=(pi/N)^2 sum_{i<j} S_i.S_j / sin^2(pi*(i-j)/N), J=1.
/// Motif positions are increasing in [1,N-1], with separation at least two.
template <uni20::Real Real> Level<Real> evaluate(std::size_t n, std::vector<std::size_t> motif)
{
  detail::validate_motif(n, motif);
  std::int64_t weight = 0;
  std::size_t momentum = 0, previous = 0;
  std::optional<std::uint64_t> dimension{1};
  auto multiply_dimension = [&](std::size_t factor) {
    if (!dimension) return;
    if (*dimension > std::numeric_limits<std::uint64_t>::max() / factor)
      dimension.reset();
    else
      *dimension *= factor;
  };
  for (std::size_t i = 0; i < motif.size(); ++i)
  {
    auto const m = motif[i];
    weight += static_cast<std::int64_t>(m) * static_cast<std::int64_t>(n - m);
    momentum = (momentum + m) % n;
    multiply_dimension(i ? m - previous - 1 : m);
    previous = m;
  }
  multiply_dimension(motif.empty() ? n + 1 : n - motif.back());
  Level<Real> out;
  out.sites = n;
  out.spinons = n - 2 * motif.size();
  out.maximum_spin = uni20::from_twice(static_cast<std::int64_t>(out.spinons));
  out.motif = std::move(motif);
  out.degeneracy = dimension;
  auto const ni = static_cast<std::int64_t>(n);
  out.energy_numerator = ni * (ni * ni - 1) - 24 * weight;
  out.momentum_index = momentum;
  Real const pi = Real{4} * std::atan(Real{1}), scale = (pi / Real(n)) * (pi / Real(n));
  out.energy = scale * Real(out.energy_numerator) / Real{24};
  // Subtract exact integers before scalar conversion, retaining small gaps.
  out.gap = scale * Real(out.energy_numerator - detail::ground_numerator(n)) / Real{24};
  out.momentum = Real{2} * pi * Real(momentum) / Real(n);
  return out;
}

/// SU(2) content of one Yangian multiplet. Consecutive unpaired sites in
/// the complement of motif U (motif+1) form symmetric spin-l/2 factors.
inline spin::Decomposition spin_decomposition(std::size_t n, std::vector<std::size_t> const& motif,
                                              spin::DecompositionOptions const& options = {})
{
  detail::validate_motif(n, motif);
  std::vector<uni20::half_int> factors;
  std::size_t first = 1;
  for (auto m : motif)
  {
    if (m > first) factors.push_back(uni20::from_twice(static_cast<std::int64_t>(m - first)));
    first = m + 2;
  }
  if (first <= n) factors.push_back(uni20::from_twice(static_cast<std::int64_t>(n - first + 1)));
  return spin::tensor_product(factors, options);
}

/// Return all minimizing motifs in a fixed Sz sector (one, or two chiral partners).
/// Degeneracy counts the whole Yangian multiplet, not its selected Sz slice.
template <uni20::Real Real> std::vector<Level<Real>> sector_ground_levels(std::size_t n, uni20::half_int sz)
{
  detail::validate_sites(n);
  auto const twice = sz.twice(), ni = static_cast<std::int64_t>(n);
  if (twice < -ni || twice > ni || (ni - twice) % 2 != 0)
    throw std::invalid_argument("Sz must lie in [-N/2,N/2] with the chain's half-integer parity");
  auto const abs_twice = static_cast<std::size_t>(twice < 0 ? -twice : twice);
  auto const m = (n - abs_twice) / 2;
  if (m == 0) return {evaluate<Real>(n, {})};
  std::vector<Level<Real>> out;
  auto const first = (n - 2 * m + 2) / 2;
  for (std::size_t shift = 0; shift <= n % 2; ++shift)
  {
    std::vector<std::size_t> motif(m);
    for (std::size_t j = 0; j < m; ++j)
      motif[j] = first + shift + 2 * j;
    out.push_back(evaluate<Real>(n, std::move(motif)));
  }
  return out;
}

template <uni20::Real Real> std::vector<Level<Real>> ground_levels(std::size_t n)
{
  return sector_ground_levels<Real>(n, uni20::from_twice(static_cast<std::int64_t>(n % 2)));
}

template <uni20::Real Real> struct Spectrum
{
    std::vector<Level<Real>> levels;
    std::size_t total_motifs = 0;
    bool complete = false; // Budget failure publishes no purported lowest levels.
};

/// Enumerate all motifs before selecting the lowest count (nullopt means all).
/// Exact integer ordering preserves degeneracies independently of precision.
template <uni20::Real Real>
Spectrum<Real> spectrum(std::size_t n, std::optional<std::size_t> count = {}, std::size_t max_motifs = 100000)
{
  detail::validate_sites(n);
  if (count && *count == 0) throw std::invalid_argument("level count must be positive");
  Spectrum<Real> out;
  auto const total = detail::motif_count(n, max_motifs);
  if (!total) return out;
  out.total_motifs = *total;
  out.levels.reserve(*total);
  std::vector<std::size_t> motif;
  // Recursion depth is bounded by the preflight Fibonacci count, not by N.
  auto visit = [&](auto&& self, std::size_t first) -> void {
    out.levels.push_back(evaluate<Real>(n, motif));
    for (std::size_t m = first; m < n; ++m)
    {
      motif.push_back(m);
      self(self, m + 2);
      motif.pop_back();
    }
  };
  visit(visit, 1);
  std::sort(out.levels.begin(), out.levels.end(), [](auto const& a, auto const& b) {
    if (a.energy_numerator != b.energy_numerator) return a.energy_numerator < b.energy_numerator;
    if (a.momentum_index != b.momentum_index) return a.momentum_index < b.momentum_index;
    return a.motif < b.motif;
  });
  if (count && *count < out.levels.size()) out.levels.resize(*count);
  out.complete = true;
  return out;
}
} // namespace bethe::haldane_shastry
