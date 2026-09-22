// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <algorithm>
#include <bethe/solver.hpp>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <uni20/common/half_int.hpp>
#include <vector>

namespace bethe::xxz
{
template <uni20::Real Real> using SolverOptions = bethe::SolverOptions<Real>;
using QuantumNumbers = std::vector<uni20::half_int>;

namespace detail
{
inline std::int64_t checked_sites(std::size_t sites)
{
  if (sites < 2 || sites > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max() / 4))
    throw std::invalid_argument("XXZ chain requires 2 <= sites <= INT64_MAX/4");
  return static_cast<std::int64_t>(sites);
}

inline std::size_t sector_roots(std::size_t sites, uni20::half_int sz)
{
  auto const n = checked_sites(sites);
  auto const twice = sz.twice();
  if (twice < -n || twice > n || (n - twice) % 2 != 0)
    throw std::invalid_argument("Sz must lie in [-N/2,N/2] with the same half-integer parity as N/2");
  return static_cast<std::size_t>((n - (twice < 0 ? -twice : twice)) / 2);
}

inline std::size_t momentum_index(std::size_t sites, std::span<uni20::half_int const> numbers)
{
  auto const n = static_cast<std::int64_t>(sites);
  auto const modulus = 2 * n;
  std::int64_t index = numbers.size() % 2 == 0 ? 0 : n;
  for (auto const number : numbers)
  {
    index = (index - number.twice()) % modulus;
    if (index < 0) index += modulus;
  }
  return static_cast<std::size_t>(index / 2);
}
} // namespace detail
} // namespace bethe::xxz
