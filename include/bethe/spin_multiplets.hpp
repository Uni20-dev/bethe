// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <uni20/common/half_int.hpp>
#include <utility>
#include <vector>

namespace bethe::spin
{
struct Multiplet
{
    uni20::half_int spin{0};
    std::uint64_t multiplicity = 0; // Number of irreps, not magnetic states.
};
enum class DecompositionStatus
{
  complete,
  work_limit,
  count_overflow
};
struct DecompositionOptions
{
    std::size_t max_updates = 1000000;
};
struct Decomposition
{
    std::vector<Multiplet> multiplets;      // Increasing spin; absent on failure.
    std::optional<std::uint64_t> dimension; // May overflow even when multiplicities fit.
    std::size_t updates = 0;
    bool complete = false;
    DecompositionStatus status = DecompositionStatus::work_limit;
};

/// Exact SU(2) Clebsch-Gordan tensor product. The empty product is a singlet.
/// One update is one addition of an input multiplicity to an output spin.
/// No truncated decomposition is published on budget or multiplicity overflow.
inline Decomposition tensor_product(std::vector<uni20::half_int> const& factors,
                                    DecompositionOptions const& options = {})
{
  std::int64_t total = 0;
  for (auto s : factors)
  {
    if (s.twice() < 0 || s.twice() > std::numeric_limits<std::int64_t>::max() - total)
      throw std::invalid_argument("spin factors must be nonnegative with representable total spin");
    total += s.twice();
  }
  Decomposition out;
  std::map<std::int64_t, std::uint64_t> current{{0, 1}};
  for (auto factor : factors)
  {
    auto const b = factor.twice();
    if (b == 0) continue;
    std::map<std::int64_t, std::uint64_t> next;
    for (auto const& [a, count] : current)
    {
      auto const end = a + b;
      for (auto j = a > b ? a - b : b - a;; j += 2)
      {
        if (out.updates == options.max_updates) return out;
        ++out.updates;
        auto& value = next[j];
        if (count > std::numeric_limits<std::uint64_t>::max() - value)
        {
          out.status = DecompositionStatus::count_overflow;
          return out;
        }
        value += count;
        if (j == end) break;
      }
    }
    current = std::move(next);
  }
  out.dimension = 0;
  for (auto const& [twice, count] : current)
  {
    out.multiplets.push_back({uni20::from_twice(twice), count});
    auto const size = static_cast<std::uint64_t>(twice) + 1;
    if (out.dimension)
    {
      if (count > (std::numeric_limits<std::uint64_t>::max() - *out.dimension) / size)
        out.dimension.reset();
      else
        *out.dimension += count * size;
    }
  }
  out.complete = true;
  out.status = DecompositionStatus::complete;
  return out;
}
} // namespace bethe::spin
