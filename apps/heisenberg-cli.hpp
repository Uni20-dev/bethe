// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#pragma once

#include "cli-common.hpp"
#include <bethe/heisenberg_excitations.hpp>

#include <limits>
#include <optional>

namespace bethe::cli
{
struct ExcitationArguments
{
    std::optional<std::size_t> count;
    std::optional<uni20::half_int> spin;
    std::optional<std::size_t> max_candidates;

    bool parse(std::string_view option, std::string_view value)
    {
      if (option == "--excitations")
        // The scan clamps this to the family size after checking max_candidates,
        // before allocating retained states. "all" does not bypass that limit.
        count = value == "all" ? std::numeric_limits<std::size_t>::max() : parse_size(value);
      else if (option == "--spin")
        spin = uni20::half_int::parse(value);
      else if (option == "--max-candidates")
        max_candidates = parse_size(value);
      else
        return false;
      return true;
    }

    void validate() const
    {
      if (!count && (spin || max_candidates))
        throw std::invalid_argument("--spin and --max-candidates require --excitations COUNT|all");
    }

    uni20::half_int selected_spin(std::size_t sites) const
    {
      return spin.value_or(uni20::from_twice(std::int64_t{sites % 2 == 0 ? 2 : 1}));
    }

    heisenberg::RealExcitationOptions options() const
    {
      return {.count = count.value_or(10), .max_candidates = max_candidates.value_or(10000)};
    }
};

inline void excitation_usage(std::ostream& out)
{
  out << "  --excitations COUNT|all             lowest COUNT, or all, multiplets in a restricted real-root family\n"
      << "  --spin S                           total spin for that scan (default: 1 even N, 1/2 odd N)\n"
      << "  --max-candidates COUNT             exhaustive scan limit (default: 10000)\n"
      << "                                     includes the sector minimum; NOT a complete spectrum\n";
}

} // namespace bethe::cli
