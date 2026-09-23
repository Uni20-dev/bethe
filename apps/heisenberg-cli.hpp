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

} // namespace bethe::cli
