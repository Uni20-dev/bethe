// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#pragma once

#include <bethe/heisenberg_excitations.hpp>
#include <uni20/core/scalar_io.hpp>

#include <charconv>
#include <ctime>
#include <fmt/format.h>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace bethe::cli
{
#if UNI20_HAS_FLOAT128
// Other modes in the pinned Uni20 scalar I/O narrow fp128 to long double.
static_assert(MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_FLOAT128,
              "fp128 CLI output requires MPLAPACK's _Float128/strfromf128 mode");
#endif

// Process CPU time, excluding report construction and printing.
class CpuTimer {
  public:
    std::string elapsed_text() const
    {
      auto const end = std::clock();
      if (start_ == std::clock_t{-1} || end == std::clock_t{-1} || end < start_) return "unavailable";
      auto const seconds = (static_cast<long double>(end) - static_cast<long double>(start_)) / CLOCKS_PER_SEC;
      return fmt::format("{:.6f} s", seconds);
    }

  private:
    std::clock_t const start_ = std::clock();
};

inline std::size_t parse_size(std::string_view text)
{
  std::size_t value = 0;
  auto const parsed = std::from_chars(text.data(), text.data() + text.size(), value);
  if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size())
    throw std::invalid_argument("invalid nonnegative integer: " + std::string(text));
  return value;
}

struct ExcitationArguments
{
    std::optional<std::size_t> count;
    std::optional<uni20::half_int> spin;
    std::optional<std::size_t> max_candidates;

    bool parse(std::string_view option, std::string_view value)
    {
      if (option == "--excitations")
        count = parse_size(value);
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
        throw std::invalid_argument("--spin and --max-candidates require --excitations COUNT");
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
  out << "  --excitations COUNT                 lowest COUNT multiplets in a restricted real-root family\n"
      << "  --spin S                           total spin for that scan (default: 1 even N, 1/2 odd N)\n"
      << "  --max-candidates COUNT             exhaustive scan limit (default: 10000)\n"
      << "                                     includes the sector minimum; NOT a complete spectrum\n";
}

inline heisenberg::QuantumNumbers parse_quantum_numbers(std::string_view text)
{
  heisenberg::QuantumNumbers numbers;
  // An explicitly empty list is the fully polarized state.
  if (text.empty()) return numbers;
  for (;;)
  {
    auto const comma = text.find(',');
    numbers.push_back(uni20::half_int::parse(text.substr(0, comma)));
    if (comma == std::string_view::npos) return numbers;
    text.remove_prefix(comma + 1);
  }
}

template <typename State> void print_roots(State const& state)
{
  std::cout << "# index rapidity I\n";
  for (std::size_t i = 0; i < state.rapidities.size(); ++i)
    std::cout << i << ' ' << uni20::format_real(state.rapidities[i]) << ' ' << state.quantum_numbers[i] << '\n';
}

inline int finish(bool converged)
{
  if (!converged)
    std::cerr << "The iteration budget was exhausted; each nonconverged energy is an unconverged estimate.\n";
  return converged ? 0 : 2;
}
} // namespace bethe::cli
