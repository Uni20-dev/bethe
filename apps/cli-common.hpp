// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#pragma once

#include "precision.hpp"
#include <uni20/common/half_int.hpp>
#include <uni20/core/scalar_io.hpp>

#include <charconv>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace bethe::cli
{
#if UNI20_HAS_FLOAT128
// Other modes in the pinned Uni20 scalar I/O narrow fp128 to long double.
static_assert(MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_FLOAT128,
              "fp128 CLI output requires MPLAPACK's _Float128/strfromf128 mode");
#endif

template <typename Function> int dispatch_precision(std::string_view precision, Function&& run)
{
  auto const request = parse_precision(precision);
  validate_precision(request, false);
  switch (request.kind)
  {
    case PrecisionKind::fp64:
      return run.template operator()<double>();
    case PrecisionKind::long_double:
    case PrecisionKind::fp80:
      return run.template operator()<long double>();
    case PrecisionKind::fp128:
#if UNI20_HAS_FLOAT128
      return run.template operator()<uni20::float128>();
#else
      break;
#endif
    case PrecisionKind::bits:
    case PrecisionKind::digits:
      break;
  }
  throw std::logic_error("unavailable native precision");
}

inline std::size_t parse_size(std::string_view text)
{
  std::size_t value = 0;
  auto const parsed = std::from_chars(text.data(), text.data() + text.size(), value);
  if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size())
    throw std::invalid_argument("invalid nonnegative integer: " + std::string(text));
  return value;
}

inline std::vector<uni20::half_int> parse_quantum_numbers(std::string_view text)
{
  std::vector<uni20::half_int> numbers;
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

inline int finish(bool converged)
{
  if (!converged)
    std::cerr << "The calculation did not converge; each nonconverged energy is an unconverged estimate.\n";
  return converged ? 0 : 2;
}
} // namespace bethe::cli
