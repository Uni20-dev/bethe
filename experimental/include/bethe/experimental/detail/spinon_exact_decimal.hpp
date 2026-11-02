// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/experimental/detail/spinon_checkpoint_codec.hpp>

namespace bethe::experimental::xxx::checkpoint_detail
{
// A round-trip decimal at geometry precision can lie on the wrong side of a
// lower-precision rounding tie. Export the exact dyadic value instead. For
// x=M*2^(e-p), the terminating decimal needs at most p+abs(e)+2 significant
// digits. Scientific formatting avoids a number of leading zeros set by e.
// Bound the requested allocation before formatting (MPFR exponents can be huge).
template <uni20::Real R> std::string exact_decimal(R const& x)
{
  auto const precision = uni20::precision_of(x);
  require(!precision.is_exact() && uni20::isfinite(x), "exact decimal requires finite approximate input");
  auto const exponent = uni20::math::frexp(x).exponent;
  require(exponent > -static_cast<std::int64_t>(max_string) && exponent < static_cast<std::int64_t>(max_string),
          "exact decimal exponent exceeds observation storage bound");
  auto const digits = precision.bit_count() + (exponent < 0 ? -exponent : exponent) + 2;
  require(digits > 0 && digits < static_cast<std::int64_t>(max_string - 64),
          "exact decimal precision exceeds observation storage bound");
  auto result = uni20::format_real(x, {.precision = static_cast<int>(digits),
                                       .notation = uni20::real_format_notation::scientific,
                                       .normalize_negative_zero = false});
  require(result.size() <= max_string, "exact decimal exceeds observation storage bound");
  return result;
}
} // namespace bethe::experimental::xxx::checkpoint_detail
