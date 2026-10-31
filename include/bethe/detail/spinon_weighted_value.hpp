// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/xxx_spinon_density.hpp>

namespace bethe::heisenberg::detail
{
// Assemble a positive density or diagnostic with its integration measure.
// A vanishing input diagnostic is allowed; a nonzero contribution lost to the
// exponent range is a failure, including native subnormal relative arithmetic.
template <uni20::Real R> R spinon_weighted_value(R const& value, R const& log_weight)
{
  bethe::detail::require_same_precision(value, log_weight);
  if (!uni20::isfinite(value) || value < R{0} || !uni20::isfinite(log_weight))
    throw SpinonDensityError(SpinonDensityFailure::exponent_range, "invalid weighted spinon contribution");
  if (value == R{0}) return value;
  R const result = uni20::math::exp(uni20::math::log(value) + log_weight);
  if (!uni20::isfinite(result) || result <= R{0})
    throw SpinonDensityError(SpinonDensityFailure::exponent_range,
                             "weighted spinon contribution exceeds exponent range");
  if constexpr (!uni20::has_runtime_precision_v<R>)
    if (result < uni20::numeric_limits<R>::min())
      throw SpinonDensityError(SpinonDensityFailure::exponent_range,
                               "subnormal spinon contribution requires a wider exponent range");
  return result;
}
} // namespace bethe::heisenberg::detail
