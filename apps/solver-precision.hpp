// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include "cli-common.hpp"
#include <bethe/solver.hpp>
#include <optional>

namespace bethe::cli
{
// Parse real-valued controls only after selecting the working precision, never
// through double or an ambient MPFR default. The compile-time capability keeps
// unavailable frontend/provider combinations from being instantiated.
template <ArbitraryPrecision Support, typename Function>
int dispatch_solver_precision(std::string_view precision, std::optional<std::string> const& tolerance,
                              std::size_t max_iterations, Function&& run)
{
  auto const request = parse_precision(precision);
  validate_precision(request, Support);
  auto configured = [&]<typename Real>(SolverOptions<Real> options) {
    options.max_iterations = max_iterations;
    if (tolerance)
    {
#if UNI20_HAS_MPREAL
      if constexpr (uni20::has_runtime_precision_v<Real>)
        options.residual_tolerance = uni20::parse_real<Real>(*tolerance, options.precision());
      else
#endif
        options.residual_tolerance = uni20::parse_real<Real>(*tolerance);
    }
    return run(std::move(options));
  };
  if (request.arbitrary())
  {
#if UNI20_HAS_MPREAL
    if constexpr (has_arbitrary_precision(Support))
      return configured(SolverOptions<uni20::mpreal>(request.working_precision()));
#endif
    throw std::logic_error("unavailable arithmetic passed validation");
  }
  return dispatch_precision(precision, [&]<typename Real> { return configured(SolverOptions<Real>{}); });
}
} // namespace bethe::cli
