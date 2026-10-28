// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include "cli-common.hpp"
#include <bethe/solver.hpp>
#include <optional>

namespace bethe::cli
{
template <uni20::Real Real, typename Options> Real parse_solver_real(std::string_view text, Options const& options)
{
  if constexpr (uni20::has_runtime_precision_v<Real>)
    return uni20::parse_real<Real>(text, options.precision());
  else
    return uni20::parse_real<Real>(text);
}

// Construct model controls at the requested precision. The compile-time
// capability keeps unavailable frontend/provider combinations from instantiating.
template <ArbitraryPrecision Support, template <typename> class Options = SolverOptions, typename Function>
int dispatch_options_precision(std::string_view precision, Function&& run)
{
  auto const request = parse_precision(precision);
  validate_precision(request, Support);
  if (request.arbitrary())
  {
#if UNI20_HAS_MPREAL
    if constexpr (has_arbitrary_precision(Support)) return run(Options<uni20::mpreal>(request.working_precision()));
#endif
    throw std::logic_error("unavailable arithmetic passed validation");
  }
  return dispatch_precision(precision, [&]<typename Real> { return run(Options<Real>{}); });
}

// Parse real-valued controls only after selecting the working precision, never
// through double or an ambient MPFR default.
template <ArbitraryPrecision Support, template <typename> class Options = SolverOptions, typename Function>
int dispatch_solver_precision(std::string_view precision, std::optional<std::string> const& tolerance,
                              std::size_t max_iterations, Function&& run)
{
  return dispatch_options_precision<Support, Options>(precision, [&]<typename Real>(Options<Real> options) {
    options.max_iterations = max_iterations;
    if (tolerance) options.residual_tolerance = parse_solver_real<Real>(*tolerance, options);
    return run(std::move(options));
  });
}
} // namespace bethe::cli
