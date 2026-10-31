// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/spinon_weighted_value.hpp>
#include <bethe/monte_carlo.hpp>
#include <bethe/xxx_spinon_density_recovery.hpp>

namespace bethe::heisenberg
{
enum class SpinonObservationFailure
{
  none,
  geometry,
  density,
  exponent_range,
  cancelled
};

struct SpinonSampleMetadata
{
    std::size_t supported_configurations = 0, density_attempts = 0;
    SpinonObservationFailure failure = SpinonObservationFailure::none;
    std::optional<SpinonRecoveryStatus> density_status;
    std::optional<SpinonDensityFailure> density_failure;
    SpinonSampleMetadata& operator+=(SpinonSampleMetadata const& other)
    {
      supported_configurations += other.supported_configurations;
      density_attempts += other.density_attempts;
      if (other.failure != SpinonObservationFailure::none)
      {
        failure = other.failure;
        density_status = other.density_status;
        density_failure = other.density_failure;
      }
      return *this;
    }
};

} // namespace bethe::heisenberg

namespace bethe::heisenberg
{
/// Apply observable log weights to one checked physical W_n. The columns are
/// all values, then their contour diagnostics, then arithmetic diagnostics.
/// No physical measure is inferred here: the adapter supplies every Jacobian,
/// moment factor and exchange scale. Density recovery always sees the same
/// geometry coordinates, imported afresh at each declared precision.
template <uni20::Real G, uni20::Real Low, uni20::Real High>
bethe::MonteCarloObservation<G, SpinonSampleMetadata>
weighted_spinon_observation(std::span<G const> beta, std::span<G const> log_weights,
                            XXXSpinonDensityRecovery<Low, High> const& recovery, std::stop_token stop = {})
{
  if (log_weights.empty() || log_weights.size() > std::numeric_limits<std::size_t>::max() / 3)
    throw std::invalid_argument("require a bounded nonempty observable weight vector");
  G const one = uni20::scalar_like(log_weights.front(), 1), zero = uni20::scalar_like(one, 0);
  for (auto const& w : log_weights)
  {
    bethe::detail::require_same_precision(w, one);
    if (!uni20::isfinite(w)) throw std::invalid_argument("observable log weight must be finite");
  }
  bethe::MonteCarloObservation<G, SpinonSampleMetadata> result;
  result.metadata.supported_configurations = 1;
  auto recovered = recovery.recover(
      [&](auto const& context) {
        std::vector<decltype(context.parse("0"))> converted;
        for (auto const& b : beta)
          converted.push_back(context.import(b));
        return converted;
      },
      stop);
  result.metadata.density_attempts = recovered.attempts.size();
  if (!recovered.accepted())
  {
    result.cancelled = recovered.status == SpinonRecoveryStatus::cancelled;
    result.metadata.failure =
        result.cancelled ? SpinonObservationFailure::cancelled : SpinonObservationFailure::density;
    result.metadata.density_status = recovered.status;
    result.metadata.density_failure = recovered.density_failure;
    result.reason = recovered.reason;
    return result;
  }
  auto const& density = *recovered.estimate;
  std::array<High, 3> factors{density.value, density.contour_tail_error, *recovered.arithmetic_screen};
  std::vector<G> values(3 * log_weights.size(), zero);
  try
  {
    for (std::size_t field = 0; field < 3; ++field)
    {
      G const factor = uni20::scalar_like(one, factors[field]);
      if (field == 0 && factor <= zero)
        throw SpinonDensityError(SpinonDensityFailure::exponent_range, "checked physical density lost positivity");
      for (std::size_t i = 0; i < log_weights.size(); ++i)
        values[field * log_weights.size() + i] = detail::spinon_weighted_value(factor, log_weights[i]);
    }
  }
  catch (SpinonDensityError const& error)
  {
    result.metadata.failure = SpinonObservationFailure::exponent_range;
    result.reason = error.what();
    return result;
  }
  result.values = std::move(values);
  return result;
}
} // namespace bethe::heisenberg
