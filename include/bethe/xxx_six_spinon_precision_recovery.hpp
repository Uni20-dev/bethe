// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/xxx_six_spinon_structure_factor.hpp>
#include <map>
#include <tuple>

namespace bethe::heisenberg
{
struct SixSpinonPrecisionAttempt
{
    std::size_t supported_roots = 0, density_attempts = 0;
    SpinonObservationFailure failure = SpinonObservationFailure::none;
    std::optional<SpinonRecoveryStatus> density_status;
    std::optional<SpinonDensityFailure> density_failure;
    std::string reason;
};

template <uni20::Real R, std::size_t Stages> struct SixSpinonPrecisionResult
{
    SpinonSpectrumObservation<R> observation;
    std::array<SixSpinonPrecisionAttempt, Stages> attempts;
    std::size_t attempted_stages = 0;
    // Includes partial work in failed stages, unlike observation's final counts.
    std::size_t evaluated_roots = 0, density_attempts = 0;
};

namespace detail
{
inline void add_precision_work(std::size_t& count, std::size_t extra)
{
  if (extra > std::numeric_limits<std::size_t>::max() - count)
    throw std::length_error("precision-recovery work counter overflow");
  count += extra;
}
} // namespace detail

struct SixSpinonPrecisionFailureWork
{
    std::size_t attempts = 0, evaluated_roots = 0, density_attempts = 0;
    SixSpinonPrecisionFailureWork& operator+=(SixSpinonPrecisionFailureWork const& other)
    {
      detail::add_precision_work(attempts, other.attempts);
      detail::add_precision_work(evaluated_roots, other.evaluated_roots);
      detail::add_precision_work(density_attempts, other.density_attempts);
      return *this;
    }
};

struct SixSpinonPrecisionStageWork
{
    using FailureKey = std::tuple<SpinonObservationFailure, std::optional<SpinonRecoveryStatus>,
                                  std::optional<SpinonDensityFailure>, std::string>;
    std::size_t attempts = 0, accepted = 0, retried = 0;
    std::size_t evaluated_roots = 0, density_attempts = 0;
    // Successful observations allocate no failure entries. A retry means that
    // a later stage actually started, not merely that one was requested.
    std::map<FailureKey, SixSpinonPrecisionFailureWork> failures;

    SixSpinonPrecisionStageWork& operator+=(SixSpinonPrecisionStageWork const& other)
    {
      detail::add_precision_work(attempts, other.attempts);
      detail::add_precision_work(accepted, other.accepted);
      detail::add_precision_work(retried, other.retried);
      detail::add_precision_work(evaluated_roots, other.evaluated_roots);
      detail::add_precision_work(density_attempts, other.density_attempts);
      for (auto const& [key, work] : other.failures)
        failures[key] += work;
      return *this;
    }
};

/// The inherited counters describe the final physical observation. Summing
/// stage work includes unsuccessful attempts as well. These records follow
/// the Monte Carlo engine's retained-sample semantics, not speculative work
/// discarded after another sample has failed.
template <std::size_t Stages> struct SixSpinonPrecisionMetadata : SpinonSampleMetadata
{
    std::array<SixSpinonPrecisionStageWork, Stages> stages;
    SixSpinonPrecisionMetadata& operator+=(SixSpinonPrecisionMetadata const& other)
    {
      detail::add_precision_work(supported_configurations, other.supported_configurations);
      detail::add_precision_work(density_attempts, other.density_attempts);
      if (other.failure != SpinonObservationFailure::none)
      {
        failure = other.failure;
        density_status = other.density_status;
        density_failure = other.density_failure;
      }
      for (std::size_t i = 0; i < Stages; ++i)
        stages[i] += other.stages[i];
      return *this;
    }
};

/// Recover one physical observation through an explicitly ordered list of
/// immutable Point objects. Each Point must describe the SAME physical query,
/// proposal and controls, differing only in its declared precision policy.
/// Every stage reconstructs geometry from the original UniformSample. This
/// adapter never redraws, selects an importance proposal or changes tolerances.
/// The last Point is the reference for continuum classification; an earlier
/// stage that disagrees is retried before evaluating its observation.
///
/// The caller retains the stage policy for metadata; the result retains all
/// attempt outcomes. A failure remains a failure when the stages are exhausted.
/// Point callables must allow concurrent const calls. Output precision must
/// retain the accepted value's precision; output conversion is checked.
template <uni20::Real R, class... Points> class SixSpinonPrecisionRecovery {
  public:
    static constexpr std::size_t stage_count = sizeof...(Points);
    using Result = SixSpinonPrecisionResult<R, stage_count>;
    using Metadata = SixSpinonPrecisionMetadata<stage_count>;

    SixSpinonPrecisionRecovery(R exemplar, Points... points)
        : exemplar_(std::move(exemplar)), points_(std::move(points)...)
    {
      static_assert(stage_count > 0, "at least one precision stage is required");
      if constexpr (uni20::has_runtime_precision_v<R>)
        if (exemplar_.precision().is_exact())
          throw std::invalid_argument("observation accumulation needs finite precision");
      if (!uni20::isfinite(exemplar_)) throw std::invalid_argument("nonfinite output exemplar");
    }

    [[nodiscard]] Result observation(bethe::UniformSample<5> const& identity, std::stop_token stop = {}) const
    {
      Result result;
      evaluate<0>(identity, stop, result);
      return result;
    }

    /// Monte Carlo callback for Point::evaluate_observations_async(). Sampling
    /// and physical result assembly remain in their existing model/engine layers.
    [[nodiscard]] bethe::MonteCarloObservation<R, Metadata> operator()(bethe::UniformSample<5> const& identity,
                                                                       std::stop_token stop = {}) const
    {
      auto recovered = observation(identity, stop);
      auto& final = recovered.observation;
      bethe::MonteCarloObservation<R, Metadata> result;
      if (final.values)
        result.values.emplace(std::make_move_iterator(final.values->begin()),
                              std::make_move_iterator(final.values->end()));
      static_cast<SpinonSampleMetadata&>(result.metadata) = {
          final.supported_roots, final.density_attempts, final.failure, final.density_status, final.density_failure};
      result.cancelled = final.failure == SpinonObservationFailure::cancelled;
      result.reason = std::move(final.reason);
      for (std::size_t i = 0; i < recovered.attempted_stages; ++i)
      {
        auto& attempt = recovered.attempts[i];
        auto& work = result.metadata.stages[i];
        work.attempts = 1;
        work.accepted = i + 1 == recovered.attempted_stages && result.values.has_value();
        work.retried = i + 1 < recovered.attempted_stages;
        work.evaluated_roots = attempt.supported_roots;
        work.density_attempts = attempt.density_attempts;
        if (!work.accepted)
          work.failures[{attempt.failure, attempt.density_status, attempt.density_failure, std::move(attempt.reason)}] =
              {1, attempt.supported_roots, attempt.density_attempts};
      }
      return result;
    }

  private:
    R exemplar_;
    std::tuple<Points...> points_;

    template <class G> static bool retryable(SpinonSpectrumObservation<G> const& value)
    {
      if (value.values) return false;
      if (value.failure == SpinonObservationFailure::geometry ||
          value.failure == SpinonObservationFailure::exponent_range)
        return true;
      // Time limits, cancellation, invalid/unsupported input and inconsistent
      // controls are terminal; a wider format is not permission to ignore them.
      return value.failure == SpinonObservationFailure::density &&
             value.density_status == SpinonRecoveryStatus::exhausted;
    }

    template <class G> R import(G const& value) const
    {
      if (uni20::numeric_limits<R>::digits(exemplar_) < uni20::numeric_limits<G>::digits(value))
        throw std::invalid_argument("output precision cannot retain the accepted observation");
      R result = uni20::scalar_like(exemplar_, value);
      auto const zero = uni20::scalar_like(exemplar_, 0);
      if (!uni20::isfinite(result) || (value != G{0} && result == zero))
        throw std::range_error("accepted observation cannot be represented in output format");
      // Nominal significand width is insufficient for a native subnormal:
      // conversion can be finite and nonzero while losing low bits.
      if (uni20::scalar_like(value, result) != value)
        throw std::range_error("accepted observation loses precision in output format");
      return result;
    }

    template <std::size_t I>
    void evaluate(bethe::UniformSample<5> const& identity, std::stop_token stop, Result& result) const
    {
      if (stop.stop_requested())
      {
        result.observation.failure = SpinonObservationFailure::cancelled;
        result.observation.reason = "observation cancelled before next precision stage";
        return;
      }
      auto value = [&] {
        auto const& point = std::get<I>(points_);
        if (point.boundary_status() != std::get<stage_count - 1>(points_).boundary_status())
        {
          // Native input rounding can create a false successful zero, or a
          // lower-threshold exception, before energy-shell conditioning runs.
          decltype(point.observation(identity, stop)) mismatch;
          mismatch.failure = SpinonObservationFailure::geometry;
          mismatch.reason = "precision stage disagrees with reference continuum classification";
          return mismatch;
        }
        return point.observation(identity, stop);
      }();
      result.attempts[I] = {value.supported_roots, value.density_attempts, value.failure,
                            value.density_status,  value.density_failure,  value.reason};
      ++result.attempted_stages;
      detail::add_precision_work(result.evaluated_roots, value.supported_roots);
      detail::add_precision_work(result.density_attempts, value.density_attempts);
      if constexpr (I + 1 < stage_count)
      {
        if (retryable(value))
        {
          evaluate<I + 1>(identity, stop, result);
          return;
        }
      }
      auto& final = result.observation;
      final.supported_roots = value.supported_roots;
      final.density_attempts = value.density_attempts;
      final.failure = value.failure;
      final.density_status = value.density_status;
      final.density_failure = value.density_failure;
      final.reason = std::move(value.reason);
      if (value.values)
        final.values =
            std::array<R, 3>{import((*value.values)[0]), import((*value.values)[1]), import((*value.values)[2])};
    }
};

template <uni20::Real R, class... Points> auto make_six_spinon_precision_recovery(R exemplar, Points... points)
{
  return SixSpinonPrecisionRecovery<R, Points...>(std::move(exemplar), std::move(points)...);
}
} // namespace bethe::heisenberg
