// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/detail/xxx_spinon_proposal.hpp>
#include <bethe/detail/xxx_spinon_surface.hpp>
#include <bethe/monte_carlo.hpp>
#include <bethe/observable_grid.hpp>
#include <bethe/xxx_spinon_observation.hpp>

namespace bethe::heisenberg
{
using SpinonSamplingOptions = bethe::MonteCarloOptions;

enum class SixSpinonSpectrumStatus
{
  estimated,
  outside_continuum,
  upper_edge,
  lower_threshold,
  no_support_sampled,
  numerical_failure,
  cancelled
};

template <uni20::Real R> struct SpinonSpectrumObservation
{
    // S6zz, absolute contour diagnostic, absolute arithmetic diagnostic.
    // No extra n! or 2pi factor. Empty means numerical failure, never zero.
    std::optional<std::array<R, 3>> values;
    std::size_t supported_roots = 0, density_attempts = 0;
    SpinonObservationFailure failure = SpinonObservationFailure::none;
    std::optional<SpinonRecoveryStatus> density_status;
    std::optional<SpinonDensityFailure> density_failure;
    std::string reason;
};

template <uni20::Real R> struct SpinonSpectrumReplicate
{
    std::uint64_t id = 0;
    std::size_t completed_samples = 0, supported_roots = 0, density_attempts = 0;
    std::optional<std::array<R, 3>> means;
    std::optional<std::array<R, 3>> maximum_absolute_observation;
    std::optional<std::size_t> failed_sample;
    SpinonObservationFailure failure = SpinonObservationFailure::none;
    std::optional<SpinonRecoveryStatus> density_status;
    std::optional<SpinonDensityFailure> density_failure;
    std::string reason;
};

template <uni20::Real R> struct SixSpinonSpectrumEstimate
{
    SixSpinonSpectrumStatus status = SixSpinonSpectrumStatus::numerical_failure;
    std::optional<R> value, sampling_error, contour_diagnostic, arithmetic_diagnostic;
    R momentum, omega, exchange, rapidity_cutoff, proposal_soft_fraction;
    std::vector<typename detail::SpinonJointProposal<R>::Component> proposal_components;
    SpinonSamplingOptions sampling;
    SpinonRecoveryPolicy recovery;
    std::vector<std::optional<SpinonSpectrumReplicate<R>>> replicates;
    std::size_t completed_samples = 0, supported_roots = 0, density_attempts = 0;
    double wall_seconds = 0;
    std::string reason;
};

/// Infinite XXX chain, longitudinal S6zz at J>0 and a FINITE rapidity cutoff.
/// Geometry/accumulation uses Geometry, independently of the density's Low/High
/// recovery contexts. The geometry precision must hold every density result.
/// Fresh ordinary Monte Carlo, not a claim of thermodynamic-cutoff convergence.
template <uni20::Real Geometry, uni20::Real Low, uni20::Real High = Low> class ThermodynamicSixSpinonStructureFactor {
  public:
    using G = Geometry;
    using Proposal = detail::SpinonJointProposal<G>;
    using Observation = SpinonSpectrumObservation<G>;
    using Replicate = SpinonSpectrumReplicate<G>;
    using Result = SixSpinonSpectrumEstimate<G>;

  private:
    struct Kernel
    {
        Proposal proposal;
        XXXSpinonDensityRecovery<Low, High> recovery;
        G exchange, zero, one, pi;
        Kernel(Proposal p, SpinonRecoveryPolicy policy, G j)
            : proposal(std::move(p)), recovery(std::move(policy)), exchange(std::move(j)),
              zero(uni20::scalar_like(proposal.cutoff(), 0)), one(uni20::scalar_like(zero, 1)),
              pi(G{4} * uni20::math::atan(one))
        {
          bethe::detail::require_same_precision(exchange, one);
          if (!uni20::isfinite(exchange) || exchange <= zero) throw std::invalid_argument("require finite J>0");
          if (uni20::numeric_limits<G>::digits(one) < recovery.policy().stages.back().guard_bits)
            throw std::invalid_argument("geometry/accumulation must retain the largest density guard precision");
        }
    };
    std::shared_ptr<Kernel const> kernel_;

  public:
    ThermodynamicSixSpinonStructureFactor(Proposal proposal, SpinonRecoveryPolicy recovery, G exchange)
        : kernel_(std::make_shared<Kernel>(std::move(proposal), std::move(recovery), std::move(exchange)))
    {}

    class Point {
      public:
        Observation observation(bethe::UniformSample<5> const& identity, std::stop_token stop = {}) const
        {
          Observation result;
          auto const& k = *kernel_;
          if (stop.stop_requested())
          {
            result.failure = SpinonObservationFailure::cancelled;
            result.reason = "observation cancelled";
            return result;
          }
          if (boundary_)
          {
            if (*boundary_ == SixSpinonSpectrumStatus::lower_threshold)
              throw std::domain_error("the unrestricted lower threshold is not a finite spectral point");
            result.values = std::array{k.zero, k.zero, k.zero};
            return result;
          }
          std::array<bethe::detail::SpinonSum<G>, 3> sums{bethe::detail::SpinonSum<G>(k.zero),
                                                          bethe::detail::SpinonSum<G>(k.zero),
                                                          bethe::detail::SpinonSum<G>(k.zero)};
          auto const uniforms = identity.at_precision(k.one);
          std::array<G, 4> free;
          try
          {
            free = k.proposal.draw(uniforms).beta;
          }
          catch (std::runtime_error const& error)
          {
            result.failure = SpinonObservationFailure::geometry;
            result.reason = error.what();
            return result;
          }
          for (auto const& surface : surfaces_)
          {
            typename detail::SpinonEnergySurface<G>::Point point;
            try
            {
              point = surface.solve(free);
            }
            catch (std::runtime_error const& error)
            {
              result.failure = SpinonObservationFailure::geometry;
              result.reason = error.what();
              return result;
            }
            if (!point.inside) continue;
            // Both labelled roots contribute. W6 owns the single factorial.
            for (auto const& beta : point.beta)
            {
              ++result.supported_roots;
              G log_jacobian;
              try
              {
                std::array<G, 6> coordinates;
                std::copy(beta.begin(), beta.end(), coordinates.begin());
                log_jacobian = k.proposal.balanced_log_jacobian(coordinates) - log_exchange_;
              }
              catch (std::runtime_error const& error)
              {
                result.failure = SpinonObservationFailure::geometry;
                result.reason = error.what();
                return result;
              }
              std::array<G, 1> const weights{log_jacobian};
              auto weighted = weighted_spinon_observation<G>(beta, weights, k.recovery, stop);
              result.density_attempts += weighted.metadata.density_attempts;
              if (!weighted.values)
              {
                result.failure = weighted.metadata.failure;
                result.density_status = weighted.metadata.density_status;
                result.density_failure = weighted.metadata.density_failure;
                result.reason = std::move(weighted.reason);
                return result;
              }
              for (std::size_t field = 0; field < 3; ++field)
                sums[field].add((*weighted.values)[field]);
            }
          }
          std::array<G, 3> values{sums[0].value(), sums[1].value(), sums[2].value()};
          for (auto const& value : values)
            if (!uni20::isfinite(value))
            {
              result.failure = SpinonObservationFailure::exponent_range;
              result.reason = "observation sum exceeded exponent range";
              return result;
            }
          result.values = std::move(values);
          return result;
        }

        Replicate replicate(SpinonSamplingOptions const& options, std::size_t index, std::stop_token stop = {}) const
        {
          validate(options);
          auto row = bethe::monte_carlo_replicate<5>(
              options, index, 3, kernel_->one,
              [point = *this](auto const& identity, std::stop_token token) { return point.sample(identity, token); },
              stop);
          return physical_record(row);
        }

        /// Owns numerical state and callables until all admitted replicates drain.
        /// The borrowed scheduler must outlive completion, even if the Async
        /// handle is dropped. No global scheduler/precision or thread pool changes.
        [[nodiscard]] uni20::async::Async<Result> evaluate_async(SpinonSamplingOptions options,
                                                                 bethe::ExecutionOptions execution = {}) const
        {
          validate(options);
          auto& scheduler = execution.continuous_scheduler();
          auto run = std::make_shared<Run>(*this, options);
          uni20::async::Async<Result> output;
          if (boundary_)
          {
            auto task = finish_boundary(std::move(run), execution.stop, output.write());
            bethe::detail::schedule_owned(scheduler, std::move(task));
          }
          else
          {
            auto work = bethe::monte_carlo_async<5>(
                options, 3, kernel_->one,
                [point = *this](auto const& identity, std::stop_token token) { return point.sample(identity, token); },
                execution);
            auto task = finish(work.read(), std::move(run), output.write());
            bethe::detail::schedule_owned(scheduler, std::move(task));
          }
          return output;
        }

        Result evaluate(SpinonSamplingOptions options, bethe::ExecutionOptions execution = {}) const
        {
          auto& scheduler = execution.continuous_scheduler();
          auto pending = evaluate_async(options, execution);
          return pending.get_wait(scheduler);
        }

      private:
        friend class ThermodynamicSixSpinonStructureFactor;
        std::shared_ptr<Kernel const> kernel_;
        G momentum_, omega_, log_exchange_;
        std::optional<SixSpinonSpectrumStatus> boundary_;
        std::vector<detail::SpinonEnergySurface<G>> surfaces_;

        Point(std::shared_ptr<Kernel const> kernel, G momentum, G omega)
            : kernel_(std::move(kernel)), momentum_(std::move(momentum)), omega_(std::move(omega)),
              log_exchange_(uni20::math::log(kernel_->exchange))
        {
          auto const& k = *kernel_;
          bethe::detail::require_same_precision(momentum_, k.one);
          bethe::detail::require_same_precision(omega_, k.one);
          if (!uni20::isfinite(momentum_) || !uni20::isfinite(omega_) || momentum_ < k.zero || momentum_ > G{2} * k.pi)
            throw std::invalid_argument("require finite frequency and momentum in [0,2*pi]");
          G const q = std::min(momentum_, G{2} * k.pi - momentum_);
          G const lower = q == k.pi ? k.zero : k.exchange * k.pi / G{2} * uni20::math::sin(q);
          G upper = k.zero;
          for (int section = 0; section < 3; ++section)
            upper = std::max(upper,
                             k.exchange * G{3} * k.pi * uni20::math::sin((momentum_ + G(2 * section) * k.pi) / G{6}));
          if (!uni20::isfinite(lower) || !uni20::isfinite(upper))
            throw std::overflow_error("spectral support exceeds range");
          if (q > k.zero && q != k.pi && lower <= k.zero)
            throw std::underflow_error("lower spectral boundary lost exchange scale");
          if (q == k.zero || omega_ < lower || omega_ > upper)
            boundary_ = SixSpinonSpectrumStatus::outside_continuum;
          else if (omega_ == lower)
            boundary_ = SixSpinonSpectrumStatus::lower_threshold;
          else if (omega_ == upper)
            boundary_ = SixSpinonSpectrumStatus::upper_edge;
          if (!boundary_)
          {
            G const normalized = omega_ / k.exchange;
            if (!uni20::isfinite(normalized) || normalized <= k.zero)
              throw std::underflow_error("frequency lost scale");
            for (int section = 0; section < 3; ++section)
              surfaces_.emplace_back(6, k.proposal.cutoff(), momentum_, section, normalized);
          }
        }

        using SamplingResult = bethe::MonteCarloEstimate<G, SpinonSampleMetadata>;
        using SamplingRow = bethe::MonteCarloReplicate<G, SpinonSampleMetadata>;

        bethe::MonteCarloObservation<G, SpinonSampleMetadata> sample(bethe::UniformSample<5> const& identity,
                                                                     std::stop_token stop) const
        {
          auto value = observation(identity, stop);
          bethe::MonteCarloObservation<G, SpinonSampleMetadata> result;
          if (value.values) result.values.emplace(value.values->begin(), value.values->end());
          result.metadata = {value.supported_roots, value.density_attempts, value.failure, value.density_status,
                             value.density_failure};
          result.cancelled = value.failure == SpinonObservationFailure::cancelled;
          result.reason = std::move(value.reason);
          return result;
        }

        static Replicate physical_record(SamplingRow const& row)
        {
          Replicate result;
          result.id = row.id;
          result.completed_samples = row.completed_samples;
          result.supported_roots = row.metadata.supported_configurations;
          result.density_attempts = row.metadata.density_attempts;
          result.failure = row.cancelled ? SpinonObservationFailure::cancelled : row.metadata.failure;
          if (!row.means && !row.cancelled && result.failure == SpinonObservationFailure::none)
            result.failure = SpinonObservationFailure::exponent_range;
          result.density_status = row.metadata.density_status;
          result.density_failure = row.metadata.density_failure;
          result.failed_sample = row.failed_sample;
          result.reason = row.reason;
          if (!row.maximum_absolute_observation.empty())
            result.maximum_absolute_observation =
                std::array{row.maximum_absolute_observation[0], row.maximum_absolute_observation[1],
                           row.maximum_absolute_observation[2]};
          if (row.means) result.means = std::array{(*row.means)[0], (*row.means)[1], (*row.means)[2]};
          return result;
        }

        void validate(SpinonSamplingOptions const& options) const
        {
          if (!options.replicates || !options.samples_per_replicate || !options.uniform_bits ||
              options.uniform_bits >= unsigned(std::numeric_limits<int>::max()) ||
              uni20::numeric_limits<G>::digits(kernel_->one) < options.uniform_bits + 1 ||
              options.replicates - 1 > std::numeric_limits<std::uint64_t>::max() - options.first_replicate ||
              options.samples_per_replicate > std::numeric_limits<std::size_t>::max() / options.replicates / 6 /
                                                  (2 * kernel_->recovery.policy().stages.size()))
            throw std::invalid_argument("invalid sampling budget/precision or overflowing work counters");
        }

        struct Run
        {
            std::shared_ptr<Kernel const> kernel;
            std::optional<SixSpinonSpectrumStatus> boundary;
            Result result;
            std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
            Run(Point const& point, SpinonSamplingOptions options) : kernel(point.kernel_), boundary(point.boundary_)
            {
              auto const& k = *kernel;
              result.momentum = point.momentum_;
              result.omega = point.omega_;
              result.exchange = k.exchange;
              result.rapidity_cutoff = k.proposal.cutoff();
              result.proposal_soft_fraction = k.proposal.soft_fraction();
              auto components = k.proposal.components();
              result.proposal_components.assign(components.begin(), components.end());
              result.sampling = options;
              result.recovery = k.recovery.policy();
              result.replicates.resize(boundary ? 0 : options.replicates);
            }
        };

        static uni20::async::AsyncTask finish_boundary(std::shared_ptr<Run> run, std::stop_token stop,
                                                       uni20::async::WriteBuffer<Result> output)
        {
          auto destination = co_await output;
          auto& r = run->result;
          if (stop.stop_requested())
          {
            r.status = SixSpinonSpectrumStatus::cancelled;
            r.reason = "cancelled before boundary evaluation";
          }
          else
          {
            r.status = *run->boundary;
            r.reason = "unrestricted continuum boundary";
            if (r.status != SixSpinonSpectrumStatus::lower_threshold)
              r.value = r.sampling_error = r.contour_diagnostic = r.arithmetic_diagnostic = run->kernel->zero;
          }
          r.wall_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - run->start).count();
          destination = std::move(r);
        }

        static uni20::async::AsyncTask finish(uni20::async::ReadBuffer<SamplingResult> work, std::shared_ptr<Run> run,
                                              uni20::async::WriteBuffer<Result> output)
        {
          auto destination = co_await output;
          auto const& sample = co_await work;
          auto& r = run->result;
          r.completed_samples = sample.completed_samples;
          r.supported_roots = sample.metadata.supported_configurations;
          r.density_attempts = sample.metadata.density_attempts;
          for (std::size_t i = 0; i < sample.replicates.size(); ++i)
            if (sample.replicates[i]) r.replicates[i] = physical_record(*sample.replicates[i]);
          if (sample.status == bethe::MonteCarloStatus::cancelled)
          {
            r.status = SixSpinonSpectrumStatus::cancelled;
            r.reason = sample.reason;
          }
          else if (sample.status != bethe::MonteCarloStatus::complete)
          {
            r.status = SixSpinonSpectrumStatus::numerical_failure;
            r.reason = sample.reason;
          }
          else if (!r.supported_roots)
          {
            r.status = SixSpinonSpectrumStatus::no_support_sampled;
            r.reason = "no supported configuration was sampled; this is not a physical zero";
          }
          else
          {
            r.status = SixSpinonSpectrumStatus::estimated;
            r.value = (*sample.mean)[0];
            r.contour_diagnostic = (*sample.mean)[1];
            r.arithmetic_diagnostic = (*sample.mean)[2];
            if (sample.covariance) r.sampling_error = uni20::math::sqrt((*sample.covariance)[0]);
            r.reason = "finite-cutoff Monte Carlo estimate; numerical diagnostics exclude sampling and cutoff errors";
          }
          r.wall_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - run->start).count();
          destination = std::move(r);
        }
    };

    Point at(G momentum, G omega) const { return Point(kernel_, std::move(momentum), std::move(omega)); }
};
} // namespace bethe::heisenberg
