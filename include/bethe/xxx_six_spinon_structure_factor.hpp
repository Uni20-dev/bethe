// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/detail/spinon_sampling_statistics.hpp>
#include <bethe/detail/xxx_spinon_proposal.hpp>
#include <bethe/detail/xxx_spinon_surface.hpp>
#include <bethe/observable_grid.hpp>
#include <bethe/xxx_spinon_density_recovery.hpp>

namespace bethe::heisenberg
{
struct SpinonSamplingOptions
{
    std::size_t replicates = 8, samples_per_replicate = 256;
    std::uint64_t seed = 0, first_replicate = 0;
    unsigned uniform_bits = 52;
};

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

enum class SpinonObservationFailure
{
  none,
  geometry,
  density,
  exponent_range,
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
        Observation observation(detail::SpinonUniforms<5> const& identity, std::stop_token stop = {}) const
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
              auto recovered = k.recovery.recover(
                  [&](auto const& context) {
                    std::vector<decltype(context.parse("0"))> converted;
                    for (auto const& b : beta)
                      converted.push_back(context.import(b));
                    return converted;
                  },
                  stop);
              result.density_attempts += recovered.attempts.size();
              if (!recovered.accepted())
              {
                result.failure = recovered.status == SpinonRecoveryStatus::cancelled
                                     ? SpinonObservationFailure::cancelled
                                     : SpinonObservationFailure::density;
                result.density_status = recovered.status;
                result.reason = recovered.reason;
                if (!recovered.attempts.empty())
                {
                  result.density_failure = recovered.attempts.back().failure;
                  result.reason += ": " + recovered.attempts.back().reason;
                }
                return result;
              }
              auto const& density = *recovered.estimate;
              std::array<High, 3> factors{density.value, density.contour_tail_error, *recovered.arithmetic_screen};
              for (std::size_t field = 0; field < factors.size(); ++field)
              {
                G const value = uni20::scalar_like(k.one, factors[field]);
                if (value == k.zero && field == 1) continue; // an exact zero tail is allowed
                G const weighted = uni20::math::exp(uni20::math::log(value) + log_jacobian);
                if (!uni20::isfinite(weighted) || weighted <= k.zero)
                {
                  result.failure = SpinonObservationFailure::exponent_range;
                  result.reason = "weighted density or diagnostic is unrepresentable";
                  return result;
                }
                sums[field].add(weighted);
              }
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
          if (index >= options.replicates) throw std::out_of_range("replicate outside requested work");
          auto const& k = *kernel_;
          Replicate result;
          result.id = options.first_replicate + index;
          detail::SpinonSampleStream<5> stream(options.seed, result.id, options.uniform_bits);
          std::array<bethe::detail::SpinonSum<G>, 3> sums{bethe::detail::SpinonSum<G>(k.zero),
                                                          bethe::detail::SpinonSum<G>(k.zero),
                                                          bethe::detail::SpinonSum<G>(k.zero)};
          for (std::size_t i = 0; i < options.samples_per_replicate; ++i)
          {
            auto value = observation(stream.next(), stop);
            result.supported_roots += value.supported_roots;
            result.density_attempts += value.density_attempts;
            if (!value.values)
            {
              result.failed_sample = i;
              result.failure = value.failure;
              result.density_status = value.density_status;
              result.density_failure = value.density_failure;
              result.reason = value.reason;
              return result;
            }
            for (std::size_t field = 0; field < sums.size(); ++field)
              sums[field].add((*value.values)[field]);
            ++result.completed_samples;
          }
          G const count = uni20::scalar_like(k.one, options.samples_per_replicate);
          std::array<G, 3> means;
          for (std::size_t field = 0; field < means.size(); ++field)
          {
            G const total = sums[field].value();
            means[field] = total / count;
            if (!uni20::isfinite(means[field]) || (total != k.zero && means[field] == k.zero))
            {
              result.failure = SpinonObservationFailure::exponent_range;
              result.reason = "replicate mean exceeded exponent range";
              return result;
            }
          }
          result.means = std::move(means);
          return result;
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
          std::size_t count = boundary_ ? 0 : options.replicates;
          auto work = bethe::evaluate_observable_grid_async(
              count,
              [point = *this, options, run, stop = execution.stop](std::size_t i) {
                // Workers own distinct preallocated slots. Retain terminal
                // records even if cancellation prevents ordered collection;
                // finish reads them only after all admitted work has drained.
                run->result.replicates[i] = point.replicate(options, i, stop);
                return i;
              },
              [](std::size_t, std::size_t) {}, execution);
          uni20::async::Async<Result> output;
          auto task = finish(work.read(), std::move(run), output.write());
          bethe::detail::schedule_owned(scheduler, std::move(task));
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

        static uni20::async::AsyncTask finish(uni20::async::ReadBuffer<std::size_t> work, std::shared_ptr<Run> run,
                                              uni20::async::WriteBuffer<Result> output)
        {
          auto destination = co_await output;
          bool cancelled = false;
          try
          {
            (void)co_await work;
          }
          catch (bethe::ExecutionCancelled const&)
          {
            cancelled = true;
          }
          auto& r = run->result;
          auto const& k = *run->kernel;
          if (run->boundary && !cancelled)
          {
            r.status = *run->boundary;
            r.reason = "unrestricted continuum boundary";
            if (r.status != SixSpinonSpectrumStatus::lower_threshold)
              r.value = r.sampling_error = r.contour_diagnostic = r.arithmetic_diagnostic = k.zero;
          }
          else
          {
            detail::SpinonReplicateStatistics<G> statistics(r.sampling.replicates, r.sampling.samples_per_replicate, 3,
                                                            k.one);
            for (std::size_t i = 0; i < r.replicates.size(); ++i)
              if (auto const& row = r.replicates[i])
              {
                r.completed_samples += row->completed_samples;
                r.supported_roots += row->supported_roots;
                r.density_attempts += row->density_attempts;
                if (row->means)
                  statistics.accept(i, row->completed_samples, {row->means->begin(), row->means->end()});
                else
                  statistics.reject(i, row->reason);
              }
            std::optional<detail::SpinonReplicateSummary<G>> summary;
            try
            {
              summary = statistics.summary();
            }
            catch (std::runtime_error const& error)
            {
              r.reason = error.what(); // e.g. an unrepresentable nonzero variance
            }
            if (cancelled)
            {
              r.status = SixSpinonSpectrumStatus::cancelled;
              r.reason = "cancelled after draining admitted work";
            }
            else if (!summary || !summary->mean)
            {
              r.status = SixSpinonSpectrumStatus::numerical_failure;
              if (summary) r.reason = "incomplete or failed replicate; no survivor-only estimate";
            }
            else if (!r.supported_roots)
            {
              r.status = SixSpinonSpectrumStatus::no_support_sampled;
              r.reason = "no supported configuration was sampled; this is not a physical zero";
            }
            else
            {
              r.status = SixSpinonSpectrumStatus::estimated;
              r.value = (*summary->mean)[0];
              r.contour_diagnostic = (*summary->mean)[1];
              r.arithmetic_diagnostic = (*summary->mean)[2];
              if (summary->covariance) r.sampling_error = uni20::math::sqrt((*summary->covariance)[0]);
              r.reason = "finite-cutoff Monte Carlo estimate; numerical diagnostics exclude sampling and cutoff errors";
            }
          }
          r.wall_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - run->start).count();
          destination = std::move(r);
        }
    };

    Point at(G momentum, G omega) const { return Point(kernel_, std::move(momentum), std::move(omega)); }
};
} // namespace bethe::heisenberg
