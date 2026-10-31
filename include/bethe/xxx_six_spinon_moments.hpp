// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/xxx_spinon_proposal.hpp>
#include <bethe/xxx_spinon_observation.hpp>

namespace bethe::heisenberg
{
template <uni20::Real R> struct SixSpinonGlobalMoments
{
    // Columns: I0, I1, contour(I0), contour(I1), arithmetic(I0), arithmetic(I1).
    // Absolute longitudinal moments. First moment includes J. No fitted scales.
    bethe::MonteCarloEstimate<R, SpinonSampleMetadata> integral;
    R exchange, rapidity_cutoff, proposal_soft_fraction;
    std::vector<typename detail::SpinonJointProposal<R>::Component> proposal_components;
    SpinonRecoveryPolicy recovery;
    std::optional<std::array<R, 2>> fraction;
};

/// Direct global moments in the same finite |beta|<R box as the spectra.
/// Physical proposal and W6 normalization remain outside the numerical engine.
/// The continuum first-moment momentum law is not imposed on these samples.
template <uni20::Real Geometry, uni20::Real Low, uni20::Real High = Low> class ThermodynamicSixSpinonMoments {
  public:
    using G = Geometry;
    using Proposal = detail::SpinonJointProposal<G>;
    using Result = SixSpinonGlobalMoments<G>;
    using Observation = bethe::MonteCarloObservation<G, SpinonSampleMetadata>;

  private:
    struct Kernel
    {
        Proposal proposal;
        XXXSpinonDensityRecovery<Low, High> recovery;
        G exchange, zero, one, pi, log_measure, log_exchange;
        Kernel(Proposal p, SpinonRecoveryPolicy policy, G j)
            : proposal(std::move(p)), recovery(std::move(policy)), exchange(std::move(j)),
              zero(uni20::scalar_like(proposal.cutoff(), 0)), one(uni20::scalar_like(zero, 1)),
              pi(G{4} * uni20::math::atan(one)), log_measure(-G{2} * uni20::math::log(G{2} * pi))
        {
          bethe::detail::require_same_precision(exchange, one);
          if (!uni20::isfinite(exchange) || exchange <= zero) throw std::invalid_argument("require finite J>0");
          if (uni20::numeric_limits<G>::digits(one) < recovery.policy().stages.back().guard_bits)
            throw std::invalid_argument("moment accumulation must retain the largest density guard precision");
          log_exchange = uni20::math::log(exchange);
        }
    };
    std::shared_ptr<Kernel const> kernel_;

    using Integral = bethe::MonteCarloEstimate<G, SpinonSampleMetadata>;
    static uni20::async::AsyncTask finish(uni20::async::ReadBuffer<Integral> input,
                                          std::shared_ptr<Kernel const> kernel,
                                          uni20::async::WriteBuffer<Result> output)
    {
      auto destination = co_await output;
      auto const& integral = co_await input;
      Result result;
      result.integral = integral;
      result.exchange = kernel->exchange;
      result.rapidity_cutoff = kernel->proposal.cutoff();
      result.proposal_soft_fraction = kernel->proposal.soft_fraction();
      auto components = kernel->proposal.components();
      result.proposal_components.assign(components.begin(), components.end());
      result.recovery = kernel->recovery.policy();
      if (integral.mean)
      {
        G const exact_first =
            G{2} * kernel->exchange / G{3} * (uni20::math::log(G{2} * kernel->one) - kernel->one / G{4});
        if (!uni20::isfinite(exact_first) || exact_first <= kernel->zero)
          throw std::runtime_error("global first-moment reference exceeds exponent range");
        result.fraction = std::array{G{4} * (*integral.mean)[0], (*integral.mean)[1] / exact_first};
        for (auto const& value : *result.fraction)
          if (!uni20::isfinite(value)) throw std::runtime_error("global moment fraction exceeds exponent range");
      }
      destination = std::move(result);
    }

  public:
    ThermodynamicSixSpinonMoments(Proposal proposal, SpinonRecoveryPolicy recovery, G exchange)
        : kernel_(std::make_shared<Kernel>(std::move(proposal), std::move(recovery), std::move(exchange)))
    {}

    Observation global_observation(bethe::UniformSample<7> const& identity, std::stop_token stop = {}) const
    {
      Observation result;
      auto const& k = *kernel_;
      if (stop.stop_requested())
      {
        result.cancelled = true;
        result.metadata.failure = SpinonObservationFailure::cancelled;
        result.reason = "global observation cancelled";
        return result;
      }
      typename Proposal::GlobalSample point;
      auto const uniforms = identity.at_precision(k.one);
      try
      {
        point = k.proposal.draw_global(uniforms);
      }
      catch (std::runtime_error const& error)
      {
        result.metadata.failure = SpinonObservationFailure::geometry;
        result.reason = error.what();
        return result;
      }
      G const base = point.log_momentum_measure + k.log_measure;
      std::array<G, 2> const weights{base, base + point.log_energy + k.log_exchange};
      return weighted_spinon_observation<G>(point.beta, weights, k.recovery, stop);
    }

    [[nodiscard]] uni20::async::Async<Result> global_async(bethe::MonteCarloOptions options,
                                                           bethe::ExecutionOptions execution = {}) const
    {
      // Metadata counters can include up to two density controls per stage.
      auto const stages = kernel_->recovery.policy().stages.size();
      bethe::detail::validate_monte_carlo(options, 6, kernel_->one);
      if (options.samples_per_replicate > std::numeric_limits<std::size_t>::max() / options.replicates / (2 * stages))
        throw std::invalid_argument("global moment work counters overflow");
      auto& scheduler = execution.continuous_scheduler();
      auto work = bethe::monte_carlo_async<7>(
          options, 6, kernel_->one,
          [self = *this](auto const& identity, std::stop_token stop) {
            return self.global_observation(identity, stop);
          },
          execution);
      uni20::async::Async<Result> output;
      auto task = finish(work.read(), kernel_, output.write());
      bethe::detail::schedule_owned(scheduler, std::move(task));
      return output;
    }
    Result global(bethe::MonteCarloOptions options, bethe::ExecutionOptions execution = {}) const
    {
      auto& scheduler = execution.continuous_scheduler();
      return global_async(options, execution).get_wait(scheduler);
    }
};
} // namespace bethe::heisenberg
