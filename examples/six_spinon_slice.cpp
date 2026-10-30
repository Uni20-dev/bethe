// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include <bethe/xxx_six_spinon_structure_factor.hpp>
#include <iostream>
#include <uni20/async/tbb_scheduler.hpp>

// A physical integrand supplies already weighted values to the generic engine.
// Sharing sample identities across frequencies retains their covariance without
// teaching the Monte Carlo controller about spinons or importance proposals.
namespace
{
constexpr std::size_t count = 3, columns = 3 * count;
using Metadata = bethe::heisenberg::SpinonSampleMetadata;
struct SliceWork
{
    std::array<Metadata, count> points;
    SliceWork& operator+=(SliceWork const& other)
    {
      for (std::size_t i = 0; i < count; ++i)
        points[i] += other.points[i];
      return *this;
    }
};
} // namespace

int main()
{
  try
  {
    namespace model = bethe::heisenberg;
    using R = uni20::mpreal;
    R const one(1, uni20::Precision::bits(384)), pi = 4 * uni20::math::atan(one);
    using Kind = model::detail::SpinonProposalKind;
    model::detail::SpinonJointProposal<R> proposal(8 * one, one / 2,
                                                   {{Kind::background, one / 2, one, one},
                                                    {Kind::one_pair, one / 4, 4 * one, one / 2},
                                                    {Kind::two_pairs, one / 4, 8 * one, one}});
    model::SpinonRecoveryPolicy recovery{
        "1e-7", "0", {{128, 192, {{{32, 12}, {64, 20}}}}, {224, 320, {{{64, 20}, {128, 24}}}}}, 60};
    model::ThermodynamicSixSpinonStructureFactor<R, R> spectrum(proposal, recovery, one);
    std::array<R, count> omega{9 * one / 5, 2 * one, 11 * one / 5};
    std::array points{spectrum.at(pi / 2, omega[0]), spectrum.at(pi / 2, omega[1]), spectrum.at(pi / 2, omega[2])};
    auto observe = [points](bethe::UniformSample<5> const& identity, std::stop_token stop) {
      bethe::MonteCarloObservation<R, SliceWork> result;
      std::vector<R> values;
      values.reserve(columns);
      for (std::size_t i = 0; i < count; ++i)
      {
        auto sample = points[i].observation(identity, stop);
        result.metadata.points[i] = {sample.supported_roots, sample.density_attempts, sample.failure,
                                     sample.density_status, sample.density_failure};
        if (!sample.values)
        {
          result.cancelled = sample.failure == model::SpinonObservationFailure::cancelled;
          result.reason = "frequency " + std::to_string(i) + ": " + sample.reason;
          return result;
        }
        values.insert(values.end(), sample.values->begin(), sample.values->end());
      }
      result.values = std::move(values);
      return result;
    };
    uni20::async::TbbScheduler scheduler(2);
    bethe::MonteCarloOptions options{8, 64, 42, 0, 52};
    auto pending =
        bethe::monte_carlo_async<5>(options, columns, one, observe, {.scheduler = &scheduler, .max_in_flight = 2});
    auto result = pending.get_wait(scheduler);
    if (result.status != bethe::MonteCarloStatus::complete)
    {
      for (auto const& row : result.replicates)
        if (row && !row->means) std::cerr << "replicate " << row->id << ": " << row->reason << '\n';
      throw std::runtime_error("incomplete slice; no survivor estimates are printed");
    }
    auto number = [](R const& x) { return uni20::format_scalar(x, {.precision = 18}); };
    auto const& mean = *result.mean;
    auto const& covariance = *result.covariance;
    std::cout << "# S6zz, J=1, k=pi/2, |beta|<8, geometry=384 bits\n"
              << "# Eight shared-stream replicates of 64 draws; seed=42, uniform_bits=52\n"
              << "# Sampling errors exclude omitted rapidities; this small pilot is not a convergence claim\n"
              << "# omega S6 sampling_SE contour_diagnostic arithmetic_diagnostic supported_roots\n";
    for (std::size_t i = 0; i < count; ++i)
    {
      std::cout << number(omega[i]) << ' ';
      if (!result.metadata.points[i].supported_configurations)
        std::cout << "NA NA NA NA 0 # no support sampled, not a zero spectrum\n";
      else
        std::cout << number(mean[3 * i]) << ' ' << number(uni20::math::sqrt(covariance[3 * i * columns + 3 * i])) << ' '
                  << number(mean[3 * i + 1]) << ' ' << number(mean[3 * i + 2]) << ' '
                  << result.metadata.points[i].supported_configurations << '\n';
    }
    std::cout << "# Covariance of the three spectral means, in frequency order\n";
    for (std::size_t i = 0; i < count; ++i)
    {
      for (std::size_t j = 0; j < count; ++j)
      {
        if (j) std::cout << ' ';
        if (result.metadata.points[i].supported_configurations && result.metadata.points[j].supported_configurations)
          std::cout << number(covariance[3 * i * columns + 3 * j]);
        else
          std::cout << "NA";
      }
      std::cout << '\n';
    }
    std::cout << "# elapsed_seconds " << result.wall_seconds << '\n';
  }
  catch (std::exception const& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
