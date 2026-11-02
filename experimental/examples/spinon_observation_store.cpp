// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include <bethe/experimental/spinon_observation_store.hpp>
#include <iostream>
#include <uni20/async/tbb_scheduler.hpp>

int main(int argc, char** argv)
try
{
  namespace x = bethe::experimental::xxx;
  namespace h = bethe::heisenberg;
  if (argc < 2 || argc > 3 || (argc == 3 && std::string_view(argv[2]) != "--stop-after-control"))
  {
    std::cerr << "Usage: bethe-experimental-observation-store-example DIRECTORY [--stop-after-control]\n";
    return 1;
  }
  int const bits = uni20::numeric_limits<long double>::digits();
  if (bits <= 53)
  {
    std::cerr << "This example requires long double wider than binary64.\n";
    return 77;
  }
  using Group = x::SixSpinonObservationGroup<long double>;
  using Store = x::SpinonObservationStore<long double, double, long double>;
  long double const pi = 4 * std::atan(1.L);
  Group::Proposal proposal(8.L, .5L,
                           {{h::detail::SpinonProposalKind::background, .5L, 1.L, 1.L},
                            {h::detail::SpinonProposalKind::one_pair, .25L, 4.L, .5L},
                            {h::detail::SpinonProposalKind::two_pairs, .25L, 8.L, 1.L}});
  std::vector<x::SpinonSectionRequest<long double>> requests;
  for (int s = 0; s < 3; ++s)
  {
    requests.push_back({"whole-" + std::to_string(s), pi / 2, 2.L, 0.L, 8.L, s});
    requests.push_back({"compact-" + std::to_string(s), pi / 2, 2.L, 0.L, 4.L, s});
    requests.push_back({"remainder-" + std::to_string(s), pi / 2, 2.L, 4.L, 8.L, s});
  }
  bethe::MonteCarloStream<5> stream(42, 0, 24);
  std::vector<x::SpinonObservationIdentity> samples;
  for (std::size_t i = 0; i < 12; ++i)
    samples.push_back({"sample-" + std::to_string(i), stream.next()});
  h::SpinonRecoveryPolicy policy{"1e-5", "0", {{53, bits, {{{32, 12}, {64, 20}}}}}, 60};
  auto store = std::make_shared<Store>(argv[1], Group(proposal, 1.L, requests), samples, policy,
                                       "bethe-observation-example-v1", "seed42 replicate0 uniform_bits24; J=1");
  if (argc == 3)
  {
    std::stop_source stop;
    for (std::size_t i = 0; i < store->size(); ++i)
    {
      auto result = store->run(i, [&](auto, auto, auto const&) { stop.request_stop(); }, stop.get_token());
      if (result.cancelled)
      {
        std::cout << "Stopped after one newly completed density control; rerun without the stop option.\n";
        return 2;
      }
    }
  }
  std::size_t replayed = 0, failed = 0;
  for (std::size_t i = 0; i < store->size(); ++i)
    replayed += store->load(i).has_value();
  std::vector<Store::Observation> observations;
  uni20::async::TbbScheduler scheduler(2);
  auto pending = x::evaluate_spinon_observation_chunk_async(store,
                                                            [&](std::size_t i, Store::Observation observation) {
                                                              if (!observation.values)
                                                              {
                                                                ++failed;
                                                                std::cout << i << ": " << observation.reason << '\n';
                                                              }
                                                              observations.push_back(std::move(observation));
                                                            },
                                                            {.scheduler = &scheduler, .max_in_flight = 2});
  auto count = pending.get_wait(scheduler);
  std::cout << count << " original samples; " << replayed << " committed observations replayed; " << failed
            << " failed observations\n";
  if (failed)
  {
    std::cout << "No mean or covariance: the original sample set is incomplete.\n";
    return 2;
  }
  bethe::ReplicateStatistics<long double> statistics(count, 1, requests.size(), 1.L);
  for (std::size_t i = 0; i < count; ++i)
    statistics.accept(i, 1, {observations[i].values->begin(), observations[i].values->begin() + requests.size()});
  auto result = statistics.summary();
  auto const whole = (*result.mean)[0] + (*result.mean)[3] + (*result.mean)[6];
  std::cout << "Whole-point sample mean: " << uni20::format_real(whole) << "; joint covariance dimension "
            << result.covariance->size() << '\n';
  std::cout << "Small restart demonstration, not a converged DSF estimate.\n";
}
catch (std::exception const& error)
{
  std::cerr << error.what() << '\n';
  return 1;
}
