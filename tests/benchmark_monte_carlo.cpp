// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
// Synthetic CPU workload: scheduling measurements, not a physical speed claim.
#include <bethe/monte_carlo.hpp>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <uni20/async/tbb_scheduler.hpp>

int main(int argc, char** argv)
{
  if (argc != 5)
  {
    std::cerr << "usage: benchmark_monte_carlo WORKERS BLOCK_SIZE REPLICATES ITERATIONS\n";
    return 2;
  }
  auto workers = std::stoul(argv[1]), grain = std::stoul(argv[2]);
  auto replicates = std::stoul(argv[3]), iterations = std::stoul(argv[4]);
  if (!workers || !grain || !replicates) return 2;
  uni20::async::TbbScheduler scheduler(workers);
  std::atomic<unsigned> active{0}, peak{0};
  auto start = std::chrono::steady_clock::now();
  auto cpu = std::clock();
  auto pending =
      bethe::monte_carlo_async<2>({.replicates = replicates, .samples_per_replicate = 768, .seed = 109}, 1, 1.,
                                  [&](auto const& identity, std::stop_token) {
                                    auto n = ++active, old = peak.load();
                                    while (n > old && !peak.compare_exchange_weak(old, n))
                                    {}
                                    auto u = identity.at_precision(1.);
                                    double x = u[0];
                                    for (std::size_t j = 0; j < iterations; ++j)
                                      x = std::sin(x) + std::sqrt(1. + x * x);
                                    --active;
                                    return bethe::MonteCarloObservation<double>{std::vector{x + u[1]}};
                                  },
                                  {.scheduler = &scheduler, .batch_size = grain, .max_in_flight = 2 * workers});
  auto result = pending.get_wait(scheduler);
  double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  double cpu_seconds = double(std::clock() - cpu) / CLOCKS_PER_SEC;
  if (result.status != bethe::MonteCarloStatus::complete) return 1;
  std::cout << std::setprecision(17) << "{\"workers\":" << workers << ",\"block\":" << grain
            << ",\"replicates\":" << replicates << ",\"samples\":768,\"iterations\":" << iterations
            << ",\"wall_seconds\":" << wall << ",\"cpu_seconds\":" << cpu_seconds << ",\"peak_observations\":" << peak
            << ",\"mean\":" << (*result.mean)[0]
            << ",\"variance\":" << (result.covariance ? (*result.covariance)[0] : 0) << "}\n";
}
