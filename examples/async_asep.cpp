// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include <bethe/asep.hpp>
#include <bethe/prepared_observable_grid.hpp>
#include <charconv>
#include <chrono>
#include <ctime>
#include <iostream>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
int count(std::string_view text)
{
  int value = -1;
  auto const [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
  if (error != std::errc{} || end != text.data() + text.size() || value < 0 || value > 10000)
    throw std::invalid_argument("counts must be integers in [0,10000]");
  return value;
}

template <uni20::Real R> int run(bool reuse, std::size_t sites, std::size_t particles, std::size_t points, int threads)
{
  using State = bethe::asep::GapState<R>;
  using Prepared = bethe::asep::PreparedGap<R>;
  bethe::asep::Options<R> options;
  auto left_rate = [points](std::size_t i) {
    return points == 1 ? R{1} / R{2} : R{9} * R(i) / (R{10} * R(points - 1));
  };
  std::vector<State> rows;
  rows.reserve(points);
  auto collect = [&rows](std::size_t, State state) { rows.push_back(std::move(state)); };
  auto solve = [=](std::size_t i) {
    return bethe::asep::relaxation_gap(sites, particles, R{1}, left_rate(i), options);
  };
  auto const start = std::chrono::steady_clock::now();
  auto const cpu_start = std::clock();
  int const concurrency = std::min(threads, static_cast<int>(std::min(points, std::size_t{32})));
  if (concurrency == 0)
  {
    if (reuse)
    {
      Prepared const background(sites, particles, options);
      for (std::size_t i = 0; i != points; ++i)
        collect(i, background.at_rates(R{1}, left_rate(i)));
    }
    else
      for (std::size_t i = 0; i != points; ++i)
        collect(i, solve(i));
  }
  else
  {
    uni20::async::TbbScheduler scheduler(concurrency);
    bethe::ExecutionOptions execution{.scheduler = &scheduler, .max_in_flight = 32};
    if (reuse)
    {
      auto background = bethe::prepare_background_async([=] { return Prepared(sites, particles, options); }, execution);
      std::exception_ptr error;
      try
      {
        // Submit the consumer before preparation is ready. The reader owns
        // the shared seed through all admitted queries and collection.
        auto done = bethe::evaluate_observable_grid_async(
            background.read(), points,
            [left_rate](Prepared const& prepared, std::size_t i) { return prepared.at_rates(R{1}, left_rate(i)); },
            collect, execution);
        (void)done.get_wait(scheduler);
      }
      catch (...)
      {
        error = std::current_exception();
      }
      // Preparation can be admitted even if submitting the consumer throws.
      try
      {
        (void)background.get_wait(scheduler);
      }
      catch (...)
      {
        if (!error) error = std::current_exception();
      }
      if (error) std::rethrow_exception(error);
    }
    else
    {
      auto done = bethe::evaluate_observable_grid_async(points, solve, collect, execution);
      (void)done.get_wait(scheduler);
    }
  }
  double const wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  double const cpu = double(std::clock() - cpu_start) / CLOCKS_PER_SEC;
  std::cout << "# shared_seed=" << reuse << " threads=" << threads << " concurrency=" << concurrency
            << " points=" << points << " wall_seconds=" << wall << " cpu_seconds=" << cpu << '\n';
  std::cout << "# Timing includes preparation, scheduler and ordered collection; excludes printing.\n"
               "# Seed work counters describe preparation; do not sum them for shared-seed CPU work.\n"
               "sites,particles,effective_particles,right_rate,left_rate,reached_ratio,lambda_real,lambda_imag,gap,"
               "residual,wave_index,wave_base_real,wave_base_imag,iterations,continuation_steps,seed_iterations,"
               "seed_newton_iterations,converged,analytic,status,scaled_roots\n";
  bool complete = true;
  for (auto const& s : rows)
  {
    std::cout << s.sites << ',' << s.particles << ',' << s.effective_particles << ','
              << uni20::format_real(s.right_rate) << ',' << uni20::format_real(s.left_rate) << ','
              << uni20::format_real(s.reached_ratio) << ',';
    if (s.eigenvalue) std::cout << uni20::format_real(s.eigenvalue->real());
    std::cout << ',';
    if (s.eigenvalue) std::cout << uni20::format_real(s.eigenvalue->imag());
    std::cout << ',';
    if (s.gap) std::cout << uni20::format_real(*s.gap);
    std::cout << ',' << uni20::format_real(s.residual_norm) << ',' << s.wave_index << ','
              << uni20::format_real(s.wave_base.real()) << ',' << uni20::format_real(s.wave_base.imag()) << ','
              << s.iterations << ',' << s.continuation_steps << ',' << s.seed_iterations << ','
              << s.seed_newton_iterations << ',' << s.converged << ',' << s.analytic << ',' << int(s.status) << ',';
    for (std::size_t j = 0; j != s.scaled_roots.size(); ++j)
    {
      if (j) std::cout << ';';
      std::cout << uni20::format_real(s.scaled_roots[j].real()) << ' ' << uni20::format_real(s.scaled_roots[j].imag());
    }
    std::cout << '\n';
    complete = complete && s.converged;
  }
  return complete ? 0 : 2;
}
} // namespace

int main(int argc, char** argv)
{
  try
  {
    if (argc != 7)
      throw std::invalid_argument("usage: bethe-async-asep-example MODE SITES PARTICLES POINTS THREADS PRECISION\n"
                                  "MODE: independent or prepared; right rate 1, left rates [0,0.9] (one point: 0.5)\n"
                                  "THREADS: 0=direct serial, positive=Async capped by points and window 32\n"
                                  "PRECISION: fp64, long-double, fp128 when configured; MPFR is not supported");
    std::string_view const mode = argv[1], precision = argv[6];
    if (mode != "independent" && mode != "prepared") throw std::invalid_argument("unknown seed mode");
    int const sites = count(argv[2]), particles = count(argv[3]), points = count(argv[4]), threads = count(argv[5]);
    if (points == 0) throw std::invalid_argument("points must be positive");
    if (precision == "fp64") return run<double>(mode == "prepared", sites, particles, points, threads);
    if (precision == "long-double") return run<long double>(mode == "prepared", sites, particles, points, threads);
#if UNI20_HAS_FLOAT128
    if (precision == "fp128") return run<uni20::float128>(mode == "prepared", sites, particles, points, threads);
#endif
    throw std::invalid_argument("unsupported precision");
  }
  catch (std::exception const& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
