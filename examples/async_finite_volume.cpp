// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include <bethe/lee_yang_excited.hpp>
#include <bethe/observable_grid.hpp>
#include <bethe/sine_gordon_excited.hpp>
#include <charconv>
#include <chrono>
#include <ctime>
#include <iostream>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
int count(std::string_view text)
{
  int n = -1;
  auto const [end, error] = std::from_chars(text.data(), text.data() + text.size(), n);
  if (error != std::errc{} || end != text.data() + text.size() || n < 0 || n > 10000)
    throw std::invalid_argument("counts must be integers in [0,10000]");
  return n;
}

template <typename Solve> int sweep(Solve solve, int threads, std::size_t points)
{
  using State = std::invoke_result_t<Solve, std::size_t>;
  std::vector<State> states;
  states.reserve(points);
  std::unique_ptr<uni20::async::TbbScheduler> scheduler;
  if (threads > 0) scheduler = std::make_unique<uni20::async::TbbScheduler>(threads);
  auto const start = std::chrono::steady_clock::now();
  auto const cpu_start = std::clock();
  if (scheduler)
  {
    auto done = bethe::evaluate_observable_grid_async(
        points, std::move(solve), [&states](std::size_t, State state) { states.push_back(std::move(state)); },
        {.scheduler = scheduler.get(), .max_in_flight = 32});
    (void)done.get_wait(*scheduler); // every admitted solve and ordered collection
  }
  else
    for (std::size_t i = 0; i != points; ++i)
      states.push_back(solve(i));
  double const wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  double const cpu = double(std::clock() - cpu_start) / CLOCKS_PER_SEC;
  std::cout << "# threads=" << threads << " points=" << points << " wall_seconds=" << wall << " cpu_seconds=" << cpu
            << '\n';
  std::cout << "# Solver and ordered collection timing; excludes scheduler construction and printing.\n"
               "length,casimir_energy,scaling_function,intervals,iterations,converged,status\n";
  bool complete = true;
  for (auto const& s : states)
  {
    std::cout << uni20::format_real(s.length) << ',';
    if (s.casimir_energy) std::cout << uni20::format_real(*s.casimir_energy);
    std::cout << ',';
    if (s.scaling_function) std::cout << uni20::format_real(*s.scaling_function);
    std::cout << ',' << s.intervals << ',' << s.iterations << ',' << s.converged << ',' << int(s.status) << '\n';
    complete = complete && s.converged;
  }
  return complete ? 0 : 2;
}

template <uni20::Real R> int run(std::string_view family, int threads, std::size_t points)
{
  // Owning callables retain all scales/options. Each state has independent
  // meshes, scratch and convergence budgets; there is no coupled continuation.
  if (family == "lee-yang-vacuum" || family == "lee-yang-excited")
  {
    bethe::lee_yang::OneParticleOptions<R> options;
    options.tolerance = R{1e-10};
    auto length = [points](std::size_t i) { return R{5} + R{5} * R(i) / R(points); };
    if (family == "lee-yang-excited")
      return sweep([=](std::size_t i) { return bethe::lee_yang::one_particle(R{1}, length(i), options); }, threads,
                   points);
    return sweep([=](std::size_t i) { return bethe::lee_yang::ground_state(R{1}, length(i), options); }, threads,
                 points);
  }
  bethe::sine_gordon::TwoSolitonOptions<R> options;
  options.tolerance = R{1e-7};
  auto length = [points](std::size_t i) { return R{1} + R{3} * R(i) / R(points); };
  if (family == "sine-gordon-excited")
    return sweep(
        [=](std::size_t i) {
          return bethe::sine_gordon::two_soliton_level(R{1}, length(i), R{2}, uni20::from_twice(1), options);
        },
        threads, points);
  if (family == "sine-gordon-vacuum")
    return sweep([=](std::size_t i) { return bethe::sine_gordon::vacuum_energy(R{1}, length(i), R{2}, options); },
                 threads, points);
  throw std::invalid_argument("unknown model/state");
}
} // namespace

int main(int argc, char** argv)
{
  try
  {
    if (argc != 5)
      throw std::invalid_argument(
          "usage: bethe-async-finite-volume-example MODEL THREADS POINTS PRECISION\n"
          "MODEL: lee-yang-vacuum, lee-yang-excited, sine-gordon-vacuum, sine-gordon-excited\n"
          "THREADS: 0=direct serial, 1..10000=Async; POINTS: 1..10000\n"
          "PRECISION: fp64, long-double, fp128 (when configured); MPFR is not supported by these solvers");
    int const threads = count(argv[2]), points = count(argv[3]);
    if (!points) throw std::invalid_argument("points must be positive");
    std::string_view const precision = argv[4];
    if (precision == "fp64") return run<double>(argv[1], threads, points);
    if (precision == "long-double") return run<long double>(argv[1], threads, points);
#if UNI20_HAS_FLOAT128
    if (precision == "fp128") return run<uni20::float128>(argv[1], threads, points);
#endif
    throw std::invalid_argument("unsupported precision");
  }
  catch (std::exception const& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
