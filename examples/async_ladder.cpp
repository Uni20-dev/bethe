// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include <bethe/ladder.hpp>
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
// Compact diagnostic cells for this example, including all nested vectors.
// Semicolons delimit sequence items; half-integer labels print as decimals.
template <typename T> void cell(T const& value)
{
  if constexpr (std::is_integral_v<T>)
    std::cout << value;
  else if constexpr (std::is_enum_v<T>)
    std::cout << static_cast<int>(value);
  else if constexpr (requires { value.has_value(); })
  {
    if (value) cell(*value);
  }
  else if constexpr (requires { value.twice(); })
    std::cout << uni20::to_string(value);
  else if constexpr (requires {
                       value.real();
                       value.imag();
                     })
  {
    cell(value.real());
    std::cout << ':';
    cell(value.imag());
  }
  else if constexpr (requires {
                       std::begin(value);
                       std::end(value);
                     })
  {
    std::cout << '[';
    bool first = true;
    for (auto const& item : value)
    {
      if (!first) std::cout << ';';
      first = false;
      cell(item);
    }
    std::cout << ']';
  }
  else
    std::cout << uni20::format_real(value);
}
template <typename... T> void row(T const&... values)
{
  bool first = true;
  auto next = [&](auto const& value) {
    if (!first) std::cout << ',';
    first = false;
    std::cout << '"';
    cell(value);
    std::cout << '"';
  };
  (next(values), ...);
  std::cout << '\n';
}

template <typename R> void print(std::vector<bethe::ladder::State<R>> const& states)
{
  std::cout << "rungs,singlets,iterations,branches,tableaux,populations,rung_coupling,field,magnetization,energy,"
               "converged,descendant,analytic,status,shape,labels,rapidities,permutation_energy,residual,"
               "momentum_index,branch_iterations,branch_converged,branch_status\n";
  for (auto const& s : states)
  {
    auto const& b = s.highest_weight;
    row(s.rungs, s.singlets, s.iterations, s.branches, s.tableaux, s.populations, s.rung_coupling, s.magnetic_field,
        s.magnetization, s.energy, s.converged, s.descendant, s.analytic, s.status, b.shape, b.labels, b.rapidities,
        b.energy, b.residual, b.momentum_index, b.iterations, b.converged, b.status);
  }
}
template <typename R> int run(bool prepared, std::size_t rungs, std::size_t points, int threads)
{
  namespace model = bethe::ladder;
  auto parameters = [points](std::size_t i) {
    R const field = points == 1 ? R{0} : R{4} * R(i) / R(points - 1) - R{2};
    return std::pair{R{2} - field / R{2}, field};
  };
  std::vector<model::State<R>> states;
  states.reserve(points);
  auto collect = [&states](std::size_t, auto state) { states.push_back(std::move(state)); };
  auto direct = [=](std::size_t i) {
    auto const [rung, field] = parameters(i);
    return model::ground_state<R>(rungs, rung, field);
  };
  auto query = [=](model::PreparedSectors<R> const& background, std::size_t i) {
    auto const [rung, field] = parameters(i);
    return background.ground_state(rung, field);
  };
  auto const start = std::chrono::steady_clock::now();
  auto const cpu_start = std::clock();
  double preparation_seconds = 0;
  auto prepare = [rungs, &preparation_seconds] {
    auto const before = std::chrono::steady_clock::now();
    model::PreparedSectors<R> background(rungs);
    preparation_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - before).count();
    return background;
  };
  int const concurrency = std::min(threads, static_cast<int>(std::min(points, std::size_t{32})));
  if (concurrency == 0)
  {
    if (prepared)
    {
      auto background = prepare();
      for (std::size_t i = 0; i != points; ++i)
        collect(i, query(background, i));
    }
    else
      for (std::size_t i = 0; i != points; ++i)
        collect(i, direct(i));
  }
  else
  {
    uni20::async::TbbScheduler scheduler(concurrency);
    bethe::ExecutionOptions execution{.scheduler = &scheduler, .max_in_flight = 32};
    if (prepared)
    {
      auto background = bethe::prepare_background_async(prepare, execution);
      try
      {
        // Submit consumers without synchronously waiting for preparation.
        auto done = bethe::evaluate_observable_grid_async(background.read(), points, query, collect, execution);
        (void)done.get_wait(scheduler);
      }
      catch (...)
      {
        auto error = std::current_exception();
        // Preparation can be accepted before consumer admission fails.
        try
        {
          (void)background.get_wait(scheduler);
        }
        catch (...)
        {}
        std::rethrow_exception(error);
      }
    }
    else
    {
      auto done = bethe::evaluate_observable_grid_async(points, direct, collect, execution);
      (void)done.get_wait(scheduler);
    }
  }
  std::cout << "# mode=" << (prepared ? "prepared" : "direct") << " threads=" << threads
            << " concurrency=" << concurrency << " points=" << points
            << " wall_seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count()
            << " cpu_seconds=" << double(std::clock() - cpu_start) / CLOCKS_PER_SEC
            << " preparation_seconds=" << preparation_seconds << '\n';
  std::cout << "# Timing includes scheduler, preparation and ordered collection; excludes printing.\n";
  print(states);
  return std::all_of(states.begin(), states.end(), [](auto const& s) { return s.converged; }) ? 0 : 2;
}
} // namespace
int main(int argc, char** argv)
{
  try
  {
    if (argc != 6)
      throw std::invalid_argument(
          "usage: bethe-async-ladder-example MODE RUNGS POINTS THREADS PRECISION\n"
          "MODE: direct, prepared; THREADS: 0=serial, 1+=Async; PRECISION: fp64, long-double, configured fp128");
    std::string_view const mode = argv[1], precision = argv[5];
    int const rungs = count(argv[2]), points = count(argv[3]), threads = count(argv[4]);
    if ((mode != "direct" && mode != "prepared") || points == 0 || rungs < 2)
      throw std::invalid_argument("invalid mode, empty sweep or fewer than two rungs");
    if (precision == "fp64") return run<double>(mode == "prepared", rungs, points, threads);
    if (precision == "long-double") return run<long double>(mode == "prepared", rungs, points, threads);
#if UNI20_HAS_FLOAT128
    if (precision == "fp128") return run<uni20::float128>(mode == "prepared", rungs, points, threads);
#endif
    throw std::invalid_argument("unsupported precision");
  }
  catch (std::exception const& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
