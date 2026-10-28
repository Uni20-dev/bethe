// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include <algorithm>
#include <bethe/observable_grid.hpp>
#include <bethe/su_fermions.hpp>
#include <bethe/takhtajan_babujian.hpp>
#include <bethe/tj.hpp>
#include <charconv>
#include <chrono>
#include <ctime>
#include <functional>
#include <iostream>
#include <iterator>
#include <string_view>
#include <tuple>
#include <type_traits>
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
template <typename R> void print(std::vector<bethe::takhtajan_babujian::State<R>> const& states)
{
  std::cout << "sites,iterations,string_quantum_numbers,centers,deviations,rapidities,energy,momentum,residual,"
               "phase_residual,modulus_residual,momentum_index,converged,status\n";
  for (auto const& s : states)
    row(s.sites, s.iterations, s.string_quantum_numbers, s.centers, s.deviations, s.rapidities, s.energy, s.momentum,
        s.residual_norm, s.phase_residual, s.modulus_residual, s.momentum_index, s.converged, s.status);
}
template <typename R> void print(std::vector<bethe::su_fermions::State<R>> const& states)
{
  std::cout << "populations,component_order,particles,iterations,stages,length,interaction,residual,target_residual,"
               "momentum,momentum_index,quantum_numbers,rapidities,level_residuals,free_modes,energy,"
               "reached_interaction,free,converged,status\n";
  for (auto const& s : states)
    row(s.populations, s.component_order, s.particles, s.iterations, s.stages, s.length, s.interaction, s.residual_norm,
        s.target_residual_norm, s.momentum, s.momentum_index, s.quantum_numbers, s.rapidities, s.level_residuals,
        s.free_modes, s.energy, s.reached_interaction, s.free, s.converged, s.status);
}
template <typename R> void print(std::vector<bethe::tj::State<R>> const& states)
{
  std::cout << "sites,num_up,num_down,holes,branch,quantum_numbers,rapidities,free_modes,energy,momentum,residual,"
               "level_residuals,momentum_index,iterations,converged,spin_reversed,status\n";
  for (auto const& s : states)
    row(s.sites, s.num_up, s.num_down, s.holes, s.branch, s.quantum_numbers, s.rapidities, s.free_modes, s.energy,
        s.momentum, s.residual_norm, s.level_residuals, s.momentum_index, s.iterations, s.converged, s.spin_reversed,
        s.status);
}
template <typename Solve> int sweep(std::size_t points, int threads, Solve solve)
{
  using State = std::invoke_result_t<Solve, std::size_t>;
  std::vector<State> states;
  states.reserve(points);
  auto collect = [&states](std::size_t, State state) { states.push_back(std::move(state)); };
  auto const start = std::chrono::steady_clock::now();
  auto const cpu_start = std::clock();
  int const concurrency = std::min(threads, static_cast<int>(std::min(points, std::size_t{32})));
  if (concurrency == 0)
    for (std::size_t i = 0; i != points; ++i)
      collect(i, solve(i));
  else
  {
    uni20::async::TbbScheduler scheduler(concurrency);
    // Captures own inputs/options, including the SU(N) population vector.
    // Every coupled solve owns its scratch. Wait on this grid, not a global drain.
    auto done = bethe::evaluate_observable_grid_async(points, std::move(solve), collect,
                                                      {.scheduler = &scheduler, .max_in_flight = 32});
    (void)done.get_wait(scheduler);
  }
  std::cout << "# threads=" << threads << " concurrency=" << concurrency << " points=" << points
            << " wall_seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count()
            << " cpu_seconds=" << double(std::clock() - cpu_start) / CLOCKS_PER_SEC << '\n';
  std::cout << "# Timing includes scheduler and ordered collection; excludes printing.\n";
  print(states);
  return std::all_of(states.begin(), states.end(), [](auto const& s) { return s.converged; }) ? 0 : 2;
}
template <typename R> int run(std::string_view model, std::size_t size, std::size_t points, int threads)
{
  if (model == "tb")
  {
    if (size < 4 || size % 2 || size + 2 * (points - 1) > 10000)
      throw std::invalid_argument("TB requires even initial size>=4 and final size<=10000");
    bethe::takhtajan_babujian::SolverOptions<R> options;
    return sweep(points, threads,
                 [=](std::size_t i) { return bethe::takhtajan_babujian::ground_state<R>(size + 2 * i, options); });
  }
  if (model == "sun")
  {
    if (size == 0 || size % 2 == 0) throw std::invalid_argument("SU(N) requires odd positive population per component");
    std::vector<std::size_t> populations(3, size);
    bethe::su_fermions::SolverOptions<R> options;
    return sweep(points, threads, [=](std::size_t i) {
      R const c = points == 1 ? R{1} : (R{1} + R{7} * R(i) / R(points - 1)) / R{2};
      return bethe::su_fermions::ground_state<R>(populations, R(3 * size), c, options);
    });
  }
  if (size < 3 || points > (size / 2 + 1) / 2)
    throw std::invalid_argument("t-J requires sites>=3 and 2*(2*points-1)<=sites");
  bethe::tj::SolverOptions<R> options;
  return sweep(points, threads,
               [=](std::size_t i) { return bethe::tj::ground_state<R>(size, 2 * i + 1, 2 * i + 1, options); });
}
} // namespace
int main(int argc, char** argv)
{
  try
  {
    if (argc != 6)
      throw std::invalid_argument("usage: bethe-async-nested-sweeps-example MODEL SIZE POINTS THREADS PRECISION\n"
                                  "MODEL: tb, sun, tj; THREADS: 0=direct, 1+=Async; "
                                  "PRECISION: fp64, long-double, configured fp128");
    std::string_view const model = argv[1], precision = argv[5];
    int const size = count(argv[2]), points = count(argv[3]), threads = count(argv[4]);
    if ((model != "tb" && model != "sun" && model != "tj") || points == 0)
      throw std::invalid_argument("invalid model or empty sweep");
    if (precision == "fp64") return run<double>(model, size, points, threads);
    if (precision == "long-double") return run<long double>(model, size, points, threads);
#if UNI20_HAS_FLOAT128
    if (precision == "fp128") return run<uni20::float128>(model, size, points, threads);
#endif
    throw std::invalid_argument("unsupported precision");
  }
  catch (std::exception const& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
