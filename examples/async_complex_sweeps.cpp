// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include <algorithm>
#include <bethe/observable_grid.hpp>
#include <bethe/tasep.hpp>
#include <bethe/xyz.hpp>
#include <charconv>
#include <chrono>
#include <ctime>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string_view>
#include <uni20/async/tbb_scheduler.hpp>
#include <vector>

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
template <typename R> void print(std::vector<bethe::tasep::GapState<R>> const& rows)
{
  std::cout << "sites,particles,effective_particles,rate,residual,iterations,seed_iterations,log_y_real,log_y_imag,"
               "lambda_real,lambda_imag,gap,converged,status,roots\n";
  for (auto const& s : rows)
  {
    std::cout << s.sites << ',' << s.particles << ',' << s.effective_particles << ',' << uni20::format_real(s.rate)
              << ',' << uni20::format_real(s.residual_norm) << ',' << s.iterations << ',' << s.seed_iterations << ','
              << uni20::format_real(s.log_y.real()) << ',' << uni20::format_real(s.log_y.imag()) << ',';
    if (s.eigenvalue) std::cout << uni20::format_real(s.eigenvalue->real());
    std::cout << ',';
    if (s.eigenvalue) std::cout << uni20::format_real(s.eigenvalue->imag());
    std::cout << ',';
    if (s.gap) std::cout << uni20::format_real(*s.gap);
    std::cout << ',' << s.converged << ',' << int(s.status) << ',';
    for (std::size_t j = 0; j != s.roots.size(); ++j)
    {
      if (j) std::cout << ';';
      std::cout << uni20::format_real(s.roots[j].real()) << ' ' << uni20::format_real(s.roots[j].imag());
    }
    std::cout << '\n';
  }
}
template <typename R> void print(std::vector<bethe::xyz::GroundState<R>> const& rows)
{
  std::cout << "sites,eta,t,Jx,Jy,Jz,energy,momentum,residual,iterations,converged,status,quantum_numbers,roots\n";
  for (auto const& s : rows)
  {
    std::cout << s.sites << ',' << uni20::format_real(s.eta) << ',' << uni20::format_real(s.t) << ','
              << uni20::format_real(s.exchange.x) << ',' << uni20::format_real(s.exchange.y) << ','
              << uni20::format_real(s.exchange.z) << ',';
    if (s.energy) std::cout << uni20::format_real(*s.energy);
    std::cout << ',' << uni20::format_real(s.momentum) << ',' << uni20::format_real(s.residual_norm) << ','
              << s.iterations << ',' << s.converged << ',' << int(s.status) << ',';
    for (std::size_t j = 0; j != s.quantum_numbers.size(); ++j)
    {
      if (j) std::cout << ';';
      std::cout << uni20::format_real(R(s.quantum_numbers[j].twice()) / R{2});
    }
    std::cout << ',';
    for (std::size_t j = 0; j != s.roots.size(); ++j)
    {
      if (j) std::cout << ';';
      std::cout << uni20::format_real(s.roots[j].real()) << ' ' << uni20::format_real(s.roots[j].imag());
    }
    std::cout << '\n';
  }
}
template <typename Solve> int sweep(std::size_t points, int threads, Solve solve)
{
  using State = std::invoke_result_t<Solve, std::size_t>;
  std::vector<State> rows;
  rows.reserve(points);
  auto collect = [&rows](std::size_t, State state) { rows.push_back(std::move(state)); };
  auto const start = std::chrono::steady_clock::now();
  auto const cpu_start = std::clock();
  int const concurrency = std::min(threads, static_cast<int>(std::min(points, std::size_t{32})));
  if (concurrency == 0)
    for (std::size_t i = 0; i != points; ++i)
      collect(i, solve(i));
  else
  {
    uni20::async::TbbScheduler scheduler(concurrency);
    // The callable owns its inputs/options. Each scalar solve owns all mutable
    // roots and scratch; collection is ordered and runs on the scheduler.
    auto done = bethe::evaluate_observable_grid_async(points, std::move(solve), collect,
                                                      {.scheduler = &scheduler, .max_in_flight = 32});
    (void)done.get_wait(scheduler);
  }
  std::cout << "# threads=" << threads << " concurrency=" << concurrency << " points=" << points
            << " wall_seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count()
            << " cpu_seconds=" << double(std::clock() - cpu_start) / CLOCKS_PER_SEC << '\n';
  std::cout << "# Timing includes scheduler and ordered collection; excludes printing.\n";
  print(rows);
  return std::all_of(rows.begin(), rows.end(), [](auto const& s) { return s.converged; }) ? 0 : 2;
}
template <typename R> int run(bool tasep, std::size_t sites, std::size_t points, int threads)
{
  if (tasep)
  {
    if (sites < 2 || sites > 256 || points > sites / 2)
      throw std::invalid_argument("TASEP requires 2<=sites<=256 and points<=sites/2 (fillings 1..points)");
    bethe::tasep::Options<R> options;
    return sweep(points, threads,
                 [=](std::size_t i) { return bethe::tasep::relaxation_gap(sites, i + 1, R{1}, options); });
  }
  if (sites < 2 || sites % 2) throw std::invalid_argument("XYZ requires even sites>=2");
  bethe::xyz::SolverOptions<R> options;
  return sweep(points, threads, [=](std::size_t i) {
    R const eta = points == 1 ? R{1} / R{2} : (R{1} + R{3} * R(i) / R(points - 1)) / R{5};
    return bethe::xyz::ground_state(sites, eta, R{3} / R{4}, options);
  });
}
} // namespace
int main(int argc, char** argv)
{
  try
  {
    if (argc != 6)
      throw std::invalid_argument("usage: bethe-async-complex-sweeps-example MODEL SITES POINTS THREADS PRECISION\n"
                                  "MODEL: tasep or xyz; THREADS: 0=direct, 1+=Async; "
                                  "PRECISION: fp64, long-double, configured fp128");
    std::string_view const model = argv[1], precision = argv[5];
    int const sites = count(argv[2]), points = count(argv[3]), threads = count(argv[4]);
    if ((model != "tasep" && model != "xyz") || points == 0)
      throw std::invalid_argument("invalid model or empty sweep");
    if (precision == "fp64") return run<double>(model == "tasep", sites, points, threads);
    if (precision == "long-double") return run<long double>(model == "tasep", sites, points, threads);
#if UNI20_HAS_FLOAT128
    if (precision == "fp128") return run<uni20::float128>(model == "tasep", sites, points, threads);
#endif
    throw std::invalid_argument("unsupported precision");
  }
  catch (std::exception const& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
