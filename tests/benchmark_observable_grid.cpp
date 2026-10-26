// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
// Opt-in matched-grid benchmark. Numerical results are checked, not timed output.
#include <bethe/hubbard_doped.hpp>
#include <bethe/lieb_liniger_thermo.hpp>
#include <bethe/observable_grid.hpp>
#include <charconv>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <uni20/async/tbb_scheduler.hpp>
#include <uni20/core/scalar_io.hpp>

using Clock = std::chrono::steady_clock;
struct Timing
{
    double wall = 0, cpu = 0;
};
struct Timer
{
    Clock::time_point start = Clock::now();
    std::clock_t cpu = std::clock();
    Timing elapsed() const noexcept
    {
      auto const end = std::clock();
      auto const seconds = cpu == std::clock_t(-1) || end == std::clock_t(-1) || end < cpu
                               ? std::numeric_limits<double>::quiet_NaN()
                               : double(end - cpu) / CLOCKS_PER_SEC;
      return {std::chrono::duration<double>(Clock::now() - start).count(), seconds};
    }
};

template <uni20::Real Real, typename Evaluate>
void measure(std::size_t count, std::size_t repeats, Timing preparation, Evaluate evaluate)
{
  std::cout << "prepare_wall=" << preparation.wall << " prepare_cpu=" << preparation.cpu << '\n';
  std::optional<Real> reference_sum;
  std::optional<std::size_t> reference_failures;
  for (std::size_t repeat = 0; repeat < repeats; ++repeat)
    for (int step = 0; step < 4; ++step)
    {
      // Alternate the sweep order, rather than always timing eight workers last.
      int const workers = 1 << (repeat % 2 ? 3 - step : step);
      uni20::async::TbbScheduler scheduler(workers);
      bethe::detail::CompensatedSum<Real> sum;
      std::size_t next = 0, failed = 0;
      Timer timer;
      bethe::evaluate_observable_grid(count, evaluate,
                                      [&](std::size_t i, auto const& point) {
                                        if (i != next++) throw std::runtime_error("reordered grid");
                                        if (!point.converged)
                                        {
                                          ++failed;
                                          return;
                                        }
                                        if constexpr (requires { *point.energy; })
                                          sum.add(*point.energy);
                                        else
                                          sum.add(point.energy);
                                      },
                                      {&scheduler});
      auto const whole = timer.elapsed();
      if (next != count || (reference_sum && *reference_sum != sum.value()) ||
          (reference_failures && *reference_failures != failed))
        throw std::runtime_error("serial/parallel grid checksum differs");
      reference_sum = sum.value();
      reference_failures = failed;
      std::cout << "repeat=" << repeat + 1 << " workers=" << workers << " grid_wall=" << whole.wall
                << " grid_cpu=" << whole.cpu << " failed=" << failed
                << " converged_energy_sum=" << uni20::format_scalar(sum.value()) << std::endl;
    }
}

template <uni20::Real Real> void benchmark(std::string_view model, std::size_t count, std::size_t repeats)
{
  Real const pi = Real{4} * std::atan(Real{1});
  Timer prepare;
  if (model == "lieb-liniger")
  {
    namespace ll = bethe::lieb_liniger::thermo;
    ll::Solver<Real> const solver(Real{4}, Real{1});
    if (!solver.background().converged) throw std::runtime_error("background did not converge");
    measure<Real>(count, repeats, prepare.elapsed(), [&](std::size_t i) {
      return solver.at_momentum(ll::Branch::type_i, Real{8} * Real(i) / Real(count - 1));
    });
  }
  else if (model == "hubbard-doped")
  {
    namespace hub = bethe::hubbard::thermo;
    hub::DopedSolver<Real> const solver(Real{4}, Real{0.5});
    if (!solver.background().converged) throw std::runtime_error("background did not converge");
    auto const [lo, hi] = solver.momentum_range(hub::DopedBranch::holon);
    measure<Real>(count, repeats, prepare.elapsed(), [&](std::size_t i) {
      return solver.at_momentum(hub::DopedBranch::holon,
                                i == count - 1 ? hi : lo + (hi - lo) * Real(i) / Real(count - 1));
    });
  }
  else if (model == "hubbard-half")
    measure<Real>(count, repeats, prepare.elapsed(), [&](std::size_t i) {
      namespace hub = bethe::hubbard::thermo;
      return hub::dispersion(hub::Branch::holon, Real{4}, pi * (Real{2} * Real(i) / Real(count - 1) - Real{1}));
    });
  else
    throw std::invalid_argument("unknown model");
}

int main(int argc, char** argv)
try
{
  if (argc < 4 || argc > 5)
    throw std::invalid_argument("usage: bethe_observable_grid_benchmark "
                                "lieb-liniger|hubbard-doped|hubbard-half fp64|long-double|fp128 points [repeats=3]");
  auto parse = [](std::string_view value) {
    std::size_t result{};
    auto const [end, error] = std::from_chars(value.data(), value.data() + value.size(), result);
    if (error != std::errc{} || end != value.data() + value.size() || result == 0)
      throw std::invalid_argument("expected a positive integer");
    return result;
  };
  auto const count = parse(argv[3]), repeats = argc == 5 ? parse(argv[4]) : 3;
  if (count < 2) throw std::invalid_argument("grid needs at least two points");
  std::string_view const model(argv[1]), precision(argv[2]);
  std::cout << std::setprecision(9) << "model=" << model << " precision=" << precision << " points=" << count
            << " repeats=" << repeats << " max_in_flight=128\n";
  if (precision == "fp64")
    benchmark<double>(model, count, repeats);
  else if (precision == "long-double")
    benchmark<long double>(model, count, repeats);
#if UNI20_HAS_FLOAT128
  else if (precision == "fp128")
    benchmark<uni20::float128>(model, count, repeats);
#endif
  else
    throw std::invalid_argument("unsupported precision");
  return 0;
}
catch (std::exception const& error)
{
  std::cerr << error.what() << '\n';
  return 1;
}
