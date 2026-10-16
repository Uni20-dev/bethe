// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
// Paired, single-thread form-factor timing on identical pre-solved states.
#include "xxx_unprepared_form_factor.hpp"
#include <bethe/heisenberg_excitations.hpp>
#include <uni20/core/scalar_io.hpp>

#include <charconv>
#include <chrono>
#include <iostream>
#include <string_view>

namespace model = bethe::heisenberg;
using Clock = std::chrono::steady_clock;

double seconds(Clock::time_point start) { return std::chrono::duration<double>(Clock::now() - start).count(); }

template <uni20::Real Real> int benchmark(std::size_t n, std::size_t repeats, std::string_view precision)
{
  using std::abs;
  auto const count = model::real_excitation_count(n, uni20::half_int{1});
  auto const scan_start = Clock::now();
  auto const scan = model::real_excitations<Real>(n, uni20::half_int{1}, {count});
  double const scan_seconds = seconds(scan_start);
  if (!scan.converged()) throw std::runtime_error("incomplete root scan");
  auto const prepare_start = Clock::now();
  model::PreparedRaisingFormFactors<Real> const prepared(n, scan.ground_state);
  double const prepare_seconds = seconds(prepare_start);
  if (!prepared.ready()) throw std::runtime_error("ground preparation failed");
  Real max_relative_change{0};
  // Untimed validation is also a warm-up; both paths use precisely the same roots.
  for (auto const& level : scan.levels)
  {
    auto const old = test_support::unprepared_raising_form_factor(n, scan.ground_state, level.state);
    auto const current = prepared(level.state);
    if (!old.converged() || !current.converged()) throw std::runtime_error("unresolved form factor");
    Real const error = abs(*current.weight / *old.weight - Real{1});
    max_relative_change = std::max(max_relative_change, error);
  }
  if (max_relative_change > Real{8192} * uni20::numeric_limits<Real>::epsilon())
    throw std::runtime_error("prepared/reference weights differ beyond the regression budget");
  auto measure = [&](auto&& evaluate) {
    auto const start = Clock::now();
    bethe::detail::CompensatedSum<Real> sum;
    for (auto const& level : scan.levels)
    {
      auto const factor = evaluate(level.state);
      if (!factor.converged()) throw std::runtime_error("unresolved timed form factor");
      sum.add(*factor.weight);
    }
    return std::pair{seconds(start), sum.value()};
  };
  std::cout << "precision=" << precision << " N=" << n << " states=" << scan.levels.size() << " repeats=" << repeats
            << " scalar_digits=" << uni20::numeric_limits<Real>::digits << "\nroot_scan_seconds=" << scan_seconds
            << " prepare_seconds=" << prepare_seconds
            << " max_relative_change=" << uni20::format_scalar(max_relative_change) << '\n';
  std::vector<double> baseline_times, prepared_times;
  for (std::size_t i = 0; i < repeats; ++i)
  {
    std::pair<double, Real> old, current;
    auto baseline = [&](auto const& state) {
      return test_support::unprepared_raising_form_factor(n, scan.ground_state, state);
    };
    // Alternate order to avoid always giving one implementation the warmest CPU.
    if (i % 2 == 0)
    {
      old = measure(baseline);
      current = measure(prepared);
    }
    else
    {
      current = measure(prepared);
      old = measure(baseline);
    }
    baseline_times.push_back(old.first);
    prepared_times.push_back(current.first);
    std::cout << "repeat=" << i + 1 << " unprepared_seconds=" << old.first << " prepared_seconds=" << current.first
              << " unprepared_weight=" << uni20::format_scalar(old.second)
              << " prepared_weight=" << uni20::format_scalar(current.second) << '\n';
  }
  std::sort(baseline_times.begin(), baseline_times.end());
  std::sort(prepared_times.begin(), prepared_times.end());
  auto median = [](auto const& values) {
    auto const size = values.size();
    return (values[size / 2] + values[(size - 1) / 2]) / 2;
  };
  std::cout << "median_unprepared_seconds=" << median(baseline_times)
            << " median_prepared_seconds=" << median(prepared_times)
            << " speedup_including_preparation=" << median(baseline_times) / (prepare_seconds + median(prepared_times))
            << '\n';
  return 0;
}

int main(int argc, char** argv)
try
{
  if (argc < 3 || argc > 4)
    throw std::invalid_argument("usage: bethe_xxx_form_factor_benchmark N fp64|long-double|fp128 [repeats=3]");
  auto parse = [](std::string_view value) {
    std::size_t result{};
    auto const [end, error] = std::from_chars(value.data(), value.data() + value.size(), result);
    if (error != std::errc{} || end != value.data() + value.size() || result == 0)
      throw std::invalid_argument("expected a positive integer");
    return result;
  };
  auto const n = parse(argv[1]), repeats = argc == 4 ? parse(argv[3]) : 3;
  if (n % 2) throw std::invalid_argument("N must be even");
  std::string_view const precision(argv[2]);
  if (precision == "fp64") return benchmark<double>(n, repeats, precision);
  if (precision == "long-double") return benchmark<long double>(n, repeats, precision);
#if UNI20_HAS_FLOAT128
  if (precision == "fp128") return benchmark<uni20::float128>(n, repeats, precision);
#endif
  throw std::invalid_argument("unsupported precision");
}
catch (std::exception const& error)
{
  std::cerr << error.what() << '\n';
  return 1;
}
