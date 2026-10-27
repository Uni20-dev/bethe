// SPDX-License-Identifier: GPL-3.0-or-later
#include "two_stage_dsf_reference.hpp"
#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <uni20/async/tbb_scheduler.hpp>
#include <uni20/core/scalar_io.hpp>

template <uni20::Real R>
int run(std::string const& mode, std::size_t n, std::string const& precision, int workers, std::size_t window,
        std::size_t budget, std::string const& path)
{
  namespace model = bethe::heisenberg;
  uni20::async::TbbScheduler scheduler(workers);
  scheduler.execute_batch(std::size_t(workers), [](auto) {});
  bethe::ExecutionOptions options{&scheduler, window, window};
  auto start = std::chrono::steady_clock::now();
  auto cpu = std::clock();
  model::TwoSpinonStructureFactor<R> result;
  if (mode == "two-stage")
    result = dsf_reference::two_spinon_structure_factor<R>(n, {}, budget, options);
  else if (mode == "sync")
    result = model::two_spinon_structure_factor<R>(n, {}, budget, options);
  else
  {
    auto pending = model::two_spinon_structure_factor_async<R>(n, {}, budget, options);
    result = bethe::detail::take_scan_result(pending, scheduler);
  }
  double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  double cpu_seconds = double(std::clock() - cpu) / CLOCKS_PER_SEC;
  // Timing excludes serialization and hashing. Record native round-trip values,
  // all diagnostics and final ordered reductions, not just a checksum of energies.
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot open result file");
  auto value = [&](R x) { out << '\t' << uni20::format_scalar(x); };
  auto optional = [&](std::optional<R> x) {
    if (x)
      value(*x);
    else
      out << "\tmissing";
  };
  auto state = [&](auto const& s) {
    out << s.converged << '\t' << s.iterations << '\t' << s.sz.twice() << '\t' << s.momentum_index << '\t'
        << s.spin_reversed;
    value(s.energy);
    value(s.residual_norm);
    value(s.momentum);
    for (auto q : s.quantum_numbers)
      out << '\t' << q.twice();
    for (R root : s.rapidities)
      value(root);
    out << '\n';
  };
  out << result.scan.spin.twice() << '\t' << result.scan.candidate_count << '\t' << result.scan.converged_count << '\n';
  state(result.scan.ground_state);
  for (auto const& level : result.scan.levels)
  {
    optional(level.gap);
    out << '\t';
    state(level.state);
  }
  out << result.scan.first_unconverged.has_value() << '\n';
  if (result.scan.first_unconverged) state(*result.scan.first_unconverged);
  for (auto const& f : result.form_factors)
  {
    out << int(f.status);
    optional(f.weight);
    optional(f.log_weight);
    value(f.pivot_margin);
    out << '\n';
  }
  for (auto const& line : result.lines)
  {
    out << line.state_id << '\t' << line.momentum_index;
    value(line.momentum);
    value(line.gap);
    value(line.weight);
    out << '\n';
  }
  for (auto const& moment : result.moments)
  {
    value(moment.weight);
    value(moment.first_moment);
    out << '\n';
  }
  value(result.integrated_weight);
  value(result.weight_fraction);
  optional(result.first_moment_fraction);
  out << '\n';
  out.close();
  if (!out) throw std::runtime_error("result write failed");
  std::cout << std::setprecision(12) << "{\"mode\":" << std::quoted(mode) << ",\"sites\":" << n
            << ",\"precision\":" << std::quoted(precision) << ",\"workers\":" << workers << ",\"window\":" << window
            << ",\"budget\":" << budget << ",\"candidates\":" << result.scan.candidate_count
            << ",\"roots_ok\":" << result.scan.converged_count << ",\"weights_ok\":" << result.lines.size()
            << ",\"wall_seconds\":" << wall << ",\"cpu_seconds\":" << cpu_seconds
            << ",\"compiler\":" << std::quoted(__VERSION__) << "}\n";
  return 0;
}
int main(int argc, char** argv)
{
  try
  {
    if (argc != 8) throw std::invalid_argument("usage: benchmark MODE N PRECISION WORKERS WINDOW BUDGET RESULTS.tsv");
    std::string mode = argv[1], precision = argv[3];
    auto n = std::stoul(argv[2]), window = std::stoul(argv[5]), budget = std::stoul(argv[6]);
    int workers = std::stoi(argv[4]);
    if (n < 2 || n > 512 || n % 2 || workers < 1 || workers > 64 || !window || window > 1024 || !budget ||
        budget > 100000)
      throw std::invalid_argument("invalid bounded experiment dimensions");
    if (mode != "two-stage" && mode != "sync" && mode != "async") throw std::invalid_argument("unknown mode");
    if (precision == "fp64") return run<double>(mode, n, precision, workers, window, budget, argv[7]);
    if (precision == "long-double") return run<long double>(mode, n, precision, workers, window, budget, argv[7]);
#if UNI20_HAS_FLOAT128
    if (precision == "fp128") return run<uni20::float128>(mode, n, precision, workers, window, budget, argv[7]);
#endif
    throw std::invalid_argument("unavailable precision");
  }
  catch (std::exception const& e)
  {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
