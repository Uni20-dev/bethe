// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#include <bethe/heisenberg.hpp>
#include <uni20/core/scalar_io.hpp>

#include "heisenberg-report.hpp"

#include <charconv>
#include <ctime>
#include <fmt/format.h>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

namespace
{
#if UNI20_HAS_FLOAT128
// The pinned Uni20 scalar I/O falls back to long double in other MPLAPACK modes.
// Reject that configuration here rather than silently narrowing fp128 output.
static_assert(MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_FLOAT128,
              "fp128 CLI output requires MPLAPACK's _Float128/strfromf128 mode");
#endif

struct Arguments
{
  std::size_t sites;
  std::string_view precision = "fp64";
  std::optional<std::string_view> tolerance = std::nullopt;
  std::size_t max_iterations = 10000;
  bool print_roots = false;
  std::optional<uni20::half_int> sz = std::nullopt;
  std::optional<std::string_view> quantum_numbers = std::nullopt;
  bool sectors = false;
  bool spinons = false;
  std::string_view format = "auto";
};

// Process CPU time, not wall time. Keep timing in the front end and stop it
// before constructing reports or printing roots, so output costs are excluded.
class CpuTimer
{
public:
  std::string elapsed_text() const
  {
    auto const end = std::clock();
    if (start_ == std::clock_t{-1} || end == std::clock_t{-1} || end < start_)
      return "unavailable";
    auto const seconds = (static_cast<long double>(end) - static_cast<long double>(start_)) / CLOCKS_PER_SEC;
    return fmt::format("{:.6f} s", seconds);
  }

private:
  std::clock_t const start_ = std::clock();
};

void usage(std::ostream& out)
{
  out << "Usage: heisenberg-energy N [options]\n"
      << "Periodic spin-1/2 Heisenberg chain, J=1, zero field.\n"
      << "Default: ground state (one representative for odd N). Modes:\n"
      << "  --sz VALUE                         lowest energy in an Sz sector, e.g. 1/2\n"
      << "  --sectors                          lowest energy in every Sz sector\n"
      << "  --spinons                          odd-N one-spinon branch (Sz=1/2)\n"
      << "  --quantum-numbers I0,I1,...         specified finite real-root state\n"
      << "  --precision fp64|long-double|fp128  (default: fp64)\n"
      << "  --tolerance VALUE                  normalized equation residual\n"
      << "  --max-iterations COUNT             update budget (default: 10000)\n"
      << "  --roots                            print the rapidities\n"
      << "  --format auto|pretty|plain         terminal report or script output (default: auto)\n"
      << "  --help                             show this help\n"
      << "fp128 requires a Uni20 build with MPLAPACK enabled.\n";
}

std::size_t parse_size(std::string_view text)
{
  std::size_t value = 0;
  auto const parsed = std::from_chars(text.data(), text.data() + text.size(), value);
  if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size())
    throw std::invalid_argument("invalid nonnegative integer: " + std::string(text));
  return value;
}

bethe::heisenberg::QuantumNumbers parse_quantum_numbers(std::string_view text)
{
  bethe::heisenberg::QuantumNumbers numbers;
  // An explicitly empty list is the fully polarized state.
  if (text.empty())
    return numbers;
  for (;;)
  {
    auto const comma = text.find(',');
    numbers.push_back(uni20::half_int::parse(text.substr(0, comma)));
    if (comma == std::string_view::npos)
      return numbers;
    text.remove_prefix(comma + 1);
  }
}

template <uni20::Real Real> void print_roots(bethe::heisenberg::RealState<Real> const& state)
{
  std::cout << "# index rapidity I\n";
  for (std::size_t i = 0; i < state.rapidities.size(); ++i)
    std::cout << i << ' ' << uni20::format_real(state.rapidities[i]) << ' ' << state.quantum_numbers[i] << '\n';
}

int finish(bool converged)
{
  if (!converged)
    std::cerr << "The iteration budget was exhausted; each nonconverged energy is an unconverged estimate.\n";
  return converged ? 0 : 2;
}

template <uni20::Real Real> int run(Arguments const& args)
{
  bethe::heisenberg::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  if (args.tolerance)
    options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  bool const pretty = args.format == "pretty" || (args.format == "auto" && terminal::is_a_terminal(stdout));
  if (!pretty)
    std::cout << "Sites: " << args.sites << '\n'
              << "Precision: " << args.precision << '\n'
              << "Residual tolerance: " << uni20::format_real(options.residual_tolerance) << '\n';
  CpuTimer const timer;
  if (args.sectors)
  {
    auto const states = bethe::heisenberg::sector_ground_states<Real>(args.sites, options);
    auto const cpu_time = timer.elapsed_text();
    if (pretty)
      return finish(bethe::cli::print_sectors(args.sites, args.precision, options, states, args.print_roots, cpu_time));
    std::cout << "# CPU time: " << cpu_time << '\n' << "# Sz momentum_index P energy residual iterations converged\n";
    bool converged = true;
    for (auto const& state : states)
    {
      std::cout << state.sz << ' ' << state.momentum_index << ' ' << uni20::format_real(state.momentum) << ' '
                << uni20::format_real(state.energy) << ' ' << uni20::format_real(state.residual_norm) << ' '
                << state.iterations << ' ' << state.converged << '\n';
      if (args.print_roots)
        print_roots(state);
      converged = converged && state.converged;
    }
    return finish(converged);
  }
  if (args.spinons)
  {
    auto const branch = bethe::heisenberg::one_spinon_branch<Real>(args.sites, options);
    auto const cpu_time = timer.elapsed_text();
    if (pretty)
      return finish(bethe::cli::print_spinons(args.sites, args.precision, options, branch, args.print_roots, cpu_time));
    std::cout << "# CPU time: " << cpu_time << '\n'
              << "# Sz=1/2; k=pi/2-2*pi*I_h/N; bulk reference e_inf=1/4-log(2)\n"
              << "# hole k momentum_index P energy E_minus_N_e_inf epsilon_inf residual iterations converged\n";
    bool converged = true;
    for (auto const& point : branch)
    {
      auto const& state = point.state;
      std::cout << point.hole << ' ' << uni20::format_real(point.spinon_momentum) << ' ' << state.momentum_index << ' '
                << uni20::format_real(state.momentum) << ' ' << uni20::format_real(state.energy) << ' '
                << uni20::format_real(point.bulk_subtracted_energy) << ' '
                << uni20::format_real(bethe::heisenberg::spinon_energy(point.spinon_momentum)) << ' '
                << uni20::format_real(state.residual_norm) << ' ' << state.iterations << ' ' << state.converged << '\n';
      if (args.print_roots)
        print_roots(state);
      converged = converged && state.converged;
    }
    return finish(converged);
  }
  auto const result =
      args.quantum_numbers
          ? bethe::heisenberg::solve_real<Real>(args.sites, parse_quantum_numbers(*args.quantum_numbers), options)
      : args.sz ? bethe::heisenberg::sector_ground_state<Real>(args.sites, *args.sz, options)
                : bethe::heisenberg::ground_state<Real>(args.sites, options);
  auto const cpu_time = timer.elapsed_text();
  if (pretty)
    return finish(bethe::cli::print_state(args.sites, args.precision, options, result,
                                          args.quantum_numbers ? "specified real-root state"
                                          : args.sz            ? "sector minimum"
                                                               : "ground state",
                                          args.print_roots, cpu_time));
  std::cout << "Sz: " << result.sz << '\n'
            << "Spin-reversed reference: " << result.spin_reversed << '\n'
            << "Momentum index: " << result.momentum_index << '\n'
            << "Momentum: " << uni20::format_real(result.momentum) << '\n'
            << "Status: " << (result.converged ? "converged" : "iteration limit reached") << '\n'
            << "Iterations: " << result.iterations << '\n'
            << "CPU time: " << cpu_time << '\n'
            << "Residual norm: " << uni20::format_real(result.residual_norm) << '\n'
            << "Total energy: " << uni20::format_real(result.energy) << '\n'
            << "Energy per site: " << uni20::format_real(result.energy / static_cast<Real>(args.sites)) << '\n';
  if (args.print_roots)
    print_roots(result);
  return finish(result.converged);
}
} // namespace

int main(int argc, char** argv)
{
  if (argc == 2 && std::string_view(argv[1]) == "--help")
  {
    usage(std::cout);
    return 0;
  }
  if (argc < 2)
  {
    usage(std::cerr);
    return 1;
  }
  try
  {
    Arguments args{.sites = parse_size(argv[1])};
    for (int i = 2; i < argc; ++i)
    {
      std::string_view const option = argv[i];
      if (option == "--roots")
        args.print_roots = true;
      else if (option == "--sectors")
        args.sectors = true;
      else if (option == "--spinons")
        args.spinons = true;
      else if (option == "--precision" || option == "--tolerance" || option == "--max-iterations" || option == "--sz" ||
               option == "--quantum-numbers" || option == "--format")
      {
        if (++i == argc)
          throw std::invalid_argument("missing value for " + std::string(option));
        if (option == "--precision")
          args.precision = argv[i];
        else if (option == "--format")
          args.format = argv[i];
        else if (option == "--tolerance")
          args.tolerance = argv[i];
        else if (option == "--sz")
          args.sz = uni20::half_int::parse(argv[i]);
        else if (option == "--quantum-numbers")
          args.quantum_numbers = argv[i];
        else
          args.max_iterations = parse_size(argv[i]);
      }
      else
        throw std::invalid_argument("unknown option: " + std::string(option));
    }
    if (static_cast<int>(args.sz.has_value()) + static_cast<int>(args.quantum_numbers.has_value()) +
            static_cast<int>(args.sectors) + static_cast<int>(args.spinons) >
        1)
      throw std::invalid_argument("--sz, --sectors, --spinons and --quantum-numbers are mutually exclusive");
    if (args.format != "auto" && args.format != "pretty" && args.format != "plain")
      throw std::invalid_argument("unknown output format: " + std::string(args.format));
    if (args.precision == "fp64")
      return run<double>(args);
    if (args.precision == "long-double")
      return run<long double>(args);
    if (args.precision == "fp128")
    {
#if UNI20_HAS_FLOAT128
      return run<uni20::float128>(args);
#else
      throw std::invalid_argument("fp128 is unavailable; configure with -DUNI20_ENABLE_MPLAPACK=ON");
#endif
    }
    throw std::invalid_argument("unknown precision: " + std::string(args.precision));
  }
  catch (std::exception const& error)
  {
    std::cerr << "heisenberg-energy: " << error.what() << '\n';
    return 1;
  }
}
