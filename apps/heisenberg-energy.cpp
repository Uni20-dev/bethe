// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#include <bethe/heisenberg.hpp>
#include <uni20/core/scalar_io.hpp>

#include <charconv>
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
  std::string_view precision = "long-double";
  std::optional<std::string_view> tolerance = std::nullopt;
  std::size_t max_iterations = 10000;
  bool print_roots = false;
};

void usage(std::ostream& out)
{
  out << "Usage: heisenberg-energy N [options]\n"
      << "Even-N periodic spin-1/2 Heisenberg ground state, J=1, zero field.\n"
      << "  --precision fp64|long-double|fp128  (default: long-double)\n"
      << "  --tolerance VALUE                  normalized equation residual\n"
      << "  --max-iterations COUNT             update budget (default: 10000)\n"
      << "  --roots                            print the rapidities\n"
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

template <uni20::Real Real> int run(Arguments const& args)
{
  bethe::heisenberg::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  if (args.tolerance)
    options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  auto const result = bethe::heisenberg::ground_state<Real>(args.sites, options);

  std::cout << "Sites: " << args.sites << '\n'
            << "Precision: " << args.precision << '\n'
            << "Status: " << (result.converged ? "converged" : "iteration limit reached") << '\n'
            << "Iterations: " << result.iterations << '\n'
            << "Residual tolerance: " << uni20::format_real(options.residual_tolerance) << '\n'
            << "Residual norm: " << uni20::format_real(result.residual_norm) << '\n'
            << "Total energy: " << uni20::format_real(result.energy) << '\n'
            << "Energy per site: " << uni20::format_real(result.energy / static_cast<Real>(args.sites)) << '\n';
  if (args.print_roots)
  {
    std::cout << "# index rapidity\n";
    for (std::size_t i = 0; i < result.rapidities.size(); ++i)
      std::cout << i << ' ' << uni20::format_real(result.rapidities[i]) << '\n';
  }
  if (!result.converged)
    std::cerr << "The iteration budget was exhausted; the reported energy is an unconverged estimate.\n";
  return result.converged ? 0 : 2;
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
      else if (option == "--precision" || option == "--tolerance" || option == "--max-iterations")
      {
        if (++i == argc)
          throw std::invalid_argument("missing value for " + std::string(option));
        if (option == "--precision")
          args.precision = argv[i];
        else if (option == "--tolerance")
          args.tolerance = argv[i];
        else
          args.max_iterations = parse_size(argv[i]);
      }
      else
        throw std::invalid_argument("unknown option: " + std::string(option));
    }
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
