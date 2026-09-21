// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include "cli-common.hpp"
#include <limits>
#include <optional>

namespace bethe::cli
{
struct XxzArguments
{
    std::size_t sites;
    std::optional<std::string_view> delta = std::nullopt;
    std::string_view precision = "fp64";
    std::optional<std::string_view> tolerance = std::nullopt;
    std::size_t max_iterations = 10000;
    std::optional<uni20::half_int> sz = std::nullopt;
    bool sectors = false;
    std::optional<std::size_t> excitation_count = std::nullopt;
    std::optional<std::size_t> max_candidates = std::nullopt;
    std::optional<std::vector<uni20::half_int>> quantum_numbers = std::nullopt;
    bool print_roots = false;
    std::string_view format = "auto";
};

inline XxzArguments parse_xxz_arguments(int argc, char** argv)
{
  XxzArguments args{.sites = parse_size(argv[1])};
  for (int i = 2; i < argc; ++i)
  {
    std::string_view const option = argv[i];
    if (option == "--roots")
      args.print_roots = true;
    else if (option == "--sectors")
      args.sectors = true;
    else if (option == "--delta" || option == "--sz" || option == "--precision" || option == "--tolerance" ||
             option == "--max-iterations" || option == "--format" || option == "--excitations" ||
             option == "--max-candidates" || option == "--quantum-numbers")
    {
      if (++i == argc) throw std::invalid_argument("missing value for " + std::string(option));
      if (option == "--delta")
        args.delta = argv[i];
      else if (option == "--sz")
        args.sz = uni20::half_int::parse(argv[i]);
      else if (option == "--precision")
        args.precision = argv[i];
      else if (option == "--tolerance")
        args.tolerance = argv[i];
      else if (option == "--format")
        args.format = argv[i];
      else if (option == "--excitations")
        args.excitation_count =
            std::string_view(argv[i]) == "all" ? std::numeric_limits<std::size_t>::max() : parse_size(argv[i]);
      else if (option == "--max-candidates")
        args.max_candidates = parse_size(argv[i]);
      else if (option == "--quantum-numbers")
        args.quantum_numbers = std::string_view(argv[i]) == "none" ? std::vector<uni20::half_int>{}
                                                                   : bethe::cli::parse_quantum_numbers(argv[i]);
      else
        args.max_iterations = parse_size(argv[i]);
    }
    else
      throw std::invalid_argument("unknown option: " + std::string(option));
  }
  if (!args.delta) throw std::invalid_argument("--delta VALUE is required");
  if (args.sectors && args.sz) throw std::invalid_argument("--sz and --sectors are mutually exclusive");
  if (args.excitation_count && (args.sectors || args.quantum_numbers))
    throw std::invalid_argument("--excitations is mutually exclusive with --sectors and --quantum-numbers");
  if (args.quantum_numbers && (args.sectors || args.sz))
    throw std::invalid_argument("--quantum-numbers is mutually exclusive with --sectors and --sz");
  if (args.max_candidates && !args.excitation_count)
    throw std::invalid_argument("--max-candidates requires --excitations COUNT|all");
  if (args.format != "auto" && args.format != "pretty" && args.format != "plain")
    throw std::invalid_argument("unknown output format: " + std::string(args.format));
  return args;
}
} // namespace bethe::cli
