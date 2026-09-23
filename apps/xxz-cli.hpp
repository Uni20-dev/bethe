// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include "program-options.hpp"
#include <limits>
#include <optional>

namespace bethe::cli
{
struct XxzArguments
{
    std::size_t sites = 0;
    std::optional<std::string> delta = std::nullopt;
    std::string precision = "fp64";
    std::optional<std::string> tolerance = std::nullopt;
    std::size_t max_iterations = 10000;
    std::optional<uni20::half_int> sz = std::nullopt;
    bool sectors = false;
    std::optional<std::size_t> excitation_count = std::nullopt;
    std::optional<std::size_t> max_candidates = std::nullopt;
    std::optional<std::vector<uni20::half_int>> quantum_numbers = std::nullopt;
    bool print_roots = false;
    std::string format = "auto";
};

inline void add_xxz_options(CLI::App& app, XxzArguments& args)
{
  count_option(app, "N", args.sites, "Number of sites")->required();
  text_option(app, "--delta", args.delta, "Anisotropy; ground states >-1 (negative odd rings excluded)")->required();
  option(app, "--sz", args.sz, "Sector magnetization; scan default: 1 even N, 1/2 odd N");
  option(app, "--sectors", args.sectors, "Lowest energy in every Sz sector")->excludes("--sz");
  all_count_option(app, "--excitations", args.excitation_count,
                   "Lowest COUNT, or all, states in the finite-real window")
      ->excludes("--sectors");
  option(app, "--max-candidates", args.max_candidates, "Exhaustive scan limit (default: 10000)")
      ->needs("--excitations");
  option(app, "--quantum-numbers", args.quantum_numbers, "Explicit real-root labels; none or empty for vacuum")
      ->excludes("--excitations")
      ->excludes("--sectors")
      ->excludes("--sz");
  option(app, "--roots", args.print_roots, "Print roots and exact labels");
  text_option(app, "--tolerance", args.tolerance, "Residual in the reported convention; default: 32 epsilon");
  count_option(app, "--max-iterations", args.max_iterations, "Update budget")->capture_default_str();
  precision_option(app, args.precision);
  app.add_option("--format", args.format, "Stdout layout")
      ->check(CLI::IsMember({"auto", "pretty", "plain"}))
      ->capture_default_str();
}
inline void validate_xxz_arguments(XxzArguments const& args)
{
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
}
} // namespace bethe::cli
