// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include "data-output-options.hpp"
#include "program-options.hpp"
#include <bethe/lee_yang.hpp>

namespace bethe::cli::lee_yang
{
struct Arguments
{
    std::optional<std::string> mass, length, tolerance, cutoff;
    std::size_t initial_intervals = 32, max_intervals = 2048, iterations = 1000, cutoffs = 3, products = 200000000;
    std::string precision = "fp64";
    DataOutputOptions output;
};
inline void add_options(CLI::App& app, Arguments& a, char const* tolerance_help)
{
  text_option(app, "--mass", a.mass, "Positive particle mass m; default 1")->type_name("REAL");
  text_option(app, "--length", a.length, "Positive circumference L")->required()->type_name("REAL");
  text_option(app, "--tolerance", a.tolerance, tolerance_help)->type_name("REAL");
  text_option(app, "--initial-cutoff", a.cutoff, "Positive initial rapidity cutoff; default automatic")
      ->type_name("REAL");
  count_option(app, "--initial-intervals", a.initial_intervals, "Initial intervals on positive half-line, >=2")
      ->capture_default_str();
  count_option(app, "--max-intervals", a.max_intervals, "Maximum intervals per grid, <=8192")->capture_default_str();
  count_option(app, "--max-iterations", a.iterations, "Fixed-point updates per source evaluation per grid")
      ->capture_default_str();
  count_option(app, "--max-cutoffs", a.cutoffs, "Total cutoff trials per state, including initial cutoff")
      ->capture_default_str();
  count_option(app, "--max-kernel-products", a.products, "Maximum folded kernel-times-logarithm terms per state")
      ->capture_default_str();
  precision_option(app, a.precision);
  add_data_output_options(app, a.output, true);
}
template <uni20::Real Real, typename Options> Options options_from(Arguments const& a, Options options)
{
  if (a.tolerance) options.tolerance = uni20::parse_real<Real>(*a.tolerance);
  if (a.cutoff) options.initial_cutoff = uni20::parse_real<Real>(*a.cutoff);
  options.initial_intervals = a.initial_intervals;
  options.max_intervals = a.max_intervals;
  options.max_iterations = a.iterations;
  options.max_cutoffs = a.cutoffs;
  options.max_kernel_products = a.products;
  return options;
}
inline char const* name(bethe::lee_yang::Status s)
{
  using bethe::lee_yang::Status;
  switch (s)
  {
    case Status::converged:
      return "converged";
    case Status::iteration_limit:
      return "iteration_limit";
    case Status::mesh_limit:
      return "mesh_limit";
    case Status::cutoff_limit:
      return "cutoff_limit";
    case Status::work_limit:
      return "work_limit";
    case Status::precision_limit:
      return "precision_limit";
  }
  return "unknown";
}
} // namespace bethe::cli::lee_yang
