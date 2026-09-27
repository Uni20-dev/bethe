// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include "program-options.hpp"
#include "result-output.hpp"
#include <bethe/sine_gordon_vacuum.hpp>

namespace bethe::cli::sine_gordon
{
struct Arguments
{
    std::optional<std::string> mass, length, coupling, tolerance, contour, cutoff;
    std::size_t initial_intervals = 64, max_intervals = 2048, iterations = 10000, cutoffs = 3;
    std::size_t kernel_evaluations = 100000, kernel_levels = 16, fourier_cutoffs = 64;
    std::string precision = "fp64";
    DataOutputOptions output;
};
inline void add_options(CLI::App& app, Arguments& a, std::string coupling_help, std::string tolerance_help)
{
  text_option(app, "--mass", a.mass, "Positive soliton mass M; default 1")->type_name("REAL");
  text_option(app, "--length", a.length, "Positive circumference L")->required()->type_name("REAL");
  text_option(app, "--p", a.coupling, std::move(coupling_help))->required()->type_name("REAL");
  text_option(app, "--tolerance", a.tolerance, std::move(tolerance_help))->type_name("REAL");
  text_option(app, "--contour-shift", a.contour, "eta in (0,pi*min(1,p)/2); default pi*min(1,p)/4")->type_name("REAL");
  text_option(app, "--initial-cutoff", a.cutoff, "Initial rapidity cutoff; default estimated from tolerance and u")
      ->type_name("REAL");
  count_option(app, "--initial-intervals", a.initial_intervals, "Initial even rapidity intervals, >=8")
      ->capture_default_str();
  count_option(app, "--max-intervals", a.max_intervals, "Maximum rapidity intervals, <=16384")->capture_default_str();
  count_option(app, "--max-iterations", a.iterations, "Total nonlinear updates per state, across both contours")
      ->capture_default_str();
  count_option(app, "--max-cutoffs", a.cutoffs, "Rapidity cutoff attempts per contour")->capture_default_str();
  count_option(app, "--max-kernel-evaluations", a.kernel_evaluations, "Fourier evaluations per kernel/source table")
      ->capture_default_str();
  count_option(app, "--max-kernel-levels", a.kernel_levels, "Fourier quadrature refinement levels per table")
      ->capture_default_str();
  count_option(app, "--max-fourier-cutoffs", a.fourier_cutoffs, "Fourier cutoff attempts per table")
      ->capture_default_str();
  precision_option(app, a.precision);
}
template <uni20::Real Real> void configure(Arguments const& a, bethe::sine_gordon::VacuumOptions<Real>& controls)
{
  if (a.tolerance) controls.tolerance = uni20::parse_real<Real>(*a.tolerance);
  if (a.contour) controls.contour_shift = uni20::parse_real<Real>(*a.contour);
  if (a.cutoff) controls.initial_cutoff = uni20::parse_real<Real>(*a.cutoff);
  controls.initial_intervals = a.initial_intervals;
  controls.max_intervals = a.max_intervals;
  controls.max_iterations = a.iterations;
  controls.max_cutoffs = a.cutoffs;
  controls.kernel.max_evaluations = a.kernel_evaluations;
  controls.kernel.max_levels = a.kernel_levels;
  controls.kernel.max_cutoffs = a.fourier_cutoffs;
}
inline char const* name(bethe::sine_gordon::VacuumStatus status)
{
  using S = bethe::sine_gordon::VacuumStatus;
  switch (status)
  {
    case S::converged:
      return "converged";
    case S::kernel_limit:
      return "kernel_limit";
    case S::iteration_limit:
      return "iteration_limit";
    case S::mesh_limit:
      return "mesh_limit";
    case S::cutoff_limit:
      return "cutoff_limit";
    case S::contour_limit:
      return "contour_limit";
    case S::precision_limit:
      return "precision_limit";
  }
  return "unknown";
}
template <uni20::Real Real>
void report_controls(RunReport& report, bethe::sine_gordon::VacuumOptions<Real> const& controls,
                     bethe::sine_gordon::VacuumState<Real> const& state)
{
  report.field("tolerance", "Absolute Y tolerance per state", controls.tolerance)
      .field("contour_shift", "Contour shift", state.contour_shift)
      .field("verification_contour_shift", "Verification contour shift", state.verification_contour_shift)
      .field("initial_cutoff", "Requested initial rapidity cutoff", controls.initial_cutoff)
      .field("initial_intervals", "Initial intervals", controls.initial_intervals)
      .field("max_intervals", "Max intervals", controls.max_intervals)
      .field("max_iterations", "Max total nonlinear updates per state", controls.max_iterations)
      .field("max_cutoffs", "Max rapidity cutoffs per contour", controls.max_cutoffs)
      .field("max_kernel_evaluations", "Max evaluations per kernel/source table", controls.kernel.max_evaluations)
      .field("max_kernel_levels", "Max kernel levels", controls.kernel.max_levels)
      .field("max_fourier_cutoffs", "Max Fourier cutoffs", controls.kernel.max_cutoffs);
}
} // namespace bethe::cli::sine_gordon
