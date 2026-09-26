// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "biquadratic/analytic.hpp"
#include "biquadratic/clusters.hpp"
#include "biquadratic/qsystem.hpp"
#include "biquadratic/real.hpp"
#include "biquadratic/singlet.hpp"

namespace bethe::apps::biquadratic
{
template <uni20::Real Real> int run(Arguments const& args, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  bethe::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  if (args.tolerance) options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  if (!uni20::isfinite(options.residual_tolerance) || options.residual_tolerance <= Real{0})
    throw std::invalid_argument("residual tolerance must be finite and positive");
  auto report = preamble(context, args, options);
  if (args.q_seed || args.q_spectrum) return run_qsystem(args, options, std::move(report));
  if (args.singlet_excitation) return run_singlet(args, options, std::move(report));
  if (args.bound_pairs)
    return run_bound_clusters<Real, ClusterFamily::pair>(args, *args.bound_pairs, options, std::move(report));
  if (args.bound_triples)
    return run_bound_clusters<Real, ClusterFamily::triple>(args, *args.bound_triples, options, std::move(report));
  if (args.bound_quartets)
    return run_bound_clusters<Real, ClusterFamily::quartet>(args, *args.bound_quartets, options, std::move(report));
  if (args.two_pairs || args.two_pair_states)
    return run_bound_clusters<Real, ClusterFamily::two_pairs>(args, args.two_pair_states.value_or(1), options,
                                                              std::move(report));
  if (args.triple_defect || args.triple_defects)
    return run_bound_clusters<Real, ClusterFamily::triple_real>(args, args.triple_defects.value_or(1), options,
                                                                std::move(report));
  if (args.pair_defect || args.pair_defects)
    return run_bound_clusters<Real, ClusterFamily::pair_real>(args, args.pair_defects.value_or(1), options,
                                                              std::move(report));
  if (args.ferromagnetic && !args.excitations && !args.numbers)
    return run_ferro_analytic<Real>(args, std::move(report));
  return run_real(args, options, std::move(report));
}
} // namespace bethe::apps::biquadratic
int main(int argc, char** argv)
{
  namespace frontend = bethe::apps::biquadratic;
  frontend::Arguments args;
  return bethe::cli::program_main(
      argc, argv, frontend::program_info(), [&](auto& app) { frontend::add_options(app, args); },
      [&](auto&) {
        frontend::validate(args);
        return bethe::cli::dispatch_precision(args.precision,
                                              [&]<uni20::Real Real> { return frontend::run<Real>(args, argc, argv); });
      });
}
