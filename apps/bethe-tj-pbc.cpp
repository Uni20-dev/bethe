// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "result-output.hpp"
#include <bethe/tj.hpp>

namespace
{
namespace model = bethe::tj;
namespace cli = bethe::cli;
struct Arguments
{
    std::size_t sites = 0, max_iterations = 10000;
    std::optional<std::size_t> particles;
    std::optional<uni20::half_int> sz;
    std::optional<std::string> tolerance;
    std::string precision = "fp64";
    cli::DataOutputOptions output;
    bool roots = false;
};
auto program_info()
{
  auto info =
      bethe::cli::program_info("bethe-tj-pbc", "Supersymmetric periodic t-J sector ground state, t=1, J=2, L>=3.",
                               bethe::citations::Tool::tj_pbc);
  info.notes = {"H=-sum(projected hopping+h.c.)+2 sum(S_i.S_j-n_i*n_j/4).",
                "No double occupancy and no chemical-potential energy shift.",
                "Doped mixed-spin sectors require odd N_up AND odd N_down.",
                "No-hole sectors and fully polarized free fermions allow all populations.",
                "Nested roots use Sutherland lambda,mu; E=2*N_h-sum 1/(lambda^2+1/4).",
                "Momentum includes fermionic translation signs, also in the no-hole limit.",
                "Other doped shell parities, excitations, open ends, and J!=2t are not implemented.",
                "See docs/tj.md for the supported branches and conventions.",
                "Use --references for literature and applicability; see CITATIONS.md."};
  return info;
}
void add_options(CLI::App& app, Arguments& args)
{
  bethe::cli::option(app, "L", args.sites, "Number of sites or rungs")->required();
  bethe::cli::option(app, "--particles", args.particles, "0 <= N <= L (default: L)");
  bethe::cli::option(app, "--sz", args.sz, "spin projection (default: 0 for even N, 1/2 for odd N)");
  bethe::cli::option(app, "--roots", args.roots, "print nested rapidities/labels or free modes");
  bethe::cli::option(app, "--tolerance", args.tolerance,
                     "max logarithmic residual divided by L (default: 32 epsilon; not an energy-error bound)");
  bethe::cli::option(app, "--max-iterations", args.max_iterations, "accepted updates (default: 10000)")
      ->capture_default_str();
  bethe::cli::precision_option(app, args.precision);
  cli::add_data_output_options(app, args.output, true);
}

void validate(Arguments const& args) { args.output.validate(); }
char const* status(model::SolveStatus value)
{
  switch (value)
  {
    case model::SolveStatus::converged:
      return "converged";
    case model::SolveStatus::iteration_limit:
      return "iteration limit reached; unconverged estimate";
    case model::SolveStatus::stalled:
      return "line search or representable precision stalled; unconverged estimate";
    case model::SolveStatus::ill_conditioned:
      return "Newton system unresolved; unconverged estimate";
  }
  return "unknown";
}
template <uni20::Real Real> int run(Arguments const& args, int argc, char** argv)
{
  model::detail::check_counts(args.sites, 0, 0);
  auto const particles = args.particles.value_or(args.sites);
  if (particles > args.sites) throw std::invalid_argument("t-J requires 0 <= particles <= L");
  auto const n = std::int64_t(particles);
  auto const sz = args.sz.value_or(uni20::from_twice(n % 2));
  auto const s = sz.twice();
  if (s < -n || s > n || (n - s) % 2)
    throw std::invalid_argument("Sz must satisfy |2*Sz|<=N and have the same parity as N");
  auto const down = std::size_t((n - s) / 2), up = particles - down;
  model::SolverOptions<Real> options;
  options.max_iterations = args.max_iterations;
  if (args.tolerance) options.residual_tolerance = uni20::parse_real<Real>(*args.tolerance);
  cli::CpuTimer const timer;
  auto const state = model::ground_state<Real>(args.sites, up, down, options);
  auto const cpu_time = timer.elapsed_text();
  char const* branch = state.branch == model::Branch::sutherland       ? "Sutherland real-root sector"
                       : state.branch == model::Branch::polarized_free ? "exact polarized free fermions"
                                                                       : "no-hole XXX reduction";
  cli::report_builder report("Supersymmetric t-J chain (periodic)");
  report.status(state.converged ? cli::semantic_glyph::success : cli::semantic_glyph::warning, status(state.status))
      .field("Hamiltonian", "t=1, J=2; projected hopping + 2*(S.S-nn/4)")
      .field("Calculation", branch)
      .field("Sites", state.sites)
      .field("Particles", particles)
      .field("N_up", up)
      .field("N_down", down)
      .field("Holes", state.holes)
      .field("Sz", uni20::to_string_fraction(sz))
      .field("Precision", args.precision)
      .field("Residual tolerance", uni20::format_real(options.residual_tolerance))
      .field("Status", status(state.status))
      .field("Total energy", uni20::format_real(state.energy))
      .field("Energy per site", uni20::format_real(state.energy / Real(state.sites)))
      .field("Momentum index", state.momentum_index)
      .field("Momentum P", uni20::format_real(state.momentum))
      .field("First-level roots", state.rapidities[0].size())
      .field("Second-level roots", state.rapidities[1].size())
      .field("Residual norm", uni20::format_real(state.residual_norm))
      .field("Iterations", state.iterations)
      .field("CPU time", cpu_time);
  std::vector<std::string> tables{"states"};
  bool const free = state.branch == model::Branch::polarized_free;
  if (args.roots)
  {
    if (free)
      tables.push_back("free_modes");
    else
      tables.insert(tables.end(), {"first_roots", "second_roots"});
  }
  cli::ResultOutput output(report, args.output, "bethe-tj-pbc", argc, argv, tables);
  output.table(
      "states", "State",
      [&](auto& table) {
        table.append(0, sz, state.energy, state.momentum_index, state.momentum, state.residual_norm, state.iterations,
                     state.converged, std::string(status(state.status)));
      },
      cli::column<std::size_t>("state_id"), cli::column<uni20::half_int>("sz", "Sz"),
      cli::column<Real>("energy", "Energy"), cli::column<std::size_t>("momentum_index"), cli::column<Real>("p", "P"),
      cli::column<Real>("residual"), cli::column<std::size_t>("iterations"), cli::column<bool>("converged"),
      cli::column<std::string>("status"));
  if (args.roots && free)
    output.table(
        "free_modes", "Occupied free modes",
        [&](auto& table) {
          for (std::size_t j = 0; j < state.free_modes.size(); ++j)
            table.append(0, j, state.free_modes[j]);
        },
        cli::column<std::size_t>("state_id"), cli::column<std::size_t>("index", "Index"),
        cli::column<std::int64_t>("mode", "j (k=2*pi*j/L)"));
  else if (args.roots)
    for (std::size_t a = 0; a < 2; ++a)
      output.table(
          a ? "second_roots" : "first_roots", a ? "Second-level rapidities" : "First-level rapidities",
          [&](auto& table) {
            for (std::size_t j = 0; j < state.rapidities[a].size(); ++j)
              table.append(0, j, state.quantum_numbers[a][j], state.rapidities[a][j]);
          },
          cli::column<std::size_t>("state_id"), cli::column<std::size_t>("index", "Index"),
          cli::column<uni20::half_int>("quantum_number", a ? "J" : "I"),
          cli::column<Real>("rapidity", a ? "mu" : "lambda"));
  output.finish();
  if (!state.converged) std::cerr << "t-J solve incomplete; consider a larger budget or higher precision.\n";
  return state.converged ? 0 : 2;
}
} // namespace
int main(int argc, char** argv)
{
  Arguments args;
  return bethe::cli::program_main(
      argc, argv, program_info(), [&](auto& app) { add_options(app, args); },
      [&](auto&) {
        validate(args);
        return bethe::cli::dispatch_precision(args.precision,
                                              [&]<uni20::Real Real> { return run<Real>(args, argc, argv); });
      });
}
