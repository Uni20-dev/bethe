// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "result-output.hpp"
#include <bethe/xxz_quantum_group_critical.hpp>
#include <bethe/xxz_quantum_group_free.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::xxz::quantum_group::critical;
struct Arguments
{
    std::size_t sites = 0;
    std::optional<std::size_t> through_lines, iterations, max_blocks, max_mode_entries;
    std::optional<std::vector<uni20::half_int>> numbers;
    std::string delta, precision = "fp64";
    std::optional<std::string> tolerance, sz;
    bool roots = false;
    cli::DataOutputOptions output;
};
auto program_info()
{
  auto info = cli::program_info("bethe-xxz-qg-obc", "Quantum-group XXZ with opposite imaginary end fields.",
                                bethe::citations::Tool::xxz_qg_obc);
  info.notes = {
      "Spin-half J=1: H=sum(SxSx+SySy+Delta SzSz)+i*sqrt(1-Delta^2)/2*(Sz_1-Sz_N).",
      "Current scope: N>=2, 0<Delta<1, positive finite-real Bethe roots only; not free-end XXZ.",
      "Delta=0: complete fixed-Sz free-fermion spectrum with Jordan-block sizes; --sz defaults to (N mod 2)/2.",
      "Default: consecutive labels with ell=N mod 2. --through-lines selects ell=N-2M, not physical SU(2) spin.",
      "--numbers selects increasing integer labels, or none for the polarized state.",
      "Regular-root branch has no completeness or degeneracy claim. Neither branch returns Jordan vectors or CFT fits.",
      "Tolerance controls max|F|/(2N), not energy error; default 32 epsilon. Failed energies are missing (exit 2).",
      "Tables: state; optional roots contains last-iterate coordinates even on failure.",
      "Delta=0 table: blocks, ordered by zero occupation then mode labels, not energy. Equal energies are not merged.",
      "See docs/xxz-nonhermitian.md; --references for literature."};
  info.examples = {
      {"bethe-xxz-qg-obc 32 --delta 0.25 --roots", "Consecutive-label sea"},
      {"bethe-xxz-qg-obc 8 --delta 0.6 --numbers 1,3 --precision fp128 --json state.json", "Selected regular state"}};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  cli::count_option(app, "N", a.sites, "Number of sites")->required();
  app.add_option("--delta", a.delta, "Anisotropy, 0<=Delta<1; zero selects the complete sector spectrum")
      ->required()
      ->type_name("REAL");
  auto* lines = cli::count_option(app, "--through-lines", a.through_lines,
                                  "Consecutive-label sea with ell=N-2M; default N mod 2");
  cli::option(app, "--numbers", a.numbers, "Explicit Bethe labels, e.g. 1,3; none selects M=0")->excludes(lines);
  cli::text_option(app, "--tolerance", a.tolerance, "Normalized logarithmic residual tolerance");
  cli::count_option(app, "--max-iterations", a.iterations, "Simultaneous root updates (regular branch)")
      ->default_str("10000");
  cli::text_option(app, "--sz", a.sz, "Delta=0 only: fixed magnetization, including negative sectors");
  cli::count_option(app, "--max-blocks", a.max_blocks, "Delta=0 only: maximum output blocks")->default_str("100000");
  cli::count_option(app, "--max-mode-entries", a.max_mode_entries, "Delta=0 only: total occupied-mode storage budget")
      ->default_str("1000000");
  app.add_flag("--roots", a.roots, "Export last-iterate roots with convergence status");
  cli::precision_option(app, a.precision);
  cli::add_data_output_options(app, a.output, true);
}
char const* name(model::Status s)
{
  switch (s)
  {
    case model::Status::converged:
      return "converged";
    case model::Status::iteration_limit:
      return "iteration_limit";
    case model::Status::precision_limit:
      return "precision_limit";
  }
  return "unknown";
}
template <uni20::Real Real> int run_free(Arguments const& a, uni20::run_context& context)
{
  namespace free = bethe::xxz::quantum_group::free;
  if (a.through_lines || a.numbers || a.tolerance || a.iterations || a.roots)
    throw std::invalid_argument("Delta=0 uses --sz; regular-root labels, roots and solver controls do not apply");
  auto const n = bethe::xxz::detail::checked_sites(a.sites);
  auto const sz = a.sz ? uni20::half_int::parse(*a.sz) : uni20::from_twice(n % 2);
  // The existing sector helper validates parity and bounds but folds Sz to
  // |Sz|. Retain the signed sector for this complete spin-basis spectrum.
  bethe::xxz::detail::sector_roots(a.sites, sz);
  auto const down = static_cast<std::size_t>((n - sz.twice()) / 2);
  free::SectorOptions controls;
  if (a.max_blocks) controls.max_blocks = *a.max_blocks;
  if (a.max_mode_entries) controls.max_mode_entries = *a.max_mode_entries;
  auto const blocks = context.measure([&] { return free::sector<Real>(a.sites, down, controls); });
  std::size_t dimension = 0, defective = 0;
  for (auto const& block : blocks)
  {
    dimension += block.block_size;
    defective += block.block_size == 2;
  }
  cli::RunReport report(context, "Non-Hermitian quantum-group XX complete sector");
  report.field("hamiltonian", "Hamiltonian", "sum(SxSx+SySy)+i/2*(Sz_1-Sz_N); spin-half J=1")
      .field("sites", "Sites", a.sites)
      .field("delta", "Delta", Real{0})
      .field("boundary_strength", "Imaginary end-field coefficient", Real{1} / Real{2})
      .field("sz", "Sz", sz)
      .field("down", "Down spins", down)
      .field("scope", "Scope", "complete fixed-Sz Hamiltonian spectrum and block sizes; no Jordan vectors or CFT fit")
      .field("energy_reference", "Energy reference", "absolute Hamiltonian energy; polarized energy is zero")
      .field("mode_labels", "Mode labels",
             "occupied k with energy cos(pi*k/N); exclude k=N/2 for even N; comma-separated")
      .field("block_semantics", "Block semantics",
             "one row per Jordan block, not distinct energy; size is algebraic dimension")
      .field("ordering", "Ordering", "zero occupation then lexicographic modes; not energy sorted")
      .field("dimension", "Sector dimension", dimension)
      .field("blocks", "Jordan blocks", blocks.size())
      .field("defective_blocks", "Size-two blocks", defective)
      .field("max_blocks", "Max blocks", controls.max_blocks)
      .field("max_mode_entries", "Max mode entries", controls.max_mode_entries)
      .field("precision", "Precision", a.precision)
      .result(true, "complete");
  cli::ResultOutput output(report, a.output, {"blocks"});
  using cli::column;
  output.table(
      "blocks", "Hamiltonian Jordan blocks",
      [&](auto& table) {
        for (std::size_t i = 0; i < blocks.size(); ++i)
        {
          auto const& b = blocks[i];
          std::string modes;
          for (auto k : b.modes)
          {
            if (!modes.empty()) modes += ',';
            modes += std::to_string(k);
          }
          table.append(i, b.energy, b.block_size, b.zero_occupation, modes);
        }
      },
      column<std::size_t>("block_id"), column<Real>("energy"), column<std::size_t>("block_size"),
      column<std::size_t>("zero_occupation"), column<std::string>("modes"));
  output.finish();
  return 0;
}
template <uni20::Real Real> int run(Arguments const& a, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  Real const delta = uni20::parse_real<Real>(a.delta);
  if (delta == Real{0}) return run_free<Real>(a, context);
  if (a.sz || a.max_blocks || a.max_mode_entries)
    throw std::invalid_argument("--sz, --max-blocks and --max-mode-entries require Delta=0");
  bethe::SolverOptions<Real> controls;
  controls.max_iterations = a.iterations.value_or(10000);
  if (a.tolerance) controls.residual_tolerance = uni20::parse_real<Real>(*a.tolerance);
  // All argument/branch validation and the solve precede opening --force targets.
  auto const s = context.measure([&] {
    return a.numbers ? model::solve_real<Real>(a.sites, delta, *a.numbers, controls)
                     : model::sea_state<Real>(a.sites, delta, a.through_lines.value_or(a.sites % 2), controls);
  });
  auto const m = s.quantum_numbers.size(), ell = a.sites - 2 * m;
  cli::RunReport report(context, "Non-Hermitian quantum-group XXZ regular state");
  report.field("hamiltonian", "Hamiltonian", "sum(SxSx+SySy+Delta SzSz)+i*sqrt(1-Delta^2)/2*(Sz_1-Sz_N); spin-half J=1")
      .field("sites", "Sites", a.sites)
      .field("delta", "Delta", delta)
      .field("boundary_strength", "Imaginary end-field coefficient",
             std::sqrt((Real{1} - delta) * (Real{1} + delta)) / Real{2})
      .field("roots", "Bethe roots", m)
      .field("through_lines", "ell=N-2M", ell)
      .field("sz", "Sz", uni20::from_twice(static_cast<std::int64_t>(ell)))
      .field("selection", "Selection", a.numbers ? "explicit labels" : "consecutive-label sea")
      .field("scope", "Scope", "regular positive-real family only; no completeness, multiplicity or Jordan claim")
      .field("energy_reference", "Energy reference",
             "absolute Hamiltonian energy; shift is relative to E_F=(N-1)*Delta/4")
      .field("root_coordinates", "Root coordinates",
             "z=tanh(lambda)/tan(gamma/2), Delta=cos(gamma); last iterate on failure")
      .field("residual", "Residual convention", "max|2N*theta1-sum(theta2_minus+theta2_plus)-2*pi*I|/(2N)")
      .field("precision", "Precision", a.precision)
      .field("tolerance", "Residual tolerance", controls.residual_tolerance)
      .field("max_iterations", "Max iterations", controls.max_iterations)
      .result(s.converged, name(s.status));
  std::vector<std::string> tables{"state"};
  if (a.roots) tables.push_back("roots");
  cli::ResultOutput output(report, a.output, std::move(tables));
  using cli::column;
  using Optional = std::optional<Real>;
  output.table(
      "state", "Regular Bethe state",
      [&](auto& table) {
        table.append(s.energy, s.energy_shift, s.residual_norm, s.iterations, s.converged, std::string(name(s.status)));
      },
      column<Optional>("energy"), column<Optional>("energy_shift"), column<Real>("residual"),
      column<std::size_t>("iterations"), column<bool>("converged"), column<std::string>("status"));
  if (a.roots)
    output.table(
        "roots", "Last-iterate Bethe roots",
        [&](auto& table) {
          for (std::size_t i = 0; i < m; ++i)
            table.append(i, s.quantum_numbers[i], s.scaled_roots[i],
                         uni20::isfinite(s.rapidities[i]) ? Optional(s.rapidities[i]) : std::nullopt, s.converged);
        },
        column<std::size_t>("root_id"), column<uni20::half_int>("I"), column<Real>("z"), column<Optional>("lambda"),
        column<bool>("converged"));
  output.finish();
  if (!s.converged)
    std::cerr << "Regular Bethe solve incomplete; energies are unavailable and roots are provisional.\n";
  return s.converged ? 0 : 2;
}
} // namespace
int main(int argc, char** argv)
{
  Arguments a;
  return cli::program_main(
      argc, argv, program_info(), [&](auto& app) { add_options(app, a); },
      [&](auto&) {
        a.output.validate();
        return cli::dispatch_precision(a.precision, [&]<typename Real>() { return run<Real>(a, argc, argv); });
      });
}
