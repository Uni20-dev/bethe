// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "citation-report.hpp"
#include "report-common.hpp"
#include <bethe/hubbard_doped.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::hubbard::thermo;
struct Arguments
{
    std::optional<std::string_view> u, momentum, tolerance;
    std::string_view density = "1", reference = "hamiltonian";
    std::string_view branch = "all", convention = "symmetric", precision = "fp64", format = "auto";
    std::size_t points = 33, max_evaluations = 1000000, max_levels = 12, max_iterations = 160;
    std::size_t initial_nodes = 16, max_nodes = 256, max_background_iterations = 64;
    bool quadrature_set = false, mesh_set = false;
    bool points_set = false;
};
void usage(std::ostream& out)
{
  out << "Usage: bethe-hubbard-dispersion --u U [options]\n"
      << "Zero-field Hubbard thermodynamic elementary lines; U>0, t=1.\n"
      << "  --density N/L                    0<n<=1 (default: 1, half filling)\n"
      << "  --branch spinon|holon|antiholon|charge-particle|all (all)\n"
      << "                                     antiholon: n=1 only; charge-particle: n<1 only\n"
      << "  --convention symmetric|unshifted  (default: symmetric, SO(4))\n"
      << "  --reference hamiltonian|fermi     DeltaE or DeltaE-mu*DeltaN (hamiltonian)\n"
      << "  --points COUNT                   uniform dressed-momentum grid, >=2 (33)\n"
      << "  --momentum P                     one dressed momentum in radians instead\n"
      << "  --precision fp64|long-double|fp128 (fp64; fp128 requires MPLAPACK)\n"
      << "  --format auto|pretty|plain|csv|tsv (auto)\n"
      << "  --tolerance VALUE                half-filled relative target (256 epsilon)\n"
      << "  --max-evaluations COUNT           quadrature samples per point (1000000)\n"
      << "  --max-levels COUNT                quadrature refinement levels, <=24 (12)\n"
      << "  --max-iterations COUNT            momentum inversion updates (160)\n"
      << "  --initial-nodes COUNT             doped positive-half quadrature order (16)\n"
      << "  --max-nodes COUNT                 doped mesh limit, <=512 (256)\n"
      << "  --max-background-iterations COUNT doped density solves across all meshes (64)\n"
      << "Doped tolerance: 4096 epsilon; absolute background / max(1,|E|) point energy.\n"
      << "--max-evaluations/--max-levels apply only at half filling.\n"
      << "  --help                           show this help and references\n"
      << "At half filling: spinon p in [0,pi], DeltaN=0, S=1/2; holon/antiholon\n"
      << "p in [-pi,pi], DeltaN=-1/+1, S=0; symmetric lines differ by a pi momentum shift.\n"
      << "H_sym=H_unshifted-U*N/2+U*L/4; excitation E_unshifted=E_sym+U*DeltaN/2.\n"
      << "These are elementary lines, not multiparticle continuum thresholds.\n"
      << "Doped unwrapped p intervals: spinon [0,pi*n], holon [-pi*n/2,3*pi*n/2],\n"
      << "charge-particle [pi*n/2,2*pi-3*pi*n/2]. Particle is a real-root addition,\n"
      << "not the gapped half-filled antiholon. All are gapless in the Fermi reference.\n"
      << "No n>1, attractive U, finite-size levels or spectral weights here.\n"
      << "CSV/TSV have # metadata and empty energies on failure (exit 2).\n"
      << "See docs/hubbard-dispersion.md for iMPS momentum conventions and errors.\n";
  cli::print_citations(out, bethe::citations::Tool::hubbard_dispersion);
}
Arguments parse(int argc, char** argv)
{
  Arguments args;
  for (int i = 1; i < argc; ++i)
  {
    std::string_view const option = argv[i];
    if (option != "--u" && option != "--branch" && option != "--convention" && option != "--points" &&
        option != "--momentum" && option != "--precision" && option != "--format" && option != "--tolerance" &&
        option != "--max-evaluations" && option != "--max-levels" && option != "--max-iterations" &&
        option != "--density" && option != "--reference" && option != "--initial-nodes" && option != "--max-nodes" &&
        option != "--max-background-iterations")
      throw std::invalid_argument("unknown option: " + std::string(option));
    if (++i == argc) throw std::invalid_argument("missing value for " + std::string(option));
    std::string_view const value = argv[i];
    if (option == "--u")
      args.u = value;
    else if (option == "--density")
      args.density = value;
    else if (option == "--reference")
      args.reference = value;
    else if (option == "--branch")
      args.branch = value;
    else if (option == "--convention")
      args.convention = value;
    else if (option == "--momentum")
      args.momentum = value;
    else if (option == "--precision")
      args.precision = value;
    else if (option == "--format")
      args.format = value;
    else if (option == "--tolerance")
      args.tolerance = value;
    else if (option == "--points")
    {
      args.points = cli::parse_size(value);
      args.points_set = true;
    }
    else if (option == "--max-evaluations")
    {
      args.max_evaluations = cli::parse_size(value);
      args.quadrature_set = true;
    }
    else if (option == "--max-levels")
    {
      args.max_levels = cli::parse_size(value);
      args.quadrature_set = true;
    }
    else if (option == "--initial-nodes")
    {
      args.initial_nodes = cli::parse_size(value);
      args.mesh_set = true;
    }
    else if (option == "--max-nodes")
    {
      args.max_nodes = cli::parse_size(value);
      args.mesh_set = true;
    }
    else if (option == "--max-background-iterations")
    {
      args.max_background_iterations = cli::parse_size(value);
      args.mesh_set = true;
    }
    else
      args.max_iterations = cli::parse_size(value);
  }
  if (!args.u) throw std::invalid_argument("--u is required");
  if (args.branch != "spinon" && args.branch != "holon" && args.branch != "antiholon" && args.branch != "all" &&
      args.branch != "charge-particle")
    throw std::invalid_argument("unknown branch: " + std::string(args.branch));
  if (args.convention != "symmetric" && args.convention != "unshifted")
    throw std::invalid_argument("unknown energy convention: " + std::string(args.convention));
  if (args.reference != "hamiltonian" && args.reference != "fermi")
    throw std::invalid_argument("unknown energy reference: " + std::string(args.reference));
  if (args.format != "auto" && args.format != "pretty" && args.format != "plain" && args.format != "csv" &&
      args.format != "tsv")
    throw std::invalid_argument("unknown output format: " + std::string(args.format));
  if (args.momentum && args.points_set) throw std::invalid_argument("--points and --momentum are mutually exclusive");
  if (!args.momentum && args.points < 2) throw std::invalid_argument("--points must be at least 2");
  return args;
}
char const* name(model::Branch b)
{
  switch (b)
  {
    case model::Branch::spinon:
      return "spinon";
    case model::Branch::holon:
      return "holon";
    case model::Branch::antiholon:
      return "antiholon";
  }
  return "unknown";
}
char const* name(model::Status s)
{
  switch (s)
  {
    case model::Status::converged:
      return "converged";
    case model::Status::quadrature_limit:
      return "quadrature_limit";
    case model::Status::momentum_limit:
      return "momentum_limit";
    case model::Status::precision_limit:
      return "precision_limit";
  }
  return "unknown";
}
char const* name(model::DopedBranch b)
{
  switch (b)
  {
    case model::DopedBranch::spinon:
      return "spinon";
    case model::DopedBranch::holon:
      return "holon";
    case model::DopedBranch::charge_particle:
      return "charge-particle";
  }
  return "unknown";
}
char const* name(model::DopedStatus s)
{
  switch (s)
  {
    case model::DopedStatus::converged:
      return "converged";
    case model::DopedStatus::mesh_limit:
      return "mesh_limit";
    case model::DopedStatus::density_limit:
      return "density_limit";
    case model::DopedStatus::linear_failure:
      return "linear_failure";
    case model::DopedStatus::momentum_limit:
      return "momentum_limit";
    case model::DopedStatus::precision_limit:
      return "precision_limit";
  }
  return "unknown";
}
template <uni20::Real Real> int run(Arguments const& args)
{
  Real const u = uni20::parse_real<Real>(*args.u), pi = Real{4} * std::atan(Real{1});
  Real const density = uni20::parse_real<Real>(args.density);
  if (!uni20::isfinite(density) || density <= Real{0} || density > Real{1})
    throw std::invalid_argument("density must be in (0,1]");
  bool const doped = density < Real{1};
  if (doped && args.branch == "antiholon")
    throw std::invalid_argument("doped real-root addition is --branch charge-particle, not the half-filled antiholon");
  if (!doped && args.branch == "charge-particle") throw std::invalid_argument("charge-particle requires density<1");
  if (doped && args.quadrature_set)
    throw std::invalid_argument("doped curves use mesh controls, not --max-evaluations/--max-levels");
  if (!doped && args.mesh_set) throw std::invalid_argument("mesh controls require density<1");
  auto const convention = args.convention == "symmetric" ? model::Convention::symmetric : model::Convention::unshifted;
  auto const reference =
      args.reference == "fermi" ? model::EnergyReference::fermi : model::EnergyReference::hamiltonian;
  model::Options<Real> options;
  options.max_evaluations = args.max_evaluations;
  options.max_levels = args.max_levels;
  options.max_iterations = args.max_iterations;
  if (args.tolerance) options.relative_tolerance = uni20::parse_real<Real>(*args.tolerance);
  std::vector<model::Point<Real>> results;
  std::vector<model::DopedPoint<Real>> doped_results;
  std::optional<model::DopedBackground<Real>> background;
  cli::CpuTimer const timer;
  bool complete = true;
  if (doped)
  {
    model::DopedOptions<Real> controls;
    if (args.tolerance) controls.tolerance = uni20::parse_real<Real>(*args.tolerance);
    controls.initial_nodes = args.initial_nodes;
    controls.max_nodes = args.max_nodes;
    controls.max_background_iterations = args.max_background_iterations;
    controls.max_iterations = args.max_iterations;
    model::DopedSolver<Real> solver(u, density, controls);
    background = solver.background();
    complete = background->converged;
    options.relative_tolerance = controls.tolerance;
    for (auto b : {model::DopedBranch::spinon, model::DopedBranch::holon, model::DopedBranch::charge_particle})
    {
      if (args.branch != "all" && args.branch != name(b)) continue;
      auto const [lo, hi] = solver.momentum_range(b);
      std::size_t const count = args.momentum ? 1 : args.points;
      for (std::size_t i = 0; i < count; ++i)
      {
        Real const p = args.momentum    ? uni20::parse_real<Real>(*args.momentum)
                       : i == count - 1 ? hi
                                        : lo + (hi - lo) * Real(i) / Real(count - 1);
        doped_results.push_back(solver.at_momentum(b, p, convention, reference));
        complete = complete && doped_results.back().converged;
      }
    }
  }
  else
    for (auto b : {model::Branch::spinon, model::Branch::holon, model::Branch::antiholon})
    {
      if (args.branch != "all" && args.branch != name(b)) continue;
      std::size_t const count = args.momentum ? 1 : args.points;
      for (std::size_t i = 0; i < count; ++i)
      {
        Real const fraction = args.momentum ? Real{0} : Real(i) / Real(count - 1);
        Real const p = args.momentum                ? uni20::parse_real<Real>(*args.momentum)
                       : b == model::Branch::spinon ? pi * fraction
                                                    : pi * (Real{2} * fraction - Real{1});
        results.push_back(model::dispersion(b, u, p, convention, options));
        if (reference == model::EnergyReference::fermi && results.back().converged)
          results.back().energy = results.back().symmetric_energy;
        complete = complete && results.back().converged;
      }
    }
  auto const cpu = timer.elapsed_text();
  cli::report_builder report(doped ? "Doped Hubbard thermodynamic dispersions"
                                   : "Half-filled Hubbard thermodynamic dispersions");
  report.field("U (t=1)", uni20::format_real(u))
      .field("Density N/L", uni20::format_real(density))
      .field("Background",
             doped ? "below half filling, zero field, infinite chain" : "half filling, zero field, infinite chain")
      .field("Energy convention", std::string(args.convention))
      .field("Energy reference", std::string(args.reference))
      .field("Energy conversion", "E_unshifted = E_symmetric + U*DeltaN/2; E_fermi = E_H - mu*DeltaN")
      .field("Momentum", doped ? "unwrapped one-site radians; hole offset pi*n/2; particle offset -pi*n/2"
                               : "one-site radians; antiholon p = holon p - pi at the same bare k")
      .field("Precision", std::string(args.precision))
      .field(doped ? "Tolerance" : "Relative tolerance", uni20::format_real(options.relative_tolerance))
      .field("Status", complete ? "converged" : "incomplete; failed energies omitted")
      .field("CPU time", cpu);
  if (background)
  {
    auto value = [](std::optional<Real> const& x) { return x ? uni20::format_real(*x) : "unavailable"; };
    report.field("Background status", name(background->status))
        .field("Fermi rapidity Q", value(background->fermi_rapidity))
        .field("Mu unshifted", value(background->mu_unshifted))
        .field("Mu symmetric", value(background->mu_symmetric))
        .field("Ground energy/site unshifted", value(background->energy_per_site_unshifted))
        .field("Nodes (positive half)", background->nodes)
        .field("Background iterations", background->iterations)
        .field("Density residual", uni20::format_real(background->density_error))
        .field("Background mesh error estimate", uni20::format_real(background->mesh_error));
  }
  else
    report.field("Mu unshifted", uni20::format_real(u / Real{2})).field("Mu symmetric", "0");
  bool const separated = args.format == "csv" || args.format == "tsv";
  char const separator = args.format == "tsv" ? '\t' : ',';
  std::vector<std::string> const columns = {
      "branch", "p",           "p_over_pi",         "energy",         "symmetric_energy", "delta_n",
      "spin",   "parameter",   "energy_quad_error", "momentum_error", "evaluations",      "iterations",
      "status", "fermi_energy"};
  auto headings = columns;
  if (doped) headings[8] = "energy_mesh_error";
  auto& table = report.table("Elementary lines (not continuum thresholds)");
  table.header_separator();
  for (auto const& c : headings)
    table.column(c);
  auto separated_row = [&](auto const& cells) {
    bool first = true;
    for (auto const& cell : cells)
    {
      if (!first) std::cout << separator;
      first = false;
      std::cout << cell;
    }
    std::cout << '\n';
  };
  if (separated)
  {
    for (auto const& [key, value] : report.fields())
      std::cout << "# " << key << ": " << value << '\n';
    separated_row(headings);
  }
  auto row = [&](std::vector<std::string> const& cells) {
    if (separated)
      separated_row(cells);
    else
      table.row(cells[0], cells[1], cells[2], cells[3], cells[4], cells[5], cells[6], cells[7], cells[8], cells[9],
                cells[10], cells[11], cells[12], cells[13]);
  };
  for (auto const& s : results)
  {
    std::vector<std::string> const cells = {name(s.branch),
                                            uni20::format_real(s.momentum),
                                            uni20::format_real(s.momentum / pi),
                                            s.energy ? uni20::format_real(*s.energy) : "",
                                            s.symmetric_energy ? uni20::format_real(*s.symmetric_energy) : "",
                                            std::to_string(s.delta_particles),
                                            uni20::to_string_fraction(s.spin),
                                            s.converged ? uni20::format_real(s.parameter) : "",
                                            s.converged ? uni20::format_real(s.energy_error) : "",
                                            s.converged ? uni20::format_real(s.momentum_error) : "",
                                            std::to_string(s.evaluations),
                                            std::to_string(s.iterations),
                                            name(s.status),
                                            s.symmetric_energy ? uni20::format_real(*s.symmetric_energy) : ""};
    row(cells);
  }
  for (auto const& s : doped_results)
    row({name(s.branch), uni20::format_real(s.momentum), uni20::format_real(s.momentum / pi),
         s.energy ? uni20::format_real(*s.energy) : "",
         s.symmetric_energy ? uni20::format_real(*s.symmetric_energy) : "", std::to_string(s.delta_particles),
         uni20::to_string_fraction(s.spin), s.converged ? uni20::format_real(s.parameter) : "",
         s.converged ? uni20::format_real(s.energy_error) : "", s.converged ? uni20::format_real(s.momentum_error) : "",
         "", std::to_string(s.iterations), name(s.status), s.fermi_energy ? uni20::format_real(*s.fermi_energy) : ""});
  if (!separated) cli::print_report(report, args.format);
  if (!complete) std::cerr << "Some dispersion points did not converge; their energies are unavailable.\n";
  return complete ? 0 : 2;
}
} // namespace

int main(int argc, char** argv)
{
  if (argc == 1)
  {
    usage(std::cerr);
    return 1;
  }
  for (int i = 1; i < argc; ++i)
    if (std::string_view(argv[i]) == "--help")
    {
      usage(std::cout);
      return 0;
    }
  try
  {
    auto const args = parse(argc, argv);
    return cli::dispatch_precision(args.precision, [&]<typename Real>() { return run<Real>(args); });
  }
  catch (std::exception const& error)
  {
    std::cerr << "Error: " << error.what() << '\n';
    return 1;
  }
}
