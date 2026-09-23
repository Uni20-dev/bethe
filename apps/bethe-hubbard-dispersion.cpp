// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "citation-report.hpp"
#include "data-output.hpp"
#include <bethe/hubbard_doped.hpp>

namespace
{
namespace cli = bethe::cli;
namespace model = bethe::hubbard::thermo;
struct Arguments
{
    std::optional<std::string_view> u, momentum, tolerance;
    std::string_view density = "1", reference = "hamiltonian";
    std::string_view branch = "all", convention = "symmetric", precision = "fp64";
    cli::DataOutputOptions output;
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
      << "CSV/TSV have # metadata and empty energies on failure (exit 2). JSON uses null.\n"
      << "See docs/hubbard-dispersion.md for iMPS momentum conventions and errors.\n";
  cli::data_output_usage(out);
  cli::print_citations(out, bethe::citations::Tool::hubbard_dispersion);
}
Arguments parse(int argc, char** argv)
{
  Arguments args;
  for (int i = 1; i < argc; ++i)
  {
    if (args.output.parse(i, argc, argv)) continue;
    std::string_view const option = argv[i];
    if (option != "--u" && option != "--branch" && option != "--convention" && option != "--points" &&
        option != "--momentum" && option != "--precision" && option != "--tolerance" && option != "--max-evaluations" &&
        option != "--max-levels" && option != "--max-iterations" && option != "--density" && option != "--reference" &&
        option != "--initial-nodes" && option != "--max-nodes" && option != "--max-background-iterations")
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
  args.output.validate();
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
template <uni20::Real Real> int run(Arguments const& args, int argc, char** argv)
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
  cli::ComputeCpuTime timer;
  std::optional<model::DopedSolver<Real>> solver;
  model::DopedOptions<Real> controls;
  bool complete = true;
  if (doped)
  {
    if (args.tolerance) controls.tolerance = uni20::parse_real<Real>(*args.tolerance);
    controls.initial_nodes = args.initial_nodes;
    controls.max_nodes = args.max_nodes;
    controls.max_background_iterations = args.max_background_iterations;
    controls.max_iterations = args.max_iterations;
    timer.measure([&] { solver.emplace(u, density, controls); });
    complete = solver->background().converged;
    options.relative_tolerance = controls.tolerance;
  }
  else
  {
    // Validate before opening any output file, even when all requested points are endpoints.
    if (!uni20::isfinite(u) || u <= Real{0}) throw std::invalid_argument("Hubbard dispersions require finite U>0");
    if (!uni20::isfinite(options.relative_tolerance) || options.relative_tolerance <= Real{0} ||
        options.relative_tolerance >= Real{1})
      throw std::invalid_argument("relative tolerance must be finite and in (0,1)");
    if (options.max_levels > 24) throw std::invalid_argument("quadrature levels must not exceed 24");
  }
  std::optional<Real> momentum;
  if (args.momentum) momentum = uni20::parse_real<Real>(*args.momentum);
  auto check_momentum = [&](Real lo, Real hi) {
    if (momentum && (!uni20::isfinite(*momentum) || *momentum < lo || *momentum > hi))
      throw std::invalid_argument("momentum is outside a selected branch's interval");
  };
  if (doped)
    for (auto b : {model::DopedBranch::spinon, model::DopedBranch::holon, model::DopedBranch::charge_particle})
    {
      if (args.branch != "all" && args.branch != name(b)) continue;
      auto const [lo, hi] = solver->momentum_range(b);
      check_momentum(lo, hi);
    }
  else
    for (auto b : {model::Branch::spinon, model::Branch::holon, model::Branch::antiholon})
      if (args.branch == "all" || args.branch == name(b))
        check_momentum(b == model::Branch::spinon ? Real{0} : -pi, pi);

  auto metadata = cli::provenance("bethe-hubbard-dispersion", argc, argv);
  auto field = [&](std::string key, std::string value) { metadata.emplace(std::move(key), std::move(value)); };
  field("U (t=1)", uni20::format_real(u));
  field("Density N/L", uni20::format_real(density));
  field("Background",
        doped ? "below half filling, zero field, infinite chain" : "half filling, zero field, infinite chain");
  field("Spectrum", "elementary lines, not multiparticle continuum thresholds");
  field("Hamiltonian", args.convention == "symmetric"
                           ? "-sum(c^dagger_i,s c_i+1,s + h.c.) + U*sum((n_i,up-1/2)*(n_i,down-1/2))"
                           : "-sum(c^dagger_i,s c_i+1,s + h.c.) + U*sum(n_i,up*n_i,down)");
  field("Energy convention", std::string(args.convention));
  field("Energy reference", std::string(args.reference));
  field("Energy conversion", "E_unshifted = E_symmetric + U*DeltaN/2; E_fermi = E_H - mu*DeltaN");
  field("Momentum", doped ? "unwrapped one-site radians; hole offset pi*n/2; particle offset -pi*n/2"
                          : "one-site radians; antiholon p = holon p - pi at the same bare k");
  field("Precision", std::string(args.precision));
  field("Branch selection", std::string(args.branch));
  field(momentum ? "Requested momentum" : "Points per branch",
        momentum ? uni20::format_real(*momentum) : std::to_string(args.points));
  field(doped ? "Tolerance" : "Relative tolerance", uni20::format_real(options.relative_tolerance));
  field("Max iterations", std::to_string(args.max_iterations));
  if (doped)
  {
    auto const& background = solver->background();
    auto value = [](std::optional<Real> const& x) { return x ? uni20::format_real(*x) : "unavailable"; };
    field("Initial nodes", std::to_string(controls.initial_nodes));
    field("Max nodes", std::to_string(controls.max_nodes));
    field("Max background iterations", std::to_string(controls.max_background_iterations));
    field("Background status", name(background.status));
    field("Fermi rapidity Q", value(background.fermi_rapidity));
    field("Mu unshifted", value(background.mu_unshifted));
    field("Mu symmetric", value(background.mu_symmetric));
    field("Ground energy/site unshifted", value(background.energy_per_site_unshifted));
    field("Nodes (positive half)", std::to_string(background.nodes));
    field("Background iterations", std::to_string(background.iterations));
    field("Density residual", uni20::format_real(background.density_error));
    field("Background mesh error estimate", uni20::format_real(background.mesh_error));
  }
  else
  {
    field("Mu unshifted", uni20::format_real(u / Real{2}));
    field("Mu symmetric", "0");
    field("Max evaluations", std::to_string(options.max_evaluations));
    field("Max levels", std::to_string(options.max_levels));
  }
  namespace data = cli::data;
  using Optional = std::optional<Real>;
  cli::DataOutput output(args.output);
  auto table = data::make_data_table(
      doped ? "Doped Hubbard thermodynamic dispersions" : "Half-filled Hubbard thermodynamic dispersions",
      {.retain = args.output.retain ? data::retention::all : data::retention::none, .metadata = std::move(metadata)},
      data::data_column<std::string>("branch"), data::data_column<Real>("p").unit("radians").round_trip(),
      data::data_column<Real>("p_over_pi").round_trip(),
      data::data_column<Optional>("energy")
          .unit("t")
          .description("Selected convention and reference; missing on failure")
          .round_trip(),
      data::data_column<Optional>("symmetric_energy")
          .unit("t")
          .description("Symmetric Hamiltonian energy difference")
          .round_trip(),
      data::data_column<int>("delta_n"), data::data_column<uni20::half_int>("spin"),
      data::data_column<Optional>("parameter")
          .description("Spin rapidity or bare charge momentum; infinity at spinon endpoints")
          .round_trip(),
      data::data_column<Optional>(doped ? "energy_mesh_error" : "energy_quad_error").unit("t").round_trip(),
      data::data_column<Optional>("momentum_error").unit("radians").round_trip(),
      data::data_column<std::optional<std::size_t>>("evaluations"), data::data_column<std::size_t>("iterations"),
      data::data_column<std::string>("status"),
      data::data_column<Optional>("fermi_energy")
          .unit("t")
          .description("DeltaE - mu*DeltaN; independent of Hamiltonian convention")
          .round_trip());
  auto summary = [&](std::string status) -> data::table_metadata {
    return {{"Status", std::move(status)}, {"CPU time", timer.text()}, {"Rows", std::to_string(table.size())}};
  };
  auto append = [&](auto const& s, std::optional<std::size_t> evaluations, Optional fermi_energy) {
    complete = complete && s.converged;
    table.append(std::string(name(s.branch)), s.momentum, s.momentum / pi, s.energy, s.symmetric_energy,
                 s.delta_particles, s.spin, s.converged ? Optional(s.parameter) : std::nullopt,
                 s.converged ? Optional(s.energy_error) : std::nullopt,
                 s.converged ? Optional(s.momentum_error) : std::nullopt, evaluations, s.iterations,
                 std::string(name(s.status)), fermi_energy);
  };
  try
  {
    output.attach(table);
    std::size_t const count = momentum ? 1 : args.points;
    if (doped)
      for (auto b : {model::DopedBranch::spinon, model::DopedBranch::holon, model::DopedBranch::charge_particle})
      {
        if (args.branch != "all" && args.branch != name(b)) continue;
        auto const [lo, hi] = solver->momentum_range(b);
        for (std::size_t i = 0; i < count; ++i)
        {
          Real const p = momentum ? *momentum : i == count - 1 ? hi : lo + (hi - lo) * Real(i) / Real(count - 1);
          auto const s = timer.measure([&] { return solver->at_momentum(b, p, convention, reference); });
          append(s, std::nullopt, s.fermi_energy);
        }
      }
    else
      for (auto b : {model::Branch::spinon, model::Branch::holon, model::Branch::antiholon})
      {
        if (args.branch != "all" && args.branch != name(b)) continue;
        for (std::size_t i = 0; i < count; ++i)
        {
          Real const fraction = momentum ? Real{0} : Real(i) / Real(count - 1);
          Real const p = momentum                     ? *momentum
                         : b == model::Branch::spinon ? pi * fraction
                                                      : pi * (Real{2} * fraction - Real{1});
          auto s = timer.measure([&] { return model::dispersion(b, u, p, convention, options); });
          if (reference == model::EnergyReference::fermi && s.converged) s.energy = s.symmetric_energy;
          append(s, s.evaluations, s.symmetric_energy);
        }
      }
    output.finish(table, summary(complete ? "converged" : "incomplete; failed energies omitted"));
  }
  catch (...)
  {
    output.abort(table, summary("aborted"));
    throw;
  }
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
    return cli::dispatch_precision(args.precision, [&]<typename Real>() { return run<Real>(args, argc, argv); });
  }
  catch (std::exception const& error)
  {
    cli::print_output_error(std::cerr, error);
    return 1;
  }
}
