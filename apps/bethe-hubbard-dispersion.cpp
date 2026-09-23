// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "data-output-options.hpp"
#include "program-options.hpp"
#include <bethe/hubbard_doped.hpp>

namespace
{
namespace cli = bethe::cli;
namespace options = uni20::cli;
namespace model = bethe::hubbard::thermo;
struct Arguments
{
    std::string u;
    std::optional<std::string> momentum, tolerance;
    std::string density = "1", reference = "hamiltonian";
    std::string branch = "all", convention = "symmetric", precision = "fp64";
    cli::DataOutputOptions output;
    std::size_t points = 33, max_evaluations = 1000000, max_levels = 12, max_iterations = 160;
    std::size_t initial_nodes = 16, max_nodes = 256, max_background_iterations = 64;
    bool quadrature_set = false, mesh_set = false;
};
auto program_info()
{
  auto info =
      cli::program_info("bethe-hubbard-dispersion", "Zero-field Hubbard thermodynamic elementary lines; U>0, t=1.",
                        bethe::citations::Tool::hubbard_dispersion);
  info.examples = {{"bethe-hubbard-dispersion --u 4 --branch spinon --points 33", "Half-filled spinon line"},
                   {"bethe-hubbard-dispersion --u 4 --density 0.75 --reference fermi --csv doped.csv",
                    "Doped lines with a screen report and CSV export"},
                   {"bethe-hubbard-dispersion --u=4.000000000000000001 --precision=long-double --format=json",
                    "Native-precision tokens and machine-readable stdout"}};
  info.notes = {
      "At half filling: spinon p in [0,pi], DeltaN=0, S=1/2; holon/antiholon p in [-pi,pi], DeltaN=-1/+1, S=0. "
      "Symmetric charge lines differ by a pi momentum shift.",
      "H_sym=H_unshifted-U*N/2+U*L/4; excitation E_unshifted=E_sym+U*DeltaN/2.",
      "These are elementary lines, not multiparticle continuum thresholds. No n>1, attractive U, finite-size levels "
      "or spectral weights here.",
      "Doped unwrapped p intervals: spinon [0,pi*n], holon [-pi*n/2,3*pi*n/2], charge-particle "
      "[pi*n/2,2*pi-3*pi*n/2]. The particle is a real-root addition, not the gapped half-filled antiholon. "
      "All doped lines are gapless in the Fermi reference.",
      "Default tolerance: 256 epsilon at half filling (relative target); 4096 epsilon when doped "
      "(absolute background / max(1,|E|) point energy), evaluated in the selected precision.",
      "fp128 requires a build with MPLAPACK. Files and machine stdout stream rows; human stdout is a final report "
      "unless --stream is set. --no-retain requires live output or --quiet.",
      "CSV/TSV have # metadata and empty energies on failure (exit 2). JSON uses null. "
      "See docs/output.md for export and overwrite rules.",
      "See docs/hubbard-dispersion.md for iMPS momentum conventions and errors. "
      "Use --references for literature and applicability; see CITATIONS.md for conventions and provenance."};
  return info;
}
void add_options(CLI::App& app, Arguments& args)
{
  auto* physics = app.add_option_group("Model");
  physics->add_option("--u", args.u, "Repulsive interaction U>0, in hopping units")->required()->type_name("REAL");
  physics->add_option("--density", args.density, "Particles per site, 0<n<=1; 1 is half filling")
      ->type_name("REAL")
      ->capture_default_str();
  physics->add_option("--branch", args.branch, "antiholon: n=1 only; charge-particle: n<1 only")
      ->check(CLI::IsMember({"spinon", "holon", "antiholon", "charge-particle", "all"}))
      ->capture_default_str();
  physics->add_option("--convention", args.convention, "Interaction convention; symmetric is SO(4)")
      ->check(CLI::IsMember({"symmetric", "unshifted"}))
      ->capture_default_str();
  physics->add_option("--reference", args.reference, "Hamiltonian DeltaE or Fermi DeltaE-mu*DeltaN")
      ->check(CLI::IsMember({"hamiltonian", "fermi"}))
      ->capture_default_str();
  auto* sampling = app.add_option_group("Momentum sampling");
  auto* points = options::add_count_option(*sampling, "--points", args.points, "Uniform dressed-momentum grid, >=2")
                     ->capture_default_str();
  cli::text_option(*sampling, "--momentum", args.momentum, "One dressed momentum in radians instead of a grid")
      ->type_name("REAL")
      ->excludes(points);
  auto* numerics = app.add_option_group("Numerics");
  numerics->add_option("--precision", args.precision, "Real scalar type; fp128 requires MPLAPACK")
      ->check(CLI::IsMember({"fp64", "long-double", "fp128"}))
      ->capture_default_str();
  cli::text_option(*numerics, "--tolerance", args.tolerance, "Target in native precision; see conventions below")
      ->type_name("REAL")
      ->default_str("256 epsilon (half-filled); 4096 epsilon (doped)");
  options::add_count_option(*numerics, "--max-iterations", args.max_iterations, "Momentum inversion updates")
      ->capture_default_str();
  auto* quadrature = app.add_option_group("Half-filled quadrature", "These controls apply only at density=1.");
  options::add_count_option(*quadrature, "--max-evaluations", args.max_evaluations, "Quadrature samples per point")
      ->capture_default_str();
  options::add_count_option(*quadrature, "--max-levels", args.max_levels, "Refinement levels, <=24")
      ->capture_default_str();
  auto* mesh = app.add_option_group("Doped mesh", "These controls require density<1.");
  options::add_count_option(*mesh, "--initial-nodes", args.initial_nodes, "Positive-half quadrature order")
      ->capture_default_str();
  options::add_count_option(*mesh, "--max-nodes", args.max_nodes, "Mesh limit, <=512")->capture_default_str();
  options::add_count_option(*mesh, "--max-background-iterations", args.max_background_iterations,
                            "Density solves across all meshes")
      ->capture_default_str();
  cli::add_data_output_options(app, args.output);
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
  Real const u = uni20::parse_real<Real>(args.u), pi = Real{4} * std::atan(Real{1});
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
  Arguments args;
  return cli::program_main(
      argc, argv, program_info(), [&](auto& app) { add_options(app, args); },
      [&](auto& app) {
        args.quadrature_set = app.count("--max-evaluations") || app.count("--max-levels");
        args.mesh_set =
            app.count("--initial-nodes") || app.count("--max-nodes") || app.count("--max-background-iterations");
        args.output.validate();
        if (!args.momentum && args.points < 2) throw std::invalid_argument("--points must be at least 2");
        return cli::dispatch_precision(args.precision, [&]<typename Real>() { return run<Real>(args, argc, argv); });
      });
}
