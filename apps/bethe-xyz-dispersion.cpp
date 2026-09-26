// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "result-output.hpp"
#include <bethe/xyz.hpp>
#include <bethe/xyz_dispersion.hpp>

namespace
{
namespace cli = bethe::cli;
struct Arguments
{
    std::string eta, t, exchange = "1", branch = "all", precision = "fp64";
    std::optional<std::string> momentum;
    std::optional<std::size_t> bound_state;
    std::size_t points = 65;
    cli::DataOutputOptions output;
};
auto program_info()
{
  auto info = cli::program_info("bethe-xyz-dispersion",
                                "XYZ thermodynamic spinons, bound branches and two-spinon continuum envelopes.",
                                bethe::citations::Tool::xyz_dispersion);
  info.examples = {
      {"bethe-xyz-dispersion --eta 0.4 --t 1 --csv xyz.csv", "Repulsive-region spinons and continuum"},
      {"bethe-xyz-dispersion --eta 0.75 --t 1 --branch bound --json bound.json", "Both allowed bound branches"},
      {"bethe-xyz-dispersion --eta 0.875 --t 0.75 --bound-state 2 --precision long-double",
       "Select a bound branch alongside the spinons and continuum"}};
  info.notes = {
      "H=J sum(Jx SxSx+Jy SySy+Jz SzSz), S=sigma/2, J>0, infinite chain, zero field and temperature. "
      "Jx,Jy,Jz are the theta ratios used by bethe-xyz-pbc; 0<eta<1, t>0. Sz is not conserved.",
      "Spinon: topological wall between x-Neel vacua, p in [0,pi]. Two-spinon rows give the union of "
      "pi-shifted continuum copies, without resolving their symmetry sectors or spectral weights.",
      "Bound s>=1 exists strictly when s*(1-eta)<eta; a merger is excluded. delta_rx=(-1)^s; both "
      "delta_rz=+1,-1 copies are output. reduced_q=p-pi*(delta_rz==-1) modulo 2*pi. "
      "p is physical momentum relative to one reference vacuum; cell_momentum=2*p modulo 2*pi.",
      "All existing bound branches are included by default (enumeration limit 4096); --bound-state selects one. "
      "The total output is limited to 1000000 rows. --momentum is in radians/site, not units of pi.",
      "These are thermodynamic lines, not exact finite-ring levels or spectral intensities. "
      "Numerical failures leave values empty (JSON null) and exit 2; an unrepresentable gap is not reported as zero.",
      "See docs/xyz-dispersion.md for the Hamiltonian map, momenta, bound-state thresholds and validation."};
  return info;
}
void add_options(CLI::App& app, Arguments& a)
{
  app.add_option("--eta", a.eta, "Elliptic parameter, 0<eta<1")->required()->type_name("REAL");
  app.add_option("--t", a.t, "Rectangular theta parameter tau=i*t, t>0")->required()->type_name("REAL");
  app.add_option("--exchange", a.exchange, "Positive overall exchange J")->capture_default_str()->type_name("REAL");
  app.add_option("--branch", a.branch, "Output branch selection")
      ->check(CLI::IsMember({"spinon", "two-spinon", "bound", "all"}))
      ->capture_default_str();
  cli::count_option(app, "--bound-state", a.bound_state, "Select one existing bound index s (default: all)");
  auto* points =
      cli::count_option(app, "--points", a.points, "Samples per branch/copy, 2..1000000")->capture_default_str();
  cli::text_option(app, "--momentum", a.momentum, "One physical momentum instead of a grid")
      ->type_name("REAL")
      ->excludes(points);
  cli::precision_option(app, a.precision);
  cli::add_data_output_options(app, a.output, true);
}
template <uni20::Real Real> int run(Arguments const& a, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  Real const eta = uni20::parse_real<Real>(a.eta), t = uni20::parse_real<Real>(a.t);
  Real const exchange = uni20::parse_real<Real>(a.exchange), pi = bethe::detail::pi<Real>();
  if (!uni20::isfinite(eta) || eta <= Real{0} || eta >= Real{1} || !uni20::isfinite(t) || t <= Real{0} ||
      !uni20::isfinite(exchange) || exchange <= Real{0})
    throw std::invalid_argument("require finite 0<eta<1, t>0 and exchange>0");
  bool const spinons = a.branch == "all" || a.branch == "spinon";
  bool const continuum = a.branch == "all" || a.branch == "two-spinon";
  bool const bounds = a.branch == "all" || a.branch == "bound";
  if (a.bound_state && (!bounds || *a.bound_state == 0 || *a.bound_state > std::numeric_limits<unsigned>::max()))
    throw std::invalid_argument("--bound-state requires branch all/bound and a positive unsigned index");
  std::optional<Real> momentum;
  if (a.momentum) momentum = uni20::parse_real<Real>(*a.momentum);
  if (momentum && (!uni20::isfinite(*momentum) || *momentum < Real{0} || *momentum > (spinons ? pi : Real{2} * pi)))
    throw std::invalid_argument("momentum outside the selected branch range");
  std::vector<unsigned> bound_indices;
  auto exists = [&](unsigned s) { return Real(s) * (Real{1} - eta) < eta; };
  if (bounds)
  {
    if (a.bound_state)
    {
      auto const s = static_cast<unsigned>(*a.bound_state);
      if (!exists(s)) throw std::invalid_argument("requested bound branch does not exist: require s*(1-eta)<eta");
      bound_indices.push_back(s);
    }
    else
      for (unsigned s = 1; exists(s); ++s)
      {
        if (s > 4096) throw std::invalid_argument("more than 4096 bound branches; select one with --bound-state");
        bound_indices.push_back(s);
      }
    if (a.branch == "bound" && bound_indices.empty())
      throw std::invalid_argument("no isolated bound branch exists at this eta");
  }
  std::size_t const count = momentum ? 1 : a.points;
  std::size_t const families = std::size_t(spinons) + std::size_t(continuum) + 2 * bound_indices.size();
  if (families > 1000000 / count)
    throw std::invalid_argument("output exceeds 1000000 rows; select fewer branches/points");

  std::optional<bethe::xyz::SpinonDispersion<Real>> band;
  std::optional<bethe::xyz::Couplings<Real>> couplings;
  std::string status = "converged", failure;
  try
  {
    context.measure([&] {
      couplings = bethe::xyz::couplings(eta, t);
      band.emplace(eta, t, exchange);
    });
  }
  catch (std::runtime_error const& e)
  {
    status = "precision_limit";
    failure = e.what();
  }
  cli::RunReport report(context, "XYZ thermodynamic excitations");
  report.field("eta", "Eta", eta)
      .field("t", "Theta t", t)
      .field("exchange", "Exchange J", exchange)
      .field("hamiltonian", "Hamiltonian", "H=J sum(Jx SxSx+Jy SySy+Jz SzSz); S=sigma/2; theta ratios as xyz-pbc")
      .field("energy_reference", "Energy reference", "excitation energy above the infinite-chain ground state")
      .field("spinon_sector", "Spinon sector", "domain wall between different x-Neel vacua; no conserved Sz")
      .field("continuum_sectors", "Two-spinon sectors",
             "envelope over pi-shifted translation copies; no symmetry resolution")
      .field("bound_sector", "Bound sectors", "delta_rx=(-1)^s; delta_rz=+1,-1 relative to one vacuum")
      .field("bound_momentum", "Bound momentum", "reduced_q=p-pi*(delta_rz==-1) modulo 2*pi")
      .field("cell_momentum", "Two-site translation phase", "2*p modulo 2*pi; radians/cell")
      .field("bound_existence", "Bound existence", "s>=1 and s*(1-eta)<eta; mergers excluded")
      .field("branches", "Branches", a.branch)
      .field("precision", "Precision", a.precision)
      .field("bound_count", "Selected bound branches", bound_indices.size())
      .field("band_status", "Spinon status", status);
  if (couplings) report.field("jx", "Jx", couplings->x).field("jy", "Jy", couplings->y).field("jz", "Jz", couplings->z);
  if (band)
    report.field("spinon_gap", "Single-spinon gap", band->gap())
        .field("band_maximum", "Spinon maximum energy", band->maximum_energy());
  if (!failure.empty()) report.field("failure", "Numerical failure", failure);
  if (a.bound_state) report.field("bound_state", "Requested bound index", *a.bound_state);
  if (momentum)
    report.field("momentum", "Requested momentum", *momentum);
  else
    report.field("points", "Points per branch/copy", a.points);

  namespace data = cli::data;
  using Optional = std::optional<Real>;
  using Index = std::optional<unsigned>;
  using Parity = std::optional<int>;
  auto table = data::make_data_table(
      "XYZ spinon and bound lines, two-spinon envelopes",
      {.retain = a.output.retain ? data::retention::all : data::retention::none, .metadata = report.metadata()},
      cli::column<std::string>("branch"), cli::column<Index>("s"), cli::column<Parity>("delta_rx"),
      cli::column<Parity>("delta_rz"), cli::column<Real>("p").unit("radians/site"), cli::column<Real>("p_over_pi"),
      cli::column<Real>("cell_momentum").unit("radians/cell"), cli::column<Optional>("reduced_q"),
      cli::column<Optional>("energy"), cli::column<Optional>("lower"), cli::column<Optional>("upper"),
      cli::column<std::string>("status"));
  cli::DataOutput output(a.output, {"dispersion"});
  bool complete = band.has_value();
  try
  {
    output.attach(table, "dispersion");
    auto emit = [&](std::string const& branch, Index s = {}, Parity rz = {}) {
      Real const end = branch == "spinon" ? pi : Real{2} * pi;
      for (std::size_t i = 0; i < count; ++i)
      {
        Real const p = momentum ? *momentum : i == count - 1 ? end : end * (Real(i) / Real(count - 1));
        Real const q = rz && *rz == -1 ? (p < pi ? p + pi : p - pi) : p;
        Optional energy, lower, upper;
        std::string row_status = status;
        if (band) try
          {
            context.measure([&] {
              if (s)
                energy = band->bound_energy(*s, q);
              else if (branch == "spinon")
                energy = band->energy(p);
              else
              {
                auto const edges = band->continuum(p, bethe::SpinonMomentum::folded);
                lower = edges.lower;
                upper = edges.upper;
              }
            });
          }
          catch (std::runtime_error const&)
          {
            row_status = "precision_limit";
            complete = false;
          }
        Real const cell = Real{2} * (p < pi ? p : p < Real{2} * pi ? p - pi : Real{0});
        table.append(branch, s, s ? Parity(*s % 2 ? -1 : 1) : Parity{}, rz, p, p / pi, cell,
                     s ? Optional(q) : Optional{}, energy, lower, upper, row_status);
      }
    };
    if (spinons) emit("spinon");
    if (continuum) emit("two-spinon");
    for (auto s : bound_indices)
      for (int rz : {1, -1})
        emit("bound", s, rz);
    report.result(complete, complete ? "converged" : "incomplete; unavailable quantities omitted");
    auto const summary = report.finish();
    output.overview(report.overview(summary));
    output.finish(table, summary);
    output.finish_document();
  }
  catch (...)
  {
    data::table_metadata aborted{{"Status", "aborted"}};
    if (!context.finished()) try
      {
        aborted = cli::run_summary(context.finish(uni20::run_outcome::failed));
      }
      catch (...)
      {}
    output.abort(table, std::move(aborted));
    throw;
  }
  if (!complete) std::cerr << "Some XYZ excitation quantities are unavailable; inspect row statuses.\n";
  return complete ? 0 : 2;
}
} // namespace
int main(int argc, char** argv)
{
  Arguments a;
  return cli::program_main(
      argc, argv, program_info(), [&](auto& app) { add_options(app, a); },
      [&](auto&) {
        a.output.validate();
        if (!a.momentum && (a.points < 2 || a.points > 1000000))
          throw std::invalid_argument("require 2<=points<=1000000");
        return cli::dispatch_precision(a.precision, [&]<typename Real>() { return run<Real>(a, argc, argv); });
      });
}
