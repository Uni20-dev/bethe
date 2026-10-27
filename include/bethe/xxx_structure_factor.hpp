// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/dynamical_structure_factor.hpp>
#include <bethe/heisenberg_excitations.hpp>
#include <bethe/xxx_form_factors.hpp>

namespace bethe::heisenberg
{
template <uni20::Real Real> struct TwoSpinonStructureFactor
{
    RealExcitationScan<RealState<Real>> scan;
    /// One diagnostic per retained converged root state; failed form factors
    /// are absent from lines, but retained here (never a silent zero weight).
    std::vector<FormFactor<Real>> form_factors;
    std::vector<LehmannLine<Real>> lines;
    std::vector<SpectralMoment<Real>> moments;
    Real integrated_weight = Real{0}; // sum_q weight/N; full raising result=1/2
    Real weight_fraction = Real{0};
    std::optional<Real> first_moment_fraction;
    bool converged() const { return scan.converged() && lines.size() == scan.candidate_count; }
};

/// Full raising-channel first-frequency sum rule for H=sum S.S, J=1.
/// E0 is the unshifted total ground energy, not E0-N/4 or an energy density.
template <uni20::Real Real> Real raising_first_moment(std::size_t sites, std::size_t momentum_index, Real ground_energy)
{
  using std::sin;
  if (sites < 2 || momentum_index >= sites || !uni20::isfinite(ground_energy))
    throw std::invalid_argument("invalid XXX first-moment parameters");
  Real const s = sin(bethe::detail::pi_like(ground_energy) * Real(momentum_index) / Real(sites));
  return -Real{8} * ground_energy / (Real{3} * Real(sites)) * s * s;
}

namespace detail
{
/// Keep the factor attached to its level through the common scanner's sorting.
/// Enumeration identities and final energy ranks are not interchangeable.
template <uni20::Real Real> struct TwoSpinonLevel
{
    RealState<Real> state;
    std::optional<Real> gap;
    FormFactor<Real> factor;
};

template <uni20::Real Real> struct TwoSpinonScan
{
    RealState<Real> ground_state;
    std::vector<TwoSpinonLevel<Real>> levels;
    std::size_t candidate_count = 0;
    std::size_t converged_count = 0;
    std::optional<RealState<Real>> first_unconverged;
};

template <uni20::Real Real>
TwoSpinonStructureFactor<Real> finish_two_spinon_structure_factor(std::size_t sites, TwoSpinonScan<Real> scan,
                                                                  Real const& zero)
{
  TwoSpinonStructureFactor<Real> result;
  result.scan.spin = uni20::half_int{1};
  result.scan.ground_state = std::move(scan.ground_state);
  result.scan.candidate_count = scan.candidate_count;
  result.scan.converged_count = scan.converged_count;
  result.scan.first_unconverged = std::move(scan.first_unconverged);
  result.scan.levels.reserve(scan.levels.size());
  result.form_factors.reserve(scan.levels.size());
  for (auto& level : scan.levels)
  {
    result.scan.levels.push_back({std::move(level.state), std::move(level.gap)});
    result.form_factors.push_back(std::move(level.factor));
  }
  auto const& ground = result.scan.ground_state;
  Real const two_pi = Real{2} * bethe::detail::pi_like(zero);
  for (std::size_t i = 0; i < result.scan.levels.size(); ++i)
  {
    auto const& level = result.scan.levels[i];
    auto& factor = result.form_factors[i];
    // Translation T moves the spin at j to j+1 and has eigenvalue exp(-iP).
    // S_q^+=sum exp(-iqj) S_j^+/sqrt(N) therefore selects q=P0-Pn.
    auto const q = (ground.momentum_index + sites - level.state.momentum_index) % sites;
    if (factor.converged() && level.gap && *level.gap >= Real{0})
      result.lines.push_back({i + 1, q, two_pi * Real(q) / Real(sites), *level.gap, *factor.weight});
    else if (factor.converged())
      factor = {}; // inconsistent gap: do not publish even a finite weight
  }
  result.moments = spectral_moments<Real>(sites, result.lines, zero);
  bethe::detail::CompensatedSum<Real> weights(zero), first(zero);
  for (auto const& moment : result.moments)
  {
    weights.add(moment.weight);
    first.add(moment.first_moment);
  }
  result.integrated_weight = weights.value() / Real(sites);
  result.weight_fraction = Real{2} * result.integrated_weight;
  if (ground.converged && ground.energy < Real{0})
    result.first_moment_fraction = first.value() / (-Real{4} * ground.energy / Real{3});
  return result;
}

template <uni20::Real Real, typename Solve, typename Ground, typename Prepare>
uni20::async::AsyncTask
two_spinon_structure_factor_task(uni20::async::IAsyncScheduler& scheduler, std::size_t sites,
                                 RealExcitationOptions options, Solve solve, Ground ground, Prepare prepare, Real zero,
                                 uni20::async::WriteBuffer<TwoSpinonStructureFactor<Real>> output)
{
  auto destination = co_await output;
  // Prepare once before starting the runners; this owned snapshot is immutable.
  using Prepared = std::invoke_result_t<Prepare&, RealState<Real> const&>;
  std::optional<Prepared> prepared;
  auto prepare_ground = [&] {
    auto reference = ground();
    if (reference.converged && !options.execution.stop.stop_requested())
      prepared.emplace(prepare(std::as_const(reference)));
    return reference;
  };
  TwoSpinonScan<Real> scan;
  scan.candidate_count = options.count;
  auto const m = sites / 2 - 1, slots = sites - m;
  bethe::detail::StateEnergy energy;
  bethe::detail::scan_real_combinations_into(
      scheduler, slots, m, -static_cast<std::int64_t>(slots - 1), options, solve, prepare_ground,
      bethe::detail::EnergyOrder::ascending, energy, scan,
      [&](RealExcitation<RealState<Real>> level, RealState<Real> const& reference) {
        FormFactor<Real> weight;
        // The common scanner only transforms converged excited roots. A failed
        // ground reference likewise must not enter the numerical factor kernel.
        if (reference.converged)
          weight = std::as_const(*prepared)(level.state);
        else
          weight.status = FormFactorStatus::roots_unconverged;
        return TwoSpinonLevel<Real>{std::move(level.state), std::move(level.gap), std::move(weight)};
      });
  auto result = finish_two_spinon_structure_factor(sites, std::move(scan), zero);
  if (options.execution.stop.stop_requested()) throw ScanCancelled{};
  destination = std::move(result);
}

/// Own callable/configuration values, but borrow the scheduler. Kept separate
/// from the public numerical wrapper so execution failures can be tested without
/// changing the scientific kernels. One submitted task owns the complete chain.
template <uni20::Real Real, typename Solve, typename Ground, typename Prepare>
uni20::async::Async<TwoSpinonStructureFactor<Real>>
submit_two_spinon_structure_factor(std::size_t sites, std::size_t max_candidates, ExecutionOptions execution,
                                   Solve solve, Ground ground, Prepare prepare, Real zero)
{
  auto& scheduler = execution.continuous_scheduler();
  checked_sites(sites);
  if (sites % 2) throw std::invalid_argument("XXX two-spinon structure factor requires even N");
  auto const candidates = real_excitation_count(sites, uni20::half_int{1}, max_candidates);
  execution.scheduler = &scheduler;
  uni20::async::Async<TwoSpinonStructureFactor<Real>> output;
  auto task = two_spinon_structure_factor_task<Real>(scheduler, sites, {candidates, max_candidates, execution},
                                                     std::move(solve), std::move(ground), std::move(prepare),
                                                     std::move(zero), output.write());
  bethe::detail::schedule_owned(scheduler, std::move(task));
  return output;
}
} // namespace detail

/// Exhaustive conventional two-spinon triplet family, not the full DSF.
/// One Async operation owns continuous solve -> factor -> collection runners.
/// Owns all options; the borrowed scheduler must outlive completion, even if the
/// handle is dropped. Cancellation finishes active kernels, then throws
/// ScanCancelled, never publishing an incomplete scan as successful execution.
template <uni20::Real Real = double>
[[nodiscard]] uni20::async::Async<TwoSpinonStructureFactor<Real>>
two_spinon_structure_factor_async(std::size_t sites, SolverOptions<Real> solver = {},
                                  std::size_t max_candidates = 10000, ExecutionOptions execution = {})
{
  auto zero = bethe::detail::solver_zero(solver);
  return detail::submit_two_spinon_structure_factor<Real>(
      sites, max_candidates, execution,
      [sites, solver](QuantumNumbers const& numbers) { return solve_real<Real>(sites, numbers, solver); },
      [sites, solver] { return ground_state<Real>(sites, solver); },
      [sites](auto const& ground) { return PreparedRaisingFormFactors<Real>(sites, ground); }, std::move(zero));
}

/// Synchronous observation of the same operation. Retains all roots and factor
/// diagnostics in energy order. The candidate budget is checked before work.
template <uni20::Real Real = double>
[[nodiscard]] TwoSpinonStructureFactor<Real>
two_spinon_structure_factor(std::size_t sites, SolverOptions<Real> const& solver = {},
                            std::size_t max_candidates = 10000, ExecutionOptions execution = {})
{
  auto& scheduler = execution.resolved_scheduler();
  execution.scheduler = &scheduler;
  auto pending = two_spinon_structure_factor_async<Real>(sites, solver, max_candidates, execution);
  return bethe::detail::take_scan_result(pending, scheduler);
}
} // namespace bethe::heisenberg
