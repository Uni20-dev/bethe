// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "batched_scan_reference.hpp"
#include <bethe/xxx_structure_factor.hpp>

// Frozen public 5909c21 two-stage control flow, for regression tests only.
// Scientific solvers and the prepared form-factor kernel remain shared.
namespace dsf_reference
{
using bethe::ExecutionOptions;
using bethe::spectral_moments;
using namespace bethe::heisenberg;
/// Exhaustive conventional two-spinon triplet family, not the full DSF.
/// Always scans and retains all candidates; the budget is checked before work.
template <uni20::Real Real = double>
[[nodiscard]] TwoSpinonStructureFactor<Real>
two_spinon_structure_factor(std::size_t sites, SolverOptions<Real> const& solver = {},
                            std::size_t max_candidates = 10000, ExecutionOptions execution = {})
{
  auto& scheduler = execution.resolved_scheduler();
  detail::checked_sites(sites);
  if (sites % 2) throw std::invalid_argument("XXX two-spinon structure factor requires even N");
  auto const candidates = real_excitation_count(sites, uni20::half_int{1}, max_candidates);
  TwoSpinonStructureFactor<Real> result;
  // Resolve once: both stages use the same scheduler even when supplied implicitly.
  execution.scheduler = &scheduler;
  auto const m = sites / 2 - 1, slots = sites - m;
  result.scan = scan_reference::scan_real_combinations<RealExcitationScan<RealState<Real>>>(
      slots, m, -static_cast<std::int64_t>(slots - 1), {candidates, max_candidates, execution},
      [&](auto const& labels) { return solve_real<Real>(sites, labels, solver); },
      [&] { return ground_state<Real>(sites, solver); });
  result.scan.spin = uni20::half_int{1};
  auto const& ground = result.scan.ground_state;
  PreparedRaisingFormFactors<Real> const prepared(sites, ground);
  auto const zero = bethe::detail::solver_zero(solver);
  Real const two_pi = Real{2} * bethe::detail::pi_like(zero);
  result.form_factors.resize(result.scan.levels.size());
  for (std::size_t offset = 0; offset < result.scan.levels.size();)
  {
    auto const size = std::min(execution.batch_size, result.scan.levels.size() - offset);
    scheduler.execute_batch(size, [&](std::size_t j) {
      auto const i = offset + j;
      result.form_factors[i] = prepared(result.scan.levels[i].state);
    });
    offset += size;
  }
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
} // namespace dsf_reference
