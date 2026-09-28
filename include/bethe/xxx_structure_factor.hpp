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
  using std::atan;
  using std::sin;
  if (sites < 2 || momentum_index >= sites || !uni20::isfinite(ground_energy))
    throw std::invalid_argument("invalid XXX first-moment parameters");
  Real const s = sin(Real{4} * atan(Real{1}) * Real(momentum_index) / Real(sites));
  return -Real{8} * ground_energy / (Real{3} * Real(sites)) * s * s;
}

/// Exhaustive conventional two-spinon triplet family, not the full DSF.
/// Always scans and retains all candidates; the budget is checked before work.
template <uni20::Real Real = double>
[[nodiscard]] TwoSpinonStructureFactor<Real> two_spinon_structure_factor(std::size_t sites,
                                                                         SolverOptions<Real> const& solver = {},
                                                                         std::size_t max_candidates = 10000)
{
  detail::checked_sites(sites);
  if (sites % 2) throw std::invalid_argument("XXX two-spinon structure factor requires even N");
  auto const candidates = real_excitation_count(sites, uni20::half_int{1}, max_candidates);
  TwoSpinonStructureFactor<Real> result;
  result.scan = real_excitations<Real>(sites, uni20::half_int{1}, {candidates, max_candidates}, solver);
  auto const& ground = result.scan.ground_state;
  using std::atan;
  Real const two_pi = Real{8} * atan(Real{1});
  for (std::size_t i = 0; i < result.scan.levels.size(); ++i)
  {
    auto const& level = result.scan.levels[i];
    auto factor = raising_form_factor(sites, ground, level.state);
    // Translation T moves the spin at j to j+1 and has eigenvalue exp(-iP).
    // S_q^+=sum exp(-iqj) S_j^+/sqrt(N) therefore selects q=P0-Pn.
    auto const q = (ground.momentum_index + sites - level.state.momentum_index) % sites;
    if (factor.converged() && level.gap && *level.gap >= Real{0})
      result.lines.push_back({i + 1, q, two_pi * Real(q) / Real(sites), *level.gap, *factor.weight});
    else if (factor.converged())
      factor = {}; // inconsistent gap: do not publish even a finite weight
    result.form_factors.push_back(std::move(factor));
  }
  result.moments = spectral_moments<Real>(sites, result.lines);
  bethe::detail::CompensatedSum<Real> weights, first;
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
} // namespace bethe::heisenberg
