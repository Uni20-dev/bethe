// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/xxz_helix_check.hpp>
#include <bethe/xxz_odd_continuation.hpp>
#include <bethe/xxz_phantom_check.hpp>
#include <bethe/xxz_regularity.hpp>
#include <bethe/xxz_wronskian.hpp>
#include <optional>

namespace bethe::xxz::detail
{
template <uni20::Real Real> struct OddSectorCandidate
{
    uni20::half_int sz;
    OddPolynomialBranch<Real> branch;
    // Absent if continuation failed. Consistency is separate from equation
    // convergence, and neither flag establishes sector minimality.
    std::optional<WronskianCheck<Real>> wronskian;
    std::optional<PolynomialRegularity<Real>> regularity;
    // Tried only when regularity is unresolved; never changes the branch.
    std::optional<PolynomialHelixCheck<Real>> helix;
    // Tried after regularity and helix checks. Keep unsuccessful attempts;
    // absence of a witness never proves the vector is zero.
    std::vector<PhantomLiftCheck<Real>> phantom_lifts;
    bool phantom_work_limited = false;
};

/// All spin-reversal-distinct odd-ring continuation branches, in increasing
/// M (decreasing positive Sz). This is an INTERNAL candidate scan, not a
/// public ground-state result or a physical-state certificate.
template <uni20::Real Real> struct OddSectorScan
{
    std::vector<OddSectorCandidate<Real>> sectors;
    bool equations_complete = false;
    bool wronskians_consistent = false;
    bool regular_states_complete = false;
    // Every converged branch passes the regular test, explicit helix match,
    // or resolved numerical phantom witness. Not a rigorous state proof or
    // a proof of sector/global minimality.
    bool state_checks_complete = false;
    std::size_t iterations = 0;
    // Populated only when EVERY sector's equations converged and its energy
    // is finite. Failed sectors may hide the true minimum: do not skip them.
    std::optional<std::size_t> lowest_index;
    // Indices numerically close to the minimum; NOT a certified degeneracy.
    std::vector<std::size_t> nearby_indices;
    Real comparison_band = Real{0};
};

/// Compare every M=0,...,floor(N/2), without assuming the minimum lies at
/// smallest |Sz| or hard-coding a finite-size phase boundary. The Newton
/// budget applies PER SECTOR, matching the other sector-scan APIs. The
/// summed count includes rejected stages. Wronskian diagnostics never
/// silently veto, approve, or relabel the candidate energies.
template <uni20::Real Real>
OddSectorScan<Real> scan_odd_polynomial_sectors(std::size_t sites, Real delta, SolverOptions<Real> const& options = {},
                                                Real wronskian_tolerance = Real{256} *
                                                                           uni20::numeric_limits<Real>::epsilon(),
                                                PhantomLiftOptions<Real> phantom_options = {})
{
  auto const n = checked_sites(sites);
  if (sites % 2 == 0) throw std::invalid_argument("odd XXZ sector scan requires odd N");
  if (!uni20::isfinite(wronskian_tolerance) || wronskian_tolerance <= Real{0})
    throw std::invalid_argument("Wronskian tolerance must be finite and positive");
  for (Real tolerance : {phantom_options.factor_tolerance, phantom_options.residual_tolerance,
                         phantom_options.amplitude_tolerance, phantom_options.root_options.tolerance})
    if (!uni20::isfinite(tolerance) || tolerance <= Real{0})
      throw std::invalid_argument("phantom lift tolerances must be finite and positive");
  auto const count = sites / 2 + 1;
  if (count > std::size_t(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(OddSectorCandidate<Real>))
    throw std::length_error("odd XXZ sector scan is too large");
  // Also validate Delta and Newton options before allocating the scan.
  auto vacuum = continue_odd_polynomial(sites, delta, uni20::from_twice(n), options);
  OddSectorScan<Real> out;
  out.sectors.reserve(count);
  out.equations_complete = true;
  out.wronskians_consistent = true;
  out.regular_states_complete = true;
  out.state_checks_complete = true;
  for (std::size_t m = 0; m < count; ++m)
  {
    OddSectorCandidate<Real> sector;
    sector.sz = uni20::from_twice(n - 2 * std::int64_t(m));
    sector.branch = m == 0 ? std::move(vacuum) : continue_odd_polynomial(sites, delta, sector.sz, options);
    auto const& branch = sector.branch;
    if (branch.iterations > std::numeric_limits<std::size_t>::max() - out.iterations)
      throw std::overflow_error("odd XXZ scan iteration total overflowed");
    out.iterations += branch.iterations;
    if (branch.equations_converged && uni20::isfinite(branch.energy))
    {
      sector.wronskian = check_odd_wronskian<Real>(sites, branch.coefficients, delta, branch.center,
                                                   branch.coordinate_scale, wronskian_tolerance);
      out.wronskians_consistent &= sector.wronskian->status == WronskianStatus::consistent;
      PolynomialBetheSystem<Real> const system(sites, m, branch.center, branch.coordinate_scale);
      sector.regularity =
          check_regular_polynomial<Real>(system, branch.coefficients, delta, options.residual_tolerance);
      out.regular_states_complete &= sector.regularity->status == RegularityStatus::regular_on_shell;
      bool resolved = sector.regularity->status == RegularityStatus::regular_on_shell;
      if (!resolved)
      {
        sector.helix = check_helix_polynomial<Real>(system, branch.coefficients, delta, (sites + 1) / 2);
        resolved = sector.helix->status == HelixMatchStatus::compatible;
      }
      if (!resolved && m > 1)
      {
        // Share amplitude budgets among all (p,chirality) hypotheses for
        // this sector. Root-iteration limits still apply per hypothesis.
        auto remaining = phantom_options;
        for (std::size_t p = 1; p < m && !resolved; ++p)
          for (int chirality : {-1, 1})
          {
            if (!remaining.max_configurations || !remaining.max_subset_updates)
            {
              sector.phantom_work_limited = true;
              break;
            }
            auto check = check_phantom_lift<Real>(system, branch.coefficients, delta, p, chirality, remaining);
            remaining.max_configurations -= check.configurations_tested;
            remaining.max_subset_updates -= check.subset_updates;
            sector.phantom_work_limited |= check.status == PhantomLiftStatus::work_limit;
            resolved = check.status == PhantomLiftStatus::nonzero_witness;
            sector.phantom_lifts.push_back(std::move(check));
            if (resolved) break;
          }
      }
      out.state_checks_complete &= resolved;
    }
    else
    {
      out.equations_complete = false;
      out.wronskians_consistent = false;
      out.regular_states_complete = false;
      out.state_checks_complete = false;
    }
    out.sectors.push_back(std::move(sector));
  }
  if (!out.equations_complete) return out;
  std::size_t lowest = 0;
  for (std::size_t i = 1; i < count; ++i)
    if (out.sectors[i].branch.energy < out.sectors[lowest].branch.energy) lowest = i;
  out.lowest_index = lowest;
  Real const minimum = out.sectors[lowest].branch.energy;
  // This scale is a reporting convention, not a solver acceptance tolerance
  // or an energy-error bound. Never change which numerical energy is lowest.
  out.comparison_band = Real{128} * uni20::numeric_limits<Real>::epsilon() * std::max(Real(sites), std::abs(minimum));
  for (std::size_t i = 0; i < count; ++i)
    if (out.sectors[i].branch.energy - minimum <= out.comparison_band) out.nearby_indices.push_back(i);
  return out;
}
} // namespace bethe::xxz::detail
