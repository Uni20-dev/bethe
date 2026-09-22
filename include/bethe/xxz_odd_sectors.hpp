// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/xxz_odd_continuation.hpp>
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
};

/// All spin-reversal-distinct odd-ring continuation branches, in increasing
/// M (decreasing positive Sz). This is an INTERNAL candidate scan, not a
/// public ground-state result or a physical-state certificate.
template <uni20::Real Real> struct OddSectorScan
{
    std::vector<OddSectorCandidate<Real>> sectors;
    bool equations_complete = false;
    bool wronskians_consistent = false;
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
                                                                           uni20::numeric_limits<Real>::epsilon())
{
  auto const n = checked_sites(sites);
  if (sites % 2 == 0) throw std::invalid_argument("odd XXZ sector scan requires odd N");
  if (!uni20::isfinite(wronskian_tolerance) || wronskian_tolerance <= Real{0})
    throw std::invalid_argument("Wronskian tolerance must be finite and positive");
  auto const count = sites / 2 + 1;
  if (count > std::size_t(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(OddSectorCandidate<Real>))
    throw std::length_error("odd XXZ sector scan is too large");
  // Also validate Delta and Newton options before allocating the scan.
  auto vacuum = continue_odd_polynomial(sites, delta, uni20::from_twice(n), options);
  OddSectorScan<Real> out;
  out.sectors.reserve(count);
  out.equations_complete = true;
  out.wronskians_consistent = true;
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
    }
    else
    {
      out.equations_complete = false;
      out.wronskians_consistent = false;
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
