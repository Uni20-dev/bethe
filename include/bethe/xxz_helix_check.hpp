// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/xxz_polynomial.hpp>
#include <bethe/xxz_spin_helix.hpp>

namespace bethe::xxz::detail
{
enum class HelixMatchStatus
{
  compatible,
  different_coupling,
  different_polynomial,
  nonfinite
};

template <uni20::Real Real> struct PolynomialHelixCheck
{
    HelixMatchStatus status = HelixMatchStatus::nonfinite;
    Real helix_delta = Real{0}, helix_energy = Real{0};
    Real coupling_error = uni20::numeric_limits<Real>::infinity();
    Real coefficient_error = uni20::numeric_limits<Real>::infinity();
    /// Bound on ||(H(input Delta)-E_helix) psi_helix||/||psi_helix||,
    /// apart from roundoff: N*|input Delta-helix Delta|/4.
    Real eigenvector_defect_bound = uni20::numeric_limits<Real>::infinity();
};

/// Compare against the repeated-root polynomial of an explicitly nonzero
/// helix eigenvector. Numerical compatibility is not exact equality, does
/// not snap Delta, and says nothing about sector minimality. In particular,
/// an endpoint root or a root-of-unity coupling ALONE is insufficient.
template <uni20::Real Real>
PolynomialHelixCheck<Real> check_helix_polynomial(PolynomialBetheSystem<Real> const& system, std::span<Real const> c,
                                                  Real delta, std::size_t winding,
                                                  Real tolerance = Real{128} * uni20::numeric_limits<Real>::epsilon())
{
  if (!uni20::isfinite(tolerance) || tolerance <= Real{0})
    throw std::invalid_argument("helix matching tolerance must be finite and positive");
  if (!uni20::isfinite(delta) || delta <= -Real{1} || delta > Real{0})
    throw std::invalid_argument("XXZ polynomial helix matching requires -1 < Delta <= 0");
  if (c.size() != system.order) throw std::invalid_argument("helix polynomial coefficient count differs from degree");
  for (Real value : c)
    if (!uni20::isfinite(value)) throw std::invalid_argument("nonfinite helix polynomial coefficient");
  auto const state = periodic_spin_helix<Real>(system.sites, winding,
                                               uni20::from_twice(std::int64_t(system.sites - 2 * system.order)));
  if (state.delta <= -Real{1} || state.delta > Real{0})
    throw std::invalid_argument("helix winding lies outside the polynomial coordinate domain");
  PolynomialHelixCheck<Real> out;
  out.helix_delta = state.delta;
  out.helix_energy = state.energy;
  out.coupling_error = std::abs(delta - state.delta);
  out.eigenvector_defect_bound = (Real(system.sites) / Real{4}) * out.coupling_error;
  if (out.coupling_error > tolerance)
  {
    out.status = HelixMatchStatus::different_coupling;
    return out;
  }
  // e^(ik)=-(1-i*z)/(1+i*z), so z=cot(k/2). This expression
  // is well-conditioned in the nonpositive-Delta domain used here.
  auto const phase = helix_phase<Real>(system.sites, winding);
  Real const z = phase.imag() / (Real{1} - phase.real());
  Real const x = (z - system.center) / system.coordinate_scale;
  if (!uni20::isfinite(x)) return out;
  std::vector<Real> expected(c.size() + 1, Real{0});
  expected[0] = Real{1};
  for (std::size_t degree = 0; degree < c.size(); ++degree)
    for (std::size_t j = degree + 2; j-- > 0;)
    {
      expected[j] = (j ? expected[j - 1] : Real{0}) - x * expected[j];
      if (!uni20::isfinite(expected[j])) return out;
    }
  out.coefficient_error = Real{0};
  for (std::size_t j = 0; j < c.size(); ++j)
  {
    Real const scale = std::max({Real{1}, std::abs(c[j]), std::abs(expected[j])});
    // Normalize before subtraction to avoid overflow of a finite pair.
    out.coefficient_error = std::max(out.coefficient_error, std::abs(c[j] / scale - expected[j] / scale));
  }
  out.status =
      out.coefficient_error <= tolerance ? HelixMatchStatus::compatible : HelixMatchStatus::different_polynomial;
  return out;
}
} // namespace bethe::xxz::detail
