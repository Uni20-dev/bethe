// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/polynomial_roots.hpp>
#include <bethe/xxz_coordinate_wave.hpp>
#include <bethe/xxz_phantom_wave.hpp>
#include <optional>

namespace bethe::xxz::detail
{
enum class PhantomLiftStatus
{
  nonzero_witness,
  reduction_unresolved,
  roots_unresolved,
  resolution_unresolved,
  work_limit,
  nonfinite
};

template <uni20::Real Real> struct PhantomLiftOptions
{
    Real factor_tolerance = Real{128} * uni20::numeric_limits<Real>::epsilon();
    Real residual_tolerance = Real{32} * uni20::numeric_limits<Real>::epsilon();
    Real amplitude_tolerance = Real{128} * uni20::numeric_limits<Real>::epsilon();
    bethe::detail::PolynomialRootOptions<Real> root_options;
    std::size_t max_configurations = 16;
    /// Sum of 2^r over attempted finite-amplitude evaluations; bounds exponential
    /// work before allocation/callbacks. Each subset costs O(r^2) arithmetic.
    std::size_t max_subset_updates = 1000000;
};

/// A resolved numerical amplitude, NOT an interval proof, eigenstate
/// certificate, or ground-state classification. Unresolved never means zero.
template <uni20::Real Real> struct PhantomLiftCheck
{
    PhantomLiftStatus status = PhantomLiftStatus::reduction_unresolved;
    PhantomReduction<Real> reduction;
    std::optional<bethe::detail::PolynomialRoots<Real>> recovered;
    std::vector<std::size_t> occupied; // Best witness among tested configurations.
    std::complex<Real> amplitude{};
    Real absolute_term_sum{}, input_variation{}, roundoff_allowance{}, resolution_ratio{};
    std::size_t configurations_tested = 0, subset_updates = 0;
    bool all_configurations_tested = false;
};

/// Try ordered configurations lexicographically until one dressed amplitude
/// exceeds its propagated root-uncertainty and arithmetic allowances. The
/// finite regularity check and the separate phantom phase condition must
/// pass first. The input polynomial and requested Delta are never changed.
template <uni20::Real Real>
PhantomLiftCheck<Real> check_phantom_lift(PolynomialBetheSystem<Real> const& system, std::span<Real const> c,
                                          Real delta, std::size_t phantom_count, int chirality,
                                          PhantomLiftOptions<Real> const& options = {})
{
  using C = std::complex<Real>;
  if (!uni20::isfinite(options.amplitude_tolerance) || options.amplitude_tolerance <= Real{0} ||
      !uni20::isfinite(options.root_options.tolerance) || options.root_options.tolerance <= Real{0})
    throw std::invalid_argument("phantom lift tolerances must be finite and positive");
  PhantomLiftCheck<Real> out;
  out.reduction = reduce_phantom_polynomial(system, c, delta, phantom_count, chirality, options.factor_tolerance,
                                            options.residual_tolerance);
  if (out.reduction.status == PhantomReductionStatus::nonfinite)
  {
    out.status = PhantomLiftStatus::nonfinite;
    return out;
  }
  if (out.reduction.status != PhantomReductionStatus::regular_reduced_equations) return out;
  auto const r = out.reduction.finite_coefficients.size(), m = system.order;
  out.status = PhantomLiftStatus::work_limit;
  if (!options.max_configurations || r > options.root_options.max_degree ||
      r >= std::numeric_limits<std::size_t>::digits)
    return out;
  auto const subsets = std::size_t{1} << r;
  std::size_t terms;
  try
  {
    terms = phantom_dressing_terms(m, phantom_count, options.max_subset_updates / subsets);
  }
  catch (std::length_error const&)
  {
    return out;
  }
  auto const cost = terms * subsets; // Bounded by max_subset_updates above.
  out.recovered =
      bethe::detail::recover_polynomial_roots<Real>(out.reduction.finite_coefficients, options.root_options);
  if (out.recovered->status != bethe::detail::PolynomialRootStatus::resolved)
  {
    out.status = PhantomLiftStatus::roots_unresolved;
    return out;
  }
  out.status = PhantomLiftStatus::resolution_unresolved;
  std::vector<C> v;
  std::vector<Real> radii;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (C x : out.recovered->roots)
  {
    C const z = system.center + system.coordinate_scale * x, denominator = Real{1} + C{0, 1} * z;
    Real const rho =
        system.coordinate_scale * out.recovered->max_root_uncertainty + Real{64} * eps * std::max(Real{1}, std::abs(z));
    Real const d = std::abs(denominator);
    if (!uni20::isfinite(d) || !uni20::isfinite(rho) || !(d > rho)) return out;
    C const momentum = -(Real{1} - C{0, 1} * z) / denominator;
    Real const radius = Real{2} * (rho / d) / (d - rho) + Real{64} * eps * std::max(Real{1}, std::abs(momentum));
    if (!uni20::isfinite(momentum.real()) || !uni20::isfinite(momentum.imag()) || !uni20::isfinite(radius) ||
        momentum == C{})
      return out;
    v.push_back(momentum);
    radii.push_back(radius);
  }
  try
  {
    CoordinateBetheWave<Real> const wave(system.sites, delta, v, subsets);
    C const phase{delta, Real(chirality) * std::sqrt(Real{1} + delta) * std::sqrt(Real{1} - delta)};
    PhantomDressing<Real> const dressing(system.sites, r, phantom_count, phase, terms);
    std::vector<std::size_t> occupied(m);
    std::iota(occupied.begin(), occupied.end(), std::size_t{0});
    while (out.configurations_tested < options.max_configurations)
    {
      if (cost > options.max_subset_updates - out.subset_updates)
      {
        out.status = PhantomLiftStatus::work_limit;
        return out;
      }
      ++out.configurations_tested;
      bethe::detail::CompensatedSum<Real> real, imag, absolute, variation;
      dressing.for_each_term(occupied, [&](auto selected, C coefficient) {
        out.subset_updates += subsets;
        auto const finite = wave.evaluate(selected, radii);
        C const term = coefficient * finite.value;
        real.add(term.real());
        imag.add(term.imag());
        absolute.add(std::abs(coefficient) * finite.absolute_term_sum);
        variation.add(std::abs(coefficient) * finite.input_variation);
      });
      C const amplitude{real.value(), imag.value()};
      Real const allowance = options.amplitude_tolerance *
                             (Real{1} + Real(system.sites) * Real(m) + Real(r) * Real(r) + Real(terms)) *
                             absolute.value();
      Real const uncertainty = variation.value() + allowance;
      if (!uni20::isfinite(std::abs(amplitude)) || !uni20::isfinite(absolute.value()) || !uni20::isfinite(uncertainty))
      {
        out.status = PhantomLiftStatus::nonfinite;
        return out;
      }
      Real const ratio = uncertainty > Real{0} && absolute.value() >= uni20::numeric_limits<Real>::min()
                             ? std::abs(amplitude) / uncertainty
                             : Real{0};
      if (out.occupied.empty() || ratio > out.resolution_ratio)
      {
        out.occupied = occupied;
        out.amplitude = amplitude;
        out.absolute_term_sum = absolute.value();
        out.input_variation = variation.value();
        out.roundoff_allowance = allowance;
        out.resolution_ratio = ratio;
      }
      // Determine exhaustion without computing the potentially huge C(N,M).
      auto j = m;
      while (j && occupied[j - 1] == system.sites - m + j - 1)
        --j;
      out.all_configurations_tested = j == 0;
      if (ratio > Real{1})
      {
        out.status = PhantomLiftStatus::nonzero_witness;
        return out;
      }
      if (!j) return out;
      ++occupied[j - 1];
      for (std::size_t k = j; k < m; ++k)
        occupied[k] = occupied[k - 1] + 1;
    }
    out.status = PhantomLiftStatus::work_limit;
  }
  catch (std::overflow_error const&)
  {
    out.status = PhantomLiftStatus::nonfinite;
  }
  return out;
}
} // namespace bethe::xxz::detail
