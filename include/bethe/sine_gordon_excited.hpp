// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/sine_gordon_holes.hpp>
#include <uni20/common/half_int.hpp>

namespace bethe::sine_gordon
{
template <uni20::Real Real> struct TwoSolitonOptions : VacuumOptions<Real>
{
    std::size_t max_root_iterations = 1024; // total across meshes/contours
    TwoSolitonOptions()
    {
      this->tolerance = Real{16777216} * uni20::numeric_limits<Real>::epsilon();
      this->max_iterations = 32768; // coupled hole trials require more updates than a vacuum solve
    }
};
template <uni20::Real Real> struct TwoSolitonState : VacuumState<Real>
{
    std::optional<Real> rapidity;
    uni20::half_int number = uni20::from_twice(1);
    Real quantization_residual = uni20::numeric_limits<Real>::infinity();
    Real hole_error = uni20::numeric_limits<Real>::infinity(); // estimate in Y and 2*u*sinh(H)
    Real source_quadrature_error = uni20::numeric_limits<Real>::infinity();
    Real source_tail_bound = uni20::numeric_limits<Real>::infinity();
    std::size_t root_iterations = 0, source_evaluations = 0;
};
namespace detail
{
template <uni20::Real Real> struct TwoSolitonMesh : VacuumMesh<Real>
{
    Real hole{}, coordinate{};
};
template <uni20::Real Real>
TwoSolitonMesh<Real> two_soliton_mesh(Real p, Real u, Real eta, Real cutoff, std::size_t n,
                                      TwoSolitonOptions<Real> const& options, TwoSolitonState<Real>& work,
                                      std::optional<Real> initial_hole = {})
{
  TwoSolitonMesh<Real> failed;
  work.quantization_residual = work.hole_error = work.source_quadrature_error = work.source_tail_bound =
      uni20::numeric_limits<Real>::infinity();
  if (work.root_iterations == options.max_root_iterations) return failed;
  Real const h = Real{2} * cutoff / Real(n);
  if (!uni20::isfinite(h) || !(h > Real{0}) || !uni20::isfinite(u * std::cosh(cutoff)))
  {
    failed.status = VacuumStatus::precision_limit;
    return failed;
  }
  Real const pi = Real{4} * std::atan(Real{1}), target = pi * Real(work.number.twice());
  auto kernel_options = options.kernel;
  kernel_options.tolerance = options.tolerance / (Real{128} * (Real{1} + cutoff));
  auto const kernel = kernel_table(p, eta, h, n, kernel_options, work.kernel_evaluations);
  if (!kernel.converged)
  {
    failed.kernel_error = kernel.error;
    failed.status = VacuumStatus::kernel_limit;
    return failed;
  }
  Real lower{}, upper = std::asinh(target / u) + Real{1};
  Real previous = upper, previous_f{}, slope = u * std::cosh(upper);
  bool have_previous = false;
  bool upper_verified = false; // Z(0)=0 by parity; verify the positive bracket endpoint numerically.
  // Coarse-grid roots are guesses only: solve the new source/counting equation
  // and check its residual before accepting any result on this grid.
  Real candidate = initial_hole.value_or(std::asinh(target / u));
  if (!uni20::isfinite(candidate) || !uni20::isfinite(slope))
  {
    failed.status = VacuumStatus::precision_limit;
    return failed;
  }
  auto inner_options = static_cast<VacuumOptions<Real> const&>(options);
  inner_options.tolerance /= Real{64};
  std::vector<std::complex<Real>> seed; // never shared between grids, cutoffs or contours
  // A safeguarded secant solves the real-axis hole equation. Every trial
  // solves the contour NLIE, not just its infinite-volume driving term.
  while (work.root_iterations < options.max_root_iterations)
  {
    ++work.root_iterations;
    Real const hole = candidate;
    auto const source = hole_table(p, eta, cutoff, n, hole, kernel_options, work.source_evaluations);
    work.quantization_residual = work.hole_error = uni20::numeric_limits<Real>::infinity();
    work.source_quadrature_error = source.quadrature_error;
    work.source_tail_bound = source.tail_bound;
    if (!source.converged)
    {
      failed.status = VacuumStatus::kernel_limit;
      return failed;
    }
    auto const inner =
        nlie_mesh<Real>(p, u, eta, cutoff, n, inner_options, work, kernel, source.source, source.counting, &seed);
    TwoSolitonMesh<Real> out;
    static_cast<VacuumMesh<Real>&>(out) = inner;
    if (!inner.converged) return out;
    out.converged = false;
    out.status = VacuumStatus::iteration_limit;
    Real const f = u * std::sinh(hole) + source.phase + inner.counting_integral - target;
    Real const uncertainty = (Real{1} + Real{8} * cutoff) * (source.error + inner.residual);
    work.quantization_residual = std::abs(f);
    out.value += Real{2} * u * std::cosh(hole);
    out.hole = hole;
    out.coordinate = Real{2} * u * std::sinh(hole);
    if (!uni20::isfinite(f) || !uni20::isfinite(out.value))
    {
      out.status = VacuumStatus::precision_limit;
      return out;
    }
    if (have_previous && hole != previous && std::abs(f - previous_f) > Real{4} * uncertainty)
    {
      Real const secant = (f - previous_f) / (hole - previous);
      if (secant > Real{0} && uni20::isfinite(secant)) slope = secant;
    }
    if (!have_previous) slope = u * std::cosh(hole) + Real{2} * source.phase_derivative;
    Real const step = have_previous ? std::abs(hole - previous) : uni20::numeric_limits<Real>::infinity();
    work.hole_error = (Real{1} + std::abs(out.coordinate)) * (step + (std::abs(f) + uncertainty) / slope);
    if (have_previous && work.hole_error <= options.tolerance / Real{8} && std::abs(f) <= options.tolerance / Real{8})
    {
      out.converged = true;
      out.status = VacuumStatus::converged;
      return out;
    }
    if (f > Real{0})
    {
      upper = hole;
      upper_verified = true;
    }
    else
      lower = hole;
    previous = hole;
    previous_f = f;
    have_previous = true;
    candidate = hole - f / slope;
    if (!upper_verified && f != Real{0})
    {
      if (!(upper > lower)) upper = Real{2} * lower + Real{1};
      candidate = upper;
    }
    // Repeating a converged numerical root permits the successive-position
    // check to pass without manufacturing a nonzero displacement.
    if (upper_verified && candidate != hole && !(candidate > lower && candidate < upper))
      candidate = lower + (upper - lower) / Real{2};
  }
  return failed;
}
template <uni20::Real Real>
std::optional<TwoSolitonMesh<Real>> resolved_two_soliton_mesh(Real p, Real u, Real eta, Real cutoff,
                                                              TwoSolitonOptions<Real> const& options,
                                                              TwoSolitonState<Real>& work)
{
  std::optional<TwoSolitonMesh<Real>> previous;
  unsigned stable = 0;
  for (std::size_t n = options.initial_intervals;;)
  {
    work.intervals = n;
    work.cutoff = cutoff;
    auto const mesh = two_soliton_mesh(p, u, eta, cutoff, n, options, work,
                                       previous ? std::optional<Real>(previous->hole) : std::nullopt);
    work.nonlinear_residual = mesh.residual;
    work.kernel_error = mesh.kernel_error;
    if (!mesh.converged)
    {
      work.status = mesh.status;
      return {};
    }
    work.mesh_error =
        previous ? std::max(std::abs(mesh.value - previous->value), std::abs(mesh.coordinate - previous->coordinate))
                 : uni20::numeric_limits<Real>::infinity();
    stable = previous && work.mesh_error <= options.tolerance / Real{4} ? stable + 1 : 0;
    if (stable == 2) return mesh;
    previous = mesh;
    if (n > options.max_intervals / 2)
    {
      work.status = VacuumStatus::mesh_limit;
      return {};
    }
    n *= 2;
  }
}
template <uni20::Real Real>
std::optional<TwoSolitonMesh<Real>> resolved_two_soliton_cutoff(Real p, Real u, Real eta, Real cutoff,
                                                                TwoSolitonOptions<Real> const& options,
                                                                TwoSolitonState<Real>& work)
{
  std::optional<TwoSolitonMesh<Real>> previous;
  for (std::size_t attempt = 0; attempt < options.max_cutoffs; ++attempt)
  {
    ++work.cutoffs;
    auto const mesh = resolved_two_soliton_mesh(p, u, eta, cutoff, options, work);
    if (!mesh) return {};
    work.cutoff_error =
        previous ? std::max(std::abs(mesh->value - previous->value), std::abs(mesh->coordinate - previous->coordinate))
                 : uni20::numeric_limits<Real>::infinity();
    if (previous && work.cutoff_error <= options.tolerance / Real{4}) return mesh;
    previous = mesh;
    cutoff += Real{1};
  }
  work.status = VacuumStatus::cutoff_limit;
  return {};
}
} // namespace detail

/// Exact continuum NLIE for the symmetric two-hole charge +/-2 family,
/// p>=1, I=1/2 or 3/2. Returns E-L*e_bulk, NOT a vacuum-relative gap.
/// Uses delta=0 (half-odd labels), no complex/special roots. Numerical
/// success additionally requires two contours, cutoffs and mesh refinement.
template <uni20::Real Real>
TwoSolitonState<Real> two_soliton_level(Real mass, Real length, Real p, uni20::half_int number = uni20::from_twice(1),
                                        TwoSolitonOptions<Real> options = {})
{
  Real const pi = Real{4} * std::atan(Real{1}), eta = options.contour_shift.value_or(pi / Real{4});
  if (!uni20::isfinite(mass) || mass <= Real{0} || !uni20::isfinite(length) || length <= Real{0} ||
      !uni20::isfinite(p) || p < Real{1} || (number.twice() != 1 && number.twice() != 3) || !uni20::isfinite(eta) ||
      eta <= Real{0} || eta >= pi / Real{2} || !uni20::isfinite(options.tolerance) || options.tolerance <= Real{0} ||
      options.initial_intervals < 8 || options.initial_intervals % 2 ||
      options.max_intervals < options.initial_intervals ||
      (options.initial_cutoff && (!uni20::isfinite(*options.initial_cutoff) || *options.initial_cutoff <= Real{0})))
    throw std::invalid_argument("invalid sine-Gordon two-hole scales, p>=1, I=1/2 or 3/2, or numerical controls");
  if (options.max_intervals > 16384) throw std::length_error("sine-Gordon two-hole mesh exceeds work ceiling");
  TwoSolitonState<Real> out;
  out.mass = mass;
  out.length = length;
  out.coupling = p;
  out.number = number;
  out.scaled_length = mass * length;
  out.contour_shift = eta;
  out.verification_contour_shift = eta * Real{0.75};
  if (!uni20::isfinite(out.scaled_length) || out.scaled_length <= Real{0})
  {
    out.status = VacuumStatus::precision_limit;
    return out;
  }
  Real const tail_scale = (-std::log(std::min(options.tolerance, Real{0.1})) + Real{8}) / out.scaled_length /
                          std::sin(out.verification_contour_shift);
  Real const cutoff = options.initial_cutoff.value_or(std::acosh(std::max(Real{2}, tail_scale)));
  if (!uni20::isfinite(cutoff))
  {
    out.status = VacuumStatus::precision_limit;
    return out;
  }
  auto const first = detail::resolved_two_soliton_cutoff(p, out.scaled_length, eta, cutoff, options, out);
  if (!first) return out;
  auto const primary = out;
  auto const second =
      detail::resolved_two_soliton_cutoff(p, out.scaled_length, out.verification_contour_shift, cutoff, options, out);
  if (!second) return out;
  out.contour_error =
      std::max(std::abs(first->value - second->value), std::abs(first->coordinate - second->coordinate));
  out.intervals = std::max(out.intervals, primary.intervals);
  out.cutoff = std::max(out.cutoff, primary.cutoff);
  out.mesh_error = std::max(out.mesh_error, primary.mesh_error);
  out.cutoff_error = std::max(out.cutoff_error, primary.cutoff_error);
  out.kernel_error = std::max(out.kernel_error, primary.kernel_error);
  out.nonlinear_residual = std::max(out.nonlinear_residual, primary.nonlinear_residual);
  out.hole_error = std::max(out.hole_error, primary.hole_error);
  out.source_quadrature_error = std::max(out.source_quadrature_error, primary.source_quadrature_error);
  out.source_tail_bound = std::max(out.source_tail_bound, primary.source_tail_bound);
  out.quantization_residual = std::max(out.quantization_residual, primary.quantization_residual);
  if (out.contour_error > options.tolerance / Real{2})
  {
    out.status = VacuumStatus::contour_limit;
    return out;
  }
  Real const energy = first->value / length;
  if (!uni20::isfinite(energy))
  {
    out.status = VacuumStatus::precision_limit;
    return out;
  }
  out.scaling_function = first->value;
  out.casimir_energy = energy;
  out.rapidity = first->hole;
  out.converged = true;
  out.status = VacuumStatus::converged;
  return out;
}
} // namespace bethe::sine_gordon
