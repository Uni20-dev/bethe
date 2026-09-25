// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/complex_math.hpp>
#include <bethe/detail/gauss_legendre.hpp>
#include <bethe/sine_gordon.hpp>

namespace bethe::sine_gordon
{
enum class VacuumStatus
{
  converged,
  kernel_limit,
  iteration_limit,
  mesh_limit,
  cutoff_limit,
  contour_limit,
  precision_limit
};
template <uni20::Real Real> struct VacuumOptions
{
    Real tolerance = Real{262144} * uni20::numeric_limits<Real>::epsilon(); // absolute error in Y=L*E_C
    std::size_t initial_intervals = 64, max_intervals = 2048, max_iterations = 10000, max_cutoffs = 3;
    std::optional<Real> contour_shift{}, initial_cutoff{};
    KernelOptions<Real> kernel{}; // tolerance is assigned from the vacuum error budget
};
template <uni20::Real Real> struct VacuumState
{
    Real mass{}, length{}, coupling{}, scaled_length{}, contour_shift{}, verification_contour_shift{}, cutoff{};
    std::optional<Real> scaling_function, casimir_energy, effective_central_charge;
    Real nonlinear_residual = uni20::numeric_limits<Real>::infinity();
    Real kernel_error = uni20::numeric_limits<Real>::infinity(), mesh_error = uni20::numeric_limits<Real>::infinity();
    Real cutoff_error = uni20::numeric_limits<Real>::infinity(),
         contour_error = uni20::numeric_limits<Real>::infinity();
    std::size_t intervals = 0, iterations = 0, kernel_evaluations = 0, cutoffs = 0;
    bool converged = false;
    VacuumStatus status = VacuumStatus::mesh_limit;
};
namespace detail
{
template <uni20::Real Real> struct KernelTable
{
    std::vector<Real> real;
    std::vector<std::complex<Real>> shifted;
    Real error = uni20::numeric_limits<Real>::infinity();
    bool converged = false;
};
// Composite native Gauss quadrature for all difference-grid entries at once.
// Expensive Fourier multipliers are shared; phase recurrences restart every
// 16 entries to limit drift. Two successive mesh agreements are required.
template <uni20::Real Real>
KernelTable<Real> kernel_table(Real p, Real eta, Real h, std::size_t n, KernelOptions<Real> options,
                               std::size_t& total_evaluations)
{
  using C = std::complex<Real>;
  using bethe::detail::CompensatedSum;
  KernelTable<Real> out{std::vector<Real>(n + 1), std::vector<C>(n + 1)};
  if (p == Real{1})
  {
    out.error = Real{0};
    out.converged = true;
    return out;
  }
  Real const pi = Real{4} * std::atan(Real{1}), width = pi * std::min(p, Real{1}) - Real{2} * eta;
  Real cutoff = Real{1} / width, tail{};
  bool tail_ok = false;
  for (std::size_t i = 0; i < options.max_cutoffs; ++i)
  {
    if (!uni20::isfinite(cutoff)) return out;
    tail = std::exp(-width * cutoff) / pi / width / (-std::expm1(-p * pi * cutoff));
    if (tail <= options.tolerance / Real{4})
    {
      tail_ok = true;
      break;
    }
    cutoff *= Real{2};
  }
  if (!tail_ok) return out;
  auto const rule = bethe::detail::gauss_legendre<Real>(16);
  unsigned stable = 0;
  std::size_t used = 0, panels = 1;
  for (std::size_t level = 0; level < options.max_levels; ++level)
  {
    if (panels > (options.max_evaluations - used) / 16) return out;
    std::vector<CompensatedSum<Real>> re(n + 1);
    std::vector<CompensatedSum<C>> shifted(n + 1);
    CompensatedSum<Real> absolute;
    Real const panel_width = cutoff / Real(panels);
    for (std::size_t panel = 0; panel < panels; ++panel)
      for (std::size_t q = 0; q < 16; ++q)
      {
        ++used;
        ++total_evaluations;
        Real const k = panel_width * (Real(panel) + (Real{1} + rule.x[q]) / Real{2});
        Real const w = panel_width * rule.w[q] / Real{2};
        auto const a = fourier_components(k, Real{0}, p), b = fourier_components(k, Real{2} * eta, p);
        absolute.add(w * (std::abs(a.real()) + std::abs(b.real()) + std::abs(b.imag())));
        C phase(1), step(std::cos(k * h), std::sin(k * h));
        for (std::size_t j = 0; j <= n; ++j)
        {
          if (j % 16 == 0) phase = C(std::cos(k * h * Real(j)), std::sin(k * h * Real(j)));
          re[j].add(w * a.real() * phase.real());
          shifted[j].add(w * C(b.real() * phase.real(), -b.imag() * phase.imag()));
          phase *= step;
        }
      }
    Real difference{};
    for (std::size_t j = 0; j <= n; ++j)
    {
      if (!uni20::isfinite(re[j].value()) || !uni20::isfinite(shifted[j].value().real()) ||
          !uni20::isfinite(shifted[j].value().imag()))
        return out;
      difference =
          std::max({difference, std::abs(re[j].value() - out.real[j]), std::abs(shifted[j].value() - out.shifted[j])});
      out.real[j] = re[j].value();
      out.shifted[j] = shifted[j].value();
    }
    out.error = tail + std::max(difference, Real{64} * uni20::numeric_limits<Real>::epsilon() * absolute.value());
    stable = level > 0 && uni20::isfinite(out.error) && out.error <= options.tolerance ? stable + 1 : 0;
    if (stable == 2)
    {
      out.converged = true;
      return out;
    }
    if (panels > options.max_evaluations / 32) return out;
    panels *= 2;
  }
  return out;
}
template <uni20::Real Real> struct VacuumMesh
{
    Real value{}, residual = uni20::numeric_limits<Real>::infinity();
    Real kernel_error = uni20::numeric_limits<Real>::infinity();
    bool converged = false;
    VacuumStatus status = VacuumStatus::iteration_limit;
};
template <uni20::Real Real>
VacuumMesh<Real> vacuum_mesh(Real p, Real u, Real eta, Real cutoff, std::size_t n, VacuumOptions<Real> const& options,
                             VacuumState<Real>& work)
{
  using C = std::complex<Real>;
  using bethe::detail::CompensatedSum;
  VacuumMesh<Real> out;
  Real const h = Real{2} * cutoff / Real(n), pi = Real{4} * std::atan(Real{1});
  if (!uni20::isfinite(h) || !(h > Real{0}))
  {
    out.status = VacuumStatus::precision_limit;
    return out;
  }
  auto kernel_options = options.kernel;
  kernel_options.tolerance = options.tolerance / (Real{64} * (Real{1} + cutoff));
  auto const kernel = kernel_table(p, eta, h, n, kernel_options, work.kernel_evaluations);
  out.kernel_error = kernel.error;
  if (!kernel.converged)
  {
    out.status = VacuumStatus::kernel_limit;
    return out;
  }
  std::vector<C> bare(n + 1), correction(n + 1), a(n + 1), next(n + 1);
  std::vector<Real> weighted_real(n + 1), weighted_imag(n + 1);
  std::vector<Real> weight(n + 1, h);
  weight.front() /= Real{2};
  weight.back() /= Real{2};
  for (std::size_t i = 0; i <= n; ++i)
  {
    Real const x = h * (Real(i) - Real(n) / Real{2});
    bare[i] = C(u * std::sin(eta) * std::cosh(x), -u * std::cos(eta) * std::sinh(x));
  }
  for (;;)
  {
    CompensatedSum<Real> energy;
    for (std::size_t j = 0; j <= n; ++j)
    {
      C const epsilon = bare[j] + correction[j];
      if (!uni20::isfinite(epsilon.real()) || !uni20::isfinite(epsilon.imag()))
      {
        out.status = VacuumStatus::precision_limit;
        return out;
      }
      C const fugacity = std::exp(-epsilon);
      // A damped transient may have Re(epsilon)<0 even on the vacuum
      // branch. Test the actual logarithm argument, not that stronger
      // sufficient condition, and stay in its open right half-plane.
      if (!uni20::isfinite(fugacity.real()) || !uni20::isfinite(fugacity.imag()) ||
          !(Real{1} + fugacity.real() > Real{0}))
      {
        out.status = VacuumStatus::precision_limit;
        return out;
      }
      a[j] = bethe::detail::complex_log1p(fugacity);
      weighted_real[j] = weight[j] * a[j].real();
      weighted_imag[j] = weight[j] * a[j].imag();
      energy.add(-weight[j] * (bare[j] * a[j]).real() / pi);
    }
    out.value = energy.value();
    out.residual = Real{0};
    if (p == Real{1} && uni20::isfinite(out.value))
    {
      out.converged = true;
      out.status = VacuumStatus::converged;
      return out;
    }
    // Vacuum parity epsilon(-theta)=conj(epsilon(theta)). Evaluate half
    // the convolution, using real components to avoid redundant products.
    for (std::size_t i = 0; i <= n / 2; ++i)
    {
      CompensatedSum<Real> real_sum, imag_sum;
      for (std::size_t j = 0; j <= n; ++j)
      {
        auto const d = i > j ? i - j : j - i;
        Real const kr = kernel.shifted[d].real(), ki = i >= j ? kernel.shifted[d].imag() : -kernel.shifted[d].imag();
        real_sum.add((kr - kernel.real[d]) * weighted_real[j] + ki * weighted_imag[j]);
        imag_sum.add(ki * weighted_real[j] - (kr + kernel.real[d]) * weighted_imag[j]);
      }
      next[i] = C(real_sum.value(), i == n / 2 ? Real{0} : imag_sum.value());
      next[n - i] = std::conj(next[i]);
      Real const residual = std::abs(next[i] - correction[i]);
      if (!uni20::isfinite(residual) || !uni20::isfinite(out.value))
      {
        out.status = VacuumStatus::precision_limit;
        return out;
      }
      out.residual = std::max(out.residual, residual);
    }
    if (out.residual <= options.tolerance / (Real{32} * (Real{1} + cutoff)))
    {
      out.converged = true;
      out.status = VacuumStatus::converged;
      return out;
    }
    if (work.iterations == options.max_iterations) return out;
    ++work.iterations;
    for (std::size_t j = 0; j <= n; ++j)
      correction[j] = (correction[j] + next[j]) / Real{2};
  }
}
template <uni20::Real Real>
std::optional<Real> resolved_mesh(Real p, Real u, Real eta, Real cutoff, VacuumOptions<Real> const& options,
                                  VacuumState<Real>& out)
{
  std::optional<Real> previous;
  unsigned stable = 0;
  for (std::size_t n = options.initial_intervals;;)
  {
    out.intervals = n;
    out.cutoff = cutoff;
    auto const mesh = vacuum_mesh(p, u, eta, cutoff, n, options, out);
    out.nonlinear_residual = mesh.residual;
    out.kernel_error = mesh.kernel_error;
    if (!mesh.converged)
    {
      out.status = mesh.status;
      return {};
    }
    out.mesh_error = previous ? std::abs(mesh.value - *previous) : uni20::numeric_limits<Real>::infinity();
    stable = previous && out.mesh_error <= options.tolerance / Real{4} ? stable + 1 : 0;
    if (stable == 2) return mesh.value;
    previous = mesh.value;
    if (n > options.max_intervals / 2)
    {
      out.status = VacuumStatus::mesh_limit;
      return {};
    }
    n *= 2;
  }
}
template <uni20::Real Real>
std::optional<Real> resolved_cutoff(Real p, Real u, Real eta, Real cutoff, VacuumOptions<Real> const& options,
                                    VacuumState<Real>& out)
{
  std::optional<Real> previous;
  for (std::size_t attempt = 0; attempt < options.max_cutoffs; ++attempt)
  {
    ++out.cutoffs;
    auto const value = resolved_mesh(p, u, eta, cutoff, options, out);
    if (!value) return {};
    out.cutoff_error = previous ? std::abs(*value - *previous) : uni20::numeric_limits<Real>::infinity();
    if (previous && out.cutoff_error <= options.tolerance / Real{4}) return value;
    previous = value;
    cutoff += Real{1};
  }
  out.status = VacuumStatus::cutoff_limit;
  return {};
}
} // namespace detail

/// Untwisted, zero-charge, bulk-subtracted vacuum energy. Numerical error
/// controls are absolute in Y=L*(E0-L*e_bulk), not relative infrared accuracy.
template <uni20::Real Real>
VacuumState<Real> vacuum_energy(Real mass, Real length, Real p, VacuumOptions<Real> options = {})
{
  Real const pi = Real{4} * std::atan(Real{1}), strip = pi * std::min(p, Real{1});
  Real const eta = options.contour_shift.value_or(strip / Real{4});
  if (!uni20::isfinite(mass) || !(mass > Real{0}) || !uni20::isfinite(length) || !(length > Real{0}) ||
      !uni20::isfinite(p) || !(p > Real{0}) || !uni20::isfinite(eta) || !(eta > Real{0} && eta < strip / Real{2}) ||
      !uni20::isfinite(options.tolerance) || !(options.tolerance > Real{0}) || options.initial_intervals < 8 ||
      options.initial_intervals % 2 || options.max_intervals < options.initial_intervals ||
      (options.initial_cutoff && (!uni20::isfinite(*options.initial_cutoff) || !(*options.initial_cutoff > Real{0}))))
    throw std::invalid_argument("invalid sine-Gordon vacuum scales, coupling, contour, tolerance or mesh controls");
  if (options.max_intervals > 16384)
    throw std::length_error("sine-Gordon vacuum mesh exceeds 16384-interval work ceiling");
  VacuumState<Real> out;
  out.mass = mass;
  out.length = length;
  out.coupling = p;
  out.scaled_length = mass * length;
  out.contour_shift = eta;
  if (!uni20::isfinite(out.scaled_length) || !(out.scaled_length > Real{0}))
  {
    out.status = VacuumStatus::precision_limit;
    return out;
  }
  // Choose the initial tail scale for the smaller of the two tested contours.
  Real const second_eta = eta * Real{0.75};
  out.verification_contour_shift = second_eta;
  Real const tail_scale =
      (-std::log(std::min(options.tolerance, Real{0.1})) + Real{8}) / out.scaled_length / std::sin(second_eta);
  Real const cutoff = options.initial_cutoff.value_or(std::acosh(std::max(Real{2}, tail_scale)));
  if (!uni20::isfinite(cutoff))
  {
    out.status = VacuumStatus::precision_limit;
    return out;
  }
  auto const first = detail::resolved_cutoff(p, out.scaled_length, eta, cutoff, options, out);
  if (!first) return out;
  auto const primary = out;
  auto const second = detail::resolved_cutoff(p, out.scaled_length, second_eta, cutoff, options, out);
  if (!second) return out;
  out.contour_error = std::abs(*first - *second);
  out.intervals = std::max(out.intervals, primary.intervals);
  out.cutoff = std::max(out.cutoff, primary.cutoff);
  out.mesh_error = std::max(out.mesh_error, primary.mesh_error);
  out.cutoff_error = std::max(out.cutoff_error, primary.cutoff_error);
  out.kernel_error = std::max(out.kernel_error, primary.kernel_error);
  out.nonlinear_residual = std::max(out.nonlinear_residual, primary.nonlinear_residual);
  if (out.contour_error > options.tolerance / Real{2})
  {
    out.status = VacuumStatus::contour_limit;
    return out;
  }
  Real const energy = *first / length, central = -Real{6} * *first / pi;
  if (!uni20::isfinite(energy) || !uni20::isfinite(central))
  {
    out.status = VacuumStatus::precision_limit;
    return out;
  }
  out.scaling_function = *first;
  out.casimir_energy = energy;
  out.effective_central_charge = central;
  out.converged = true;
  out.status = VacuumStatus::converged;
  return out;
}
} // namespace bethe::sine_gordon
