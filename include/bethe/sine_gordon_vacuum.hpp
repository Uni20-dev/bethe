// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/complex_math.hpp>
#include <bethe/detail/gauss_legendre.hpp>
#include <bethe/sine_gordon.hpp>
#include <span>

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
    uni20::Precision precision() const { return uni20::native_precision<Real>(); }
};
#if UNI20_HAS_MPREAL
template <> struct VacuumOptions<uni20::mpreal>
{
    explicit VacuumOptions(uni20::Precision p) : tolerance(default_tolerance(p)), kernel(p), precision_(p) {}
    uni20::mpreal tolerance;
    std::size_t initial_intervals = 64, max_intervals = 2048, max_iterations = 10000, max_cutoffs = 3;
    std::optional<uni20::mpreal> contour_shift{}, initial_cutoff{};
    KernelOptions<uni20::mpreal> kernel;
    uni20::Precision precision() const { return precision_; }

  private:
    static uni20::mpreal default_tolerance(uni20::Precision p)
    {
      if (p.is_exact()) throw std::invalid_argument("sine-Gordon vacuum requires finite working precision");
      return 262144 * uni20::numeric_limits<uni20::mpreal>::epsilon(p);
    }
    uni20::Precision precision_;
};
#endif
template <uni20::Real Real> struct VacuumState
{
    VacuumState()
      requires(!uni20::has_runtime_precision_v<Real>)
        : VacuumState(Real{})
    {}
    explicit VacuumState(Real const& zero)
        : mass(zero), length(zero), coupling(zero), scaled_length(zero), contour_shift(zero),
          verification_contour_shift(zero), cutoff(zero),
          nonlinear_residual(uni20::scalar_like(zero, std::numeric_limits<double>::infinity())),
          kernel_error(nonlinear_residual), mesh_error(nonlinear_residual), cutoff_error(nonlinear_residual),
          contour_error(nonlinear_residual)
    {}
    Real mass, length, coupling, scaled_length, contour_shift, verification_contour_shift, cutoff;
    std::optional<Real> scaling_function, casimir_energy, effective_central_charge;
    Real nonlinear_residual, kernel_error, mesh_error, cutoff_error, contour_error;
    std::size_t intervals = 0, iterations = 0, kernel_evaluations = 0, cutoffs = 0;
    bool converged = false;
    VacuumStatus status = VacuumStatus::mesh_limit;
};
namespace detail
{
template <uni20::Real Real> struct KernelTable
{
    std::vector<Real> real;
    std::vector<uni20::complex<Real>> shifted;
    Real error;
    bool converged = false;
};
template <uni20::Real Real> std::size_t batch_fourier_order(uni20::Precision precision)
{
  if constexpr (uni20::has_runtime_precision_v<Real>)
    return std::clamp<std::size_t>(precision.bit_count() / 4, 16, 128);
  else
    return 16;
}
// Composite Gauss quadrature at the requested precision for all difference-grid entries at once.
// Expensive Fourier multipliers are shared; phase recurrences restart every
// 16 entries to limit drift. Two successive mesh agreements are required.
template <uni20::Real Real>
KernelTable<Real> kernel_table(Real p, Real eta, Real h, std::size_t n, KernelOptions<Real> options,
                               std::size_t& total_evaluations)
{
  using C = uni20::complex<Real>;
  using bethe::detail::CompensatedSum;
  Real const zero = uni20::scalar_at_precision<Real>(0, options.precision()), one = uni20::scalar_like(zero, 1);
  C const czero(zero, zero);
  KernelTable<Real> out{std::vector<Real>(n + 1, zero), std::vector<C>(n + 1, czero),
                        uni20::scalar_like(zero, std::numeric_limits<double>::infinity())};
  if (p == Real{1})
  {
    out.error = zero;
    out.converged = true;
    return out;
  }
  Real const pi = Real{4} * uni20::math::atan(one), width = pi * std::min(p, Real{1}) - Real{2} * eta;
  Real cutoff = one / width, tail = zero;
  bool tail_ok = false;
  for (std::size_t i = 0; i < options.max_cutoffs; ++i)
  {
    if (!uni20::isfinite(cutoff)) return out;
    tail = uni20::math::exp(-width * cutoff) / pi / width / (-uni20::math::expm1(-p * pi * cutoff));
    if (tail <= options.tolerance / Real{4})
    {
      tail_ok = true;
      break;
    }
    cutoff *= Real{2};
  }
  if (!tail_ok) return out;
  // Higher arithmetic precision benefits from increasing the existing Gauss
  // order before refining panels. Native rules and all work-budget semantics
  // are unchanged; every sampled node still counts as one Fourier evaluation.
  std::size_t const order = batch_fourier_order<Real>(options.precision());
  if (order > options.max_evaluations || options.max_levels == 0) return out;
  auto const rule = bethe::detail::gauss_legendre<Real>(order, options.precision());
  unsigned stable = 0;
  std::size_t used = 0, panels = 1;
  for (std::size_t level = 0; level < options.max_levels; ++level)
  {
    if (panels > (options.max_evaluations - used) / order) return out;
    std::vector<CompensatedSum<Real>> re(n + 1, CompensatedSum<Real>(zero));
    std::vector<CompensatedSum<C>> shifted(n + 1, CompensatedSum<C>(czero));
    CompensatedSum<Real> absolute(zero);
    Real const panel_width = cutoff / Real(panels);
    for (std::size_t panel = 0; panel < panels; ++panel)
      for (std::size_t q = 0; q < order; ++q)
      {
        ++used;
        ++total_evaluations;
        Real const k = panel_width * (Real(panel) + (Real{1} + rule.x[q]) / Real{2});
        Real const w = panel_width * rule.w[q] / Real{2};
        auto const a = fourier_components(k, zero, p), b = fourier_components(k, Real{2} * eta, p);
        absolute.add(w * (uni20::math::abs(a.real()) + uni20::math::abs(b.real()) + uni20::math::abs(b.imag())));
        C phase(one, zero), step(uni20::math::cos(k * h), uni20::math::sin(k * h));
        for (std::size_t j = 0; j <= n; ++j)
        {
          if (j % 16 == 0) phase = C(uni20::math::cos(k * h * Real(j)), uni20::math::sin(k * h * Real(j)));
          re[j].add(w * a.real() * phase.real());
          shifted[j].add(w * C(b.real() * phase.real(), -b.imag() * phase.imag()));
          phase *= step;
        }
      }
    Real difference = zero;
    for (std::size_t j = 0; j <= n; ++j)
    {
      if (!uni20::isfinite(re[j].value()) || !uni20::isfinite(shifted[j].value().real()) ||
          !uni20::isfinite(shifted[j].value().imag()))
        return out;
      difference = std::max({difference, uni20::math::abs(re[j].value() - out.real[j]),
                             uni20::math::abs(shifted[j].value() - out.shifted[j])});
      out.real[j] = re[j].value();
      out.shifted[j] = shifted[j].value();
    }
    out.error = tail + std::max(difference, Real{64} * uni20::numeric_limits<Real>::epsilon(options.precision()) *
                                                absolute.value());
    stable = level > 0 && uni20::isfinite(out.error) && out.error <= options.tolerance ? stable + 1 : 0;
    if (stable == 2)
    {
      out.converged = true;
      return out;
    }
    if (panels > options.max_evaluations / (2 * order)) return out;
    panels *= 2;
  }
  return out;
}
template <uni20::Real Real> struct VacuumMesh
{
    VacuumMesh()
      requires(!uni20::has_runtime_precision_v<Real>)
        : VacuumMesh(Real{})
    {}
    explicit VacuumMesh(Real const& zero)
        : value(zero), residual(uni20::scalar_like(zero, std::numeric_limits<double>::infinity())),
          counting_integral(zero), kernel_error(residual)
    {}
    Real value, residual, counting_integral, kernel_error;
    bool converged = false;
    VacuumStatus status = VacuumStatus::iteration_limit;
};
template <uni20::Real Real>
VacuumMesh<Real>
nlie_mesh(Real p, Real u, Real eta, Real cutoff, std::size_t n, VacuumOptions<Real> const& options,
          VacuumState<Real>& work, KernelTable<Real> const& kernel, std::span<uni20::complex<Real> const> source = {},
          std::span<uni20::complex<Real> const> counting_kernel = {}, std::vector<uni20::complex<Real>>* seed = nullptr)
{
  using C = uni20::complex<Real>;
  using bethe::detail::CompensatedSum;
  Real const zero = uni20::scalar_at_precision<Real>(0, options.precision()), one = uni20::scalar_like(zero, 1);
  C const czero(zero, zero);
  VacuumMesh<Real> out(zero);
  Real const h = Real{2} * cutoff / Real(n), pi = Real{4} * uni20::math::atan(one);
  if (!uni20::isfinite(h) || !(h > Real{0}))
  {
    out.status = VacuumStatus::precision_limit;
    return out;
  }
  out.kernel_error = kernel.error;
  if (!kernel.converged)
  {
    out.status = VacuumStatus::kernel_limit;
    return out;
  }
  std::vector<C> bare(n + 1, czero), correction(n + 1, czero), a(n + 1, czero), next(n + 1, czero);
  // A previous converged solution on THIS mesh may seed a nearby hole trial.
  // It is still checked against the new source and the full residual below.
  if (seed && seed->size() == n + 1) correction = *seed;
  std::vector<Real> weighted_real(n + 1, zero), weighted_imag(n + 1, zero);
  std::vector<Real> weight(n + 1, h);
  weight.front() /= Real{2};
  weight.back() /= Real{2};
  for (std::size_t i = 0; i <= n; ++i)
  {
    Real const x = h * (Real(i) - Real(n) / Real{2});
    bare[i] = C(u * uni20::math::sin(eta) * uni20::math::cosh(x), -u * uni20::math::cos(eta) * uni20::math::sinh(x));
  }
  for (;;)
  {
    CompensatedSum<Real> energy(zero), counting(zero);
    for (std::size_t j = 0; j <= n; ++j)
    {
      C const epsilon = bare[j] + correction[j] + (source.empty() ? czero : source[j]);
      if (!uni20::isfinite(epsilon.real()) || !uni20::isfinite(epsilon.imag()))
      {
        out.status = VacuumStatus::precision_limit;
        return out;
      }
      C const fugacity = uni20::math::exp(-epsilon);
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
      if (!counting_kernel.empty()) counting.add(Real{2} * weight[j] * (counting_kernel[j] * a[j]).imag());
    }
    out.value = energy.value();
    out.counting_integral = counting.value();
    out.residual = zero;
    if (p == Real{1} && uni20::isfinite(out.value))
    {
      out.converged = true;
      out.status = VacuumStatus::converged;
      return out;
    }
    // Even-state parity epsilon(-theta)=conj(epsilon(theta)). Evaluate half
    // the convolution, using real components to avoid redundant products.
    for (std::size_t i = 0; i <= n / 2; ++i)
    {
      CompensatedSum<Real> real_sum(zero), imag_sum(zero);
      for (std::size_t j = 0; j <= n; ++j)
      {
        auto const d = i > j ? i - j : j - i;
        Real const kr = kernel.shifted[d].real(), ki = i >= j ? kernel.shifted[d].imag() : -kernel.shifted[d].imag();
        real_sum.add((kr - kernel.real[d]) * weighted_real[j] + ki * weighted_imag[j]);
        imag_sum.add(ki * weighted_real[j] - (kr + kernel.real[d]) * weighted_imag[j]);
      }
      next[i] = C(real_sum.value(), i == n / 2 ? zero : imag_sum.value());
      next[n - i] = uni20::math::conj(next[i]);
      Real const residual = uni20::math::abs(next[i] - correction[i]);
      if (!uni20::isfinite(residual) || !uni20::isfinite(out.value))
      {
        out.status = VacuumStatus::precision_limit;
        return out;
      }
      out.residual = std::max(out.residual, residual);
    }
    if (out.residual <= options.tolerance / (Real{32} * (Real{1} + cutoff)))
    {
      if (seed) *seed = correction;
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
VacuumMesh<Real> vacuum_mesh(Real p, Real u, Real eta, Real cutoff, std::size_t n, VacuumOptions<Real> const& options,
                             VacuumState<Real>& work)
{
  Real const h = Real{2} * cutoff / Real(n);
  if (!uni20::isfinite(h) || !(h > Real{0}) || !uni20::isfinite(u * uni20::math::cosh(cutoff)))
  {
    VacuumMesh<Real> failed(uni20::scalar_at_precision<Real>(0, options.precision()));
    failed.status = VacuumStatus::precision_limit;
    return failed;
  }
  auto kernel_options = options.kernel;
  kernel_options.tolerance = options.tolerance / (Real{64} * (Real{1} + cutoff));
  auto const kernel = kernel_table(p, eta, h, n, kernel_options, work.kernel_evaluations);
  return nlie_mesh(p, u, eta, cutoff, n, options, work, kernel);
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
    out.mesh_error =
        previous ? uni20::math::abs(mesh.value - *previous)
                 : uni20::scalar_at_precision<Real>(std::numeric_limits<double>::infinity(), options.precision());
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
    out.cutoff_error =
        previous ? uni20::math::abs(*value - *previous)
                 : uni20::scalar_at_precision<Real>(std::numeric_limits<double>::infinity(), options.precision());
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
  Real const zero = uni20::scalar_at_precision<Real>(0, options.precision()), one = uni20::scalar_like(zero, 1);
  for (auto const& value : {mass, length, p, options.tolerance})
    bethe::detail::require_same_precision(value, zero);
  if (options.contour_shift) bethe::detail::require_same_precision(*options.contour_shift, zero);
  if (options.initial_cutoff) bethe::detail::require_same_precision(*options.initial_cutoff, zero);
  if (options.kernel.precision() != options.precision())
    throw std::invalid_argument("sine-Gordon vacuum kernel controls require matching precision");
  bethe::detail::require_same_precision(options.kernel.tolerance, zero);
  Real const pi = Real{4} * uni20::math::atan(one), strip = pi * std::min(p, Real{1});
  Real const eta = options.contour_shift.value_or(strip / Real{4});
  if (!uni20::isfinite(mass) || !(mass > Real{0}) || !uni20::isfinite(length) || !(length > Real{0}) ||
      !uni20::isfinite(p) || !(p > Real{0}) || !uni20::isfinite(eta) || !(eta > Real{0} && eta < strip / Real{2}) ||
      !uni20::isfinite(options.tolerance) || !(options.tolerance > Real{0}) || options.initial_intervals < 8 ||
      options.initial_intervals % 2 || options.max_intervals < options.initial_intervals ||
      (options.initial_cutoff && (!uni20::isfinite(*options.initial_cutoff) || !(*options.initial_cutoff > Real{0}))))
    throw std::invalid_argument("invalid sine-Gordon vacuum scales, coupling, contour, tolerance or mesh controls");
  if (options.max_intervals > 16384)
    throw std::length_error("sine-Gordon vacuum mesh exceeds 16384-interval work ceiling");
  VacuumState<Real> out(zero);
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
  Real const second_eta = eta * (Real{3} * one / Real{4});
  out.verification_contour_shift = second_eta;
  Real const tail_scale = (-uni20::math::log(std::min(options.tolerance,
                                                      [&] {
                                                        if constexpr (uni20::has_runtime_precision_v<Real>)
                                                          return one / Real{10};
                                                        else
                                                          return Real{0.1};
                                                      }())) +
                           Real{8}) /
                          out.scaled_length / uni20::math::sin(second_eta);
  Real const cutoff = options.initial_cutoff.value_or(uni20::math::acosh(std::max(Real{2} * one, tail_scale)));
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
  out.contour_error = uni20::math::abs(*first - *second);
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
