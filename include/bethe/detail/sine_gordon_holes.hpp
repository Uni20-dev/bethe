// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/sine_gordon_vacuum.hpp>

namespace bethe::sine_gordon::detail
{
template <uni20::Real Real> struct HoleTable
{
    std::vector<uni20::complex<Real>> source, counting;
    HoleTable(std::size_t size, Real const& zero)
        : source(size, uni20::complex<Real>(zero, zero)), counting(size, uni20::complex<Real>(zero, zero)), phase(zero),
          phase_derivative(zero), quadrature_error(uni20::scalar_like(zero, std::numeric_limits<double>::infinity())),
          tail_bound(quadrature_error), error(quadrature_error)
    {}
    Real phase, phase_derivative, quadrature_error, tail_bound, error;
    bool converged = false;
};

// Fixed real holes +/-H. On x+i*eta, source=-i[chi(x-H+i*eta)+chi(x+H+i*eta)].
// Counting kernel is G(H-x-i*eta), to reconstruct the REAL-axis Z(H).
// Batch Fourier quadrature shares expensive factors and trigonometric steps;
// this is not interpolation of a precomputed, lower-precision phase.
template <uni20::Real Real>
HoleTable<Real> hole_table(Real p, Real eta, Real cutoff, std::size_t n, Real hole, KernelOptions<Real> options,
                           std::size_t& total_evaluations)
{
  using C = uni20::complex<Real>;
  using bethe::detail::CompensatedSum;
  Real const zero = uni20::scalar_at_precision<Real>(0, options.precision()), one = uni20::scalar_like(zero, 1);
  C const czero(zero, zero);
  HoleTable<Real> out(n + 1, zero);
  if (p == Real{1})
  {
    out.error = out.quadrature_error = out.tail_bound = zero;
    out.converged = true;
    return out;
  }
  Real const pi = Real{4} * uni20::math::atan(one), width = pi * std::min(p, Real{1}) - eta;
  Real const h = Real{2} * cutoff / Real(n);
  Real bound = zero, kmax = one / width;
  bool tail_ok = false;
  for (std::size_t attempt = 0; attempt < options.max_cutoffs; ++attempt)
  {
    if (!uni20::isfinite(kmax)) return out;
    Real const tail = uni20::math::exp(-width * kmax) / (pi * width * (-uni20::math::expm1(-p * pi * kmax)));
    bound = tail * std::max({one, Real{4} * pi / kmax, Real{2} * pi});
    out.tail_bound = bound;
    if (bound <= options.tolerance / Real{4})
    {
      tail_ok = true;
      break;
    }
    kmax *= Real{2};
  }
  if (!tail_ok) return out;
  std::size_t const order = batch_fourier_order<Real>(options.precision());
  if (order > options.max_evaluations || options.max_levels == 0) return out;
  auto const rule = bethe::detail::gauss_legendre<Real>(order, options.precision());
  unsigned stable = 0;
  std::size_t panels = 1, used = 0;
  for (std::size_t level = 0; level < options.max_levels; ++level)
  {
    if (panels > (options.max_evaluations - used) / order) return out;
    std::vector<CompensatedSum<C>> sources(n + 1, CompensatedSum<C>(czero)), kernels(n + 1, CompensatedSum<C>(czero));
    CompensatedSum<Real> phase(zero), derivative(zero), absolute(zero);
    Real const panel_width = kmax / Real(panels);
    for (std::size_t panel = 0; panel < panels; ++panel)
      for (std::size_t q = 0; q < order; ++q)
      {
        ++used;
        ++total_evaluations;
        Real const k = panel_width * (Real(panel) + (Real{1} + rule.x[q]) / Real{2});
        Real const w = panel_width * rule.w[q] / Real{2};
        auto const a = fourier_components(k, zero, p), b = fourier_components(k, eta, p);
        Real const ch = uni20::math::cos(k * hole), sh = uni20::math::sin(k * hole);
        phase.add(w * Real{2} * pi * a.real() * uni20::math::sin(Real{2} * k * hole) / k);
        derivative.add(w * Real{2} * pi * a.real() * uni20::math::cos(Real{2} * k * hole));
        absolute.add(w * (uni20::math::abs(a.real()) + uni20::math::abs(b.real()) + uni20::math::abs(b.imag())) *
                     (Real{1} + Real{4} * pi * (Real{1} + cutoff + hole + eta)));
        C angle = czero, step(uni20::math::cos(k * h), uni20::math::sin(k * h));
        for (std::size_t j = 0; j <= n; ++j)
        {
          if (j % 16 == 0)
          {
            Real const x = h * (Real(j) - Real(n) / Real{2});
            angle = C(uni20::math::cos(k * x), uni20::math::sin(k * x));
          }
          sources[j].add(w * Real{4} * pi * ch / k * C(b.imag() * angle.real(), -b.real() * angle.imag()));
          kernels[j].add(w * C(b.real() * (ch * angle.real() + sh * angle.imag()),
                               b.imag() * (sh * angle.real() - ch * angle.imag())));
          angle *= step;
        }
      }
    if (!uni20::isfinite(phase.value()) || !uni20::isfinite(derivative.value())) return out;
    Real difference = std::max(uni20::math::abs(phase.value() - out.phase),
                               uni20::math::abs(derivative.value() - out.phase_derivative));
    for (std::size_t j = 0; j <= n; ++j)
    {
      if (!uni20::isfinite(sources[j].value().real()) || !uni20::isfinite(sources[j].value().imag()) ||
          !uni20::isfinite(kernels[j].value().real()) || !uni20::isfinite(kernels[j].value().imag()))
        return out;
      difference = std::max({difference, uni20::math::abs(sources[j].value() - out.source[j]),
                             uni20::math::abs(kernels[j].value() - out.counting[j])});
      out.source[j] = sources[j].value();
      out.counting[j] = kernels[j].value();
    }
    out.phase = phase.value();
    out.phase_derivative = derivative.value();
    out.quadrature_error =
        std::max(difference, Real{64} * uni20::numeric_limits<Real>::epsilon(options.precision()) * absolute.value());
    out.error = bound + out.quadrature_error;
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
} // namespace bethe::sine_gordon::detail
