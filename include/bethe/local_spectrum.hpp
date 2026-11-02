// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/solver.hpp>
#include <span>

namespace bethe
{
/// Discrete local autocorrelation weight: no momentum or Fourier size factor.
template <uni20::Real Real> struct LocalSpectralLine
{
    std::size_t state_id;
    Real gap;
    Real weight;
};

template <uni20::Real Real> struct LocalSpectralMoments
{
    Real weight;
    Real first_moment;
};

/// <0|A|n><n|B|0>, which may be signed or complex. It is not a probability.
template <uni20::Real Real> struct CrossSpectralLine
{
    std::size_t state_id;
    Real gap;
    uni20::complex<Real> coefficient;
};

template <uni20::Real Real> struct CrossSpectralMoments
{
    uni20::complex<Real> value;
    uni20::complex<Real> first_moment;
};

namespace detail
{
template <uni20::Real Real> void validate_spectral_gap(Real const& gap, Real const& zero)
{
  require_same_precision(gap, zero);
  if (!uni20::isfinite(gap) || gap < 0) throw std::invalid_argument("a spectral gap must be finite and nonnegative");
}

template <uni20::Real Real> void validate_spectral_weight(Real const& weight, Real const& zero)
{
  require_same_precision(weight, zero);
  if (!uni20::isfinite(weight) || weight < 0)
    throw std::invalid_argument("an autocorrelation weight must be finite and nonnegative");
}

template <uni20::Real Real> class SpectralComplexSum {
  public:
    explicit SpectralComplexSum(Real const& zero) : real_(zero), imag_(zero) {}
    void add(Real real, Real imag)
    {
      if (!uni20::isfinite(real) || !uni20::isfinite(imag)) throw std::overflow_error("nonfinite spectral term");
      real_.add(real);
      imag_.add(imag);
    }
    uni20::complex<Real> value() const
    {
      Real const real = real_.value(), imag = imag_.value();
      if (!uni20::isfinite(real) || !uni20::isfinite(imag)) throw std::overflow_error("nonfinite spectral sum");
      return {real, imag};
    }

  private:
    CompensatedSum<Real> real_, imag_;
};

template <uni20::Real Real> void validate_cross_coefficient(uni20::complex<Real> const& c, Real const& zero)
{
  // Check the complex container before accessing a possibly unset MPC scalar.
#if UNI20_HAS_MPREAL
  if constexpr (uni20::has_runtime_precision_v<Real>)
    if (!c.initialized() || c.is_exact() || uni20::precision_of(c) != uni20::precision_of(zero))
      throw std::invalid_argument("cross coefficients require the working precision");
#endif
  require_same_precision(c.real(), zero);
  require_same_precision(c.imag(), zero);
  if (!uni20::isfinite(c.real()) || !uni20::isfinite(c.imag()))
    throw std::invalid_argument("a cross coefficient must be finite");
}

template <uni20::Real Real>
void add_spectral_phase(SpectralComplexSum<Real>& sum, Real const& gap, Real const& real, Real const& imag,
                        Real const& time)
{
  Real const phase = gap * time;
  if (!uni20::isfinite(phase)) throw std::overflow_error("nonfinite spectral phase");
  Real const cosine = uni20::math::cos(phase), sine = uni20::math::sin(phase);
  sum.add(real * cosine + imag * sine, imag * cosine - real * sine);
}

template <uni20::Real Real> Real spectral_time_zero(Real const& time)
{
  auto zero = precision_zero(time);
  if (!uni20::isfinite(time)) throw std::invalid_argument("correlation time must be finite");
  return zero;
}
} // namespace detail

/// Sum supplied lines only. Runtime scalars require an explicit finite-precision
/// exemplar, including for empty spectra; exact/unset/mixed inputs are rejected.
template <uni20::Real Real>
[[nodiscard]] LocalSpectralMoments<Real> local_spectral_moments(std::span<LocalSpectralLine<Real> const> lines,
                                                                Real const& exemplar = Real{0})
{
  Real const zero = detail::precision_zero(exemplar);
  detail::CompensatedSum<Real> weight(zero), first(zero);
  for (auto const& line : lines)
  {
    detail::validate_spectral_gap(line.gap, zero);
    detail::validate_spectral_weight(line.weight, zero);
    Real const term = line.gap * line.weight;
    if (!uni20::isfinite(term)) throw std::overflow_error("nonfinite local spectral moment");
    weight.add(line.weight);
    first.add(term);
  }
  if (!uni20::isfinite(weight.value()) || !uni20::isfinite(first.value()))
    throw std::overflow_error("nonfinite local spectral moment sum");
  return {weight.value(), first.value()};
}

template <uni20::Real Real>
[[nodiscard]] CrossSpectralMoments<Real> cross_spectral_moments(std::span<CrossSpectralLine<Real> const> lines,
                                                                Real const& exemplar = Real{0})
{
  Real const zero = detail::precision_zero(exemplar);
  detail::SpectralComplexSum<Real> value(zero), first(zero);
  for (auto const& line : lines)
  {
    detail::validate_spectral_gap(line.gap, zero);
    detail::validate_cross_coefficient<Real>(line.coefficient, zero);
    value.add(line.coefficient.real(), line.coefficient.imag());
    first.add(line.gap * line.coefficient.real(), line.gap * line.coefficient.imag());
  }
  return {value.value(), first.value()};
}

/// C(t)=sum w*exp(-i*gap*t). No broadening, time stepping, 2*pi or rescaling.
/// The time scalar supplies the working precision, even for an empty spectrum.
template <uni20::Real Real>
[[nodiscard]] uni20::complex<Real> local_time_correlation(std::span<LocalSpectralLine<Real> const> lines,
                                                          Real const& time)
{
  Real const zero = detail::spectral_time_zero(time);
  detail::SpectralComplexSum<Real> sum(zero);
  for (auto const& line : lines)
  {
    detail::validate_spectral_gap(line.gap, zero);
    detail::validate_spectral_weight(line.weight, zero);
    detail::add_spectral_phase(sum, line.gap, line.weight, zero, time);
  }
  return sum.value();
}

/// Complex cross products need not obey the autocorrelation identity
/// C(-t)=conj(C(t)). For Hermitian A and B, conjugating the coefficients
/// reverses the operator pair.
template <uni20::Real Real>
[[nodiscard]] uni20::complex<Real> cross_time_correlation(std::span<CrossSpectralLine<Real> const> lines,
                                                          Real const& time)
{
  Real const zero = detail::spectral_time_zero(time);
  detail::SpectralComplexSum<Real> sum(zero);
  for (auto const& line : lines)
  {
    detail::validate_spectral_gap(line.gap, zero);
    detail::validate_cross_coefficient<Real>(line.coefficient, zero);
    detail::add_spectral_phase(sum, line.gap, line.coefficient.real(), line.coefficient.imag(), time);
  }
  return sum.value();
}
} // namespace bethe
