// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/local_spectrum.hpp>
#include <optional>
#include <span>
#include <stdexcept>
#include <vector>

namespace bethe
{
/// A discrete contribution to S(q,w)=2*pi*sum weight*delta(w-gap).
/// Broadening and any operator/channel normalization belong to the caller.
template <uni20::Real Real> struct LehmannLine
{
    std::size_t state_id;
    std::size_t momentum_index;
    Real momentum;
    Real gap;
    Real weight;
};

template <uni20::Real Real> struct SpectralMoment
{
    Real weight = Real{0};
    Real first_moment = Real{0};
};

/// Sum accepted lines only. A partial family must not be rescaled to a full sum
/// rule. Model-specific full sum rules and scan failures live outside this layer.
/// Supply a finite-precision zero for runtime scalars, including empty spectra.
template <uni20::Real Real>
std::vector<SpectralMoment<Real>> spectral_moments(std::size_t momenta, std::span<LehmannLine<Real> const> lines,
                                                   Real const& zero = Real{0})
{
  detail::CompensatedSum<Real> const initial(detail::precision_zero(zero));
  std::vector<detail::CompensatedSum<Real>> weights(momenta, initial), first(momenta, initial);
  for (auto const& line : lines)
  {
    detail::require_same_precision(line.gap, zero);
    detail::require_same_precision(line.weight, zero);
    if (line.momentum_index >= momenta || !uni20::isfinite(line.gap) || line.gap < Real{0} ||
        !uni20::isfinite(line.weight) || line.weight < Real{0})
      throw std::invalid_argument("invalid Lehmann line");
    Real const first_moment = line.gap * line.weight;
    if (!uni20::isfinite(first_moment)) throw std::overflow_error("nonfinite Lehmann first moment");
    weights[line.momentum_index].add(line.weight);
    first[line.momentum_index].add(first_moment);
  }
  std::vector<SpectralMoment<Real>> result(momenta);
  for (std::size_t q = 0; q < momenta; ++q)
  {
    if (!uni20::isfinite(weights[q].value()) || !uni20::isfinite(first[q].value()))
      throw std::overflow_error("nonfinite Lehmann moment sum");
    result[q] = {weights[q].value(), first[q].value()};
  }
  return result;
}

/// One finite Lehmann time sum per momentum index, without broadening or
/// rescaling an incomplete family. The time scalar supplies working precision.
template <uni20::Real Real>
[[nodiscard]] std::vector<uni20::complex<Real>>
spectral_time_correlation(std::size_t momenta, std::span<LehmannLine<Real> const> lines, Real const& time)
{
  Real const zero = detail::spectral_time_zero(time);
  std::vector<detail::SpectralComplexSum<Real>> sums(momenta, detail::SpectralComplexSum<Real>(zero));
  for (auto const& line : lines)
  {
    if (line.momentum_index >= momenta) throw std::invalid_argument("invalid Lehmann momentum index");
    detail::validate_spectral_gap(line.gap, zero);
    detail::validate_spectral_weight(line.weight, zero);
    detail::add_spectral_phase(sums[line.momentum_index], line.gap, line.weight, zero, time);
  }
  std::vector<uni20::complex<Real>> result;
  result.reserve(momenta);
  for (auto const& sum : sums)
    result.push_back(sum.value());
  return result;
}
} // namespace bethe
