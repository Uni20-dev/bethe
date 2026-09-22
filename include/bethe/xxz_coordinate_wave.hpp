// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/xxz_common.hpp>
#include <bit>
#include <complex>

namespace bethe::xxz::detail
{
/// A coordinate Bethe numerator and the sum of absolute permutation terms.
/// The latter is a cancellation diagnostic, NOT a norm or an error bound.
template <uni20::Real Real> struct CoordinateAmplitude
{
    std::complex<Real> value{};
    Real absolute_term_sum{};
};

/// Coordinate Bethe wavefunction from finite, nonzero v_j=exp(i*k_j).
/// No Bethe-equation, nonzero-state, or boundary-condition claim is made.
/// Pair normalization is global (independent of the occupied sites).
/// See docs/xxz-coordinate-wave.md for the subset recurrence and limitations.
template <uni20::Real Real> class CoordinateBetheWave {
  public:
    using Complex = std::complex<Real>;

    CoordinateBetheWave(std::size_t sites, Real delta, std::span<Complex const> momentum_factors,
                        std::size_t max_subsets = 65536)
        : sites_(sites)
    {
      checked_sites(sites);
      auto const r = momentum_factors.size();
      if (r > sites || !uni20::isfinite(delta))
        throw std::invalid_argument("coordinate wave requires r <= N and finite Delta");
      // Check the shift and budget before copying roots or allocating scratch.
      if (r >= std::numeric_limits<std::size_t>::digits || (std::size_t{1} << r) > max_subsets)
        throw std::length_error("coordinate wave exceeds the subset budget");
      subsets_ = std::size_t{1} << r;
      for (auto v : momentum_factors)
        if (!finite(v) || v == Complex{})
          throw std::invalid_argument("coordinate wave momentum factors must be finite and nonzero");
      factors_.assign(momentum_factors.begin(), momentum_factors.end());
      pairs_.resize(r * r);
      for (std::size_t i = 0; i < r; ++i)
        for (std::size_t j = i + 1; j < r; ++j)
        {
          Complex const a = Real{1} + factors_[i] * factors_[j] - Real{2} * delta * factors_[i];
          Complex const b = Real{1} + factors_[i] * factors_[j] - Real{2} * delta * factors_[j];
          Real const scale = std::max({Real{1}, std::abs(a), std::abs(b)});
          if (!finite(a) || !finite(b) || !uni20::isfinite(scale))
            throw std::overflow_error("nonfinite coordinate-wave pair factor");
          pairs_[i * r + j] = a / scale;
          pairs_[j * r + i] = b / scale;
        }
    }

    std::size_t subsets_per_amplitude() const { return subsets_; }

    CoordinateAmplitude<Real> evaluate(std::span<std::size_t const> occupied) const
    {
      auto const r = factors_.size();
      if (occupied.size() != r) throw std::invalid_argument("coordinate wave configuration has the wrong sector size");
      for (std::size_t j = 0; j < r; ++j)
        if (occupied[j] >= sites_ || (j && occupied[j - 1] >= occupied[j]))
          throw std::invalid_argument("coordinate wave sites must be strictly increasing in [0,N)");
      std::vector<Complex> planes(r * r);
      for (std::size_t a = 0; a < r; ++a)
        for (std::size_t j = 0; j < r; ++j)
          planes[a * r + j] = power(factors_[j], occupied[a]);
      std::vector<CoordinateAmplitude<Real>> sums(subsets_);
      sums[0] = {Complex{1}, Real{1}};
      for (std::size_t mask = 1; mask < subsets_; ++mask)
      {
        auto const a = std::popcount(mask) - 1;
        bethe::detail::CompensatedSum<Real> real, imag, absolute;
        for (auto remaining = mask; remaining; remaining &= remaining - 1)
        {
          auto const j = std::countr_zero(remaining);
          auto const previous = mask ^ (std::size_t{1} << j);
          Complex weight = planes[a * r + j];
          Real magnitude = std::abs(weight);
          for (auto others = previous; others; others &= others - 1)
          {
            auto const i = std::countr_zero(others);
            weight *= pairs_[i * r + j];
            magnitude *= std::abs(pairs_[i * r + j]);
          }
          if (std::popcount(previous >> j) % 2) weight = -weight;
          Complex const term = weight * sums[previous].value;
          Real const absolute_term = magnitude * sums[previous].absolute_term_sum;
          if (!finite(term) || !uni20::isfinite(absolute_term))
            throw std::overflow_error("nonfinite coordinate-wave term");
          real.add(term.real());
          imag.add(term.imag());
          absolute.add(absolute_term);
        }
        sums[mask] = {{real.value(), imag.value()}, absolute.value()};
        if (!finite(sums[mask].value) || !uni20::isfinite(sums[mask].absolute_term_sum))
          throw std::overflow_error("nonfinite coordinate-wave sum");
      }
      return sums.back();
    }

  private:
    static bool finite(Complex z) { return uni20::isfinite(z.real()) && uni20::isfinite(z.imag()); }

    static Complex power(Complex z, std::size_t exponent)
    {
      Complex value{1};
      while (exponent)
      {
        if (exponent & 1U) value *= z;
        exponent >>= 1;
        if (exponent) z *= z;
        if (!finite(value) || !finite(z)) throw std::overflow_error("nonfinite coordinate-wave plane factor");
      }
      return value;
    }

    std::size_t sites_, subsets_;
    std::vector<Complex> factors_, pairs_;
};
} // namespace bethe::xxz::detail
