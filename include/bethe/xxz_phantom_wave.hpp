// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/xxz_phantom.hpp>
#include <numeric>

namespace bethe::xxz::detail
{
// Compute C(n,k) without overflowing an intermediate product. The caller's
// work limit is checked BEFORE any wavefunction callback is invoked.
inline std::size_t phantom_dressing_terms(std::size_t n, std::size_t k, std::size_t limit)
{
  if (k > n) throw std::invalid_argument("phantom term count requires k <= n");
  k = std::min(k, n - k);
  std::size_t count = 1;
  if (count > limit) throw std::length_error("phantom dressing exceeds the term budget");
  for (std::size_t j = 1; j <= k; ++j)
  {
    auto const divisor = std::gcd(count, j);
    auto const reduced = count / divisor;
    auto const factor = (n - k + j) / (j / divisor);
    if (reduced > limit / factor) throw std::length_error("phantom dressing exceeds the term budget");
    count = reduced * factor;
  }
  return count;
}

/// Coordinate-space map adding p one-sided phantom roots to an r-spin
/// amplitude callback. It can have a kernel: a nonzero input does NOT imply
/// a nonzero output. No eigenstate or ground-state claim is made here.
/// See docs/xxz-phantom-wave.md for the intertwining identity and cost.
template <uni20::Real Real> class PhantomDressing {
  public:
    using Complex = std::complex<Real>;

    PhantomDressing(std::size_t sites, std::size_t finite_count, std::size_t phantom_count, Complex phase,
                    std::size_t max_terms = 100000)
        : sites_(sites), finite_count_(finite_count), phantom_count_(phantom_count), phase_(phase)
    {
      checked_sites(sites);
      if (finite_count > sites || phantom_count > sites - finite_count)
        throw std::invalid_argument("phantom dressing requires r+p <= N");
      if (!finite(phase) || std::abs(std::abs(phase) - Real{1}) > Real{64} * uni20::numeric_limits<Real>::epsilon())
        throw std::invalid_argument("phantom dressing phase must be a finite unit complex number");
      terms_ = phantom_dressing_terms(finite_count + phantom_count, phantom_count, max_terms);
    }

    std::size_t terms_per_amplitude() const { return terms_; }

    /// Zero is required for the periodic/twisted intertwining identity.
    /// The map itself remains defined away from this condition.
    Real commensurability_error() const
    {
      auto const twice = 2 * finite_count_;
      auto const exponent = twice <= sites_ ? sites_ - twice : twice - sites_;
      return std::abs(phantom_phase_power(phase_, exponent) - Complex{1});
    }

    /// Visit each finite subset and its coefficient. The span is temporary
    /// scratch storage: the visitor must not retain it. Invalid coordinates
    /// are rejected before the first visit. This does not allocate a spin
    /// basis, matrix, or full wavefunction.
    template <typename Visitor> void for_each_term(std::span<std::size_t const> occupied, Visitor&& visitor) const
    {
      auto const m = finite_count_ + phantom_count_;
      if (occupied.size() != m) throw std::invalid_argument("phantom dressing configuration has the wrong sector size");
      for (std::size_t j = 0; j < m; ++j)
        if (occupied[j] >= sites_ || (j && occupied[j - 1] >= occupied[j]))
          throw std::invalid_argument("phantom dressing sites must be strictly increasing in [0,N)");
      std::vector<std::size_t> phantom(phantom_count_), selected;
      std::iota(phantom.begin(), phantom.end(), std::size_t{0});
      selected.reserve(finite_count_);
      for (std::size_t term = 0; term < terms_; ++term)
      {
        Complex coefficient{1};
        for (std::size_t h = 0; h < phantom_count_; ++h)
        {
          auto const rank = phantom[h];
          // Each phantom crosses r+h-rank finite particles to its right.
          // The exponent is <3*N, so checked_sites prevents overflow.
          auto const exponent = occupied[rank] + 2 * (finite_count_ + h - rank);
          coefficient *= phantom_phase_power(phase_, exponent);
        }
        if (!finite(coefficient)) throw std::overflow_error("nonfinite phantom dressing coefficient");
        selected.clear();
        for (std::size_t j = 0, h = 0; j < m; ++j)
          if (h < phantom_count_ && phantom[h] == j)
            ++h;
          else
            selected.push_back(occupied[j]);
        visitor(std::span<std::size_t const>(selected), coefficient);
        if (term + 1 == terms_) break;
        // Lexicographic combinations, with no recursion proportional to N.
        auto j = phantom_count_;
        while (j && phantom[j - 1] == finite_count_ + j - 1)
          --j;
        ++phantom[j - 1]; // Another term exists, hence j>0.
        for (std::size_t h = j; h < phantom_count_; ++h)
          phantom[h] = phantom[h - 1] + 1;
      }
    }

    template <typename Amplitude>
    Complex amplitude(std::span<std::size_t const> occupied, Amplitude&& finite_amplitude) const
    {
      bethe::detail::CompensatedSum<Real> real, imag;
      for_each_term(occupied, [&](std::span<std::size_t const> selected, Complex coefficient) {
        Complex const value = finite_amplitude(selected);
        if (!finite(value)) throw std::overflow_error("nonfinite reduced-state amplitude");
        Complex const term = coefficient * value;
        if (!finite(term)) throw std::overflow_error("nonfinite dressed-state amplitude term");
        real.add(term.real());
        imag.add(term.imag());
      });
      Complex const result{real.value(), imag.value()};
      if (!finite(result)) throw std::overflow_error("nonfinite dressed-state amplitude sum");
      return result;
    }

  private:
    static bool finite(Complex value) { return uni20::isfinite(value.real()) && uni20::isfinite(value.imag()); }
    std::size_t sites_, finite_count_, phantom_count_, terms_;
    Complex phase_;
};
} // namespace bethe::xxz::detail
