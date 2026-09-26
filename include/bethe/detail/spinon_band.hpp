// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <algorithm>
#include <bethe/spinon.hpp>

namespace bethe
{
/// Fixed sum p1+p2=Q (mod 2*pi), or the union with its pi-shifted image.
/// The latter is the kinematic continuum seen after two-site-cell folding;
/// neither choice makes a statement about spectral weights.
enum class SpinonMomentum
{
  unfolded,
  folded
};

template <uni20::Real Real> struct SpinonContinuum
{
    Real lower{}, upper{};
};

namespace detail
{
/// Common kinematics for e(p)=hypot(m*cos(p), I*sin(p)), 0<=p<=pi.
/// Useful also for elliptic branches in other spin chains. Parameters must
/// satisfy 0<=m<=I, I>0 and be representable in Real.
template <uni20::Real Real> class SpinonBand {
  public:
    SpinonBand(Real maximum, Real gap) : maximum_(maximum), gap_(gap)
    {
      if (!uni20::isfinite(maximum) || !uni20::isfinite(gap) || maximum <= Real{0} || gap < Real{0} || gap > maximum)
        throw std::invalid_argument("spinon band requires finite 0<=gap<=maximum, maximum>0");
    }
    Real gap() const { return gap_; }
    Real maximum_energy() const { return maximum_; }
    Real energy(Real p) const
    {
      validate_spinon_parameters(p, Real{1});
      // Subtract from pi BEFORE taking the sine: a tiny but nonzero sin(pi)
      // would completely hide the exponentially small massive endpoint gap.
      p = std::min(p, pi<Real>() - p);
      if (p == Real{0}) return gap_;
      if (p == pi<Real>() / Real{2}) return maximum_;
      Real const value = std::hypot(gap_ * std::cos(p), maximum_ * std::sin(p));
      if (value == Real{0}) throw std::underflow_error("positive spinon energy underflows scalar range");
      return value;
    }
    SpinonContinuum<Real> continuum(Real q, SpinonMomentum convention = SpinonMomentum::unfolded) const
    {
      Real const pi_value = pi<Real>();
      if (!uni20::isfinite(q) || q < Real{0} || q > Real{2} * pi_value)
        throw std::invalid_argument("two-spinon momentum must be finite and in [0, 2*pi]");
      if (convention != SpinonMomentum::unfolded && convention != SpinonMomentum::folded)
        throw std::invalid_argument("unknown two-spinon momentum convention");
      q = std::min(q, Real{2} * pi_value - q);
      auto result = unfolded(q);
      if (convention == SpinonMomentum::folded)
      {
        auto const shifted = unfolded(pi_value - q);
        // In the gapless case both images have the same lower edge.
        // Do not lose a tiny positive Q when pi-Q rounds to pi.
        if (gap_ != Real{0}) result.lower = std::min(result.lower, shifted.lower);
        result.upper = std::max(result.upper, shifted.upper);
      }
      if (!uni20::isfinite(result.lower) || !uni20::isfinite(result.upper))
        throw std::overflow_error("two-spinon energy exceeds scalar range");
      return result;
    }

  private:
    SpinonContinuum<Real> unfolded(Real q) const
    {
      if (gap_ == Real{0})
      {
        // Keep a positive subnormal Q intact: Q/2 may round to zero even
        // when the energy is representable after multiplication by I.
        Real const upper_factor = q < uni20::numeric_limits<Real>::min() ? q : Real{2} * std::sin(q / Real{2});
        return {energy(q), maximum_ * upper_factor};
      }
      Real const equal = Real{2} * energy(q / Real{2}), endpoint = gap_ + energy(q);
      Real lower = endpoint;
      if (q <= pi<Real>() / Real{2})
      {
        // cos(Q_kappa)=(I-m)/(I+m). Comparing half-angle squares avoids
        // rounding Q_kappa to zero when m/I is much smaller than epsilon.
        Real const s = std::sin(q / Real{2}), ratio = gap_ / maximum_;
        lower = s * s <= ratio / (Real{1} + ratio) ? equal : (maximum_ + gap_) * std::sin(q);
      }
      // The only maximum candidates are equal momenta and an endpoint.
      // The additional unequal-momentum stationary point is a minimum.
      return {lower, std::max(equal, endpoint)};
    }
    Real maximum_, gap_;
};

/// Four identical gapless spinons, e(p)=I*sin(p), 0<=p<=pi, fixed total
/// Q modulo 2*pi. Shared energy kinematics only: a model must separately
/// specify its allowed spin/fusion sectors and cannot infer spectral weights.
template <uni20::Real Real>
SpinonContinuum<Real> gapless_four_spinon_continuum(Real maximum, Real q,
                                                    SpinonMomentum convention = SpinonMomentum::unfolded)
{
  SpinonBand<Real> const band(maximum, Real{0});
  // Reuse momentum validation and the exact two-spinon lower boundary.
  auto const two = band.continuum(q, convention);
  q = std::min(q, Real{2} * pi<Real>() - q);
  if (convention == SpinonMomentum::folded) q = std::min(q, pi<Real>() - q);
  // Concavity at fixed sum puts the maximum at four equal momenta, with
  // the total-momentum alias nearest 2*pi. Soft spectators give the minimum.
  Real const upper = maximum * (Real{4} * std::cos(q / Real{4}));
  if (!uni20::isfinite(upper)) throw std::overflow_error("four-spinon energy exceeds scalar range");
  return {two.lower, upper};
}
} // namespace detail
} // namespace bethe
