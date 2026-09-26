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
      return std::hypot(gap_ * std::cos(p), maximum_ * std::sin(p));
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
        result.lower = std::min(result.lower, shifted.lower);
        result.upper = std::max(result.upper, shifted.upper);
      }
      if (!uni20::isfinite(result.lower) || !uni20::isfinite(result.upper))
        throw std::overflow_error("two-spinon energy exceeds scalar range");
      return result;
    }

  private:
    SpinonContinuum<Real> unfolded(Real q) const
    {
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
} // namespace detail
} // namespace bethe
