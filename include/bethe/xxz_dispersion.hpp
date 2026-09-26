// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/elliptic_theta.hpp>
#include <bethe/detail/spinon_band.hpp>

namespace bethe::xxz
{
/// Zero-field, infinite-chain H=J sum (SxSx+SySy+Delta SzSz), J>0,
/// -1<Delta<infinity. The massive branch is a spin-1/2 domain wall between
/// the two Neel vacua, not a local excitation within one vacuum sector.
/// Caux--Mossel--Perez Castillo, JSTAT (2008) P08006, Eqs. (21)--(23).
template <uni20::Real Real> class SpinonDispersion {
  public:
    explicit SpinonDispersion(Real delta, Real exchange = Real{1}) : band_(make_band(delta, exchange)) {}
    Real gap() const { return band_.gap(); }
    Real maximum_energy() const { return band_.maximum_energy(); }
    Real energy(Real p) const { return band_.energy(p); }
    SpinonContinuum<Real> continuum(Real q, SpinonMomentum convention = SpinonMomentum::unfolded) const
    {
      return band_.continuum(q, convention);
    }

  private:
    static bethe::detail::SpinonBand<Real> make_band(Real delta, Real exchange)
    {
      bethe::detail::validate_spinon_parameters(Real{0}, exchange);
      if (!uni20::isfinite(delta) || delta <= -Real{1})
        throw std::invalid_argument("XXZ spinons require finite Delta>-1");
      Real const pi = bethe::detail::pi<Real>();
      if (delta <= Real{1}) return {spinon_energy(pi / Real{2}, delta, exchange), Real{0}};

      Real const eta = std::acosh(delta), t = eta / pi;
      // q=exp(-eta), K=pi*theta3(q)^2/2, k'=theta4(q)^2/theta3(q)^2.
      // For eta<pi use the modular transform. In particular theta4 ->
      // theta2 at exp(-pi^2/eta): never compute k' by subtracting k^2 from 1.
      bool const dual = t < Real{1};
      Real const theta_t = dual ? Real{1} / t : t;
      auto const top = bethe::detail::elliptic_theta(3, std::complex<Real>{}, theta_t);
      auto const bottom = bethe::detail::elliptic_theta(dual ? 2 : 4, std::complex<Real>{}, theta_t);
      // Reconstruct sinh(eta) algebraically: sinh(acosh(Delta)) loses
      // O(log(Delta)*epsilon) relative accuracy at very large anisotropy.
      Real const sh = delta >= Real{2} ? delta * std::sqrt((Real{1} - Real{1} / delta) * (Real{1} + Real{1} / delta))
                                       : std::sqrt(delta - Real{1}) * std::sqrt(delta + Real{1});
      Real const factor = dual ? (pi / Real{2}) * (sh / eta) : sh / Real{2};
      Real const scale = exchange * factor;
      Real const maximum = scale * top.derivative[0].real() * top.derivative[0].real();
      Real const log_gap = std::log(scale) + Real{2} * (bottom.log_scale + std::log(bottom.derivative[0].real()));
      // Reject subnormal masses: they cannot retain the advertised relative
      // precision. An underflow must never masquerade as a gapless result.
      if (log_gap < std::log(uni20::numeric_limits<Real>::min()))
        throw std::underflow_error("massive XXZ spinon gap is below the normal scalar range; use higher precision");
      Real const gap = bottom.log_scale == Real{0} ? scale * bottom.derivative[0].real() * bottom.derivative[0].real()
                                                   : std::exp(log_gap);
      if (!uni20::isfinite(maximum) || !uni20::isfinite(gap) || !(gap > Real{0}))
        throw std::overflow_error("XXZ spinon energy exceeds scalar range");
      // At enormous anisotropy the bandwidth is below relative precision.
      // Independent rounding of exp(log(gap)) must not invert the band.
      return {maximum, std::min(gap, maximum)};
    }
    bethe::detail::SpinonBand<Real> band_;
};
} // namespace bethe::xxz
