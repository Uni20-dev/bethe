// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/spinon_band.hpp>
#include <uni20/core/elliptic_theta.hpp>

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
      Real const zero = bethe::detail::precision_zero(delta);
#if UNI20_HAS_MPREAL
      if constexpr (uni20::has_runtime_precision_v<Real>)
        if (!exchange.is_exact()) bethe::detail::require_same_precision(exchange, delta);
#endif
      exchange = uni20::scalar_like(delta, exchange);
      bethe::detail::validate_spinon_parameters(zero, exchange);
      if (!uni20::isfinite(delta) || delta <= -Real{1})
        throw std::invalid_argument("XXZ spinons require finite Delta>-1");
      Real const pi = bethe::detail::spinon_pi(delta);
      if (delta <= Real{1}) return {spinon_energy(pi / Real{2}, delta, exchange), zero};

      Real const eta = uni20::math::acosh(delta), t = eta / pi;
      // q=exp(-eta), K=pi*theta3(q)^2/2, k'=theta4(q)^2/theta3(q)^2.
      // For eta<pi use the modular transform. In particular theta4 ->
      // theta2 at exp(-pi^2/eta): never compute k' by subtracting k^2 from 1.
      bool const dual = t < Real{1};
      Real const theta_t = dual ? Real{1} / t : t;
      auto const top = uni20::elliptic_theta_constant(3, theta_t);
      auto const bottom = uni20::elliptic_theta_constant(dual ? 2 : 4, theta_t);
      // Reconstruct sinh(eta) algebraically: sinh(acosh(Delta)) loses
      // O(log(Delta)*epsilon) relative accuracy at very large anisotropy.
      Real const sh = delta >= Real{2}
                          ? delta * uni20::math::sqrt((Real{1} - Real{1} / delta) * (Real{1} + Real{1} / delta))
                          : uni20::math::sqrt(delta - Real{1}) * uni20::math::sqrt(delta + Real{1});
      Real const factor = dual ? (pi / Real{2}) * (sh / eta) : sh / Real{2};
      Real const scale = exchange * factor;
      Real const maximum = scale * top.scaled_value * top.scaled_value;
      Real const log_gap =
          uni20::math::log(scale) + Real{2} * (bottom.log_scale + uni20::math::log(bottom.scaled_value));
      // Reject subnormal masses: they cannot retain the advertised relative
      // precision. An underflow must never masquerade as a gapless result.
      Real const minimum = [&] {
        if constexpr (uni20::has_runtime_precision_v<Real>)
          return uni20::math::next_up(zero);
        else
          return uni20::numeric_limits<Real>::min();
      }();
      if (log_gap < uni20::math::log(minimum))
        throw std::underflow_error("massive XXZ spinon gap is below the normal scalar range; use higher precision");
      Real const gap =
          bottom.log_scale == Real{0} ? scale * bottom.scaled_value * bottom.scaled_value : uni20::math::exp(log_gap);
      if (!uni20::isfinite(maximum) || !uni20::isfinite(gap) || !(gap > Real{0}))
        throw std::overflow_error("XXZ spinon energy exceeds scalar range");
      // At enormous anisotropy the bandwidth is below relative precision.
      // Independent rounding of exp(log(gap)) must not invert the band.
      return {maximum, std::min(gap, maximum)};
    }
    bethe::detail::SpinonBand<Real> band_;
};
} // namespace bethe::xxz
