// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/carlson_rf.hpp>
#include <bethe/xxz_dispersion.hpp>
#include <uni20/core/elliptic_theta_jet.hpp>

namespace bethe::experimental::xxz
{
/// Complementary-modulus kinematics. Do NOT form 1-k'^2 then subtract it
/// from one: the small mass is essential when resolving soft spinons.
template <uni20::Real R> class SpectralKinematics {
  public:
    explicit SpectralKinematics(R delta, R exchange = R{1})
        : zero_(bethe::detail::precision_zero(delta)), pi_(bethe::detail::spinon_pi(delta)), band_(delta, exchange),
          eta_(uni20::math::acosh(delta))
    {
      if (!uni20::isfinite(delta) || delta <= R{1}) throw std::invalid_argument("spectral kinematics requires Delta>1");
      kp_ = band_.gap() / band_.maximum_energy();
      if (!(kp_ > R{0} && kp_ < R{1}))
        throw std::underflow_error("XXZ mass or bandwidth unresolved; increase precision");
      auto const theta = uni20::elliptic_theta_constant(3, eta_ / pi());
      K_ = pi() / R{2} * uni20::math::exp(R{2} * theta.log_scale) * theta.scaled_value * theta.scaled_value;
    }
    R gap() const { return band_.gap(); }
    R maximum_energy() const { return band_.maximum_energy(); }
    R complementary_modulus() const { return kp_; }
    R energy(R p) const
    {
      validate_momentum(p);
      return band_.energy(p);
    }
    R zero() const { return zero_; }
    uni20::Precision precision() const { return uni20::precision_of(zero_); }
    R jacobian(R p) const { return band_.maximum_energy() / (R{2} * K_ * energy(p)); }
    R angle(R p) const
    {
      validate_momentum(p);
      bool const reflect = p > pi() / R{2};
      R const d = reflect ? pi() - p : p;
      if (d == R{0}) return reflect ? -pi() / R{2} : pi() / R{2};
      if (d == pi() / R{2}) return zero_;
      R const s = uni20::math::sin(d), c = uni20::math::cos(d), dn = uni20::math::hypot(s, kp_ * c);
      R const value = pi() / (R{2} * K_) * c * bethe::detail::carlson_rf_unit(s * s, dn * dn);
      return reflect ? -value : value;
    }
    R momentum(R phi) const
    {
      bethe::detail::require_same_precision(phi, zero_);
      if (!uni20::isfinite(phi) || uni20::math::abs(phi) > pi() / R{2})
        throw std::invalid_argument("spectral angle outside [-pi/2,pi/2]");
      if (phi == R{0}) return pi() / R{2};
      if (uni20::math::abs(phi) == pi() / R{2}) return phi > R{0} ? zero_ : pi();
      // DLMF 22.2.4-5: cn/sn=sqrt(k') theta2(phi)/theta1(phi).
      // Work in logarithmic ratios and reflect AFTER finding the small
      // momentum, so near-XXX theta prefactors never under/overflow.
      auto const u = uni20::complex<R>{uni20::math::abs(phi) / pi(), R{0}};
      auto const a = uni20::elliptic_theta_jet(2, u, eta_ / pi());
      auto const b = uni20::elliptic_theta_jet(1, u, eta_ / pi());
      R const log_ratio = uni20::math::log(kp_) / R{2} + a.log_scale - b.log_scale +
                          uni20::math::log(uni20::math::abs(a.derivative[0])) -
                          uni20::math::log(uni20::math::abs(b.derivative[0]));
      R const small = log_ratio > R{0} ? pi() / R{2} - uni20::math::atan(uni20::math::exp(-log_ratio))
                                       : uni20::math::atan(uni20::math::exp(log_ratio));
      R const p = phi > R{0} ? small : pi() - small;
      if (!(p > R{0} && p < pi())) throw std::underflow_error("soft momentum unresolved; increase precision");
      return p;
    }

  private:
    R pi() const { return pi_; }
    void validate_momentum(R p) const
    {
      bethe::detail::require_same_precision(p, zero_);
      if (!uni20::isfinite(p) || p < R{0} || p > pi()) throw std::invalid_argument("spinon momentum outside [0,pi]");
    }
    R zero_, pi_;
    bethe::xxz::SpinonDispersion<R> band_;
    R eta_, kp_, K_;
};
} // namespace bethe::experimental::xxz
