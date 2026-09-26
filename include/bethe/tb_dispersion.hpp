// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/spinon_band.hpp>
#include <uni20/common/half_int.hpp>

namespace bethe::takhtajan_babujian
{
/// Zero-field infinite spin-1 chain H=J sum[S.S-(S.S)^2], J>0.
/// Our J equals J_paper/4 in Vlijm--Caux (2014), Eqs. (1.2),(1.4), Sec. 2.
/// A spinon has spin 1/2, e(p)=2*pi*J*sin(p), 0<=p<=pi. Integer-spin
/// periodic chains require even spinon number; no state-counting or weights
/// follow from the energy envelope. TB has SU(2)_2, not XXX's SU(2)_1 content.
template <uni20::Real Real> class SpinonDispersion {
  public:
    explicit SpinonDispersion(Real exchange = Real{1}) : band_(scale(exchange), Real{0}) {}
    static constexpr uni20::half_int spin() { return uni20::from_twice(1); }
    Real gap() const { return Real{0}; }
    Real velocity() const { return band_.maximum_energy(); }
    Real maximum_energy() const { return band_.maximum_energy(); }
    Real energy(Real p) const { return band_.energy(p); }
    /// Two-spinon spin sectors S=0,1. Local spin response from a singlet
    /// selects S=1; a scalar bond operator can select S=0.
    SpinonContinuum<Real> continuum(Real q, SpinonMomentum convention = SpinonMomentum::unfolded) const
    {
      return band_.continuum(q, convention);
    }
    /// Four-spinon sectors include S=0,1,2. A local quadrupole from a singlet
    /// selects S=2 and cannot create a two-spinon S=2 state.
    SpinonContinuum<Real> four_spinon_continuum(Real q, SpinonMomentum convention = SpinonMomentum::unfolded) const
    {
      return bethe::detail::gapless_four_spinon_continuum(maximum_energy(), q, convention);
    }

  private:
    static Real scale(Real exchange)
    {
      if (!uni20::isfinite(exchange) || exchange <= Real{0})
        throw std::invalid_argument("TB exchange must be finite and positive");
      Real const value = exchange * (Real{2} * bethe::detail::pi<Real>());
      if (!uni20::isfinite(value)) throw std::overflow_error("TB spinon band exceeds scalar range");
      return value;
    }
    bethe::detail::SpinonBand<Real> band_;
};
} // namespace bethe::takhtajan_babujian
