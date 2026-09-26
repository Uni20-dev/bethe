// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/elliptic_theta.hpp>
#include <bethe/detail/spinon_band.hpp>

namespace bethe::xyz
{
/// Infinite-chain excitations for the same theta-ratio Hamiltonian as xyz.hpp:
/// H=J sum(Jx SxSx+Jy SySy+Jz SzSz), S=sigma/2, 0<eta<1, t>0, J>0.
/// Johnson--Krinsky--McCoy (1973), Eqs. (7.8), (7.11), (7.12).
/// A spinon is a domain wall between the two x-Neel vacua; Sz is not conserved.
/// See docs/xyz-dispersion.md for the parity-dependent physical momentum map.
template <uni20::Real Real> class SpinonDispersion {
  public:
    explicit SpinonDispersion(Real eta, Real t, Real exchange = Real{1})
        : eta_(eta), t_(t), band_(make_band(eta, t, exchange))
    {}
    Real gap() const { return band_.gap(); }
    Real maximum_energy() const { return band_.maximum_energy(); }
    Real energy(Real p) const { return band_.energy(p); }

    /// Kinematic edges at fixed p1+p2=q. folded includes the pi-shifted copy.
    /// This does not assign matrix elements or finite-ring multiplicities.
    SpinonContinuum<Real> continuum(Real q, SpinonMomentum convention = SpinonMomentum::unfolded) const
    {
      return band_.continuum(q, convention);
    }

    /// Strict existence: the endpoint s*(1-eta)=eta is a continuum merger,
    /// not an isolated bound branch. No fuzzy tolerance changes this condition.
    bool bound_exists(unsigned s) const { return s != 0 && Real(s) * (Real{1} - eta_) < eta_; }

    /// Change of eigenvalue of Rx=product(sigma_x), relative to the vacuum.
    int bound_x_parity(unsigned s) const
    {
      require_bound(s);
      return s % 2 ? -1 : 1;
    }

    Real bound_gap(unsigned s) const { return normal_energy((Real{2} * bound_parameter(s)) * gap()); }

    /// q in [0,2*pi] is Q=P-P0-pi*delta_nu_z (mod 2*pi).
    /// Both delta_nu_z=0,1 copies are present in the thermodynamic description.
    /// A lattice bound curve need not lie below the full two-spinon continuum
    /// at every q; do not interpret its name as a spectral-weight assertion.
    Real bound_energy(unsigned s, Real q) const
    {
      Real const pi = bethe::detail::pi<Real>();
      if (!uni20::isfinite(q) || q < Real{0} || q > Real{2} * pi)
        throw std::invalid_argument("XYZ bound momentum must be finite and in [0, 2*pi]");
      Real const a = bound_parameter(s);
      q = std::min(q, Real{2} * pi - q);
      if (q == Real{0}) return normal_energy((Real{2} * a) * gap());
      if (q == pi) return normal_energy((Real{2} / a) * maximum_energy());
      Real const sine = std::sin(q / Real{2}), cosine = std::cos(q / Real{2});
      // Retain the small mass without forming sqrt(1-k1*k1).
      Real const first = std::hypot(maximum_energy() * sine, gap() * a * cosine);
      Real const second = std::hypot(sine, a * cosine);
      return normal_energy((Real{2} * second / a) * first);
    }

  private:
    void require_bound(unsigned s) const
    {
      if (!bound_exists(s)) throw std::invalid_argument("XYZ bound branch requires s>=1 and s*(1-eta)<eta");
    }

    static Real normal_energy(Real energy)
    {
      if (!uni20::isfinite(energy)) throw std::overflow_error("XYZ excitation energy exceeds scalar range");
      if (energy < uni20::numeric_limits<Real>::min())
        throw std::underflow_error("XYZ excitation energy is below the normal scalar range; use higher precision");
      return energy;
    }

    static Real log_theta(unsigned kind, Real u, Real t)
    {
      auto const jet = bethe::detail::elliptic_theta(kind, std::complex<Real>(u, 0), t);
      Real const value = jet.derivative[0].real();
      if (!(value > Real{0})) throw std::runtime_error("XYZ positive theta value lost at working precision");
      return jet.log_scale + std::log(value);
    }

    Real bound_parameter(unsigned s) const
    {
      require_bound(s);
      // a_s=sn(s*K1*(1-eta)/t, k1'). In the complementary nome its
      // normalized argument is s*(1-eta)/(2*eta), strictly between 0 and 1/2.
      Real const u = (Real(s) * (Real{1} - eta_) / eta_) / Real{2}, tc = t_ / eta_;
      Real const log_a =
          log_theta(3, Real{0}, tc) - log_theta(2, Real{0}, tc) + log_theta(1, u, tc) - log_theta(4, u, tc);
      Real const a = std::exp(log_a);
      if (!uni20::isfinite(a) || !(a > Real{0}) || a > Real{1} + Real{128} * uni20::numeric_limits<Real>::epsilon())
        throw std::runtime_error("XYZ bound parameter lost at working precision");
      return std::min(a, Real{1});
    }

    static bethe::detail::SpinonBand<Real> make_band(Real eta, Real t, Real exchange)
    {
      bethe::detail::validate_spinon_parameters(Real{0}, exchange);
      if (!uni20::isfinite(eta) || eta <= Real{0} || eta >= Real{1} || !uni20::isfinite(t) || t <= Real{0})
        throw std::invalid_argument("XYZ spinons require finite 0<eta<1 and t>0");
      auto const numerator = bethe::detail::elliptic_theta(1, std::complex<Real>(eta, 0), t);
      auto const denominator = bethe::detail::elliptic_theta(1, std::complex<Real>{}, t);
      Real const ratio = bethe::detail::theta_ratio(numerator, 0, denominator, 1).real();
      Real const nome_t = eta / t;
      if (!(nome_t > Real{0}) || !uni20::isfinite(nome_t))
        throw std::overflow_error("XYZ excitation nome exceeds scalar range");
      bool const dual = nome_t < Real{1};
      Real const theta_t = dual ? t / eta : nome_t;
      auto const top = bethe::detail::elliptic_theta(3, std::complex<Real>{}, theta_t);
      auto const bottom = bethe::detail::elliptic_theta(dual ? 2 : 4, std::complex<Real>{}, theta_t);
      // I=pi*theta1(eta|it)*theta3(0|i*eta/t)^2/(2*t*theta1'(0|it)).
      // After the modular transform t cancels, avoiding an intermediate tiny
      // prefactor and evaluating the mass from theta2, never from 1-k1^2.
      Real const factor = (bethe::detail::pi<Real>() / Real{2}) * (ratio / (dual ? eta : t));
      Real const scale = exchange * factor;
      Real const maximum = normal_energy(scale * top.derivative[0].real() * top.derivative[0].real());
      Real const log_gap = std::log(scale) + Real{2} * (bottom.log_scale + std::log(bottom.derivative[0].real()));
      if (log_gap < std::log(uni20::numeric_limits<Real>::min()))
        throw std::underflow_error("XYZ spinon gap is below the normal scalar range; use higher precision");
      Real const gap =
          normal_energy(bottom.log_scale == Real{0} ? scale * bottom.derivative[0].real() * bottom.derivative[0].real()
                                                    : std::exp(log_gap));
      return {maximum, std::min(gap, maximum)};
    }
    Real eta_, t_;
    bethe::detail::SpinonBand<Real> band_;
};
} // namespace bethe::xyz
