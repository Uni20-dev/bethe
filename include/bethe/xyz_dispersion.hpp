// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/spinon_band.hpp>
#include <uni20/core/elliptic_theta_jet.hpp>

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

    /// An owning immutable branch: theta functions are evaluated at preparation,
    /// while each energy query uses only elementary functions. Safe to retain
    /// after the parent dispersion is destroyed and share among Async readers.
    class BoundBranch {
      public:
        unsigned index() const { return index_; }
        int x_parity() const { return index_ % 2 ? -1 : 1; }
        Real gap() const { return normal_energy((Real{2} * a_) * spinon_gap_); }

        /// q is reduced Q in [0,2*pi], not a parity-shifted physical momentum.
        Real energy(Real q) const
        {
          bethe::detail::require_same_precision(q, maximum_);
          Real const pi = bethe::detail::spinon_pi(maximum_);
          if (!uni20::isfinite(q) || q < Real{0} || q > Real{2} * pi)
            throw std::invalid_argument("XYZ bound momentum must be finite and in [0, 2*pi]");
          q = std::min(q, Real{2} * pi - q);
          if (q == Real{0}) return gap();
          if (q == pi) return normal_energy((Real{2} / a_) * maximum_);
          Real const sine = uni20::math::sin(q / Real{2}), cosine = uni20::math::cos(q / Real{2});
          // Retain the small mass without forming sqrt(1-k1*k1).
          Real const first = uni20::math::hypot(maximum_ * sine, spinon_gap_ * a_ * cosine);
          Real const second = uni20::math::hypot(sine, a_ * cosine);
          return normal_energy((Real{2} * second / a_) * first);
        }

      private:
        friend class SpinonDispersion;
        BoundBranch(unsigned index, Real a, Real maximum, Real gap)
            : index_(index), a_(a), maximum_(maximum), spinon_gap_(gap)
        {}
        unsigned index_;
        Real a_, maximum_, spinon_gap_;
    };

    BoundBranch prepare_bound(unsigned s) const { return {s, bound_parameter(s), maximum_energy(), gap()}; }
    Real bound_gap(unsigned s) const { return prepare_bound(s).gap(); }

    /// q in [0,2*pi] is Q=P-P0-pi*delta_nu_z (mod 2*pi).
    /// Both delta_nu_z=0,1 copies are present in the thermodynamic description.
    /// A lattice bound curve need not lie below the full two-spinon continuum
    /// at every q; do not interpret its name as a spectral-weight assertion.
    /// Prepare a branch once with prepare_bound(s) when evaluating a curve.
    Real bound_energy(unsigned s, Real q) const { return prepare_bound(s).energy(q); }

  private:
    void require_bound(unsigned s) const
    {
      if (!bound_exists(s)) throw std::invalid_argument("XYZ bound branch requires s>=1 and s*(1-eta)<eta");
    }

    static Real minimum_energy(Real const& exemplar)
    {
      if constexpr (uni20::has_runtime_precision_v<Real>)
        return uni20::math::next_up(bethe::detail::precision_zero(exemplar));
      else
        return uni20::numeric_limits<Real>::min();
    }

    static Real normal_energy(Real energy)
    {
      if (!uni20::isfinite(energy)) throw std::overflow_error("XYZ excitation energy exceeds scalar range");
      if (energy < minimum_energy(energy))
        throw std::underflow_error("XYZ excitation energy is below the normal scalar range; use higher precision");
      return energy;
    }

    static Real log_theta(unsigned kind, Real u, Real t)
    {
      auto const jet = uni20::elliptic_theta_jet(kind, uni20::complex<Real>(u, bethe::detail::precision_zero(t)), t);
      Real const value = jet.derivative[0].real();
      if (!(value > Real{0})) throw std::runtime_error("XYZ positive theta value lost at working precision");
      return jet.log_scale + uni20::math::log(value);
    }

    Real bound_parameter(unsigned s) const
    {
      require_bound(s);
      // a_s=sn(s*K1*(1-eta)/t, k1'). In the complementary nome its
      // normalized argument is s*(1-eta)/(2*eta), strictly between 0 and 1/2.
      Real const u = (Real(s) * (Real{1} - eta_) / eta_) / Real{2}, tc = t_ / eta_;
      Real const log_a = log_theta(3, bethe::detail::precision_zero(tc), tc) -
                         log_theta(2, bethe::detail::precision_zero(tc), tc) + log_theta(1, u, tc) -
                         log_theta(4, u, tc);
      Real const a = uni20::math::exp(log_a);
      if (!uni20::isfinite(a) || !(a > Real{0}) || a > Real{1} + Real{128} * uni20::numeric_limits<Real>::epsilon(eta_))
        throw std::runtime_error("XYZ bound parameter lost at working precision");
      return std::min(a, uni20::scalar_like(a, 1));
    }

    static bethe::detail::SpinonBand<Real> make_band(Real eta, Real t, Real exchange)
    {
      Real const zero = bethe::detail::precision_zero(eta);
      bethe::detail::require_same_precision(t, eta);
#if UNI20_HAS_MPREAL
      if constexpr (uni20::has_runtime_precision_v<Real>)
        if (!exchange.is_exact()) bethe::detail::require_same_precision(exchange, eta);
#endif
      exchange = uni20::scalar_like(eta, exchange);
      bethe::detail::validate_spinon_parameters(zero, exchange);
      if (!uni20::isfinite(eta) || eta <= Real{0} || eta >= Real{1} || !uni20::isfinite(t) || t <= Real{0})
        throw std::invalid_argument("XYZ spinons require finite 0<eta<1 and t>0");
      auto const numerator = uni20::elliptic_theta_jet(1, uni20::complex<Real>(eta, zero), t);
      auto const denominator = uni20::elliptic_theta_jet(1, uni20::complex<Real>(zero, zero), t);
      Real const ratio = uni20::elliptic_theta_ratio(numerator, 0, denominator, 1).real();
      Real const nome_t = eta / t;
      if (!(nome_t > Real{0}) || !uni20::isfinite(nome_t))
        throw std::overflow_error("XYZ excitation nome exceeds scalar range");
      bool const dual = nome_t < Real{1};
      Real const theta_t = dual ? t / eta : nome_t;
      auto const top = uni20::elliptic_theta_jet(3, uni20::complex<Real>(zero, zero), theta_t);
      auto const bottom = uni20::elliptic_theta_jet(dual ? 2 : 4, uni20::complex<Real>(zero, zero), theta_t);
      // I=pi*theta1(eta|it)*theta3(0|i*eta/t)^2/(2*t*theta1'(0|it)).
      // After the modular transform t cancels, avoiding an intermediate tiny
      // prefactor and evaluating the mass from theta2, never from 1-k1^2.
      Real const factor = (bethe::detail::spinon_pi(eta) / Real{2}) * (ratio / (dual ? eta : t));
      Real const scale = exchange * factor;
      Real const maximum = normal_energy(scale * top.derivative[0].real() * top.derivative[0].real());
      Real const log_gap =
          uni20::math::log(scale) + Real{2} * (bottom.log_scale + uni20::math::log(bottom.derivative[0].real()));
      if (log_gap < uni20::math::log(minimum_energy(scale)))
        throw std::underflow_error("XYZ spinon gap is below the normal scalar range; use higher precision");
      Real const gap =
          normal_energy(bottom.log_scale == Real{0} ? scale * bottom.derivative[0].real() * bottom.derivative[0].real()
                                                    : uni20::math::exp(log_gap));
      return {maximum, std::min(gap, maximum)};
    }
    Real eta_, t_;
    bethe::detail::SpinonBand<Real> band_;
};
} // namespace bethe::xyz
