// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/spinon_band.hpp>

namespace bethe::su3
{
enum class Particle
{
  fundamental,
  antifundamental
};
enum class Continuum
{
  two_soliton, ///< 3 x bar3 = 1 + 8.
  four_soliton ///< (3 x bar3)^2; includes further adjoint channels.
};
enum class Momentum
{
  unfolded,
  three_site_folded
};

/// Infinite, zero-field SU(3) chain H=J sum P, J>0. Elementary holes are
/// 3 (first nesting sea) and bar3 (second sea). They are not single-particle
/// states of a balanced L=0 mod 3 periodic chain. No form factors are implied.
/// Doikou--Nepomechie (1998), Sec. 2.3; Voros--Penc (2021), Eq. (64), with
/// twice that paper's energies. Use Eq. (64)'s ranges, not the swapped ranges
/// in the prose below its Eq. (65).
template <uni20::Real Real> class ExcitationDispersion {
  public:
    explicit ExcitationDispersion(Real exchange = Real{1}) : exchange_(exchange)
    {
      if (!uni20::isfinite(exchange) || exchange <= Real{0})
        throw std::invalid_argument("SU(3) exchange must be finite and positive");
      amplitude_ = exchange * (Real{4} * pi_ / (Real{3} * std::sqrt(Real{3})));
      if (!uni20::isfinite(amplitude_) || amplitude_ > uni20::numeric_limits<Real>::max() / Real{1.5})
        throw std::overflow_error("SU(3) elementary band exceeds scalar range");
    }
    Real exchange() const { return exchange_; }
    Real velocity() const { return exchange_ * (Real{2} * pi_ / Real{3}); }
    Real momentum_max(Particle particle) const
    {
      check_particle(particle);
      return particle == Particle::fundamental ? Real{4} * pi_ / Real{3} : Real{2} * pi_ / Real{3};
    }
    Real maximum_energy(Particle particle) const
    {
      check_particle(particle);
      return amplitude_ * (particle == Particle::fundamental ? Real{1.5} : Real{0.5});
    }
    Real energy(Particle particle, Real p) const
    {
      Real const end = momentum_max(particle);
      if (!uni20::isfinite(p) || p < Real{0} || p > end)
        throw std::invalid_argument("momentum outside SU(3) elementary branch");
      // cos(p-centre)-cos(centre), evaluated without cancellation at either
      // gapless endpoint. Both representations have the same endpoint speed.
      Real const distance = std::min(p, end - p);
      // Halving a subnormal momentum can erase it even when v*p is nonzero.
      // The curvature correction there is far below relative roundoff.
      Real const value =
          distance < uni20::numeric_limits<Real>::min()
              ? velocity() * distance
              : amplitude_ * (Real{2} * std::sin(distance / Real{2}) * std::sin((end - distance) / Real{2}));
      if (distance > Real{0} && value == Real{0})
        throw std::underflow_error("positive SU(3) excitation energy underflows scalar range");
      return value;
    }
    /// Bounds for a specified particle content, NOT the full spectral support.
    /// Q is in [0,2*pi]. Folding takes the union at Q, Q+2*pi/3, Q+4*pi/3.
    /// A union's envelope need not describe its internal gaps or weights.
    SpinonContinuum<Real> continuum(Real q, Continuum content = Continuum::two_soliton,
                                    Momentum convention = Momentum::unfolded) const
    {
      if (!uni20::isfinite(q) || q < Real{0} || q > Real{2} * pi_)
        throw std::invalid_argument("SU(3) continuum momentum must lie in [0,2*pi]");
      if (content != Continuum::two_soliton && content != Continuum::four_soliton)
        throw std::invalid_argument("unknown SU(3) continuum content");
      if (convention != Momentum::unfolded && convention != Momentum::three_site_folded)
        throw std::invalid_argument("unknown SU(3) momentum convention");
      q = std::min(q, Real{2} * pi_ - q);
      Real const a = pi_ / Real{3}, soft = momentum_max(Particle::antifundamental);
      SpinonContinuum<Real> result;
      if (convention == Momentum::three_site_folded)
      {
        Real const r = q <= a ? q : q <= soft ? soft - q : q - soft;
        result.lower = energy(Particle::antifundamental, r);
        result.upper = content == Continuum::two_soliton ? amplitude_ * (Real{2} * std::cos((a - r) / Real{2}))
                                                         : amplitude_ * (Real{4} * std::cos(r / Real{4}));
      }
      else if (content == Continuum::two_soliton)
      {
        result.lower = q <= soft ? energy(Particle::antifundamental, q) : energy(Particle::fundamental, q - soft);
        result.upper = q < a ? energy(Particle::fundamental, q) : amplitude_ * (Real{2} * std::sin(q / Real{2}));
      }
      else
      {
        result.lower = energy(Particle::antifundamental, q <= soft ? q : q - soft);
        result.upper = amplitude_ * (Real{4} * std::cos(q / Real{4}));
      }
      if (!uni20::isfinite(result.lower) || !uni20::isfinite(result.upper))
        throw std::overflow_error("SU(3) continuum exceeds scalar range");
      return result;
    }

  private:
    static void check_particle(Particle particle)
    {
      if (particle != Particle::fundamental && particle != Particle::antifundamental)
        throw std::invalid_argument("unknown SU(3) particle representation");
    }
    Real pi_ = bethe::detail::pi<Real>();
    Real exchange_, amplitude_;
};
} // namespace bethe::su3
