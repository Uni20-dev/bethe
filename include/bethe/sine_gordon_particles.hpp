// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/solver.hpp>
#include <span>
#include <stdexcept>
#include <vector>

namespace bethe::sine_gordon
{
enum class ParticleKind
{
  soliton,
  antisoliton,
  breather
};
struct Particle
{
    ParticleKind kind = ParticleKind::soliton;
    unsigned index = 0; // positive only for B_index
    int charge() const { return kind == ParticleKind::soliton ? 1 : kind == ParticleKind::antisoliton ? -1 : 0; }
};

/// Infinite-volume particle spectrum. M is the physical soliton mass;
/// p=beta^2/(8*pi-beta^2), velocity=hbar=1. No lattice momentum folding.
template <uni20::Real Real> class ParticleSpectrum {
  public:
    ParticleSpectrum(Real mass, Real p) : mass_(mass), p_(p)
    {
      if (!uni20::isfinite(mass) || mass <= Real{0} || !uni20::isfinite(p) || p <= Real{0})
        throw std::invalid_argument("sine-Gordon particles require finite M>0 and p>0");
    }
    bool exists(Particle particle) const
    {
      switch (particle.kind)
      {
        case ParticleKind::soliton:
        case ParticleKind::antisoliton:
          return particle.index == 0;
        case ParticleKind::breather:
          // Strict pole condition: do not count n*p=1 threshold states.
          return particle.index > 0 && Real(particle.index) * p_ < Real{1};
      }
      return false;
    }
    Real mass(Particle particle) const
    {
      if (!exists(particle)) throw std::invalid_argument("particle absent from stable sine-Gordon spectrum");
      if (particle.kind != ParticleKind::breather) return mass_;
      Real const pi = Real{4} * std::atan(Real{1});
      Real const ratio = Real{2} * std::sin((Real(particle.index) * p_) * (pi / Real{2}));
      Real const value = mass_ * ratio;
      if (!uni20::isfinite(value)) throw std::overflow_error("breather mass exceeds scalar range");
      if (value == Real{0}) throw std::underflow_error("positive breather mass underflows scalar range");
      return value;
    }
    Real energy(Particle particle, Real momentum) const { return relativistic_energy(mass(particle), momentum); }
    /// Infimum for this particle content at fixed total momentum. There is
    /// no finite upper edge; this does not assert any observable's weight.
    Real threshold(std::span<Particle const> particles, Real momentum) const
    {
      if (particles.size() < 2) throw std::invalid_argument("a continuum needs at least two particles");
      bethe::detail::CompensatedSum<Real> sum;
      for (auto particle : particles)
        sum.add(mass(particle));
      return relativistic_energy(sum.value(), momentum);
    }
    /// Bounded enumeration: never silently truncate a large weak-coupling spectrum.
    std::vector<Particle> particles(unsigned max_breathers = 1024) const
    {
      std::vector<Particle> out{{ParticleKind::soliton, 0}, {ParticleKind::antisoliton, 0}};
      for (unsigned n = 1; Real(n) * p_ < Real{1}; ++n)
      {
        if (n > max_breathers || n == std::numeric_limits<unsigned>::max())
          throw std::length_error("stable breather enumeration exceeds budget");
        out.push_back({ParticleKind::breather, n});
      }
      return out;
    }

  private:
    static Real relativistic_energy(Real mass, Real momentum)
    {
      if (!uni20::isfinite(momentum)) throw std::invalid_argument("momentum must be finite");
      Real const value = std::hypot(mass, momentum);
      if (!uni20::isfinite(value)) throw std::overflow_error("relativistic energy exceeds scalar range");
      return value;
    }
    Real mass_, p_;
};
} // namespace bethe::sine_gordon
