// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/xxz_common.hpp>
#include <cmath>

namespace bethe::xxz::detail
{
/// Rayleigh quotients of two explicit fixed-sector trial states, independent
/// of the polynomial equations. These are upper bounds on the sector minimum
/// in exact arithmetic, NOT lower bounds or certificates of a candidate.
template <uni20::Real Real> struct OddSectorTrialEnergies
{
    Real free_sea{}, helix{};
    Real upper_bound() const { return std::min(free_sea, helix); }
};

/// Odd N, folded sector M<=N/2. Neither trial state depends on Delta, so
/// both expectations are affine functions of the requested coupling.
template <uni20::Real Real> class OddSectorVariationalBounds {
  public:
    OddSectorVariationalBounds(std::size_t sites, std::size_t roots)
    {
      checked_sites(sites);
      if (sites % 2 == 0 || roots > sites / 2)
        throw std::invalid_argument("odd XXZ variational bounds require odd N and M <= N/2");
      Real const n = Real(sites), m = Real(roots), pi = Real{4} * std::atan(Real{1});
      Real const angle = pi / n, cosine = -std::cos(angle);
      Real const ratio = roots == 0 ? Real{0} : std::sin(pi * (m / n)) / std::sin(angle);
      free_xy_ = cosine * ratio;
      // Wick's theorem for the free sea: <n_j n_(j+1)>=rho^2-|C_1|^2,
      // with |C_1|=sin(pi*M/N)/(N*sin(pi/N)).
      free_zz_ = n / Real{4} - m + (m * m - ratio * ratio) / n;
      // Uniform-magnitude projected helix, q=pi+pi/N. Counting unlike
      // neighbors gives <sum SzSz>=N/4-M*(N-M)/(N-1).
      Real const pairs = m * (Real(sites - roots) / Real(sites - 1));
      helix_xy_ = cosine * pairs;
      helix_zz_ = n / Real{4} - pairs;
    }

    OddSectorTrialEnergies<Real> evaluate(Real delta) const
    {
      if (!uni20::isfinite(delta) || delta <= -Real{1} || delta > Real{0})
        throw std::invalid_argument("odd XXZ variational bounds require -1 < Delta <= 0");
      return {free_xy_ + delta * free_zz_, helix_xy_ + delta * helix_zz_};
    }

  private:
    Real free_xy_, free_zz_, helix_xy_, helix_zz_;
};
} // namespace bethe::xxz::detail
