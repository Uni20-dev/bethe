// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/xxz_common.hpp>
#include <complex>

namespace bethe::xxz
{
namespace detail
{
// Inputs reduced modulo n; checked_sites leaves room for each addition.
inline std::size_t helix_product_mod(std::size_t a, std::size_t b, std::size_t n)
{
  std::size_t result = 0;
  a %= n;
  for (; b; b >>= 1, a = (a + a) % n)
    if (b & 1) result = (result + a) % n;
  return result;
}

// Reduce integer phases BEFORE floating conversion. Reflection about pi
// also avoids subtracting nearly equal angles for the odd-ring collision.
template <uni20::Real Real> std::complex<Real> helix_phase(std::size_t n, std::size_t index)
{
  if (!index) return {Real{1}, Real{0}};
  bool const negative = index > n / 2;
  auto const k = negative ? n - index : index;
  if (2 * k == n) return {-Real{1}, Real{0}};
  if (4 * k == n) return {Real{0}, negative ? -Real{1} : Real{1}};
  Real const pi = Real{4} * std::atan(Real{1});
  bool const near_pi = k > n / 4;
  Real const angle = pi * (Real(near_pi ? n - 2 * k : 2 * k) / Real(n));
  Real const cosine = (near_pi ? -Real{1} : Real{1}) * std::cos(angle);
  Real const sine = (negative ? -Real{1} : Real{1}) * std::sin(angle);
  return {cosine, sine};
}
} // namespace detail

/// Fixed-Sz projection of a periodic spin helix. This is an explicit
/// eigenstate at its own commensurate Delta, NOT a ground-state result.
/// See Popkov, Zhang and Kluemper (2021), and docs/xxz-spin-helix.md.
template <uni20::Real Real> struct SpinHelixState
{
    std::size_t sites = 0, down_spins = 0, winding = 0;
    uni20::half_int sz;
    Real delta = Real{0}, energy = Real{0};
    std::size_t momentum_index = 0;
    Real momentum = Real{0};

    /// Unit-magnitude amplitude on an ordered list of occupied sites.
    /// The full squared norm is binomial(N,M); no exponential-size vector
    /// or potentially underflowing normalization factor is constructed.
    std::complex<Real> unnormalized_amplitude(std::span<std::size_t const> occupied) const
    {
      detail::checked_sites(sites);
      if (winding >= sites || down_spins > sites || occupied.size() != down_spins)
        throw std::invalid_argument("spin-helix amplitude has invalid winding or sector size");
      std::size_t sum = 0;
      for (std::size_t j = 0; j < occupied.size(); ++j)
      {
        if (occupied[j] >= sites || (j && occupied[j - 1] >= occupied[j]))
          throw std::invalid_argument("spin-helix occupied sites must be strictly increasing in [0,N)");
        sum = (sum + occupied[j]) % sites;
      }
      return detail::helix_phase<Real>(sites, detail::helix_product_mod(winding, sum, sites));
    }
};

/// Construct the eigenstate with pitch 2*pi*winding/N and Delta=cos(pitch).
/// Every physical Sz and winding=0,...,N-1 is supported, on either parity.
/// No input coupling is rounded to a commensurate value.
template <uni20::Real Real = double>
SpinHelixState<Real> periodic_spin_helix(std::size_t sites, std::size_t winding, uni20::half_int sz)
{
  auto const n = detail::checked_sites(sites);
  detail::sector_roots(sites, sz); // Validate without folding the stored sector.
  if (winding >= sites) throw std::invalid_argument("spin-helix winding must lie in [0,N)");
  SpinHelixState<Real> state;
  state.sites = sites;
  state.down_spins = std::size_t((n - sz.twice()) / 2);
  state.winding = winding;
  state.sz = sz;
  state.delta = detail::helix_phase<Real>(sites, winding).real();
  state.energy = (Real(sites) / Real{4}) * state.delta;
  state.momentum_index = detail::helix_product_mod(winding, state.down_spins, sites);
  state.momentum = Real{8} * std::atan(Real{1}) * (Real(state.momentum_index) / Real(sites));
  return state;
}
} // namespace bethe::xxz
