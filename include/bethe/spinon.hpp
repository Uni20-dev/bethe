// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#pragma once

#include <bethe/solver.hpp>
#include <uni20/core/types.hpp>
#if UNI20_HAS_FLOAT128
#include <mplapack_binary128.h>
#endif
#include <uni20/core/math.hpp>

#include <cmath>
#include <stdexcept>

namespace bethe::detail
{
template <uni20::Real Real> Real pi()
{
  using std::atan;
  return Real{4} * atan(Real{1});
}

template <uni20::Real Real> Real spinon_pi(Real const& exemplar)
{
#if UNI20_HAS_MPREAL
  if constexpr (uni20::has_runtime_precision_v<Real>)
    return pi_like(precision_zero(exemplar));
  else
#endif
    return pi<Real>();
}

template <uni20::Real Real> void validate_spinon_parameters(Real k, Real exchange)
{
  if (!uni20::isfinite(k) || k < Real{0} || k > spinon_pi(k))
    throw std::invalid_argument("spinon momentum must be finite and in [0, pi]");
  if (!uni20::isfinite(exchange) || exchange <= Real{0})
    throw std::invalid_argument("antiferromagnetic exchange must be finite and positive");
}
} // namespace bethe::detail

namespace bethe::heisenberg
{
/// Zero-field thermodynamic single-spinon energy, H=J sum S_i.S_(i+1).
/// Momentum convention: 0 <= k <= pi. This is not a finite-chain energy.
/// Runtime precision is carried by k; an exact exchange (including the default
/// J=1) is allowed, but a finite-precision exchange must match k.
template <uni20::Real Real> [[nodiscard]] Real spinon_energy(Real k, Real exchange = Real{1})
{
  using std::sin;
  bethe::detail::validate_spinon_parameters(k, exchange);
#if UNI20_HAS_MPREAL
  if constexpr (uni20::has_runtime_precision_v<Real>)
    if (!exchange.is_exact()) bethe::detail::require_same_precision(exchange, k);
#endif
  Real const pi = bethe::detail::spinon_pi(k);
  if (k == Real{0} || k == pi) return uni20::scalar_like(k, 0);
  return exchange * (pi / Real{2}) * sin(k);
}

/// Zero-field thermodynamic ground-state energy per site, including J/4.
/// Runtime precision requires an explicit finite-precision exchange argument.
template <uni20::Real Real = double> [[nodiscard]] Real bulk_energy_density(Real exchange = Real{1})
{
  using std::log;
  auto const zero = bethe::detail::precision_zero(exchange);
  bethe::detail::validate_spinon_parameters(zero, exchange);
  return exchange * (Real{1} / Real{4} - log(uni20::scalar_like(exchange, 2)));
}
} // namespace bethe::heisenberg

namespace bethe::xxz
{
/// Zero-field thermodynamic spinon dispersion for -1 < Delta <= 1.
/// H=J sum (SxSx + SySy + Delta SzSz), k in [0, pi]. No finite-size XXZ
/// solver is implied by this function. Delta=1 uses the explicit XXX limit.
template <uni20::Real Real> [[nodiscard]] Real spinon_energy(Real k, Real delta, Real exchange = Real{1})
{
  using std::acos;
  using std::sin;
  bethe::detail::validate_spinon_parameters(k, exchange);
  if (!uni20::isfinite(delta) || delta <= -Real{1} || delta > Real{1})
    throw std::invalid_argument("gapless XXZ dispersion requires -1 < Delta <= 1");
  if (delta == Real{1}) return heisenberg::spinon_energy(k, exchange);
  Real const gamma = acos(delta);
  return heisenberg::spinon_energy(k, exchange) * (sin(gamma) / gamma);
}
} // namespace bethe::xxz
