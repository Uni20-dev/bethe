// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/elliptic_theta.hpp>
#include <span>

namespace bethe::xyz
{
template <uni20::Real Real> struct Couplings
{
    Real x{}, y{}, z{};
};

/// Rectangular Hermitian parametrization: theta_j(pi*eta | i*t)/theta_j(0 | i*t).
/// H=sum(Jx*Sx*Sx+Jy*Sy*Sy+Jz*Sz*Sz), S=sigma/2. This is only the
/// coupling parametrization; it does not select roots or compute an eigenstate.
template <uni20::Real Real> Couplings<Real> couplings(Real eta, Real t)
{
  if (!uni20::isfinite(eta)) throw std::invalid_argument("XYZ eta must be finite and real");
  auto ratio = [&](unsigned kind) {
    auto const numerator = bethe::detail::elliptic_theta(kind, std::complex<Real>(eta, 0), t);
    auto const denominator = bethe::detail::elliptic_theta(kind, std::complex<Real>{}, t);
    return bethe::detail::theta_ratio(numerator, 0, denominator, 0).real();
  };
  return {ratio(4), ratio(3), ratio(2)};
}
/// Evaluate the regular-root energy expression for even N and M=N/2.
/// This does NOT solve/validate Bethe equations or certify a physical state.
/// A complex result is retained; poles are rejected rather than regularized.
template <uni20::Real Real>
std::complex<Real> candidate_energy(std::size_t sites, std::span<std::complex<Real> const> roots, Real eta, Real t)
{
  using Complex = std::complex<Real>;
  if (sites < 2 || sites % 2 || roots.size() != sites / 2 || !uni20::isfinite(eta) || !(eta > Real{0} && eta < Real{1}))
    throw std::invalid_argument("XYZ regular energy requires even N>=2, N/2 roots, and real 0<eta<1");
  auto const at_eta = bethe::detail::elliptic_theta(1, Complex(eta, 0), t);
  auto const at_zero = bethe::detail::elliptic_theta(1, Complex{}, t);
  auto const factor = bethe::detail::theta_ratio(at_eta, 0, at_zero, 1);
  auto g = [&](Complex u) {
    auto const jet = bethe::detail::elliptic_theta(1, u, t);
    return factor * bethe::detail::theta_ratio(jet, 1, jet, 0);
  };
  bethe::detail::CompensatedSum<Complex> energy;
  energy.add(Real(sites) * g(Complex(eta, 0)) / Real{4});
  for (auto root : roots)
    energy.add((g(root - eta / Real{2}) - g(root + eta / Real{2})) / Real{2});
  if (!bethe::detail::theta_detail::finite(energy.value()))
    throw std::overflow_error("XYZ candidate energy exceeds range");
  return energy.value();
}
} // namespace bethe::xyz
