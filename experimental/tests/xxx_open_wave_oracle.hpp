// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/solver.hpp>
#include <complex>
#include <span>
#include <stdexcept>
#include <vector>

namespace test_support
{
/// Independent finite-Hilbert-space ABA oracle, NOT a scalable state solver.
/// R(v)=v*I+i*P; free-end double row T(v)*T_hat(v), v=u-i/2.
/// Fixed imaginary u=-i/2 also tests the exact local inverse problem.
template <uni20::Real Real>
std::vector<uni20::complex<Real>> open_aba_create(unsigned sites, std::span<uni20::complex<Real> const> state,
                                                  uni20::complex<Real> u, Real boundary_shift = Real{0})
{
  using C = uni20::complex<Real>;
  if (sites < 2 || sites > 10 || state.size() != (std::size_t{1} << sites))
    throw std::invalid_argument("small-chain ABA oracle requires 2<=N<=10");
  auto const dim = state.size();
  Real const zero = bethe::detail::precision_zero(u.real()), one = uni20::scalar_like(zero, 1);
  std::vector<C> vector(2 * dim, C(zero, zero)), next(2 * dim, C(zero, zero));
  std::copy(state.begin(), state.end(), vector.begin() + dim); // Auxiliary spin down.
  auto r = [&](unsigned site, Real shift) {
    for (std::size_t bits = 0; bits < vector.size(); ++bits)
      next[bits] = (u + shift - C(zero, one / 2)) * vector[bits];
    for (std::size_t bits = 0; bits < vector.size(); ++bits)
    {
      bool const auxiliary = (bits & dim) != 0;
      bool const physical = ((bits >> site) & 1U) != 0;
      auto const swapped = auxiliary == physical ? bits : bits ^ dim ^ (std::size_t{1} << site);
      next[swapped] += C(zero, one) * vector[bits];
    }
    vector.swap(next);
  };
  for (unsigned j = sites; j > 0; --j)
    r(j - 1, j == 1 ? boundary_shift : zero);
  for (unsigned j = 0; j < sites; ++j)
    r(j, j == 0 ? -boundary_shift : zero);
  vector.resize(dim); // Auxiliary spin up: B lowers the physical spin.
  return vector;
}

template <uni20::Real Real>
std::vector<uni20::complex<Real>> open_aba_wave(unsigned sites, std::span<Real const> rapidities,
                                                Real boundary_shift = Real{0})
{
  using C = uni20::complex<Real>;
  using std::sqrt;
  if (sites < 2 || sites > 10) throw std::invalid_argument("small-chain ABA oracle requires 2<=N<=10");
  Real const zero = bethe::detail::precision_zero(boundary_shift), one = uni20::scalar_like(zero, 1);
  std::vector<C> wave(std::size_t{1} << sites, C(zero, zero));
  wave[0] = C(one, zero);
  for (Real z : rapidities)
  {
    wave = open_aba_create<Real>(sites, wave, C(z / 2, zero), boundary_shift);
    bethe::detail::CompensatedSum<Real> norm(zero);
    for (C value : wave)
      norm.add(value.real() * value.real() + value.imag() * value.imag());
    Real const scale = sqrt(norm.value());
    if (!uni20::isfinite(scale) || scale <= 0) throw std::runtime_error("invalid ABA oracle norm");
    for (C& value : wave)
      value /= scale;
  }
  return wave;
}
} // namespace test_support
