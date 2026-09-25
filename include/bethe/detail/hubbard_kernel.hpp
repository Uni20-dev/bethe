// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <array>
#include <bethe/detail/gauss_legendre.hpp>
#include <bethe/solver.hpp>
#include <cstdint>
#include <utility>
#include <vector>

namespace bethe::detail
{
// R(x) = integral cos(w*x)/(1+exp(2*u*w)) dw/pi.
// Euler's beta-integral expansion of sum (-1)^m/(m+1+i*z),
// z=x/(2u), avoids numerical differences of almost equal digammas.
template <uni20::Real Real> Real hubbard_r(Real x, Real u)
{
  Real const pi = Real{4} * std::atan(Real{1}), z = x / (Real{2} * u);
  uni20::complex<Real> const a{Real{1}, z};
  uni20::complex<Real> term = Real{1} / (Real{2} * a);
  CompensatedSum<Real> sum;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (int j = 0; j < 2 * uni20::numeric_limits<Real>::digits + 32; ++j)
  {
    sum.add(term.real());
    term *= Real(j + 1) / (Real{2} * (a + Real(j + 1)));
    if (std::abs(term.real()) + std::abs(term.imag()) <= eps * std::abs(sum.value()) / Real{16}) break;
  }
  return sum.value() / (Real{2} * pi * u);
}

// Integral_0^x R(t)dt = arg[Gamma(1+i*x/(4u))/Gamma(1/2+i*x/(4u))]/pi.
// Evaluate the log-gamma RATIO: paired recurrence and paired leading Stirling
// terms avoid subtracting two large phases. Re(z)>=32 makes 16 Bernoulli
// terms sufficient through binary128. No double-precision coefficients.
template <uni20::Real Real> Real hubbard_r_primitive(Real x, Real u)
{
  if (x == Real{0}) return Real{0};
  Real const b = std::abs(x) / (Real{4} * u), pi = Real{4} * std::atan(Real{1});
  // Remaining asymptotic correction is O(1/b); also avoid squaring huge b.
  if (b > Real{1} / uni20::numeric_limits<Real>::epsilon()) return std::copysign(Real{.25}, x);
  Real const a{32}, b2 = b * b;
  CompensatedSum<Real> phase;
  phase.add(Real{.5} * std::atan2(b, a + Real{1}));
  phase.add(a * std::atan2(-b / Real{2}, (a + Real{1}) * (a + Real{.5}) + b2));
  phase.add(b / Real{2} * std::log1p((a + Real{.75}) / ((a + Real{.5}) * (a + Real{.5}) + b2)));
  for (int j = 0; j < 32; ++j)
    phase.add(std::atan2(b / Real{2}, (Real(j) + Real{1}) * (Real(j) + Real{.5}) + b2));
  constexpr std::array<std::pair<std::int64_t, std::int64_t>, 16> bernoulli = {{{1, 12},
                                                                                {-1, 360},
                                                                                {1, 1260},
                                                                                {-1, 1680},
                                                                                {1, 1188},
                                                                                {-691, 360360},
                                                                                {1, 156},
                                                                                {-3617, 122400},
                                                                                {43867, 244188},
                                                                                {-174611, 125400},
                                                                                {77683, 5796},
                                                                                {-236364091, 1506960},
                                                                                {657931, 300},
                                                                                {-3392780147, 93960},
                                                                                {1723168255201, 2492028},
                                                                                {-7709321041217, 505920}}};
  auto const inv1 = Real{1} / uni20::complex<Real>{a + Real{1}, b};
  auto const inv0 = Real{1} / uni20::complex<Real>{a + Real{.5}, b};
  auto power1 = inv1, power0 = inv0;
  for (auto const [numerator, denominator] : bernoulli)
  {
    phase.add(Real(numerator) / Real(denominator) * (power1.imag() - power0.imag()));
    power1 *= inv1 * inv1;
    power0 *= inv0 * inv0;
  }
  return std::copysign(phase.value() / pi, x);
}

} // namespace bethe::detail
