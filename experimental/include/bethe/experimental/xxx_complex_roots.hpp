// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

// Regular, unreduced finite-root reference machinery. All roots are standard
// ABA lambda, not z. Unresolved string poles fail; this is not a string limit.
#include <bethe/experimental/xxx_singlet.hpp>
#include <limits>

namespace bethe::experimental::xxx::complex_roots
{
template <uni20::Real R> using C = uni20::complex<R>;
template <typename T> using Sum = bethe::detail::CompensatedSum<T>;

namespace detail
{
template <uni20::Real R> void validate_scalar(C<R> const& value, R const& zero)
{
#if UNI20_HAS_MPREAL
  if constexpr (uni20::has_runtime_precision_v<R>)
    if (!value.initialized() || value.is_exact() || uni20::precision_of(value) != uni20::precision_of(zero))
      throw std::invalid_argument("complex roots require one finite working precision");
#endif
  bethe::detail::require_same_precision(value.real(), zero);
  bethe::detail::require_same_precision(value.imag(), zero);
  if (!uni20::isfinite(value.real()) || !uni20::isfinite(value.imag()))
    throw std::invalid_argument("complex roots must be finite");
}

template <uni20::Real R> R validate(std::span<C<R> const> roots, R const& exemplar)
{
  R const zero = bethe::detail::precision_zero(exemplar);
  for (auto const& root : roots)
    validate_scalar<R>(root, zero);
  return zero;
}

inline void validate_size(int n, std::size_t roots)
{
  if (n < 2 || roots > static_cast<std::size_t>(n))
    throw std::invalid_argument("invalid finite-chain size or root count");
}

template <uni20::Real R> C<R> log_equation_unchecked(int n, std::span<C<R> const> roots, std::size_t i, R const& zero)
{
  using uni20::math::log;
  R const one = uni20::scalar_like(zero, 1);
  C<R> const imaginary(zero, one);
  Sum<C<R>> sum(C<R>(zero, zero));
  sum.add((n * one) * (log(roots[i] + imaginary / (2 * one)) - log(roots[i] - imaginary / (2 * one))));
  for (std::size_t j = 0; j < roots.size(); ++j)
    if (i != j) sum.add(log(roots[i] - roots[j] - imaginary) - log(roots[i] - roots[j] + imaginary));
  return sum.value();
}
} // namespace detail

template <uni20::Real R> C<R> log_equation(int n, std::span<C<R> const> roots, std::size_t i, R const& exemplar = R{0})
{
  detail::validate_size(n, roots.size());
  R const zero = detail::validate<R>(roots, exemplar);
  if (i >= roots.size()) throw std::invalid_argument("invalid root index");
  return detail::log_equation_unchecked<R>(n, roots, i, zero);
}

template <uni20::Real R> R phase(R const& angle)
{
  R const zero = bethe::detail::precision_zero(angle);
  return uni20::math::remainder(angle, 2 * bethe::detail::pi_like(zero));
}

// Fix atan on the imaginary-axis cuts: the quantum-number convention is odd
// across the two cuts, independently of the sign of the real signed zero.
template <uni20::Real R> C<R> bethe_atan(C<R> const& z)
{
#if UNI20_HAS_MPREAL
  if constexpr (uni20::has_runtime_precision_v<R>)
    if (!z.initialized() || z.is_exact()) throw std::invalid_argument("Bethe atan requires finite working precision");
#endif
  R const zero = bethe::detail::precision_zero(z.real());
  detail::validate_scalar<R>(z, zero);
  R const one = uni20::scalar_like(zero, 1);
  if (z.real() == zero && uni20::math::abs(z.imag()) > one)
  {
    R const sign = z.imag() > zero ? one : -one;
    return {sign * bethe::detail::pi_like(zero) / 2,
            sign * uni20::math::log1p(2 * one / (uni20::math::abs(z.imag()) - one)) / 2};
  }
  return uni20::math::atan(z);
}

template <uni20::Real R> R rational_residual(int n, std::span<C<R> const> roots, R const& exemplar = R{0})
{
  detail::validate_size(n, roots.size());
  R const zero = detail::validate<R>(roots, exemplar);
  R result = zero;
  for (std::size_t i = 0; i < roots.size(); ++i)
  {
    auto const l = detail::log_equation_unchecked<R>(n, roots, i, zero);
    R const a = phase(l.imag()), sine = uni20::math::sin(a / 2), magnitude = uni20::math::exp(l.real());
    R const error = uni20::math::abs(
        C<R>(uni20::math::expm1(l.real()) - 2 * magnitude * sine * sine, magnitude * uni20::math::sin(a)));
    if (!uni20::isfinite(error)) return uni20::scalar_like(zero, std::numeric_limits<double>::infinity());
    result = std::max(result, error);
  }
  return result;
}

template <uni20::Real R> std::vector<C<R>> gaudin(int n, std::span<C<R> const> roots, R const& exemplar = R{0})
{
  detail::validate_size(n, roots.size());
  R const zero = detail::validate<R>(roots, exemplar), one = uni20::scalar_like(zero, 1);
  auto const m = roots.size();
  std::vector<C<R>> matrix(m * m, C<R>(zero, zero));
  for (std::size_t i = 0; i < m; ++i)
  {
    Sum<C<R>> diagonal(C<R>(zero, zero));
    diagonal.add(n * one / (roots[i] * roots[i] + one / 4));
    for (std::size_t j = 0; j < m; ++j)
      if (i != j)
      {
        auto const d = roots[i] - roots[j];
        matrix[i * m + j] = 2 * one / (one + d * d);
        diagonal.add(-matrix[i * m + j]);
      }
    matrix[i * m + i] = diagonal.value();
  }
  return matrix;
}

// Arithmetic rejection diagnostic, not an error bound. Includes norm poles,
// coincident roots and internal string poles before any determinant.
template <uni20::Real R> R separation(std::span<C<R> const> roots, R const& exemplar = R{0})
{
  R const zero = detail::validate<R>(roots, exemplar), one = uni20::scalar_like(zero, 1);
  R result = one;
  for (std::size_t i = 0; i < roots.size(); ++i)
  {
    result = std::min(result, uni20::math::abs(roots[i] * roots[i] + one / 4));
    for (std::size_t j = 0; j < i; ++j)
    {
      auto const d = roots[i] - roots[j];
      result = std::min({result, uni20::math::abs(d), uni20::math::abs(one + d * d)});
    }
  }
  return result;
}

template <uni20::Real R> R residual_tolerance(int n, std::span<C<R> const> roots, R const& exemplar = R{0})
{
  detail::validate_size(n, roots.size());
  R const zero = bethe::detail::precision_zero(exemplar);
  R const s = separation<R>(roots, zero), eps = uni20::numeric_limits<R>::epsilon(zero);
  if (!(s > 64 * uni20::math::sqrt(eps))) throw std::runtime_error("unresolved complex roots; increase precision");
  // Coincident roots are inadmissible, but are not poles of the Bethe
  // equations: their scattering factors have denominators d +/- i, not d.
  // Letting |d| inflate this allowance can accept a Newton iterate collapsing
  // onto a duplicate-root solution. Keep that geometry gate above separate
  // from the roundoff amplification at actual equation poles.
  R const one = uni20::scalar_like(zero, 1);
  R poles = one;
  for (std::size_t i = 0; i < roots.size(); ++i)
  {
    poles = std::min(poles, uni20::math::abs(roots[i] * roots[i] + one / 4));
    for (std::size_t j = 0; j < i; ++j)
    {
      auto const d = roots[i] - roots[j];
      poles = std::min(poles, uni20::math::abs(one + d * d));
    }
  }
  return 256 * eps * n / poles;
}

template <uni20::Real R> R energy(int n, std::span<C<R> const> roots, R const& exemplar = R{0})
{
  detail::validate_size(n, roots.size());
  R const zero = detail::validate<R>(roots, exemplar), one = uni20::scalar_like(zero, 1);
  Sum<C<R>> sum(C<R>(zero, zero));
  sum.add(C<R>(n * one / 4, zero));
  for (auto const& z : roots)
    sum.add(-one / ((2 * one) * (z * z + one / 4)));
  auto const value = sum.value();
  if (!uni20::isfinite(value.real()) || !uni20::isfinite(value.imag()) ||
      uni20::math::abs(value.imag()) > 512 * uni20::numeric_limits<R>::epsilon(zero) * n)
    throw std::runtime_error("nonreal or unresolved energy");
  return value.real();
}

// CHM unreduced finite-root raising determinant. The two column products are
// independent for complex roots, NOT conjugates. Fourier raising weight has
// no 2*pi or thermodynamic coefficient. This is not the stabilized pair formula.
// An optional absolute excited momentum index must come from validated family
// labels. It is cross-checked against these roots before applying the exact
// singlet selection rule; unlabeled calls evaluate the determinant only.
template <uni20::Real R>
heisenberg::FormFactor<R> raising_weight(PreparedSinglet<R> const& ground, std::span<C<R> const> lambda,
                                         std::optional<int> validated_momentum_index = std::nullopt)
{
  if (ground.sites() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    throw std::invalid_argument("complex-root helpers require a site count representable as int");
  if (lambda.size() + 1 != ground.sites() / 2) throw std::invalid_argument("invalid singlet/triplet root counts");
  if (validated_momentum_index &&
      (*validated_momentum_index < 0 || static_cast<std::size_t>(*validated_momentum_index) >= ground.sites()))
    throw std::invalid_argument("excited momentum index must be in [0,N)");
  heisenberg::FormFactor<R> result;
  result.status = ground.status();
  if (!ground.ready()) return result;
  int const n = static_cast<int>(ground.sites());
  R const zero = detail::validate<R>(lambda, *ground.working_zero()), one = uni20::scalar_like(zero, 1);
  result.pivot_margin = zero;
  result.status = heisenberg::FormFactorStatus::precision_limit;
  R const eps = uni20::numeric_limits<R>::epsilon(zero), s = separation<R>(lambda, zero);
  if (!(s > 64 * uni20::math::sqrt(eps))) return result;
  R const measured = rational_residual<R>(n, lambda, zero);
  if (!uni20::isfinite(measured)) return result;
  R const tolerance = 4 * residual_tolerance<R>(n, lambda, zero);
  if (measured > tolerance)
  {
    result.status = heisenberg::FormFactorStatus::roots_unconverged;
    return result;
  }
  auto const mu = ground.aba_rapidities();
  auto const m = mu.size();
  auto const norm = bethe::detail::log_absolute_determinant<R>(gaudin<R>(n, lambda, zero), lambda.size(), zero);
  result.pivot_margin = std::min(*ground.ground_pivot_margin(), norm.pivot_margin);
  if (!norm.log_absolute) return result;
  if (validated_momentum_index)
  {
    R const pi = bethe::detail::pi_like(zero);
    // Do not infer an integer sector by rounding a noisy phase. Require the
    // supplied label to agree within the original precision-scaled tolerance,
    // with a strict separation from neighbouring lattice momenta.
    if (!uni20::isfinite(tolerance) || tolerance >= pi / (2 * one * n)) return result;
    C<R> const imaginary(zero, one);
    Sum<C<R>> momentum(C<R>(zero, zero));
    for (auto const& v : lambda)
      momentum.add(uni20::math::log(v + imaginary / (2 * one)) - uni20::math::log(v - imaginary / (2 * one)));
    R const angular_error = uni20::math::abs(phase(momentum.value().imag() - 2 * pi * *validated_momentum_index / n));
    R const magnitude_error = uni20::math::abs(momentum.value().real());
    if (!uni20::isfinite(angular_error) || !uni20::isfinite(magnitude_error) || angular_error > tolerance ||
        magnitude_error > tolerance)
    {
      result.status = heisenberg::FormFactorStatus::roots_unconverged;
      return result;
    }
    // PreparedSinglet validates the centred even-ring ground quantum numbers.
    // Its momentum is zero for N=0 mod 4 and pi for N=2 mod 4. Only zero
    // transfer, not necessarily zero absolute excited momentum, is forbidden.
    int const ground_momentum_index = n % 4 == 0 ? 0 : n / 2;
    if (*validated_momentum_index == ground_momentum_index)
    {
      result.weight = zero;
      result.status = heisenberg::FormFactorStatus::converged;
      return result; // log_weight stays absent: log(0) is not finite.
    }
  }
  std::vector<C<R>> h(m * m, C<R>(zero, zero));
  Sum<R> logarithm(zero);
  logarithm.add(*ground.ground_log_prefactor() - *norm.log_absolute);
  C<R> const imaginary(zero, one);
  for (std::size_t b = 0; b < lambda.size(); ++b)
  {
    auto const v = lambda[b];
    Sum<C<R>> minus(C<R>(zero, zero)), plus(C<R>(zero, zero));
    plus.add((n * one) * (uni20::math::log(v + imaginary / (2 * one)) - uni20::math::log(v - imaginary / (2 * one))));
    for (R const& u : mu)
    {
      minus.add(uni20::math::log(u - v - imaginary));
      plus.add(uni20::math::log(u - v + imaginary));
    }
    R const scale = std::max(minus.value().real(), plus.value().real());
    auto const a = uni20::math::exp(minus.value() - scale), p = uni20::math::exp(plus.value() - scale);
    logarithm.add(2 * scale);
    for (std::size_t row = 0; row < m; ++row)
    {
      auto const d = mu[row] - v;
      if (uni20::math::abs(d) <= 128 * eps) return result;
      h[row * m + b] = (a / (d - imaginary) - p / (d + imaginary)) / d;
    }
  }
  for (std::size_t row = 0; row < m; ++row)
    h[row * m + m - 1] = one / (mu[row] * mu[row] + one / 4);
  auto const determinant = bethe::detail::log_absolute_determinant<R>(std::move(h), m, zero);
  result.pivot_margin = std::min(result.pivot_margin, determinant.pivot_margin);
  if (!determinant.log_absolute) return result;
  logarithm.add(2 * *determinant.log_absolute);
  for (auto const& v : lambda)
    logarithm.add(-uni20::math::log(uni20::math::abs(v * v + one / 4)));
  for (std::size_t i = 0; i < lambda.size(); ++i)
    for (std::size_t j = 0; j < i; ++j)
    {
      auto const d = lambda[i] - lambda[j];
      logarithm.add(-uni20::math::log(uni20::math::abs(one + d * d)));
    }
  R const weight = uni20::math::exp(logarithm.value());
  if (!uni20::isfinite(weight) || weight <= zero || weight > n * one) return result;
  result.weight = weight;
  result.log_weight = logarithm.value();
  result.status = heisenberg::FormFactorStatus::converged;
  return result;
}
} // namespace bethe::experimental::xxx::complex_roots
