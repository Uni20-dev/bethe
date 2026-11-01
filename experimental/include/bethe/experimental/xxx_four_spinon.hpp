// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/experimental/xxx_four_spinon_state.hpp>

namespace bethe::experimental::xxx
{
namespace four_detail
{
// Direct real nonsymmetric representation. The pair self factor times the
// original Gaudin determinant is this determinant, including finite delta.
template <uni20::Real R> std::vector<R> gaudin_matrix(FourSpinonState<R> const& state)
{
  auto const reals = state.real_rapidities();
  auto const r = reals.size(), s = r + 2;
  R const zero = state.working_zero(), one = uni20::scalar_like(zero, 1);
  int const n = state.request().sites;
  C<R> const u(state.centre(), one / 2 + state.deviation());
  R const d = 2 * state.deviation() * (one + state.deviation());
  std::vector<R> matrix(s * s, zero);
  Sum<C<R>> upper(C<R>(zero, zero));
  upper.add(n * one / (u * u + one / 4));
  for (std::size_t i = 0; i < r; ++i)
  {
    R const x = reals[i];
    C<R> const plus = 2 * one / (one + (x - u) * (x - u));
    upper.add(-plus);
    Sum<R> diagonal(zero);
    diagonal.add(n * one / (x * x + one / 4) - 2 * plus.real());
    for (std::size_t j = 0; j < r; ++j)
      if (i != j)
      {
        matrix[i * s + j] = kernel(2, x - reals[j]);
        diagonal.add(-matrix[i * s + j]);
      }
    matrix[i * s + i] = diagonal.value();
    matrix[i * s + r] = matrix[r * s + i] = 2 * plus.real();
    matrix[i * s + r + 1] = -2 * plus.imag();
    matrix[(r + 1) * s + i] = d * plus.imag();
  }
  matrix[r * s + r] = 2 * upper.value().real();
  matrix[r * s + r + 1] = -2 * upper.value().imag();
  matrix[(r + 1) * s + r] = d * upper.value().imag();
  matrix[(r + 1) * s + r + 1] = 2 * one + d * upper.value().real();
  return matrix;
}

// Unit-modulus phases and a determinant-preserving pair-column addition make
// this matrix real. Assemble the analytic result, not rounded differences.
template <uni20::Real R>
std::vector<R> raising_matrix(std::span<R const> mu, FourSpinonState<R> const& state, Sum<R>& logarithm)
{
  R const zero = state.working_zero(), one = uni20::scalar_like(zero, 1), delta = state.deviation();
  auto const reals = state.real_rapidities();
  auto const m = mu.size(), r = reals.size();
  int const n = state.request().sites;
  std::vector<R> h(m * m, zero);
  for (std::size_t b = 0; b < r; ++b)
  {
    R const x = reals[b];
    Sum<C<R>> log_product(C<R>(zero, zero));
    for (auto const& v : mu)
      log_product.add(uni20::math::log(C<R>(v - x, -one)));
    auto const logp = log_product.value();
    logarithm.add(2 * logp.real());
    R const angle = logp.imag() - n * uni20::math::atan2(one / 2, x);
    R const sine = uni20::math::sin(angle), cosine = uni20::math::cos(angle);
    for (std::size_t a = 0; a < m; ++a)
    {
      R const d = mu[a] - x;
      if (uni20::math::abs(d) <= 128 * uni20::numeric_limits<R>::epsilon(zero))
        throw Failure(FourSpinonStatus::precision_limit, "coincident roots in raising matrix");
      h[a * m + b] = 2 * (d * sine + cosine) / (d * (d * d + one));
    }
  }
  C<R> const u(state.centre(), one / 2 + delta), imaginary(zero, one);
  Sum<C<R>> loga(C<R>(zero, zero)), plus(C<R>(zero, zero)), minus(C<R>(zero, zero));
  for (auto const& x : reals)
    loga.add(uni20::math::log(u - x + imaginary) - uni20::math::log(u - x - imaginary));
  for (auto const& x : mu)
  {
    plus.add(uni20::math::log(x - u + imaginary));
    minus.add(uni20::math::log(x - u - imaginary));
  }
  C<R> const logb = minus.value() - loga.value() - plus.value();
  R const log_one_plus = uni20::math::log1p(delta), log_abs_delta = state.log_abs_deviation();
  R const scale_upper = std::max(zero, logb.real());
  R const scale_lower = std::max(zero, log_abs_delta + logb.real() - log_one_plus);
  C<R> const b = uni20::math::exp(logb - scale_upper);
  C<R> const small_b = (delta < zero ? -one : one) *
                       uni20::math::exp(uni20::math::conj(logb) + (log_abs_delta - log_one_plus - scale_lower));
  logarithm.add(2 * (log_one_plus + loga.value().real() + 2 * plus.value().real() + scale_upper + scale_lower));
  for (std::size_t a = 0; a < m; ++a)
  {
    R const d = mu[a] - state.centre(), base = d * d + one / 4 - delta * delta;
    R const denominator = base * base + 4 * delta * delta * d * d;
    C<R> const divided(zero, 4 * d / denominator);
    C<R> const kminus = one / ((mu[a] - u) * (mu[a] - u - imaginary)), kplus = uni20::math::conj(kminus);
    C<R> const upper =
        uni20::math::exp(-scale_upper) * divided + (-b * kminus + uni20::math::conj(b) * kplus) / (one + delta);
    h[a * m + r] = upper.imag();
    // H_l + delta/2*exp(scale_upper-scale_lower)*H_u is real. Its real
    // part needs neither that potentially overflowing ratio nor subtraction.
    h[a * m + r + 1] = uni20::math::exp(-scale_lower) * base / denominator - (small_b * kplus).real();
    h[a * m + m - 1] = one / (mu[a] * mu[a] + one / 4);
  }
  return h;
}
} // namespace four_detail

template <uni20::Real R>
heisenberg::FormFactor<R> four_spinon_raising_weight(PreparedSinglet<R> const& ground, FourSpinonState<R> const& state)
{
  if (!state.ready() || ground.sites() != state.request().sites)
    throw std::invalid_argument("finite-four form factor needs matching ready states");
  heisenberg::FormFactor<R> result;
  result.status = ground.status();
  if (!ground.ready()) return result;
  R const zero = *ground.working_zero(), one = uni20::scalar_like(zero, 1);
  bethe::detail::require_same_precision(state.working_zero(), zero);
  result.pivot_margin = zero;
  result.status = heisenberg::FormFactorStatus::precision_limit;
  auto const r = state.real_rapidities().size(), m = ground.aba_rapidities().size();
  if (r + 3 != m) throw std::invalid_argument("finite-four root counts do not match singlet");
  auto const gaudin = bethe::detail::log_absolute_determinant<R>(four_detail::gaudin_matrix(state), r + 2, zero);
  result.pivot_margin = std::min(*ground.ground_pivot_margin(), gaudin.pivot_margin);
  if (!gaudin.log_absolute) return result;
  four_detail::Sum<R> logarithm(zero);
  logarithm.add(*ground.ground_log_prefactor() - *gaudin.log_absolute);
  std::vector<R> h;
  try
  {
    h = four_detail::raising_matrix(ground.aba_rapidities(), state, logarithm);
  }
  catch (four_detail::Failure const&)
  {
    return result;
  }
  auto const determinant = bethe::detail::log_absolute_determinant<R>(std::move(h), m, zero);
  result.pivot_margin = std::min(result.pivot_margin, determinant.pivot_margin);
  if (!determinant.log_absolute) return result;
  logarithm.add(2 * *determinant.log_absolute);
  auto const lambda = state.roots();
  for (auto const& z : lambda)
    logarithm.add(-uni20::math::log(uni20::math::abs(z * z + one / 4)));
  for (std::size_t i = 0; i < lambda.size(); ++i)
    for (std::size_t j = 0; j < i; ++j)
      if (!(i == r + 1 && j == r))
      {
        auto const d = lambda[i] - lambda[j];
        logarithm.add(-uni20::math::log(uni20::math::abs(one + d * d)));
      }
  R const weight = uni20::math::exp(logarithm.value());
  if (!uni20::isfinite(weight) || weight <= zero || weight > ground.sites() * one) return result;
  result.weight = weight;
  result.log_weight = logarithm.value();
  result.status = heisenberg::FormFactorStatus::converged;
  return result;
}
} // namespace bethe::experimental::xxx
