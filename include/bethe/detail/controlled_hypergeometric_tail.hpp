// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/spinon_numerics.hpp>
#include <limits>
#include <span>
#include <vector>

namespace bethe::detail
{
struct TailDomainError : std::invalid_argument
{
    using std::invalid_argument::invalid_argument;
};

template <uni20::Real R> struct ControlledTail
{
    uni20::complex<R> multiplier;
    R error_multiplier, arithmetic_indicator;
};

/// Geometry-independent combinatorics at one explicit precision/order.
/// Immutable after construction; share the plan, never the scratch workspace.
/// Retain the original recurrences/rounding, including computed zeros.
template <uni20::Real R> class ControlledTailPlan {
  public:
    ControlledTailPlan(R const& precision_anchor, int order) : order_(order)
    {
      if (order < 1 || order > std::numeric_limits<int>::max() - 2)
        throw std::invalid_argument("invalid tail plan order");
      if (!uni20::isfinite(precision_anchor)) throw std::invalid_argument("nonfinite tail precision anchor");
      if constexpr (uni20::has_runtime_precision_v<R>)
        if (precision_anchor.precision().is_exact())
          throw std::invalid_argument("tail plan requires finite explicit precision");
      one_ = (precision_anchor - precision_anchor) + R{1};
      for (int k = 0; k < order; ++k)
      {
        auto& row = generalized_.emplace_back(order + 2 - k, one_);
        for (int j = 1; j < int(row.size()); ++j)
          row[j] = row[j - 1] * (R(2 - k - j) / (one_ * R(j)));
      }
      for (int degree = 0; degree <= order; ++degree)
      {
        auto& row = binomial_.emplace_back(degree + 1, one_);
        for (int k = 1; k <= degree; ++k)
          row[k] = row[k - 1] * R(degree - k + 1) / R(k);
      }
    }
    int order() const { return order_; }
    R const& one() const { return one_; }
    std::span<R const> generalized(int k) const { return generalized_.at(k); }
    std::span<R const> binomial(int degree) const { return binomial_.at(degree); }

  private:
    int order_;
    R one_;
    std::vector<std::vector<R>> generalized_, binomial_;
};

/// Caller-owned buffers, reusable sequentially even after a rejected tail.
/// Every active element is reset at the input precision on every invocation.
template <uni20::Real R> struct ControlledTailWorkspace
{
    using C = uni20::complex<R>;
    std::vector<C> ratio, coefficients, polynomial_a, polynomial_b, base, left, translated, right, defect;
};

/// General rational-ratio tail. Analytic remainder inequality evaluated in
/// ordinary floating point, NOT an interval enclosure. See the independent
/// Python reference_hypergeometric_tail.py for the derivation and oracle.
template <uni20::Real R>
ControlledTail<R>
controlled_hypergeometric_tail(std::span<uni20::complex<R> const> a, std::span<uni20::complex<R> const> b, int power,
                               int start, ControlledTailPlan<R> const& plan, ControlledTailWorkspace<R>& workspace)
{
  using C = uni20::complex<R>;
  using bethe::detail::SpinonSum;
  using uni20::math::abs;
  using uni20::math::exp;
  int const order = plan.order();
  if (a.empty() || a.size() != b.size() || power < 3 || start < 2 || order < 1)
    throw std::invalid_argument("require equal nonempty parameters, power>=3, start>=2, order>=1");
  if (order > std::numeric_limits<int>::max() - 2 || power > std::numeric_limits<int>::max() - order)
    throw std::invalid_argument("tail order/power exceed index representation");
  C zero = a.front() - a.front();
  R z = zero.real(), one = z + R{1}, eps = uni20::numeric_limits<R>::epsilon(one), n = one * R(start);
  // Approximate MPFR operands of different precisions must not share a plan.
  (void)(one * plan.one());
  SpinonSum<C> exponent(zero);
  SpinonSum<R> moduli(z), constant_sum(z);
  R denominator = one;
  for (std::size_t k = 0; k < a.size(); ++k)
  {
    if (!uni20::isfinite(a[k]) || !uni20::isfinite(b[k])) throw std::invalid_argument("tail parameters must be finite");
    exponent.add(b[k] - a[k]);
    for (auto x : {a[k], b[k]})
    {
      R m = abs(x);
      if (m >= n) throw TailDomainError("tail start must exceed every parameter modulus");
      moduli.add(m);
      constant_sum.add(m * m / (one - m / n));
    }
    denominator *= one - abs(b[k]) / n;
  }
  R decay = exponent.value().real();
  if (decay <= R{2} || abs(exponent.value() - R(power)) > R{64} * eps * std::max(one, moduli.value()))
    throw std::invalid_argument("power must equal sum(b-a)");

  auto& ratio = workspace.ratio;
  ratio.assign(order + 2, zero);
  ratio[0] = one;
  for (auto const& x : a)
    for (int k = order + 1; k > 0; --k)
      ratio[k] += x * ratio[k - 1];
  for (auto const& x : b)
    for (int k = 1; k <= order + 1; ++k)
      ratio[k] -= x * ratio[k - 1];
  auto& c = workspace.coefficients;
  c.assign(order + 1, zero);
  c[0] = C(one / R(power - 1), z);
  for (int degree = 1; degree <= order; ++degree)
  {
    SpinonSum<C> terms(zero);
    for (int k = 0; k < degree; ++k)
    {
      auto choose = plan.generalized(k);
      C value = zero;
      for (int j = 0; j < degree + 2 - k; ++j)
      {
        value += ratio[degree + 1 - k - j] * choose[j];
      }
      terms.add(c[k] * value);
    }
    c[degree] = terms.value() / R(power + degree - 1);
  }
  auto multiply = [&](auto const& x, auto const& y, std::vector<C>& result) {
    result.assign(x.size() + y.size() - 1, zero);
    for (std::size_t j = 0; j < x.size(); ++j)
      for (std::size_t k = 0; k < y.size(); ++k)
        result[j + k] += x[j] * y[k];
  };
  auto polynomial = [&](auto parameters, std::vector<C>& result) {
    result.assign(parameters.size() + 1, zero);
    result[0] = one;
    std::size_t degree = 0;
    for (auto const& x : parameters)
    {
      ++degree;
      for (std::size_t k = degree; k > 0; --k)
        // Preserve the convolution's addition order: old[k-1]*x + old[k].
        result[k] = result[k - 1] * x + result[k];
    }
  };
  auto& base = workspace.base;
  auto& left = workspace.left;
  auto& translated = workspace.translated;
  auto& right = workspace.right;
  auto& defect = workspace.defect;
  polynomial(a, workspace.polynomial_a);
  polynomial(b, workspace.polynomial_b);
  multiply(workspace.polynomial_b, plan.binomial(order - 1), base);
  multiply(base, c, left);
  translated.assign(order + 1, zero);
  for (int k = 0; k <= order; ++k)
  {
    auto choose = plan.binomial(order - k);
    for (int j = 0; j <= order - k; ++j)
      translated[k + j] += c[k] * choose[j];
  }
  multiply(workspace.polynomial_a, translated, right);
  defect.assign(std::max({left.size(), right.size(), base.size() + 1}), zero);
  for (std::size_t k = 0; k < defect.size(); ++k)
    defect[k] = (k < left.size() ? left[k] : zero) - (k < right.size() ? right[k] : zero) -
                (k > 0 && k <= base.size() ? base[k - 1] : zero);
  // Do not erase the mathematically vanishing low-order coefficients: their
  // finite-precision residuals are part of this actual computed defect.
  auto weighted = [&](std::vector<C> const& values) {
    SpinonSum<R> sum(z);
    R factor = n;
    for (std::size_t k = 0; k < values.size(); ++k)
    {
      sum.add(abs(values[k]) * factor * (one + n / (decay + R(k) - R{2})));
      factor /= n;
    }
    return sum.value();
  };
  SpinonSum<C> sum(zero);
  R factor = n;
  for (auto const& value : c)
  {
    sum.add(value * factor);
    factor /= n;
  }
  C q = sum.value();
  R constant = (constant_sum.value() + decay / (one - one / n)) / R{2};
  R envelope = exp(constant * (one / n + one / (n * n))) / denominator;
  R bound = envelope * weighted(defect);
  base.insert(base.begin(), zero);
  R cancellation = weighted(left) + weighted(right) + weighted(base);
  R cost = one * R(a.size() + order);
  R arithmetic = R{64} * eps * cost * cost * (abs(q) + envelope * cancellation);
  if (!uni20::isfinite(q) || !uni20::isfinite(bound) || !uni20::isfinite(arithmetic))
    throw SpinonFunctionError(SpinonFunctionFailure::arithmetic, "nonfinite controlled tail");
  return {q, bound, arithmetic};
}

/// Convenience path for one tail; repeated callers should share a plan/workspace.
template <uni20::Real R>
ControlledTail<R> controlled_hypergeometric_tail(std::span<uni20::complex<R> const> a,
                                                 std::span<uni20::complex<R> const> b, int power, int start, int order)
{
  if (a.empty()) throw std::invalid_argument("tail parameters must not be empty");
  ControlledTailPlan<R> plan(a.front().real(), order);
  ControlledTailWorkspace<R> workspace;
  return controlled_hypergeometric_tail<R>(a, b, power, start, plan, workspace);
}
} // namespace bethe::detail
