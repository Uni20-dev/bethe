// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <array>
#include <bethe/detail/gauss_legendre.hpp>
#include <bethe/detail/spinon_numerics.hpp>
#include <bethe/solver.hpp>

namespace bethe::heisenberg::detail
{
template <uni20::Real R> struct FourSpinonShiftIntegral
{
    R value, arithmetic_indicator, quadrature_difference, omitted_tail_bound;
    bool used_quadrature = false;
};

/// Common-shift geometry only. W4*d^4p = F(gaps)*sin^2(P/2)*d^4beta.
/// At the representative P=pi, F=W4*exp(log_rapidity_measure). The shift
/// integrals are A0 below and A1=pi^2. Neither F nor a physical prefactor is
/// supplied by this class.
template <uni20::Real R> class FourSpinonCommonShift {
  public:
    struct Representative
    {
        std::array<R, 4> beta;
        R log_rapidity_measure;
    };

    explicit FourSpinonCommonShift(R exemplar)
        : zero_(bethe::detail::precision_zero(exemplar)), one_(zero_ + R{1}), pi_(R{4} * uni20::math::atan(one_)),
          eps_(uni20::numeric_limits<R>::epsilon(one_))
    {}

    Representative representative(std::array<R, 3> const& gaps) const
    {
      check(gaps);
      std::array<R, 4> delta{-gaps[0] - gaps[1] - gaps[2], -gaps[1] - gaps[2], -gaps[2], zero_};
      // Re prod(1+i*exp(delta_j+eta))=1-e2*exp(2eta)+e4*exp(4eta).
      // The smaller positive root has P=pi. Use log-sum-exp and a rationalized
      // root, so neither a large common shift nor close roots are subtracted.
      R sum = zero_, log_e4 = zero_;
      R const largest_pair = delta[2];
      for (int j = 0; j < 4; ++j)
      {
        log_e4 += delta[j];
        for (int i = 0; i < j; ++i)
          sum += uni20::math::exp(delta[i] + delta[j] - largest_pair);
      }
      R const log_e2 = largest_pair + uni20::math::log(sum);
      R const discriminant = one_ - R{4} * uni20::math::exp(log_e4 - R{2} * log_e2);
      R const shift =
          (uni20::math::log(R{2} * one_) - log_e2 - uni20::math::log1p(uni20::math::sqrt(discriminant))) / R{2};
      Representative result{delta, zero_};
      for (auto& beta : result.beta)
      {
        beta += shift;
        R const a = uni20::math::abs(beta);
        result.log_rapidity_measure -=
            a + uni20::math::log1p(uni20::math::exp(-R{2} * a)) - uni20::math::log(R{2} * one_);
      }
      if (!uni20::isfinite(result.log_rapidity_measure)) throw std::runtime_error("common-shift measure overflow");
      return result;
    }

    FourSpinonShiftIntegral<R> zeroth(std::array<R, 3> const& gaps) const
    {
      check(gaps);
      auto const& [a, b, c] = gaps;
      R const span = a + b + c;
      if (span == zero_) return {R{8} * one_ / R{3}, R{8} * eps_, zero_, zero_};
      if (a + b > zero_ && b + c > zero_)
      {
        // Regroup partial fractions: the negative middle term decays as
        // b*exp(-b), rather than subtracting terms growing with the middle gap.
        R const d = one_ / uni20::math::tanh((a + b) / R{2});
        R const e = one_ / uni20::math::tanh((b + c) / R{2});
        R const f = one_ / uni20::math::tanh(span / R{2});
        R const x = pair_shift(a) * d * f, z = pair_shift(c) * f * e;
        R const be = b == zero_ ? one_ : b * uni20::math::exp(-b) / (-uni20::math::expm1(-b));
        R const y = R{2} * be * (one_ + uni20::math::exp(-a)) * (one_ + uni20::math::exp(-c)) *
                    (-uni20::math::expm1(-span - b)) /
                    ((-uni20::math::expm1(-a - b)) * (-uni20::math::expm1(-span)) * (-uni20::math::expm1(-b - c)));
        bethe::detail::SpinonSum<R> sum(zero_);
        sum.add(x);
        sum.add(-y);
        sum.add(z);
        R const value = sum.value(), arithmetic = R{128} * eps_ * (x + y + z);
        if (uni20::isfinite(value) && uni20::isfinite(arithmetic) && value > zero_ && arithmetic <= tolerance() * value)
          return {value, arithmetic, zero_, zero_};
      }
      return quadrature(gaps);
    }

    R first() const { return pi_ * pi_; }

    /// Positive quadrature is also an independent reference for the closed
    /// expression. Its finite shift tails have an explicit exponential bound;
    /// order differences and arithmetic indicators are separate diagnostics.
    FourSpinonShiftIntegral<R> quadrature(std::array<R, 3> const& gaps) const
    {
      check(gaps);
      std::array<R, 4> const delta{zero_, gaps[0], gaps[0] + gaps[1], gaps[0] + gaps[1] + gaps[2]};
      R const tol = tolerance();
      R const tail = -uni20::math::log(tol / (R{64} * one_)) / R{2};
      R const left = -delta.back() - tail, length = delta.back() + R{2} * tail;
      // Construct panels without converting a runtime-precision real to a
      // native integer. The explicit cap also bounds exceptional fallback work.
      std::size_t panels = 1;
      while (uni20::scalar_like(one_, panels) < length)
      {
        if (panels == 16384) throw std::runtime_error("common-shift quadrature panel budget exhausted");
        ++panels;
      }
      R const width = length / uni20::scalar_like(one_, panels);
      R const tail_bound = R{16} * uni20::math::exp(-R{2} * tail);
      R previous = zero_;
      for (std::size_t order : {16, 32, 64, 128})
      {
        auto const rule = bethe::detail::gauss_legendre<R>(order, uni20::precision_of(one_));
        bethe::detail::SpinonSum<R> sum(zero_);
        for (std::size_t panel = 0; panel < panels; ++panel)
          for (std::size_t j = 0; j < order; ++j)
          {
            R const eta = left + width * (uni20::scalar_like(one_, panel) + (one_ + rule.x[j]) / R{2});
            R reduced = zero_;
            unsigned positive = 0;
            for (auto const& d : delta)
            {
              R const beta = eta + d;
              if (beta > zero_)
              {
                ++positive;
                reduced -= uni20::math::atan(uni20::math::exp(-beta));
              }
              else
                reduced += uni20::math::atan(uni20::math::exp(beta));
            }
            // sin(P/2) without subtracting a rounded multiple of pi in a tail.
            R const amplitude = positive % 2 ? uni20::math::cos(reduced) : uni20::math::sin(reduced);
            sum.add(width / R{2} * rule.w[j] * amplitude * amplitude);
          }
        R const value = sum.value(), difference = uni20::math::abs(value - previous);
        R const arithmetic = R{64} * eps_ * value;
        if (order > 16 && difference + arithmetic + tail_bound <= tol * std::max(one_, value))
          return {value, arithmetic, difference, tail_bound, true};
        previous = value;
      }
      throw std::runtime_error("common-shift quadrature did not converge within its order budget");
    }

  private:
    void check(std::array<R, 3> const& gaps) const
    {
      for (auto const& gap : gaps)
      {
        bethe::detail::require_same_precision(gap, one_);
        if (!uni20::isfinite(gap) || gap < zero_)
          throw std::invalid_argument("require finite nonnegative rapidity gaps");
      }
      R const span = gaps[0] + gaps[1] + gaps[2];
      if (!uni20::isfinite(span) || span > R{8192})
        throw std::invalid_argument("common-shift span exceeds budget 8192");
    }
    R pair_shift(R const& gap) const { return gap == zero_ ? R{2} * one_ : gap / uni20::math::tanh(gap / R{2}); }
    R tolerance() const { return std::max(R{4096} * eps_, uni20::math::sqrt(eps_) / R{100}); }
    R zero_, one_, pi_, eps_;
};

/// Three-gap cubature after integrating the common shift. Cutoffs restrict
/// the gaps individually, not the absolute rapidities. Emitted measures
/// include the ordered-domain 4!, but not the observable's (2*pi)^-2.
template <uni20::Real R> class FourSpinonGapGeometry {
  public:
    struct Node
    {
        std::array<R, 4> beta;
        R log_measure;
        FourSpinonShiftIntegral<R> zeroth_shift;
        R first_shift;
    };

    FourSpinonGapGeometry(std::size_t order, std::array<R, 3> cutoffs, std::size_t max_nodes = 1000000)
        : cutoffs_(std::move(cutoffs)), zero_(bethe::detail::precision_zero(cutoffs_[0])), one_(zero_ + R{1}),
          shift_(one_)
    {
      auto parameters = checked_parameters(order, cutoffs_, max_nodes);
      count_ = parameters.count;
      auto const& cuts = parameters.cuts;
      auto const rule = bethe::detail::gauss_legendre<R>(order, uni20::precision_of(one_));
      for (int axis = 0; axis < 3; ++axis)
        for (std::size_t panel = 1; panel < cuts[axis].size(); ++panel)
        {
          R const left = cuts[axis][panel - 1], width = cuts[axis][panel] - left;
          for (std::size_t j = 0; j < order; ++j)
          {
            R const gap = left + width * (one_ + rule.x[j]) / R{2};
            R const weight = width * rule.w[j] / R{2};
            if (!(gap > left && gap < cuts[axis][panel] && weight > zero_))
              throw std::runtime_error("gap quadrature node is unresolved at working precision");
            nodes_[axis].push_back({gap, uni20::math::log(weight)});
          }
        }
    }

    // Validate a deferred query without constructing quadrature nodes or evaluating shifts.
    static void validate_parameters(std::size_t order, std::array<R, 3> const& cutoffs, std::size_t max_nodes)
    {
      (void)checked_parameters(order, cutoffs, max_nodes);
    }

    std::size_t size() const { return count_; }
    std::array<R, 3> const& cutoffs() const { return cutoffs_; }

    Node node(std::size_t index) const
    {
      if (index >= count_) throw std::out_of_range("gap quadrature index outside grid");
      auto const& c = nodes_[2][index % nodes_[2].size()];
      index /= nodes_[2].size();
      auto const& b = nodes_[1][index % nodes_[1].size()];
      auto const& a = nodes_[0][index / nodes_[1].size()];
      std::array<R, 3> const gaps{a.value, b.value, c.value};
      auto point = shift_.representative(gaps);
      R const log_measure =
          uni20::math::log(R{24} * one_) + a.log_weight + b.log_weight + c.log_weight + point.log_rapidity_measure;
      return {std::move(point.beta), log_measure, shift_.zeroth(gaps), shift_.first()};
    }

    template <class F> void integrate(F&& sample) const
    {
      for (std::size_t i = 0; i < count_; ++i)
        sample(node(i));
    }

  private:
    struct Parameters
    {
        std::size_t count = 1;
        std::array<std::vector<R>, 3> cuts;
    };
    static Parameters checked_parameters(std::size_t order, std::array<R, 3> const& cutoffs, std::size_t max_nodes)
    {
      R const zero = bethe::detail::precision_zero(cutoffs[0]), one = zero + R{1};
      if (order < 2 || order > 256 || max_nodes == 0)
        throw std::invalid_argument("require gap order 2..256 and positive node budget");
      R span = zero;
      Parameters result;
      for (int axis = 0; axis < 3; ++axis)
      {
        R const& cutoff = cutoffs[axis];
        bethe::detail::require_same_precision(cutoff, one);
        if (!uni20::isfinite(cutoff) || cutoff <= zero || cutoff > R{8192})
          throw std::invalid_argument("require positive finite gap cutoffs with total span at most 8192");
        span += cutoff;
        result.cuts[axis].push_back(zero);
        for (R end = R{2} * one; end < cutoff; end *= R{2})
          result.cuts[axis].push_back(end);
        result.cuts[axis].push_back(cutoff);
        std::size_t const n = order * (result.cuts[axis].size() - 1);
        if (n > max_nodes / result.count) throw std::invalid_argument("gap quadrature exceeds node budget");
        result.count *= n;
      }
      if (span > R{8192}) throw std::invalid_argument("common-shift span exceeds budget 8192");
      return result;
    }
    struct GapNode
    {
        R value, log_weight;
    };
    std::array<R, 3> cutoffs_;
    R zero_, one_;
    FourSpinonCommonShift<R> shift_;
    std::size_t count_;
    std::array<std::vector<GapNode>, 3> nodes_;
};
} // namespace bethe::heisenberg::detail
