// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <algorithm>
#include <array>
#include <limits>
#include <span>
#include <stdexcept>
#include <uni20/core/math.hpp>
#include <uni20/core/runtime_precision.hpp>
#include <utility>
#include <vector>

namespace bethe::heisenberg::detail
{
enum class SpinonProposalKind
{
  background,
  cauchy,
  one_pair,
  two_pairs
};

/// Normalized proposal for the FOUR free rapidities in a six-spinon coarea
/// chart. Its density is against d^4 beta, not momenta or physical Szz.
/// Immutable after construction. Callers supply five open uniforms, so random
/// stream ownership and correlations across frequencies stay outside this class.
template <uni20::Real R> class SpinonJointProposal {
  public:
    struct Component
    {
        SpinonProposalKind kind;
        R weight, scale, gap_scale; // scales are used only by the relevant kind
    };
    struct Sample
    {
        std::array<R, 4> beta;
        R log_density; // complete mixture, including overlapping components
        std::size_t component, pairing;
    };

    SpinonJointProposal(R cutoff, R soft_fraction, std::vector<Component> components)
        : cutoff_(std::move(cutoff)), soft_(std::move(soft_fraction)), components_(std::move(components)),
          one_(uni20::scalar_like(cutoff_, 1)), zero_(uni20::scalar_like(cutoff_, 0))
    {
      if constexpr (uni20::has_runtime_precision_v<R>)
        if (one_.precision().is_exact()) throw std::invalid_argument("proposal requires explicit finite precision");
      check_scalar(cutoff_);
      check_scalar(soft_);
      if (!(cutoff_ > zero_ && soft_ > zero_ && soft_ < one_) || components_.empty() ||
          components_.front().kind != SpinonProposalKind::background)
        throw std::invalid_argument("require positive cutoff, 0<soft<1, and a background component");
      pi_ = R{4} * uni20::math::atan(one_);
      eps_ = uni20::numeric_limits<R>::epsilon(one_);
      log_two_ = uni20::math::log(R{2} * one_);
      lower_ = R{2} * uni20::math::atan(uni20::math::exp(-cutoff_));
      width_ = pi_ - R{2} * lower_;
      soft_norm_ = R{2} * uni20::math::atan(cutoff_ / R{2});
      if (!(width_ > zero_ && soft_norm_ > zero_))
        throw std::runtime_error("proposal cutoff is unresolved at this precision");
      total_weight_ = zero_;
      for (auto const& c : components_)
      {
        check_scalar(c.weight);
        if (!(c.weight > zero_)) throw std::invalid_argument("proposal weights must be positive");
        if (c.kind != SpinonProposalKind::background)
        {
          if (c.kind != SpinonProposalKind::cauchy && c.kind != SpinonProposalKind::one_pair &&
              c.kind != SpinonProposalKind::two_pairs)
            throw std::invalid_argument("unknown proposal kind");
          check_scalar(c.scale);
          if (!(c.scale > zero_)) throw std::invalid_argument("proposal scale must be positive");
          if (c.kind != SpinonProposalKind::cauchy)
          {
            check_scalar(c.gap_scale);
            if (!(c.gap_scale > zero_)) throw std::invalid_argument("proposal gap scale must be positive");
          }
        }
        total_weight_ += c.weight;
      }
      if (!uni20::isfinite(total_weight_) || uni20::math::abs(total_weight_ - one_) > R{16} * eps_)
        throw std::invalid_argument("proposal weights must sum to one");
    }

    R const& cutoff() const { return cutoff_; }

    R native_log_density(R const& beta) const
    {
      check_beta(beta);
      std::array<R, 2> terms{uni20::math::log(one_ - soft_) - uni20::math::log(width_) - log_cosh(beta),
                             uni20::math::log(soft_) + cauchy_log_density(beta, R{2} * one_, cutoff_)};
      return log_sum(terms);
    }

    R native_cdf(R const& beta) const
    {
      check_beta(beta);
      R const depth = uni20::math::abs(beta);
      // Evaluate both lower-tail probabilities without subtracting close
      // angles or 1/2. The Cauchy difference uses atan2 and scaled arguments
      // so cutoff*depth cannot overflow. Reflection is applied only once.
      R const difference = uni20::math::exp(-depth) * (-uni20::math::expm1(depth - cutoff_));
      R const flat = R{2} * uni20::math::atan(difference / (one_ + uni20::math::exp(-depth - cutoff_))) / width_;
      R const cauchy = uni20::math::atan2((cutoff_ - depth) / cutoff_, R{2} / cutoff_ + depth / R{2}) / soft_norm_;
      R const lower = (one_ - soft_) * flat + soft_ * cauchy;
      return beta <= zero_ ? lower : one_ - lower;
    }

    R native_inverse(R const& uniform) const
    {
      check_uniform(uniform);
      bool const reflected = uniform > one_ / R{2};
      R const u = reflected ? one_ - uniform : uniform;
      if (u == one_ / R{2}) return zero_;
      R lo = -cutoff_, hi = zero_, x = lo + R{2} * u * cutoff_;
      auto const digits = uni20::numeric_limits<R>::digits(one_);
      if (digits > (std::numeric_limits<int>::max() - 32) / 2)
        throw std::invalid_argument("proposal inverse precision exceeds iteration index range");
      int const bits = int(digits);
      for (int iteration = 0; iteration < 2 * bits + 32; ++iteration)
      {
        R const residual = native_cdf(x) - u;
        R const density = uni20::math::exp(native_log_density(x));
        // A representable beta step moves the CDF by PDF*ulp(beta). A solely
        // relative residual in u can be smaller than that attainable change,
        // even with a stable CDF and an ordinary open midpoint input.
        R const tolerance = R{32} * eps_ * (u + density * std::max(one_, uni20::math::abs(x)));
        if (uni20::math::abs(residual) <= tolerance) return reflected ? -x : x;
        if (residual < zero_)
          lo = x;
        else
          hi = x;
        R const newton = x - residual / density;
        R const trial = iteration % 4 != 3 && lo < newton && newton < hi ? newton : (lo + hi) / R{2};
        if (!(lo < trial && trial < hi)) throw std::runtime_error("proposal inverse exhausted precision");
        x = trial;
      }
      throw std::runtime_error("proposal inverse did not converge");
    }

    R log_density(std::array<R, 4> const& beta) const
    {
      std::array<R, 4> native;
      R native_product = zero_;
      for (std::size_t i = 0; i < beta.size(); ++i)
      {
        native[i] = native_log_density(beta[i]);
        native_product += native[i];
      }
      std::vector<R> terms;
      for (auto const& c : components_)
      {
        R value = native_product;
        if (c.kind == SpinonProposalKind::cauchy)
        {
          value = zero_;
          for (auto const& b : beta)
            value += cauchy_log_density(b, c.scale, cutoff_);
        }
        else if (c.kind == SpinonProposalKind::one_pair || c.kind == SpinonProposalKind::two_pairs)
        {
          std::array<R, 6> pairs;
          for (std::size_t a = 0; a < pairs.size(); ++a)
          {
            auto const [i, j] = pairs_[a];
            pairs[a] = pair_log_density(beta[i], beta[j], c.scale, c.gap_scale);
          }
          if (c.kind == SpinonProposalKind::one_pair)
          {
            for (std::size_t a = 0; a < pairs.size(); ++a)
              for (int k = 0; k < 4; ++k)
                if (k != pairs_[a][0] && k != pairs_[a][1]) pairs[a] += native[k];
            value = log_sum(pairs) - uni20::math::log(R{6} * one_);
          }
          else
          {
            std::array<R, 3> matchings{pairs[0] + pairs[5], pairs[1] + pairs[4], pairs[2] + pairs[3]};
            value = log_sum(matchings) - uni20::math::log(R{3} * one_);
          }
        }
        terms.push_back(uni20::math::log(c.weight / total_weight_) + value);
      }
      return log_sum(terms);
    }

    Sample draw(std::array<R, 5> const& u) const
    {
      for (auto const& x : u)
        check_uniform(x);
      R const selector = u[0] * total_weight_;
      R cumulative = zero_, previous = zero_;
      std::size_t selected = 0;
      for (; selected < components_.size(); ++selected)
      {
        previous = cumulative;
        cumulative += components_[selected].weight;
        if (selector < cumulative) break;
      }
      if (selected == components_.size()) throw std::runtime_error("proposal selector exhausted precision");
      auto const& c = components_[selected];
      std::array<R, 4> beta;
      std::size_t pairing = 0;
      if (c.kind == SpinonProposalKind::background || c.kind == SpinonProposalKind::cauchy)
        for (int i = 0; i < 4; ++i)
          beta[i] = c.kind == SpinonProposalKind::background ? native_inverse(u[i + 1])
                                                             : cauchy_draw(u[i + 1], c.scale, cutoff_);
      else
      {
        std::size_t const count = c.kind == SpinonProposalKind::one_pair ? 6 : 3;
        R const position = R(count) * (selector - previous) / c.weight;
        for (; pairing + 1 < count && position >= R(pairing + 1); ++pairing)
        {}
        if (!(zero_ <= position && position < R(count))) throw std::runtime_error("pair selector exhausted precision");
        auto const first = pairs_[pairing];
        auto const pair = pair_draw(u[1], u[2], c.scale, c.gap_scale);
        beta[first[0]] = pair[0];
        beta[first[1]] = pair[1];
        std::array<int, 2> rest;
        int j = 0;
        for (int i = 0; i < 4; ++i)
          if (i != first[0] && i != first[1]) rest[j++] = i;
        if (c.kind == SpinonProposalKind::one_pair)
        {
          beta[rest[0]] = native_inverse(u[3]);
          beta[rest[1]] = native_inverse(u[4]);
        }
        else
        {
          auto const second = pair_draw(u[3], u[4], c.scale, c.gap_scale);
          beta[rest[0]] = second[0];
          beta[rest[1]] = second[1];
        }
      }
      return {beta, log_density(beta), selected, pairing};
    }

    /// Pushforward intensities for all 15 solved pairs. Each chart emits BOTH
    /// labelled roots. These are intensities, not normalized probabilities.
    std::array<R, 15> chart_log_intensities(std::array<R, 6> const& beta) const
    {
      std::array<R, 6> cosh_logs;
      for (int i = 0; i < 6; ++i)
      {
        check_beta(beta[i]);
        cosh_logs[i] = log_cosh(beta[i]);
      }
      std::array<R, 15> result;
      int a = 0;
      for (int i = 0; i < 6; ++i)
        for (int j = i + 1; j < 6; ++j)
        {
          if (beta[i] == beta[j])
          {
            // Only a symbolic infinity is converted here, never finite data.
            result[a++] = uni20::scalar_like(one_, -std::numeric_limits<double>::infinity());
            continue;
          }
          std::array<R, 4> free;
          R jacobian = zero_;
          int k = 0;
          for (int m = 0; m < 6; ++m)
            if (m != i && m != j)
            {
              free[k++] = beta[m];
              jacobian += cosh_logs[m];
            }
          R const d = uni20::math::abs(beta[i] - beta[j]);
          R const log_sinh = d + uni20::math::log(-uni20::math::expm1(-R{2} * d)) - log_two_;
          result[a++] =
              log_density(free) + jacobian + uni20::math::log(pi_ / R{2}) + log_sinh - cosh_logs[i] - cosh_logs[j];
        }
      return result;
    }

    /// log[1/mean(h_ij)] per labelled root. Includes proposal and coarea
    /// Jacobians, no physical density or extra state-counting factorial.
    R balanced_log_jacobian(std::array<R, 6> const& beta) const
    {
      auto const intensities = chart_log_intensities(beta);
      return uni20::math::log(R{15} * one_) - log_sum(intensities);
    }

  private:
    void check_scalar(R const& value) const
    {
      if (!uni20::isfinite(value)) throw std::invalid_argument("nonfinite proposal coordinate/control");
      (void)(value * one_); // runtime precision compatibility, no ambient default
    }
    void check_beta(R const& beta) const
    {
      check_scalar(beta);
      if (!(uni20::math::abs(beta) < cutoff_)) throw std::invalid_argument("rapidity outside open cutoff box");
    }
    void check_uniform(R const& u) const
    {
      check_scalar(u);
      if (!(zero_ < u && u < one_)) throw std::invalid_argument("require open uniform coordinate");
    }
    R log_cosh(R const& x) const
    {
      R const a = uni20::math::abs(x);
      return a + uni20::math::log1p(uni20::math::exp(-R{2} * a)) - log_two_;
    }
    R log_sum(std::span<R const> values) const
    {
      R const largest = *std::max_element(values.begin(), values.end());
      if (!uni20::isfinite(largest)) throw std::runtime_error("degenerate or unresolved proposal intensity");
      R sum = zero_;
      for (auto const& x : values)
        sum += uni20::math::exp(x - largest);
      R const result = largest + uni20::math::log(sum);
      if (!uni20::isfinite(result)) throw std::runtime_error("proposal density exhausted exponent range");
      return result;
    }
    R cauchy_log_density(R const& x, R const& scale, R const& limit) const
    {
      R const a = uni20::math::abs(x), largest = std::max(a, scale), smallest = std::min(a, scale);
      R const ratio = smallest / largest;
      return uni20::math::log(scale) - R{2} * uni20::math::log(largest) - uni20::math::log1p(ratio * ratio) -
             uni20::math::log(R{2} * uni20::math::atan(limit / scale));
    }
    R pair_log_density(R const& x, R const& y, R const& scale, R const& gap) const
    {
      R const center = (x + y) / R{2}, width = R{2} * (cutoff_ - uni20::math::abs(center));
      if (!(width > zero_)) throw std::runtime_error("soft pair lost cutoff resolution");
      return cauchy_log_density(center, scale, cutoff_) + cauchy_log_density(x - y, gap, width);
    }
    R cauchy_draw(R const& u, R const& scale, R const& limit) const
    {
      R const result = scale * uni20::math::tan((R{2} * u - one_) * uni20::math::atan(limit / scale));
      if (!uni20::isfinite(result) || !(uni20::math::abs(result) < limit))
        throw std::runtime_error("Cauchy draw lost interior support");
      return result;
    }
    std::array<R, 2> pair_draw(R const& u, R const& v, R const& scale, R const& gap) const
    {
      R const center = cauchy_draw(u, scale, cutoff_);
      R const width = R{2} * (cutoff_ - uni20::math::abs(center));
      R const delta = cauchy_draw(v, gap, width);
      return {center + delta / R{2}, center - delta / R{2}};
    }

    static constexpr std::array<std::array<int, 2>, 6> pairs_{{{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}};
    R cutoff_, soft_;
    std::vector<Component> components_;
    R one_, zero_, pi_, eps_, log_two_, lower_, width_, soft_norm_, total_weight_;
};
} // namespace bethe::heisenberg::detail
