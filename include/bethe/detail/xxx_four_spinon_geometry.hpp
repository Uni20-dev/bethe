// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <algorithm>
#include <array>
#include <bethe/detail/gauss_legendre.hpp>
#include <bethe/solver.hpp>
#include <optional>

namespace bethe::heisenberg::detail
{
/// Ordered four-momentum patch with total Q <= 2*pi. The two closest endpoint
/// distances are a=r*v, b=r*(1-v): either both at zero or at opposite ends.
/// The remaining ordered pair is eliminated by energy conservation for a
/// spectrum, or integrated directly for moments. No physical W4 coefficient.
template <uni20::Real R> class FourSpinonPatch {
  public:
    struct Bounds
    {
        R lower_energy, upper_energy, hard_sum, hard_lower;
    };
    struct Point
    {
        std::array<R, 4> beta;
        R energy, log_jacobian;
    };

    FourSpinonPatch(R total, R fraction, bool mixed)
        : total_(std::move(total)), fraction_(std::move(fraction)), mixed_(mixed),
          one_(bethe::detail::precision_zero(total_) + R{1}), zero_(one_ - one_), pi_(R{4} * uni20::math::atan(one_)),
          eps_(uni20::numeric_limits<R>::epsilon(one_))
    {
      check(total_);
      check(fraction_);
      if (!(zero_ < total_ && total_ < R{2} * pi_ && zero_ < fraction_ && fraction_ < (mixed_ ? one_ : one_ / R{2})) ||
          (mixed_ && total_ <= pi_))
        throw std::invalid_argument("invalid ordered four-spinon patch");
    }

    R radius() const
    {
      R const d = one_ - R{2} * fraction_, m = std::max(fraction_, one_ - fraction_);
      R result = mixed_ ? std::min((total_ - pi_) / (R{2} * m - d), (R{3} * pi_ - total_) / (R{2} * m + d))
                        : std::min(total_ / (R{3} - R{2} * fraction_), (R{2} * pi_ - total_) / d);
      if (!(result > zero_) || !uni20::isfinite(result)) throw std::runtime_error("unresolved four-spinon radius");
      return result;
    }

    R switch_radius() const
    {
      R const d = one_ - R{2} * fraction_;
      return mixed_ ? (d > zero_ ? R((R{2} * pi_ - total_) / d) : zero_) : R(total_ - pi_);
    }

    Bounds bounds(R const& r) const
    {
      check(r);
      if (r < zero_ || r > radius()) throw std::invalid_argument("radius outside ordered patch");
      R const a = r * fraction_, b = r * (one_ - fraction_);
      R const hard_sum = mixed_ ? total_ - pi_ + r * (one_ - R{2} * fraction_) : total_ - r;
      R const lo = mixed_ ? std::max(a, b) + std::max(zero_, R(hard_sum - pi_)) : std::max(b, R(total_ - pi_ - a));
      R const soft = pi_ / R{2} * (uni20::math::sin(a) + uni20::math::sin(b));
      R const upper = pi_ * uni20::math::sin(hard_sum / R{2});
      return {soft + upper * uni20::math::cos(hard_sum / R{2} - lo), soft + upper, hard_sum, lo};
    }

    /// Actual radial support, located before quadrature. The two energy bounds
    /// are monotone on this patch; empty support is not a failed kernel sample.
    std::optional<std::pair<R, R>> support(R const& omega) const
    {
      check(omega);
      R const maximum = radius();
      auto const first = bounds(zero_), last = bounds(maximum);
      R const margin =
          R{64} * eps_ * std::max(uni20::math::abs(first.lower_energy), uni20::math::abs(last.upper_energy));
      if (omega < first.lower_energy - margin || omega > last.upper_energy + margin) return {};
      if (omega <= first.lower_energy + margin || omega >= last.upper_energy - margin)
        throw std::runtime_error("patch support boundary requires more precision or a limiting calculation");
      auto inverse = [&](bool upper) {
        R left = zero_, right = maximum;
        auto const digits = uni20::numeric_limits<R>::digits(one_);
        for (decltype(digits + 1) i = 0; i < digits + 8; ++i)
        {
          R const mid = (left + right) / R{2};
          if (mid == left || mid == right) break;
          auto const b = bounds(mid);
          if ((upper ? b.upper_energy : b.lower_energy) < omega)
            left = mid;
          else
            right = mid;
        }
        return R((left + right) / R{2});
      };
      R const left = omega > first.upper_energy ? inverse(true) : zero_;
      R const right = omega < last.lower_energy ? inverse(false) : maximum;
      if (!(right > left)) throw std::runtime_error("radial support lost precision");
      return std::pair{left, right};
    }

    Point spectral_point(R const& log_radius, R const& omega) const
    {
      check(omega);
      R const r = radial_value(log_radius);
      auto const bd = bounds(r);
      R const soft = pi_ / R{2} * (uni20::math::sin(r * fraction_) + uni20::math::sin(r * (one_ - fraction_)));
      R const upper = pi_ * uni20::math::sin(bd.hard_sum / R{2}), hard = omega - soft;
      R const squared = (upper - hard) * (upper + hard);
      if (!(hard > zero_ && squared > zero_ && omega > bd.lower_energy && omega < bd.upper_energy))
        throw std::runtime_error("unresolved or excluded four-spinon energy root");
      R const jacobian = uni20::math::sqrt(squared);
      R const theta = uni20::math::atan2(jacobian, hard);
      R const first = bd.hard_sum / R{2} - theta, second = bd.hard_sum / R{2} + theta;
      return {rapidities(log_radius, first, second), omega, -uni20::math::log(jacobian)};
    }

    Point moment_point(R const& log_radius, R const& unit) const
    {
      check(unit);
      if (!(zero_ < unit && unit < one_)) throw std::invalid_argument("hard-pair coordinate must be open");
      R const r = radial_value(log_radius);
      auto const bd = bounds(r);
      R const width = bd.hard_sum / R{2} - bd.hard_lower;
      if (!(width > zero_)) throw std::runtime_error("hard-pair interval lost precision");
      R const first = bd.hard_lower + width * unit, second = bd.hard_sum - first;
      R const energy = pi_ / R{2} *
                       (uni20::math::sin(r * fraction_) + uni20::math::sin(r * (one_ - fraction_)) +
                        uni20::math::sin(first) + uni20::math::sin(second));
      return {rapidities(log_radius, first, second), energy, uni20::math::log(width)};
    }

  private:
    void check(R const& x) const
    {
      bethe::detail::require_same_precision(x, one_);
      if (!uni20::isfinite(x)) throw std::invalid_argument("nonfinite four-spinon coordinate");
    }
    R radial_value(R const& log_radius) const
    {
      check(log_radius);
      if (log_radius >= uni20::math::log(radius())) throw std::invalid_argument("log radius must be inside patch");
      // r may underflow; the soft rapidities still retain log_radius below.
      return uni20::math::exp(log_radius);
    }
    R from_log_distance(R const& log_distance) const
    {
      R const leading = log_distance - uni20::math::log(R{2} * one_);
      // log[tan(x/2)/(x/2)]=O(x^2). Retain log(x) when this correction
      // is smaller than epsilon, including x below the exponent range.
      if (R{2} * log_distance < uni20::math::log(eps_)) return leading;
      return uni20::math::log(uni20::math::tan(uni20::math::exp(log_distance) / R{2}));
    }
    std::array<R, 4> rapidities(R const& log_radius, R const& first, R const& second) const
    {
      if (!(zero_ < first && first < second && second < pi_))
        throw std::runtime_error("hard momenta lost interior ordering");
      auto hard_beta = [&](R const& p) {
        return p <= pi_ / R{2} ? from_log_distance(uni20::math::log(p)) : -from_log_distance(uni20::math::log(pi_ - p));
      };
      R const a = from_log_distance(log_radius + uni20::math::log(fraction_));
      R const b = from_log_distance(log_radius + uni20::math::log1p(-fraction_));
      return mixed_ ? std::array{a, hard_beta(first), hard_beta(second), -b}
                    : std::array{a, b, hard_beta(first), hard_beta(second)};
    }
    R total_, fraction_;
    bool mixed_;
    R one_, zero_, pi_, eps_;
};

/// Finite log-radius cutoff t: r/r_max >= exp(-t) on patches reaching r=0.
/// Positive radial supports are integrated in full, with a log/sine-squared
/// map. No omitted-tail estimate or convergence verdict is implied.
template <uni20::Real R> class FourSpinonGeometry {
  public:
    struct Node
    {
        std::array<R, 4> beta;
        R energy, log_measure; // includes ordered-domain 4!, but no W4 or 2*pi moment measure
    };

    FourSpinonGeometry(std::size_t order, R soft_cutoff)
        : cutoff_(std::move(soft_cutoff)), one_(bethe::detail::precision_zero(cutoff_) + R{1}), zero_(one_ - one_),
          pi_(R{4} * uni20::math::atan(one_)),
          rule_(bethe::detail::gauss_legendre<R>(checked_order(order), uni20::precision_of(one_)))
    {
      if (!uni20::isfinite(cutoff_) || cutoff_ <= zero_ || order < 2 || order > 256)
        throw std::invalid_argument("require quadrature order 2..256 and finite positive soft cutoff");
    }

    template <class F> std::size_t spectrum(R const& momentum, R const& omega, F&& sample) const
    {
      bethe::detail::require_same_precision(omega, one_);
      if (!uni20::isfinite(omega) || omega <= zero_) throw std::invalid_argument("require finite positive frequency");
      return patches(
          momentum,
          [&](auto const& patch, R const& log_v_weight, bool reflected) {
            auto const support = patch.support(omega);
            if (!support) return std::size_t{0};
            return radial(support->first, support->second, [&](R const& log_r, R const& log_r_weight) {
              auto const point = patch.spectral_point(log_r, omega);
              emit(point, log_v_weight + log_r_weight, reflected, sample);
              return std::size_t{1};
            });
          },
          omega);
    }

    template <class F> std::size_t moments(R const& momentum, F&& sample) const
    {
      return patches(momentum, [&](auto const& patch, R const& log_v_weight, bool reflected) {
        return radial(
            zero_, patch.radius(),
            [&](R const& log_r, R const& log_r_weight) {
              for (std::size_t i = 0; i < rule_.x.size(); ++i)
              {
                auto const [unit, weight] = sine_node(i);
                auto const point = patch.moment_point(log_r, unit);
                emit(point, log_v_weight + log_r_weight + uni20::math::log(weight), reflected, sample);
              }
              return rule_.x.size();
            },
            patch.switch_radius());
      });
    }

  private:
    static std::size_t checked_order(std::size_t order)
    {
      if (order < 2 || order > 256) throw std::invalid_argument("quadrature order must be in 2..256");
      return order;
    }
    std::pair<R, R> sine_node(std::size_t i) const
    {
      R const angle = pi_ * (one_ + rule_.x[i]) / R{2};
      R const sine = uni20::math::sin(angle / R{2});
      return {sine * sine, pi_ / R{4} * uni20::math::sin(angle) * rule_.w[i]};
    }
    R maximum_patch_energy(R const& total, R const& v, bool mixed) const
    {
      // Also defined at the fraction endpoints, where constructing a density
      // point would be invalid. These are limits of the geometric bounds.
      R const d = one_ - R{2} * v, m = std::max(v, R(one_ - v));
      R const r = mixed
                      ? std::min((total - pi_) / (R{2} * m - d), (R{3} * pi_ - total) / (R{2} * m + d))
                      : (d == zero_ ? R(total / R{2}) : std::min(total / (R{3} - R{2} * v), (R{2} * pi_ - total) / d));
      R const hard_sum = mixed ? total - pi_ + r * d : total - r;
      return pi_ / R{2} * (uni20::math::sin(r * v) + uni20::math::sin(r * (one_ - v))) +
             pi_ * uni20::math::sin(hard_sum / R{2});
    }
    std::optional<std::pair<R, R>> fraction_support(R const& total, bool mixed, R const& omega) const
    {
      R const end = mixed ? one_ : one_ / R{2};
      R peak = end;
      if (mixed)
      {
        // Below v=1/2 the stationary point is cos(A*v)=3*cos(A*(1-v)),
        // A=Q-pi. The alternate radius branch is increasing up to this
        // branch, and the whole v>1/2 branch is decreasing. Thus the maximum
        // is unique (possibly at v=0); its two sides are safe to bisect.
        R const a = total - pi_;
        peak =
            std::max(zero_, R(uni20::math::atan2(one_ - R{3} * uni20::math::cos(a), R{3} * uni20::math::sin(a)) / a));
      }
      R const maximum = maximum_patch_energy(total, peak, mixed);
      R const margin = R{64} * uni20::numeric_limits<R>::epsilon(one_) * maximum;
      if (omega > maximum + margin) return {};
      if (omega >= maximum - margin)
        throw std::runtime_error("fraction support edge requires more precision or a limiting calculation");
      auto inverse = [&](R left, R right, bool increasing) {
        auto const digits = uni20::numeric_limits<R>::digits(one_);
        for (decltype(digits + 1) i = 0; i < digits + 8; ++i)
        {
          R const mid = (left + right) / R{2};
          if (mid == left || mid == right) break;
          R const energy = maximum_patch_energy(total, mid, mixed);
          if (increasing ? energy < omega : energy > omega)
            left = mid;
          else
            right = mid;
        }
        return R((left + right) / R{2});
      };
      R const left = omega > maximum_patch_energy(total, zero_, mixed) ? inverse(zero_, peak, true) : zero_;
      R const right = omega > maximum_patch_energy(total, end, mixed) ? inverse(peak, end, false) : end;
      if (!(left < right)) throw std::runtime_error("fraction support interval lost precision");
      return std::pair{left, right};
    }
    template <class F> std::size_t patches(R momentum, F&& sample, std::optional<R> omega = {}) const
    {
      bethe::detail::require_same_precision(momentum, one_);
      if (!uni20::isfinite(momentum) || !(zero_ < momentum && momentum < R{2} * pi_))
        throw std::invalid_argument("four-spinon geometry requires 0<k<2*pi");
      bool const fold = momentum > pi_;
      if (fold) momentum = R{2} * pi_ - momentum;
      std::size_t count = 0;
      // The second physical section k+2*pi is reflected into Q=2*pi-k.
      // Reflect the actual beta values back before calling the density.
      for (int section = 0; section < 2; ++section)
      {
        R const total = section ? R{2} * pi_ - momentum : momentum;
        if (!(zero_ < total && total < R{2} * pi_))
          throw std::runtime_error("momentum section requires more precision");
        for (bool mixed : {false, true})
        {
          if (mixed && total <= pi_) continue;
          std::pair<R, R> interval{zero_, mixed ? one_ : one_ / R{2}};
          if (omega)
          {
            auto const allowed = fraction_support(total, mixed, *omega);
            if (!allowed) continue;
            interval = *allowed;
          }
          std::vector<R> cuts{zero_, one_ / R{2}};
          if (mixed) cuts.push_back(one_);
          if (total > R{3} * pi_ / R{2}) cuts.push_back((R{2} * total - R{3} * pi_) / (R{2} * (total - pi_)));
          std::erase_if(cuts, [&](R const& v) { return v <= interval.first || v >= interval.second; });
          cuts.push_back(interval.first);
          cuts.push_back(interval.second);
          std::sort(cuts.begin(), cuts.end());
          for (std::size_t panel = 1; panel < cuts.size(); ++panel)
            for (std::size_t i = 0; i < rule_.x.size(); ++i)
            {
              auto const [unit, weight] = sine_node(i);
              R const width = cuts[panel] - cuts[panel - 1];
              R const fraction = cuts[panel - 1] + width * unit;
              FourSpinonPatch<R> const patch(total, fraction, mixed);
              count += sample(patch, uni20::math::log(width * weight), bool(section) != fold);
            }
        }
      }
      return count;
    }
    template <class F>
    std::size_t radial(R const& left, R const& right, F&& sample, std::optional<R> switch_radius = {}) const
    {
      std::size_t count = 0;
      if (left > zero_)
      {
        R const lo = uni20::math::log(left), width = uni20::math::log(right) - lo;
        for (std::size_t i = 0; i < rule_.x.size(); ++i)
        {
          auto const [unit, weight] = sine_node(i);
          R const log_r = lo + width * unit;
          count += sample(log_r, R{2} * log_r + uni20::math::log(width * weight));
        }
      }
      else
      {
        std::vector<R> cuts{zero_, cutoff_};
        for (R t = R{2} * one_; t < cutoff_; t *= R{2})
          cuts.push_back(t);
        if (switch_radius && zero_ < *switch_radius && *switch_radius < right)
        {
          R const t = uni20::math::log(right) - uni20::math::log(*switch_radius);
          if (t < cutoff_) cuts.push_back(t);
        }
        std::sort(cuts.begin(), cuts.end());
        cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
        for (std::size_t panel = 1; panel < cuts.size(); ++panel)
          for (std::size_t i = 0; i < rule_.x.size(); ++i)
          {
            R const width = cuts[panel] - cuts[panel - 1];
            R const t = cuts[panel - 1] + width * (one_ + rule_.x[i]) / R{2};
            R const log_r = uni20::math::log(right) - t;
            count += sample(log_r, R{2} * log_r + uni20::math::log(width * rule_.w[i] / R{2}));
          }
      }
      return count;
    }
    template <class F>
    void emit(typename FourSpinonPatch<R>::Point const& point, R const& log_weight, bool reflected, F& sample) const
    {
      Node node{point.beta, point.energy, uni20::math::log(R{24} * one_) + log_weight + point.log_jacobian};
      if (reflected)
      {
        std::reverse(node.beta.begin(), node.beta.end());
        for (auto& b : node.beta)
          b = -b;
      }
      sample(node);
    }
    R cutoff_, one_, zero_, pi_;
    bethe::detail::GaussLegendreRule<R> rule_;
};
} // namespace bethe::heisenberg::detail
