// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <array>
#include <bethe/solver.hpp>
#include <limits>
#include <span>
#include <vector>

namespace bethe::heisenberg::detail
{
/// Solve the fixed-momentum/fixed-energy constraints in one labelled chart,
/// at J=1. The n-2 input rapidities are independent of the proposal law.
/// Numerical boundary ambiguity throws; only resolved support exclusions are
/// returned as outside. This is not an interval-certified geometry calculation.
/// Free rapidities retain their endpoint distances, without forming pi-distance.
/// A conditioning screen targets sqrt(epsilon) in the solved rapidities and log
/// Jacobian. Runtime precision allows two local retries at twice/four times the
/// input precision; results retain the input precision. Input rounding, density
/// accuracy and integration error are not covered by this arithmetic screen.
template <uni20::Real R> class SpinonEnergySurface {
  public:
    struct Point
    {
        bool inside;
        std::array<std::vector<R>, 2> beta; // both labelled solved-pair roots
        R log_coarea_jacobian;              // per root, no proposal or dp/d beta factors
        uni20::Precision::bit_count_type working_bits = 0;
        unsigned refinements = 0;
    };

    SpinonEnergySurface(int spinons, R cutoff, R momentum, int section, R omega)
        : n_(spinons), section_(section), cutoff_(std::move(cutoff)), momentum_(std::move(momentum)),
          omega_(std::move(omega)), zero_(uni20::scalar_like(cutoff_, 0)), one_(uni20::scalar_like(cutoff_, 1))
    {
      if constexpr (uni20::has_runtime_precision_v<R>)
        if (one_.precision().is_exact()) throw std::invalid_argument("surface requires explicit finite precision");
      check_scalar(cutoff_);
      check_scalar(momentum_);
      check_scalar(omega_);
      pi_ = R{4} * uni20::math::atan(one_);
      tolerance_ = R{64} * uni20::numeric_limits<R>::epsilon(one_) * pi_;
      accuracy_ = uni20::math::sqrt(uni20::numeric_limits<R>::epsilon(one_));
      if ((n_ != 2 && n_ != 4 && n_ != 6) || section < 0 || section >= n_ / 2 || cutoff_ <= zero_ ||
          momentum_ < zero_ || momentum_ >= R{2} * pi_ || omega_ <= zero_)
        throw std::invalid_argument("invalid spinon surface/section/cutoff");
      total_ = momentum_ + R(2 * section) * pi_;
      lower_ = R{2} * uni20::math::atan(uni20::math::exp(-cutoff_));
      if (lower_ == zero_) throw std::runtime_error("surface cutoff exceeds exponent range");
      if (!(zero_ < lower_ && lower_ < pi_ / R{2}))
        throw std::runtime_error("surface cutoff unresolved; increase geometry precision");
    }

    Point solve(std::span<R const> free_beta) const
    {
      if (free_beta.size() != std::size_t(n_ - 2)) throw std::invalid_argument("wrong number of free rapidities");
      for (auto const& b : free_beta)
      {
        check_scalar(b);
        if (!(uni20::math::abs(b) < cutoff_)) throw std::invalid_argument("free rapidity outside open cutoff box");
      }
      try
      {
        return solve_once(free_beta);
      }
      catch (Unresolved const&)
      {
        if constexpr (!uni20::has_runtime_precision_v<R>) throw;
      }
      if constexpr (uni20::has_runtime_precision_v<R>)
      {
        auto bits = one_.precision().bit_count();
        for (unsigned attempt = 1; attempt <= 2; ++attempt)
        {
          if (bits > std::numeric_limits<decltype(bits)>::max() / 2) break;
          bits *= 2;
          auto const precision = uni20::Precision::bits(bits);
          SpinonEnergySurface refined(n_, R(cutoff_, precision), R(momentum_, precision), section_,
                                      R(omega_, precision));
          refined.accuracy_ = R(accuracy_, precision);
          std::vector<R> free;
          for (auto const& b : free_beta)
            free.emplace_back(b, precision);
          try
          {
            auto point = refined.solve_once(free);
            if (point.inside)
            {
              // The integration weights recompute velocity differences from
              // the returned beta values. A good high-precision Jacobian alone
              // does not suffice if rounding those values destroys their gap.
              R const a = point.beta[0][n_ - 2], b = point.beta[0][n_ - 1];
              R const gap = uni20::math::abs(a - b);
              R const rounding = uni20::math::abs(R(R(a, one_.precision()), precision) - a) +
                                 uni20::math::abs(R(R(b, one_.precision()), precision) - b);
              // |d log sinh(gap)| <= (1+2/gap)*|d gap| when
              // |d gap| < gap/2; the two log cosh derivatives are bounded by 1.
              if (rounding >= gap / R{2} || rounding * (R{3} + R{2} / gap) > refined.accuracy_)
                throw std::runtime_error("surface pair gap unresolved at output precision");
            }
            for (auto& branch : point.beta)
            {
              for (auto& b : branch)
                b = R(b, one_.precision());
              if (point.inside && (branch[n_ - 2] == branch[n_ - 1] || uni20::math::abs(branch[n_ - 2]) >= cutoff_ ||
                                   uni20::math::abs(branch[n_ - 1]) >= cutoff_))
                throw std::runtime_error("surface result unresolved at output precision");
            }
            point.log_coarea_jacobian = R(point.log_coarea_jacobian, one_.precision());
            point.refinements = attempt;
            return point;
          }
          catch (Unresolved const&)
          {}
        }
      }
      throw std::runtime_error("surface geometry refinement exhausted; increase geometry precision");
    }

  private:
    class Unresolved : public std::runtime_error {
      public:
        using std::runtime_error::runtime_error;
    };

    Point solve_once(std::span<R const> free_beta) const
    {
      auto const bits = uni20::numeric_limits<R>::digits(one_);
      auto outside = [&] { return Point{false, {}, zero_, bits}; };
      if (outside_interval(total_, R(n_) * lower_, R(n_) * pi_ - R(n_) * lower_)) return outside();
      R const maximum = R(n_) * pi_ / R{2} * uni20::math::sin(total_ / R(n_));
      if (omega_ > maximum + tolerance_) return outside();
      if (uni20::math::abs(omega_ - maximum) <= tolerance_)
        throw std::runtime_error("upper surface edge unresolved; requires refinement or a limiting calculation");

      bethe::detail::CompensatedSum<R> distances(zero_), energies(zero_);
      int positive = 0;
      for (auto const& b : free_beta)
      {
        R const e = uni20::math::exp(-uni20::math::abs(b));
        R const distance = R{2} * uni20::math::atan(e);
        // The original beta already establishes open cutoff support. Keep the
        // endpoint's integer pi contribution separate from its tiny distance.
        if (b > zero_) ++positive;
        distances.add(b > zero_ ? -distance : distance);
        energies.add(pi_ * e / (one_ + e * e));
      }
      R const half = (momentum_ + R(2 * section_ - positive) * pi_ - distances.value()) / R{2};
      if (outside_interval(half, lower_, pi_ - lower_)) return outside();
      R const residual = omega_ - energies.value();
      bool const reflected = half > pi_ / R{2};
      R const center = reflected ? pi_ - half : half;
      R const amplitude = pi_ * uni20::math::sin(center);
      R const max_angle = center - lower_;
      R const minimum = amplitude * uni20::math::cos(max_angle);
      if (outside_interval(residual, minimum, amplitude)) return outside();

      R const discriminant = (amplitude - residual) * (amplitude + residual);
      if (!(discriminant > zero_)) throw Unresolved("unresolved surface discriminant");
      R const root = uni20::math::sqrt(discriminant);
      R const theta = uni20::math::atan2(root, residual);
      R const small = center - theta, large = center + theta;
      bool const large_reflected = large > pi_ / R{2};
      // pi-(center+theta) = (pi-2*center)+(center-theta).
      R const large_distance = large_reflected ? (pi_ - R{2} * center) + small : large;
      if (!(small > lower_ && large_distance > lower_)) throw Unresolved("solved momenta lost cutoff resolution");

      // Boundary tests alone do not detect a finite but inaccurate Jacobian.
      // These conservative forward-error indicators include conditioning of
      // the discriminant, atan2 and log(tan(p/2)); they are not interval bounds.
      R const log_error = tolerance_ * (one_ / (amplitude - residual) + one_ / (amplitude + residual));
      R const angle_error = tolerance_ * (one_ + (amplitude + uni20::math::abs(residual)) / root);
      R const beta_error = angle_error / std::min(uni20::math::sin(small), uni20::math::sin(large_distance));
      if (!uni20::isfinite(log_error) || !uni20::isfinite(beta_error) || log_error > accuracy_ ||
          beta_error > accuracy_)
        throw Unresolved("surface conditioning unresolved; increase geometry precision");

      R const small_beta = uni20::math::log(uni20::math::tan(small / R{2}));
      R const large_log = uni20::math::log(uni20::math::tan(large_distance / R{2}));
      R const large_beta = large_reflected ? -large_log : large_log;
      R const a = reflected ? -small_beta : large_beta, b = reflected ? -large_beta : small_beta;
      if (!(uni20::math::abs(a) < cutoff_ && uni20::math::abs(b) < cutoff_))
        throw Unresolved("solved rapidities lost cutoff resolution");
      std::vector<R> beta(free_beta.begin(), free_beta.end());
      beta.push_back(a);
      beta.push_back(b);
      auto swapped = beta;
      std::swap(swapped[n_ - 2], swapped[n_ - 1]);
      R const log_jacobian = -uni20::math::log(root);
      if (!uni20::isfinite(log_jacobian)) throw std::runtime_error("surface Jacobian exhausted exponent range");
      return {true, {std::move(beta), std::move(swapped)}, log_jacobian, bits};
    }

  private:
    void check_scalar(R const& value) const
    {
      if (!uni20::isfinite(value)) throw std::invalid_argument("nonfinite surface coordinate/control");
      (void)(value * one_);
    }
    bool outside_interval(R const& x, R const& lo, R const& hi) const
    {
      if (x < lo - tolerance_ || x > hi + tolerance_) return true;
      if (x <= lo + tolerance_ || x >= hi - tolerance_)
        throw Unresolved("surface boundary unresolved; increase geometry precision");
      return false;
    }
    int n_, section_;
    R cutoff_, momentum_, omega_, zero_, one_, pi_, tolerance_, accuracy_, total_, lower_;
};
} // namespace bethe::heisenberg::detail
