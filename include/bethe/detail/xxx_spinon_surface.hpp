// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <array>
#include <bethe/solver.hpp>
#include <span>
#include <vector>

namespace bethe::heisenberg::detail
{
/// Solve the fixed-momentum/fixed-energy constraints in one labelled chart,
/// at J=1. The n-2 input rapidities are independent of the proposal law.
/// Numerical boundary ambiguity throws; only resolved support exclusions are
/// returned as outside. This is not an interval-certified geometry calculation.
template <uni20::Real R> class SpinonEnergySurface {
  public:
    struct Point
    {
        bool inside;
        std::array<std::vector<R>, 2> beta; // both labelled solved-pair roots
        R log_coarea_jacobian;              // per root, no proposal or dp/d beta factors
    };

    SpinonEnergySurface(int spinons, R cutoff, R momentum, int section, R omega)
        : n_(spinons), cutoff_(std::move(cutoff)), omega_(std::move(omega)), zero_(uni20::scalar_like(cutoff_, 0)),
          one_(uni20::scalar_like(cutoff_, 1))
    {
      if constexpr (uni20::has_runtime_precision_v<R>)
        if (one_.precision().is_exact()) throw std::invalid_argument("surface requires explicit finite precision");
      check_scalar(cutoff_);
      check_scalar(momentum);
      check_scalar(omega_);
      pi_ = R{4} * uni20::math::atan(one_);
      tolerance_ = R{64} * uni20::numeric_limits<R>::epsilon(one_) * pi_;
      if ((n_ != 2 && n_ != 4 && n_ != 6) || section < 0 || section >= n_ / 2 || cutoff_ <= zero_ || momentum < zero_ ||
          momentum >= R{2} * pi_ || omega_ <= zero_)
        throw std::invalid_argument("invalid spinon surface/section/cutoff");
      total_ = momentum + R(2 * section) * pi_;
      lower_ = R{2} * uni20::math::atan(uni20::math::exp(-cutoff_));
      upper_ = pi_ - lower_;
      if (!(zero_ < lower_ && lower_ < upper_ && upper_ < pi_))
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
      auto outside = [&] { return Point{false, {}, zero_}; };
      if (outside_interval(total_, R(n_) * lower_, R(n_) * upper_)) return outside();
      R const maximum = R(n_) * pi_ / R{2} * uni20::math::sin(total_ / R(n_));
      if (omega_ > maximum + tolerance_) return outside();
      if (uni20::math::abs(omega_ - maximum) <= tolerance_)
        throw std::runtime_error("upper surface edge unresolved; requires refinement or a limiting calculation");

      bethe::detail::CompensatedSum<R> momenta(zero_), energies(zero_);
      for (auto const& b : free_beta)
      {
        R const e = uni20::math::exp(-uni20::math::abs(b));
        R const distance = R{2} * uni20::math::atan(e);
        R const p = b <= zero_ ? distance : pi_ - distance;
        if (!(lower_ < p && p < upper_)) throw std::runtime_error("free momentum lost cutoff resolution");
        momenta.add(p);
        energies.add(pi_ * e / (one_ + e * e));
      }
      R const half = (total_ - momenta.value()) / R{2};
      if (outside_interval(half, lower_, upper_)) return outside();
      R const residual = omega_ - energies.value();
      R const amplitude = pi_ * uni20::math::sin(half);
      R const max_angle = std::min(half - lower_, upper_ - half);
      R const minimum = amplitude * uni20::math::cos(max_angle);
      if (outside_interval(residual, minimum, amplitude)) return outside();

      R const discriminant = (amplitude - residual) * (amplitude + residual);
      if (!(discriminant > zero_)) throw std::runtime_error("unresolved surface discriminant");
      R const root = uni20::math::sqrt(discriminant);
      R const theta = uni20::math::atan2(root, residual);
      R const first = half + theta, second = half - theta;
      if (!(lower_ < first && first < upper_ && lower_ < second && second < upper_))
        throw std::runtime_error("solved momenta lost cutoff resolution");
      R const a = rapidity(first), b = rapidity(second);
      if (!(uni20::math::abs(a) < cutoff_ && uni20::math::abs(b) < cutoff_))
        throw std::runtime_error("solved rapidities lost cutoff resolution");
      std::vector<R> beta(free_beta.begin(), free_beta.end());
      beta.push_back(a);
      beta.push_back(b);
      auto swapped = beta;
      std::swap(swapped[n_ - 2], swapped[n_ - 1]);
      R const log_jacobian = -uni20::math::log(root);
      if (!uni20::isfinite(log_jacobian)) throw std::runtime_error("surface Jacobian exhausted exponent range");
      return {true, {std::move(beta), std::move(swapped)}, log_jacobian};
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
        throw std::runtime_error("surface boundary unresolved; increase geometry precision");
      return false;
    }
    R rapidity(R const& p) const
    {
      // Reflect near pi so the tangent does not approach its pole.
      return p <= pi_ / R{2} ? uni20::math::log(uni20::math::tan(p / R{2}))
                             : -uni20::math::log(uni20::math::tan((pi_ - p) / R{2}));
    }

    int n_;
    R cutoff_, omega_, zero_, one_, pi_, tolerance_, total_, lower_, upper_;
};
} // namespace bethe::heisenberg::detail
