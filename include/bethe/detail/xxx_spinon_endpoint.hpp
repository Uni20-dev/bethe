// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <uni20/core/math.hpp>
#include <uni20/core/runtime_precision.hpp>

namespace bethe::heisenberg::detail
{
/// Positive endpoint distance and energy at J=1, with their common exp(-|beta|)
/// factor retained as a logarithm. Projection into ordinary arithmetic carries
/// an absolute omission bound; a projected zero is not a physical zero.
template <uni20::Real R> struct SpinonScaledEndpoint
{
    R log_scale, distance_factor, energy_factor;

    struct Projection
    {
        R distance, energy, distance_error, energy_error;
    };

    static SpinonScaledEndpoint from_depth(R const& depth, R const& pi, R const& log_floor)
    {
      R const one = uni20::scalar_like(pi, 1), zero = one - one;
      // Below the projection floor the corrections to both factors are O(t^2).
      // No exp is evaluated outside the ordinary arithmetic range.
      R const t = -depth < log_floor ? zero : uni20::math::exp(-depth);
      R const distance = t == zero ? R{2} * one : R{2} * uni20::math::atan(t) / t;
      return {-depth, distance, pi / (one + t * t)};
    }

    R log_distance() const { return log_scale + uni20::math::log(distance_factor); }

    Projection project(R const& floor, R const& log_floor) const
    {
      R const zero = floor - floor;
      if (log_scale < log_floor)
        // atan(t) <= t and pi*t/(1+t^2) <= pi*t. The factors in this
        // branch are exactly 2 and pi, so these bound the omitted quantities.
        return {zero, zero, distance_factor * floor, energy_factor * floor};
      R const t = uni20::math::exp(log_scale);
      return {distance_factor * t, energy_factor * t, zero, zero};
    }
};
} // namespace bethe::heisenberg::detail
