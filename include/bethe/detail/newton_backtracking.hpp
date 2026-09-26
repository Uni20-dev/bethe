// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <uni20/core/math.hpp>
#include <utility>
#include <vector>

namespace bethe::detail
{
struct UnchangedNewtonCoordinates
{
    template <typename Vector> void operator()(Vector&) const {}
};

template <uni20::Real Real> bool newton_decreases(Real norm, Real previous, Real damping, Real tolerance)
{
  return norm <= tolerance || norm < (Real{1} - damping / Real{10000}) * previous;
}

/// Try a full Newton correction followed by at most digits halvings. Only an
/// accepted trial replaces x; every trial starts from the original iterate.
/// The caller supplies correction components (vector or dense-matrix storage),
/// optional coordinate normalization, and domain/residual acceptance. Linear
/// solves, convergence decisions and iteration counters remain with the caller.
template <uni20::Real Real, typename Correction, typename Accept, typename Normalize = UnchangedNewtonCoordinates>
bool backtrack_newton(std::vector<Real>& x, Correction&& correction, Accept&& accept, Normalize normalize = {})
{
  auto trial = x;
  Real damping{1};
  for (int backtrack = 0; backtrack <= uni20::numeric_limits<Real>::digits; ++backtrack)
  {
    for (std::size_t i = 0; i < x.size(); ++i)
      trial[i] = x[i] + damping * correction(i);
    normalize(trial);
    if (accept(trial, damping))
    {
      x = std::move(trial);
      return true;
    }
    damping /= Real{2};
  }
  return false;
}
} // namespace bethe::detail
