// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <algorithm>
#include <array>
#include <bethe/solver.hpp>
#include <cmath>
#include <map>
#include <optional>
#include <stdexcept>
#include <vector>

namespace bethe::detail
{
enum class ExtremaStatus
{
  converged,
  objective_failure,
  evaluation_limit,
  iteration_limit,
  mesh_limit,
  precision_limit
};
template <uni20::Real Real> struct ObjectiveSample
{
    Real value{}, error{};
};
template <uni20::Real Real> struct Extremum
{
    Real position{}, value{}, error{};
};
template <uni20::Real Real> struct ExtremaOptions
{
    Real tolerance = Real{65536} * uni20::numeric_limits<Real>::epsilon();
    Real position_tolerance = std::sqrt(uni20::numeric_limits<Real>::epsilon());
    std::size_t initial_intervals = 16, max_intervals = 128;
    std::size_t max_evaluations = 20000, max_iterations = 256;
};
template <uni20::Real Real> struct ExtremaResult
{
    std::optional<Extremum<Real>> minimum, maximum;
    std::size_t evaluations = 0, iterations = 0, intervals = 0, meshes = 0;
    bool converged = false;
    ExtremaStatus status = ExtremaStatus::mesh_limit;
};

// Heuristic global search on a closed finite interval, NOT a certificate.
// Enumerate extrema on nested meshes, refine candidates and both endpoint
// neighborhoods, and require two successive agreements of both edge values.
// Function returns optional<ObjectiveSample<Real>>; any failed sample aborts.
// No output extremum is published unless the entire search succeeds.
template <uni20::Real Real, typename Function>
ExtremaResult<Real> bounded_extrema(Function const& function, Real lo, Real hi, ExtremaOptions<Real> const& o = {})
{
  if (!uni20::isfinite(lo) || !uni20::isfinite(hi) || !(lo < hi) || !uni20::isfinite(hi - lo) ||
      !uni20::isfinite(o.tolerance) || !(o.tolerance > Real{0}) || !uni20::isfinite(o.position_tolerance) ||
      !(o.position_tolerance > Real{0}))
    throw std::invalid_argument("extrema require finite lo<hi and positive finite tolerances");
  if (o.initial_intervals < 4 || o.initial_intervals > o.max_intervals || o.max_intervals > 4096)
    throw std::invalid_argument("extrema require 4<=initial_intervals<=max_intervals<=4096");
  ExtremaResult<Real> out;
  std::map<Real, ObjectiveSample<Real>> cache;
  auto sample = [&](Real x) -> std::optional<Extremum<Real>> {
    auto it = cache.find(x);
    if (it == cache.end())
    {
      if (out.evaluations == o.max_evaluations)
      {
        out.status = ExtremaStatus::evaluation_limit;
        return {};
      }
      ++out.evaluations;
      auto const value = function(x);
      if (!value || !uni20::isfinite(value->value) || !uni20::isfinite(value->error) || value->error < Real{0})
      {
        out.status = ExtremaStatus::objective_failure;
        return {};
      }
      it = cache.emplace(x, *value).first;
    }
    return Extremum<Real>{x, it->second.value, it->second.error};
  };
  auto better = [](Extremum<Real> const& a, Extremum<Real> const& b, bool maximum) {
    return maximum ? a.value > b.value : a.value < b.value;
  };
  Real const fraction = (std::sqrt(Real{5}) - Real{1}) / Real{2};
  Real const xtol = o.position_tolerance * std::max({Real{1}, std::abs(lo), std::abs(hi)});
  if (!uni20::isfinite(xtol)) throw std::invalid_argument("extrema position tolerance overflows");
  auto local = [&](Real left, Real right, bool maximum) -> std::optional<Extremum<Real>> {
    auto l = sample(left), r = sample(right);
    if (!l || !r) return {};
    auto best = better(*l, *r, maximum) ? *l : *r;
    auto a = sample(right - fraction * (right - left)), b = sample(left + fraction * (right - left));
    if (!a || !b) return {};
    for (std::size_t iteration = 0;; ++iteration)
    {
      for (auto const& point : {*a, *b})
        if (better(point, best, maximum)) best = point;
      Real const uncertainty = std::max({best.error, a->error, b->error});
      if (uncertainty > o.tolerance / Real{16})
      {
        out.status = ExtremaStatus::precision_limit;
        return {};
      }
      Real const spread =
          std::max(std::abs(a->value - best.value) + a->error, std::abs(b->value - best.value) + b->error);
      if (right - left <= xtol && spread <= o.tolerance / Real{8})
      {
        // Value spread at the resolved local bracket, plus objective uncertainty.
        best.error += spread;
        return best;
      }
      if (iteration == o.max_iterations)
      {
        out.status = ExtremaStatus::iteration_limit;
        return {};
      }
      ++out.iterations;
      if (better(*a, *b, maximum))
      {
        right = b->position;
        b = a;
        Real const x = right - fraction * (right - left);
        if (!(left < x && x < b->position))
        {
          out.status = ExtremaStatus::precision_limit;
          return {};
        }
        a = sample(x);
      }
      else
      {
        left = a->position;
        a = b;
        Real const x = left + fraction * (right - left);
        if (!(a->position < x && x < right))
        {
          out.status = ExtremaStatus::precision_limit;
          return {};
        }
        b = sample(x);
      }
      if (!a || !b) return {};
    }
  };
  std::array<Extremum<Real>, 2> previous{};
  unsigned stable = 0;
  for (std::size_t n = o.initial_intervals;;)
  {
    out.intervals = n;
    ++out.meshes;
    std::vector<Extremum<Real>> grid;
    for (std::size_t i = 0; i <= n; ++i)
    {
      Real const x = i == n ? hi : lo + (hi - lo) * (Real(i) / Real(n));
      auto const point = sample(x);
      if (!point) return out;
      grid.push_back(*point);
    }
    std::array<Extremum<Real>, 2> current{grid.front(), grid.front()};
    for (std::size_t edge = 0; edge < 2; ++edge)
    {
      bool const maximum = edge == 1;
      for (auto const& point : grid)
        if (better(point, current[edge], maximum)) current[edge] = point;
      std::vector<std::pair<Real, Real>> brackets{{grid[0].position, grid[1].position},
                                                  {grid[n - 1].position, grid[n].position}};
      for (std::size_t i = 1; i < n; ++i)
        if (!better(grid[i - 1], grid[i], maximum) && !better(grid[i + 1], grid[i], maximum) &&
            (better(grid[i], grid[i - 1], maximum) || better(grid[i], grid[i + 1], maximum)))
          brackets.emplace_back(grid[i - 1].position, grid[i + 1].position);
      for (auto const& [a, b] : brackets)
      {
        auto const candidate = local(a, b, maximum);
        if (!candidate) return out;
        // Retain the refined uncertainty even when the exact mesh point wins.
        if (better(*candidate, current[edge], maximum) || candidate->value == current[edge].value)
          current[edge] = *candidate;
      }
    }
    bool good = out.meshes > 1;
    for (std::size_t edge = 0; edge < 2; ++edge)
    {
      current[edge].error += std::abs(current[edge].value - previous[edge].value);
      good = good && uni20::isfinite(current[edge].error) && current[edge].error <= o.tolerance;
    }
    stable = good ? stable + 1 : 0;
    previous = current;
    if (stable == 2)
    {
      out.minimum = current[0];
      out.maximum = current[1];
      out.converged = true;
      out.status = ExtremaStatus::converged;
      return out;
    }
    if (n > o.max_intervals / 2) return out;
    n *= 2;
  }
}
} // namespace bethe::detail
