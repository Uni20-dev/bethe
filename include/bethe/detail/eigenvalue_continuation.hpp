// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/newton.hpp>
#include <bethe/detail/newton_backtracking.hpp>
#include <optional>

namespace bethe::detail
{
template <uni20::Real Real> struct ContinuationPrediction
{
    Real trust, max_slope;
};

/// Correct one eigenvalue-variable predictor, staying in its expanded-coordinate
/// trust box. Count attempted least-squares corrections, including failed solves.
/// Deliberately retain the 12-attempt stage limit and pre-correction convergence
/// checks of the Richardson/Gaudin solvers (no extra check after attempt 12).
template <uni20::Real Real, typename Options, typename Evaluate, typename Expand>
bool correct_eigenvalue_stage(std::vector<Real>& candidate, Real trust, Options const& options, std::size_t& iterations,
                              Evaluate&& evaluate, Expand&& expand)
{
  auto const predictor = expand(candidate);
  for (unsigned update = 0; update < 12; ++update)
  {
    std::vector<Real> jac;
    auto const evaluation = evaluate(candidate, &jac);
    if (evaluation.norm <= options.residual_tolerance) return true;
    if (iterations == options.max_iterations) return false;
    ++iterations;
    auto correction = evaluation.residual;
    for (auto& v : correction)
      v = -v;
    if (!least_squares_step(std::move(jac), correction, candidate.size())) return false;
    bool const improved = backtrack_newton(
        candidate, [&](std::size_t i) { return correction[i]; },
        [&](auto const& trial, Real damping) {
          auto const values = expand(trial);
          for (std::size_t i = 0; i < values.size(); ++i)
            if (!uni20::isfinite(values[i]) || std::abs(values[i] - predictor[i]) > trust) return false;
          return newton_decreases(evaluate(trial, nullptr).norm, evaluation.norm, damping, options.residual_tolerance);
        });
    if (!improved) return false;
  }
  return false;
}

struct NoContinuationCommit
{
    template <typename Real> void operator()(Real) const {}
};

/// Common adaptive controller; the model supplies tangent/trust information and
/// correction plus physical acceptance. Rejected stages never replace x/reached
/// or call commit. Budgets are shared across all stages, not reset on retries.
/// Model work records provide iterations, stages, rejected_stages and status.
template <uni20::Real Real, typename Options, typename Work, typename Predict, typename Correct,
          typename Commit = NoContinuationCommit>
void continue_eigenvalues(std::vector<Real>& x, Real& reached, Real step, Real target, Options const& options,
                          Work& work, Predict&& predict, Correct&& correct, Commit commit = {})
{
  using Status = decltype(work.status);
  while (reached < target && work.iterations < options.max_iterations && work.stages < options.max_stages)
  {
    std::vector<Real> slope;
    auto const prediction = predict(slope);
    if (!prediction)
    {
      work.status = Status::ill_conditioned;
      break;
    }
    if (prediction->max_slope > Real{0}) step = std::min(step, prediction->trust / prediction->max_slope);
    step = std::min(step, target - reached);
    Real const next = step == target - reached ? target : reached + step;
    if (!(next > reached))
    {
      work.status = Status::stalled;
      break;
    }
    ++work.stages;
    auto candidate = x;
    for (std::size_t i = 0; i < x.size(); ++i)
      candidate[i] += step * slope[i];
    if (correct(candidate, next, step, slope, prediction->trust))
    {
      x = std::move(candidate);
      reached = next;
      commit(next);
      step *= Real{1.5};
    }
    else
    {
      ++work.rejected_stages;
      step /= Real{2};
    }
  }
}
} // namespace bethe::detail
