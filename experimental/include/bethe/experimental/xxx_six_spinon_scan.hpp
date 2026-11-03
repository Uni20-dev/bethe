// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/detail/continuous_work.hpp>
#include <bethe/dynamical_structure_factor.hpp>
#include <bethe/experimental/xxx_six_spinon_catalogue.hpp>

namespace bethe::experimental::xxx::six_spinon
{
template <uni20::Real R> struct ScanRecord
{
    std::string candidate_id;
    SolveResult<R> solved;
    std::optional<heisenberg::FormFactor<R>> factor;
    /// Failed/rejected candidates have no line; they are never zero-filled.
    std::optional<LehmannLine<R>> line;
    std::string diagnostic;
};

template <uni20::Real R> struct ScanResult
{
    int sites;
    R ground_energy;
    std::vector<ScanRecord<R>> records;
    std::vector<LehmannLine<R>> lines;
    std::vector<SpectralMoment<R>> moments;
    /// Success of these requests only, not completeness of any physical sector.
    bool converged() const { return lines.size() == records.size(); }
};

namespace detail
{
template <uni20::Real R>
ScanRecord<R> evaluate(Configuration const& request, PreparedSinglet<R> const& ground, std::stop_token stop)
{
  R const zero = *ground.working_zero();
  ScanRecord<R> row;
  row.solved = solve<R>(request, zero);
  row.diagnostic = row.solved.diagnostic;
  if (!row.solved.converged() || stop.stop_requested()) return row;
  auto const& state = *row.solved.state;
  int const q = transfer_index(request);
  if (modulo(-std::int64_t(state.momentum_index), request.sites) != q)
  {
    row.diagnostic = "state momentum differs from request";
    return row;
  }
  R const gap = state.energy - *ground.energy();
  if (!uni20::isfinite(gap) || gap <= zero)
  {
    row.diagnostic = "nonpositive or nonfinite excitation gap";
    return row;
  }
  row.factor = raising_weight(ground, row.solved);
  if (!row.factor->converged())
  {
    row.diagnostic = "six-spinon form factor rejected at the working precision";
    return row;
  }
  R const momentum = 2 * bethe::detail::pi_like(zero) * q / request.sites;
  row.line = LehmannLine<R>{0, static_cast<std::size_t>(q), momentum, gap, *row.factor->weight};
  return row;
}

template <uni20::Real R, typename Evaluate>
uni20::async::AsyncTask scan_task(uni20::async::IAsyncScheduler& scheduler, PreparedSinglet<R> ground,
                                  std::vector<Configuration> requests, ExecutionOptions execution, Evaluate evaluate,
                                  uni20::async::WriteBuffer<ScanResult<R>> output)
{
  auto destination = co_await output;
  if (execution.stop.stop_requested()) throw ExecutionCancelled("finite six-spinon scan cancelled");
  ScanResult<R> result{static_cast<int>(ground.sites()), *ground.energy(), {}, {}, {}};
  result.records.resize(requests.size());
  auto const completed = bethe::detail::continuous_work(
      scheduler, execution.max_in_flight, requests.size(),
      [&](std::size_t i, std::stop_token stop) {
        result.records[i] = evaluate(requests[i], std::as_const(ground), stop);
        result.records[i].candidate_id = candidate_id(requests[i]);
      },
      execution.stop);
  if (completed != requests.size() || execution.stop.stop_requested())
    throw ExecutionCancelled("finite six-spinon scan cancelled");
  for (std::size_t i = 0; i < result.records.size(); ++i)
    if (auto& line = result.records[i].line)
    {
      line->state_id = i + 1;
      result.lines.push_back(*line);
    }
  result.moments = spectral_moments<R>(ground.sites(), result.lines, *ground.working_zero());
  destination = std::move(result);
}

/// Owning execution seam for numerical evaluation and controlled lifetime tests.
template <uni20::Real R, typename Evaluate>
uni20::async::Async<ScanResult<R>> submit_scan(PreparedSinglet<R> ground, std::vector<Configuration> requests,
                                               ExecutionOptions execution, std::size_t max_candidates,
                                               Evaluate evaluate)
{
  auto& scheduler = execution.continuous_scheduler();
  if (!ground.ready()) throw std::invalid_argument("six-spinon scan requires a ready prepared singlet");
  if (ground.sites() < 12 || ground.sites() % 4 || ground.sites() > std::numeric_limits<int>::max())
    throw std::invalid_argument("six-spinon scan requires N>=12 divisible by four");
  if (!max_candidates) throw std::invalid_argument("max_candidates must be positive");
  if (requests.size() > max_candidates) throw std::length_error("six-spinon scan exceeds max_candidates");
  for (auto& request : requests)
  {
    if (request.sites != static_cast<int>(ground.sites()))
      throw std::invalid_argument("six-spinon request and ground have different sizes");
    request = canonical_request(std::move(request));
  }
  std::sort(requests.begin(), requests.end(), request_less);
  if (std::adjacent_find(requests.begin(), requests.end(), [](auto const& a, auto const& b) {
        return !request_less(a, b) && !request_less(b, a);
      }) != requests.end())
    throw std::invalid_argument("duplicate six-spinon request");
  uni20::async::Async<ScanResult<R>> output;
  auto task =
      scan_task(scheduler, std::move(ground), std::move(requests), execution, std::move(evaluate), output.write());
  bethe::detail::schedule_owned(scheduler, std::move(task));
  return output;
}
} // namespace detail

/// Owns the immutable prepared-ground snapshot and canonical requests. The scheduler
/// is borrowed until completion. Active kernels drain before cancellation/errors
/// propagate; dropping the handle does not cancel. Numerical failures remain records.
template <uni20::Real R>
[[nodiscard]] uni20::async::Async<ScanResult<R>>
scan_async(PreparedSinglet<R> ground, std::vector<Configuration> requests, ExecutionOptions execution = {},
           std::size_t max_candidates = 10000)
{
  return detail::submit_scan(std::move(ground), std::move(requests), execution, max_candidates, detail::evaluate<R>);
}

template <uni20::Real R>
[[nodiscard]] ScanResult<R> scan(PreparedSinglet<R> ground, std::vector<Configuration> requests,
                                 ExecutionOptions execution = {}, std::size_t max_candidates = 10000)
{
  auto& scheduler = execution.continuous_scheduler();
  execution.scheduler = &scheduler;
  auto pending = scan_async(std::move(ground), std::move(requests), execution, max_candidates);
  return bethe::detail::take_result(pending, scheduler);
}
} // namespace bethe::experimental::xxx::six_spinon
