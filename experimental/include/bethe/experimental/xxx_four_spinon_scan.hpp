// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/detail/continuous_work.hpp>
#include <bethe/dynamical_structure_factor.hpp>
#include <bethe/experimental/xxx_four_spinon.hpp>
#include <bethe/real_excitations.hpp>
#include <cstdint>
#include <tuple>

namespace bethe::experimental::xxx
{
namespace four_detail
{
inline int modulo(std::int64_t value, int sites)
{
  value %= sites;
  return static_cast<int>(value < 0 ? value + sites : value);
}

inline bool request_less(FourSpinonRequest const& a, FourSpinonRequest const& b)
{
  return std::tie(a.sites, a.holes, a.branch) < std::tie(b.sites, b.holes, b.branch);
}
} // namespace four_detail

/// Stable label identity, independent of scan selection, energy or completion order.
inline std::string four_spinon_id(FourSpinonRequest request)
{
  request = four_detail::canonical_request(request);
  std::string id = "N" + std::to_string(request.sites) + "-s4-h";
  for (std::size_t i = 0; i < request.holes.size(); ++i)
    id += (i ? "," : "") + std::to_string(request.holes[i]);
  return id + "-b" + std::to_string(request.branch);
}

/// Fourier transfer q=P0-Pn, as an index in [0,N). N divisible by four has P0=0.
inline int four_spinon_transfer_index(FourSpinonRequest request)
{
  request = four_detail::canonical_request(request);
  std::int64_t transfer = request.branch;
  for (int hole : request.holes)
    transfer -= hole;
  return four_detail::modulo(transfer, request.sites);
}

/// Enumerate quantum-number requests, not guaranteed resolved states. No parity
/// folding: every request has multiplicity one. A fixed-q catalogue solves the
/// last-hole constraint in O(N^3), rather than enumerating all momenta first.
/// The budget bounds storage and rejects before ground/state solving.
inline std::vector<FourSpinonRequest> four_spinon_catalogue(int sites, std::optional<int> transfer_index = {},
                                                            std::size_t max_candidates = 10000)
{
  (void)four_detail::canonical_request({sites, {-2, -1, 0, 1}, 0});
  if (!max_candidates) throw std::invalid_argument("max_candidates must be positive");
  if (transfer_index && (*transfer_index < 0 || *transfer_index >= sites))
    throw std::invalid_argument("momentum-transfer index must lie in [0,N)");
  std::vector<FourSpinonRequest> result;
  if (!transfer_index)
    result.reserve(3 * bethe::detail::bounded_binomial(static_cast<std::size_t>(sites / 2) + 1, 4, max_candidates / 3));
  auto append = [&](int a, int b, int c, int d, int branch) {
    if (result.size() == max_candidates) throw std::length_error("four-spinon catalogue exceeds max_candidates");
    result.push_back({sites, {a, b, c, d}, branch});
  };
  int const h = sites / 4;
  for (int a = -h; a < h - 2; ++a)
    for (int b = a + 1; b < h - 1; ++b)
      for (int c = b + 1; c < h; ++c)
        if (transfer_index)
        {
          for (int branch : {-1, 0, 1})
          {
            int const d = four_detail::modulo(std::int64_t(branch) - *transfer_index - a - b - c + h, sites) - h;
            if (d > c && d <= h) append(a, b, c, d, branch);
          }
        }
        else
          for (int d = c + 1; d <= h; ++d)
            for (int branch : {-1, 0, 1})
              append(a, b, c, d, branch);
  std::sort(result.begin(), result.end(), four_detail::request_less);
  return result;
}

template <uni20::Real R> struct FourSpinonRecord
{
    std::string candidate_id;
    FourSpinonSolve<R> solved;
    /// Absent when the solve failed or no factor was attempted.
    std::optional<heisenberg::FormFactor<R>> factor;
    /// A rejected state has no line, never a zero-weight placeholder.
    std::optional<LehmannLine<R>> line;
    std::string diagnostic;
};

template <uni20::Real R> struct FourSpinonScan
{
    int sites;
    R ground_energy;
    /// Canonical label order, including all requested failed/unsupported states.
    std::vector<FourSpinonRecord<R>> records;
    std::vector<LehmannLine<R>> lines;
    /// N momentum bins of accepted weights; omitted lines are not estimated.
    std::vector<SpectralMoment<R>> moments;
    bool converged() const { return lines.size() == records.size(); }
};

namespace four_detail
{
template <uni20::Real R>
FourSpinonRecord<R> evaluate(FourSpinonRequest const& request, PreparedSinglet<R> const& ground, std::stop_token stop)
{
  R const zero = *ground.working_zero();
  FourSpinonRecord<R> row;
  row.solved = solve_four_spinon<R>(request, zero);
  row.diagnostic = row.solved.diagnostic;
  if (!row.solved.converged() || stop.stop_requested()) return row;
  auto const& state = *row.solved.state;
  int const q = four_spinon_transfer_index(request);
  if (modulo(-std::int64_t(state.momentum_index()), request.sites) != q)
  {
    row.diagnostic = "state momentum differs from request";
    return row;
  }
  R const gap = state.energy() - *ground.energy();
  if (!uni20::isfinite(gap) || gap <= zero)
  {
    row.diagnostic = "nonpositive or nonfinite excitation gap";
    return row;
  }
  row.factor = four_spinon_raising_weight(ground, state);
  if (!row.factor->converged())
  {
    row.diagnostic = "four-spinon form factor rejected at the working precision";
    return row;
  }
  R const momentum = 2 * bethe::detail::pi_like(zero) * q / request.sites;
  row.line = LehmannLine<R>{0, static_cast<std::size_t>(q), momentum, gap, *row.factor->weight};
  return row;
}

template <uni20::Real R, typename Evaluate>
uni20::async::AsyncTask scan_task(uni20::async::IAsyncScheduler& scheduler, PreparedSinglet<R> ground,
                                  std::vector<FourSpinonRequest> requests, ExecutionOptions execution,
                                  Evaluate evaluate, uni20::async::WriteBuffer<FourSpinonScan<R>> output)
{
  auto destination = co_await output;
  if (execution.stop.stop_requested()) throw ExecutionCancelled("finite four-spinon scan cancelled");
  FourSpinonScan<R> result{static_cast<int>(ground.sites()), *ground.energy(), {}, {}, {}};
  result.records.resize(requests.size());
  auto const completed = bethe::detail::continuous_work(
      scheduler, execution.max_in_flight, requests.size(),
      [&](std::size_t i, std::stop_token stop) {
        result.records[i] = evaluate(requests[i], std::as_const(ground), stop);
        result.records[i].candidate_id = four_spinon_id(requests[i]);
      },
      execution.stop);
  if (completed != requests.size() || execution.stop.stop_requested())
    throw ExecutionCancelled("finite four-spinon scan cancelled");
  for (std::size_t i = 0; i < result.records.size(); ++i)
    if (auto& line = result.records[i].line)
    {
      line->state_id = i + 1;
      result.lines.push_back(*line);
    }
  result.moments = spectral_moments<R>(ground.sites(), result.lines, *ground.working_zero());
  destination = std::move(result);
}

/// Owning execution seam, also used by controlled failure/lifetime tests.
template <uni20::Real R, typename Evaluate>
uni20::async::Async<FourSpinonScan<R>> submit_scan(PreparedSinglet<R> ground, std::vector<FourSpinonRequest> requests,
                                                   ExecutionOptions execution, std::size_t max_candidates,
                                                   Evaluate evaluate)
{
  auto& scheduler = execution.continuous_scheduler();
  if (!ground.ready()) throw std::invalid_argument("four-spinon scan requires a ready prepared singlet");
  if (ground.sites() < 8 || ground.sites() % 4 || ground.sites() > std::numeric_limits<int>::max())
    throw std::invalid_argument("four-spinon scan requires N>=8 divisible by four");
  if (!max_candidates) throw std::invalid_argument("max_candidates must be positive");
  if (requests.size() > max_candidates) throw std::length_error("four-spinon scan exceeds max_candidates");
  for (auto& request : requests)
  {
    request = canonical_request(request);
    if (request.sites != static_cast<int>(ground.sites()))
      throw std::invalid_argument("four-spinon request and ground have different sizes");
  }
  std::sort(requests.begin(), requests.end(), request_less);
  if (std::adjacent_find(requests.begin(), requests.end()) != requests.end())
    throw std::invalid_argument("duplicate four-spinon request");
  uni20::async::Async<FourSpinonScan<R>> output;
  auto task =
      scan_task(scheduler, std::move(ground), std::move(requests), execution, std::move(evaluate), output.write());
  bethe::detail::schedule_owned(scheduler, std::move(task));
  return output;
}
} // namespace four_detail

/// Owns its immutable ground snapshot and canonicalized requests. The scheduler
/// is borrowed until completion. Active kernels finish before cancellation or
/// execution exceptions propagate; dropping the handle does not cancel work.
/// A completed operation may still contain explicit numerical failures.
template <uni20::Real R>
[[nodiscard]] uni20::async::Async<FourSpinonScan<R>>
scan_four_spinons_async(PreparedSinglet<R> ground, std::vector<FourSpinonRequest> requests,
                        ExecutionOptions execution = {}, std::size_t max_candidates = 10000)
{
  return four_detail::submit_scan(std::move(ground), std::move(requests), execution, max_candidates,
                                  four_detail::evaluate<R>);
}

template <uni20::Real R>
[[nodiscard]] FourSpinonScan<R> scan_four_spinons(PreparedSinglet<R> ground, std::vector<FourSpinonRequest> requests,
                                                  ExecutionOptions execution = {}, std::size_t max_candidates = 10000)
{
  auto& scheduler = execution.continuous_scheduler();
  execution.scheduler = &scheduler;
  auto result = scan_four_spinons_async(std::move(ground), std::move(requests), execution, max_candidates);
  return bethe::detail::take_result(result, scheduler);
}
} // namespace bethe::experimental::xxx
