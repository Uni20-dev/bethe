// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/experimental/xxx_boundary_correlations.hpp>
#include <bethe/experimental/xxx_boundary_form_factors.hpp>
#include <bethe/heisenberg_excitations.hpp>
#include <bethe/local_spectrum.hpp>

namespace bethe::experimental::xxx::open
{
template <uni20::Real R> struct BoundaryRecord
{
    std::string candidate_id;
    heisenberg::QuantumNumbers quantum_numbers;
    int reflection_parity = 1;
    /// Absent if the ground preparation failed before attempting excited roots.
    std::optional<model::RealState<R>> state;
    std::optional<BoundaryFormFactor<R>> factor;
    std::optional<LocalSpectralLine<R>> line;
    std::string diagnostic;
};

template <uni20::Real R> struct BoundarySpectrum
{
    explicit BoundarySpectrum(std::size_t n, R const& zero)
        : sites(n), moments{zero, zero}, end_to_end_moments{{zero, zero}, {zero, zero}}
    {}
    std::size_t sites;
    model::RealState<R> ground_state;
    BoundaryFormFactorStatus ground_status = BoundaryFormFactorStatus::precision_limit;
    std::optional<R> ground_residual;
    BoundaryTwoSiteDensity<R> density;
    int ground_reflection_parity = 1;
    /// Every requested label tuple, in lexicographic order, including failures.
    std::vector<BoundaryRecord<R>> records;
    /// Local longitudinal weights; no Fourier N or 2*pi factor.
    std::vector<LocalSpectralLine<R>> lines;
    LocalSpectralMoments<R> moments;
    std::vector<CrossSpectralLine<R>> end_to_end_lines;
    CrossSpectralMoments<R> end_to_end_moments;

    bool family_converged() const
    {
      return ground_status == BoundaryFormFactorStatus::converged && records.size() == lines.size();
    }
    /// Working-precision acceptance, not a certified truncation error bar.
    bool sum_rules_consistent() const
    {
      auto const zz = density.zz();
      if (!retained_moments() || !zz || !uni20::isfinite(moments.weight) || !uni20::isfinite(moments.first_moment)) return false;
      R const zero = uni20::scalar_like(moments.weight, 0), quarter = uni20::scalar_like(zero, 1) / 4;
      R const gate = 4096 * uni20::scalar_like(zero, sites) * uni20::numeric_limits<R>::epsilon(zero);
      return moments.weight >= zero && moments.first_moment >= zero && moments.weight <= quarter + gate &&
             moments.first_moment <= -*zz + gate;
    }
    bool converged() const { return family_converged() && density.converged() && sum_rules_consistent(); }
    /// Absent if spectral preparation failed; an unattempted sum is not zero.
    std::optional<LocalSpectralMoments<R>> retained_moments() const
    {
      if (ground_status != BoundaryFormFactorStatus::converged) return {};
      return moments;
    }
    std::optional<CrossSpectralMoments<R>> retained_end_to_end_moments() const
    {
      if (ground_status != BoundaryFormFactorStatus::converged) return {};
      return end_to_end_moments;
    }
    /// Uniform absolute truncation bound for exact retained weights/energies.
    /// Numerical errors are separate. An exceeded norm is absent, never clamped.
    std::optional<R> omitted_weight() const
    {
      if (ground_status != BoundaryFormFactorStatus::converged || !uni20::isfinite(moments.weight)) return {};
      R const zero = uni20::scalar_like(moments.weight, 0), quarter = uni20::scalar_like(zero, 1) / 4;
      if (moments.weight < zero || moments.weight > quarter) return {};
      return quarter - moments.weight;
    }
    /// Missing local first moment, from the full ground-root boundary bond.
    std::optional<R> omitted_first_moment() const
    {
      auto const zz = density.zz();
      if (!retained_moments() || !zz || !uni20::isfinite(moments.first_moment)) return {};
      R const missing = -*zz - moments.first_moment;
      if (!uni20::isfinite(missing) || missing < uni20::scalar_like(missing, 0)) return {};
      return missing;
    }
};

namespace boundary_detail
{
inline std::string candidate_id(std::size_t sites, std::span<uni20::half_int const> numbers)
{
  std::string id = "N" + std::to_string(sites) + "-obc-I";
  for (std::size_t i = 0; i < numbers.size(); ++i)
    id += (i ? "," : "") + std::to_string(numbers[i].twice() / 2);
  return id;
}

template <uni20::Real R>
BoundaryRecord<R> evaluate(std::size_t sites, heisenberg::QuantumNumbers numbers,
                           BoundaryFormFactors<R> const& prepared, SolverOptions<R> const& solver, std::stop_token stop)
{
  BoundaryRecord<R> row;
  row.candidate_id = candidate_id(sites, numbers);
  row.reflection_parity = real_root_reflection_parity(sites, numbers);
  row.quantum_numbers = std::move(numbers);
  row.state = model::solve_real<R>(sites, row.quantum_numbers, solver);
  if (!row.state->converged)
  {
    row.diagnostic = "excited roots did not converge";
    return row;
  }
  if (stop.stop_requested()) return row;
  row.factor = prepared(*row.state);
  if (!row.factor->converged())
  {
    row.diagnostic = "boundary form factor rejected at working precision";
    return row;
  }
  R const gap = row.state->energy - *prepared.ground_energy();
  R const zero = *prepared.working_zero();
  if (!uni20::isfinite(gap) || gap <= zero)
  {
    row.diagnostic = "nonpositive or nonfinite excitation gap";
    return row;
  }
  row.line = LocalSpectralLine<R>{0, gap, *row.factor->weight / 2};
  return row;
}

template <uni20::Real R> struct PrepareBoundary
{
    BoundaryFormFactors<R> operator()(std::size_t sites, model::RealState<R> const& ground) const
    {
      return BoundaryFormFactors<R>(sites, ground);
    }
};

template <uni20::Real R, typename Evaluate, typename Prepare>
uni20::async::AsyncTask spectrum_task(uni20::async::IAsyncScheduler& scheduler, std::size_t sites, std::size_t count,
                                      SolverOptions<R> solver, ExecutionOptions execution, Evaluate evaluate,
                                      Prepare prepare, uni20::async::WriteBuffer<BoundarySpectrum<R>> output)
{
  auto destination = co_await output;
  if (execution.stop.stop_requested()) throw ExecutionCancelled("boundary spectrum cancelled");
  R const zero = bethe::detail::solver_zero(solver);
  BoundarySpectrum<R> result(sites, zero);
  result.ground_state = model::ground_state<R>(sites, solver);
  if (execution.stop.stop_requested()) throw ExecutionCancelled("boundary spectrum cancelled");
  // Independent ground-root observable: spectral normalization must not gate it.
  result.density = boundary_two_site_density(sites, result.ground_state);
  if (execution.stop.stop_requested()) throw ExecutionCancelled("boundary spectrum cancelled");
  auto prepared = prepare(sites, result.ground_state);
  result.ground_status = prepared.status();
  result.ground_residual = result.density.ground_residual;
  result.ground_reflection_parity = real_root_reflection_parity(sites, result.ground_state.quantum_numbers);
  auto const m = sites / 2 - 1;
  heisenberg::QuantumNumbers labels(m);
  std::iota(labels.begin(), labels.end(), uni20::half_int{1});
  result.records.resize(count);
  for (std::size_t i = 0; i < count; ++i)
  {
    auto& row = result.records[i];
    row.quantum_numbers = labels;
    row.candidate_id = candidate_id(sites, labels);
    row.reflection_parity = real_root_reflection_parity(sites, labels);
    if (!prepared.ready()) row.diagnostic = "ground preparation failed; excited state not attempted";
    if (i + 1 < count) bethe::detail::advance_combination<uni20::half_int>(labels, uni20::half_int(sites - m));
  }
  if (prepared.ready())
  {
    auto const completed = bethe::detail::continuous_work(
        scheduler, execution.max_in_flight, count,
        [&](std::size_t i, std::stop_token stop) {
          result.records[i] =
              evaluate(sites, result.records[i].quantum_numbers, std::as_const(prepared), std::as_const(solver), stop);
        },
        execution.stop);
    if (completed != count || execution.stop.stop_requested()) throw ExecutionCancelled("boundary spectrum cancelled");
    for (std::size_t i = 0; i < count; ++i)
      if (auto& line = result.records[i].line)
      {
        line->state_id = i + 1;
        result.lines.push_back(*line);
        R const coefficient = result.ground_reflection_parity * result.records[i].reflection_parity * line->weight;
        result.end_to_end_lines.push_back({i + 1, line->gap, {coefficient, zero}});
      }
    result.moments = local_spectral_moments<R>(result.lines, zero);
    result.end_to_end_moments = cross_spectral_moments<R>(result.end_to_end_lines, zero);
  }
  if (execution.stop.stop_requested()) throw ExecutionCancelled("boundary spectrum cancelled");
  destination = std::move(result);
}

/// Owning numerical seam, also used by admitted-work lifetime tests.
template <uni20::Real R, typename Evaluate, typename Prepare = PrepareBoundary<R>>
uni20::async::Async<BoundarySpectrum<R>> submit_spectrum(std::size_t sites, SolverOptions<R> solver,
                                                         std::size_t max_candidates, ExecutionOptions execution,
                                                         Evaluate evaluate, Prepare prepare = {})
{
  (void)heisenberg::detail::checked_sites(sites);
  if (sites % 2) throw std::invalid_argument("boundary singlet dynamics requires even N");
  auto const count = heisenberg::real_excitation_count(sites, uni20::half_int{1}, max_candidates);
  auto& scheduler = execution.continuous_scheduler();
  R const zero = bethe::detail::solver_zero(solver);
  if (!uni20::isfinite(solver.residual_tolerance) || solver.residual_tolerance <= zero)
    throw std::invalid_argument("residual tolerance must be finite and positive");
  uni20::async::Async<BoundarySpectrum<R>> output;
  auto task = spectrum_task(scheduler, sites, count, std::move(solver), execution, std::move(evaluate),
                            std::move(prepare), output.write());
  bethe::detail::schedule_owned(scheduler, std::move(task));
  return output;
}
} // namespace boundary_detail

/// Owned solve+weight work per candidate; no intermediate root/weight batches.
/// The caller must keep the borrowed scheduler alive until the final result.
/// Ordinary numerical failures retain records. Cancellation/errors drain active
/// kernels and publish no partial spectrum. Dropping the handle does not cancel.
template <uni20::Real R = double>
[[nodiscard]] uni20::async::Async<BoundarySpectrum<R>>
boundary_real_triplet_spectrum_async(std::size_t sites, SolverOptions<R> solver = {},
                                     std::size_t max_candidates = 10000, ExecutionOptions execution = {})
{
  return boundary_detail::submit_spectrum(sites, std::move(solver), max_candidates, execution,
                                          boundary_detail::evaluate<R>);
}

template <uni20::Real R = double>
[[nodiscard]] BoundarySpectrum<R> boundary_real_triplet_spectrum(std::size_t sites, SolverOptions<R> solver = {},
                                                                 std::size_t max_candidates = 10000,
                                                                 ExecutionOptions execution = {})
{
  auto& scheduler = execution.continuous_scheduler();
  execution.scheduler = &scheduler;
  auto pending = boundary_real_triplet_spectrum_async<R>(sites, std::move(solver), max_candidates, execution);
  return bethe::detail::take_result(pending, scheduler);
}
} // namespace bethe::experimental::xxx::open
