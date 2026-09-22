// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/su4_seas.hpp>
#include <functional>
#include <optional>

namespace bethe::ladder
{
template <uni20::Real Real> struct SolverOptions : bethe::SolverOptions<Real>
{
    // max_iterations counts attempted Newton corrections across the whole
    // scan. Already-exact seas need none; it is not reset for each branch.
    std::size_t max_branches = 10000; // Whole scan, not per singlet sector.
};
enum class SolveStatus
{
  converged,
  iteration_limit,
  branch_limit,
  stalled,
  ill_conditioned
};

/// Spin-1/2 ladder, L periodic rungs, leg coefficient 1, four-spin coefficient 4:
/// H=sum [S.S_next+T.T_next+4(S.S_next)(T.T_next)] + J_r*sum S.T.
/// Equivalently sum(P_rung-1/4)+J_r*(L/4-N_s). No magnetic field.
template <uni20::Real Real> struct State
{
    std::size_t rungs = 0, singlets = 0, iterations = 0, branches = 0, tableaux = 0;
    detail::Shape populations{}; // singlet,t+,t0,t-; a balanced triplet representative.
    detail::Branch<Real> highest_weight;
    Real rung_coupling{};
    // Only solved branches contribute. On an incomplete scan this is a
    // candidate energy, not a minimum; check converged before using it as one.
    std::optional<Real> energy;
    bool converged = false, descendant = false, analytic = false;
    SolveStatus status = SolveStatus::iteration_limit;
};
template <uni20::Real Real> struct SectorScan
{
    std::vector<State<Real>> sectors;
    std::size_t iterations = 0, branches = 0, tableaux = 0;
    bool complete = false;
    SolveStatus status = SolveStatus::iteration_limit;
};

namespace detail
{
inline Shape populations(std::size_t l, std::size_t singlets)
{
  auto const t = l - singlets, q = t / 3, r = t % 3;
  return {singlets, q + (r > 0), q + (r > 1), q};
}
template <uni20::Real Real> void validate(std::size_t l, Real rung, SolverOptions<Real> const& options)
{
  if (l < 2 || l > std::size_t(std::numeric_limits<std::int64_t>::max() / 8))
    throw std::invalid_argument("integrable ladder requires 2<=rungs<=INT64_MAX/8");
  if (!uni20::isfinite(rung)) throw std::invalid_argument("rung coupling must be finite");
  if (!uni20::isfinite(options.residual_tolerance) || options.residual_tolerance <= Real{0})
    throw std::invalid_argument("ladder residual tolerance must be finite and positive");
}
template <uni20::Real Real> void set_energy(State<Real>& s, Branch<Real> const& branch)
{
  s.highest_weight = branch;
  auto sorted = s.populations;
  std::sort(sorted.begin(), sorted.end(), std::greater<>{});
  s.descendant = sorted != branch.shape;
  bethe::detail::CompensatedSum<Real> sum;
  sum.add(branch.energy);
  sum.add(-Real(s.rungs) / Real{4});
  sum.add(s.rung_coupling * (Real(s.rungs) / Real{4} - Real(s.singlets)));
  s.energy = sum.value();
  if (!uni20::isfinite(*s.energy)) throw std::runtime_error("ladder energy is not representable at this precision");
}
template <uni20::Real Real> State<Real> singlet_product(std::size_t l, Real rung)
{
  State<Real> s;
  s.rungs = s.singlets = l;
  s.rung_coupling = rung;
  s.populations = {l, 0, 0, 0};
  Branch<Real> b;
  b.shape = s.populations;
  b.energy = Real(l);
  b.converged = true;
  b.status = BranchStatus::converged;
  set_energy(s, b);
  s.converged = s.analytic = true;
  s.status = SolveStatus::converged;
  return s;
}
template <uni20::Real Real>
SectorScan<Real> scan(std::size_t l, Real rung, std::optional<std::size_t> only, SolverOptions<Real> const& options)
{
  validate(l, rung, options);
  if (only && *only > l) throw std::invalid_argument("singlet count must not exceed the number of rungs");
  if (only && *only == l)
    return {.sectors = {singlet_product(l, rung)}, .complete = true, .status = SolveStatus::converged};
  auto const order = l + l / 2, elements = std::size_t(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(Real);
  if (order > elements / order) throw std::length_error("ladder Newton matrix is too large");
  SectorScan<Real> out;
  auto const begin = only ? *only : 0, end = only ? *only : l;
  for (std::size_t ns = begin; ns <= end; ++ns)
  {
    State<Real> s;
    s.rungs = l;
    s.singlets = ns;
    s.rung_coupling = rung;
    s.populations = populations(l, ns);
    out.sectors.push_back(std::move(s));
  }
  Shape shape{};
  // Enumerate all highest weights, including those more dominant than a
  // requested population. Their descendants can be lower than that weight's
  // own filled-sea state on a ring; omitting them gives wrong sector minima.
  auto visit = [&]() {
    bool relevant = false;
    for (auto const& s : out.sectors)
      if (dominates(shape, s.populations)) relevant = true;
    if (!relevant) return true;
    ++out.tableaux;
    for (auto const& shift : shifts(shape))
    {
      if (out.branches == options.max_branches)
      {
        out.status = SolveStatus::branch_limit;
        return false;
      }
      ++out.branches;
      auto const branch =
          solve<Real>(shape, shift, options.residual_tolerance, options.max_iterations - out.iterations);
      out.iterations += branch.iterations;
      if (!branch.converged)
      {
        out.status = branch.status == BranchStatus::iteration_limit   ? SolveStatus::iteration_limit
                     : branch.status == BranchStatus::ill_conditioned ? SolveStatus::ill_conditioned
                                                                      : SolveStatus::stalled;
        return false;
      }
      for (auto& s : out.sectors)
        if (dominates(shape, s.populations) && (!s.energy || branch.energy < s.highest_weight.energy))
          set_energy(s, branch);
    }
    return true;
  };
  std::function<bool(std::size_t, std::size_t, std::size_t)> enumerate = [&](std::size_t row, std::size_t remaining,
                                                                             std::size_t maximum) {
    if (row == 3)
    {
      if (remaining > maximum) return true;
      shape[row] = remaining;
      return visit();
    }
    auto const minimum = (remaining + (3 - row)) / (4 - row);
    for (std::size_t value = std::min(maximum, remaining);; --value)
    {
      if (value < minimum) break;
      shape[row] = value;
      if (!enumerate(row + 1, remaining - value, value)) return false;
      if (value == 0) break;
    }
    return true;
  };
  out.complete = enumerate(0, l, l);
  if (out.complete) out.status = SolveStatus::converged;
  for (auto& s : out.sectors)
  {
    s.converged = out.complete;
    s.status = out.status;
    s.iterations = out.iterations;
    s.branches = out.branches;
    s.tableaux = out.tableaux;
  }
  return out;
}
} // namespace detail

/// Minimize the triplet populations too; each SU(3) multiplet has a weight in
/// the balanced triplet sector, so one representative suffices for each N_s.
template <uni20::Real Real = double>
[[nodiscard]] SectorScan<Real> sector_ground_states(std::size_t rungs, Real rung_coupling,
                                                    SolverOptions<Real> const& options = {})
{
  return detail::scan(rungs, rung_coupling, std::nullopt, options);
}

template <uni20::Real Real = double>
[[nodiscard]] State<Real> sector_ground_state(std::size_t rungs, std::size_t singlets, Real rung_coupling,
                                              SolverOptions<Real> const& options = {})
{
  return detail::scan(rungs, rung_coupling, singlets, options).sectors.front();
}

template <uni20::Real Real = double>
[[nodiscard]] State<Real> ground_state(std::size_t rungs, Real rung_coupling, SolverOptions<Real> const& options = {})
{
  detail::validate(rungs, rung_coupling, options);
  // Operator bound sum(P-1)>=-4*N_triplet proves the rung product is a
  // ground state for J_r>=4, independently of parity and iteration budgets.
  if (rung_coupling >= Real{4}) return detail::singlet_product(rungs, rung_coupling);
  auto const scan = sector_ground_states(rungs, rung_coupling, options);
  auto best = scan.sectors.front();
  for (auto const& s : scan.sectors)
    if (s.energy && (!best.energy || *s.energy < *best.energy)) best = s;
  return best;
}
} // namespace bethe::ladder
