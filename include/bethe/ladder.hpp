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
/// H=sum [S.S_next+T.T_next+4(S.S_next)(T.T_next)] + J_r*sum S.T - h*sum(Sz+Tz).
/// Equivalently sum(P_rung-1/4)+J_r*(L/4-N_s)-h*(N_+-N_-).
template <uni20::Real Real> struct State
{
    std::size_t rungs = 0, singlets = 0, iterations = 0, branches = 0, tableaux = 0;
    detail::Shape populations{}; // singlet,t+,t0,t-; minimizing representative, not a multiplet average.
    detail::Branch<Real> highest_weight;
    Real rung_coupling{}, magnetic_field{};
    std::int64_t magnetization = 0; // total physical Sz=N_+-N_-, not Sz per rung
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
// GL(4)->GL(3) interlacing at fixed singlet count. The triplet diagram mu
// satisfies lambda_i >= mu_i >= lambda_(i+1), sum(mu)=L-N_s.
// Its largest physical spin projection is mu_0-mu_2. Maximizing mu_0
// and then minimizing mu_2 gives the extremal weight in constant time.
inline Shape polarized_populations(Shape const& shape, std::size_t singlets, bool reverse)
{
  if (singlets < shape[3] || singlets > shape[0])
    throw std::invalid_argument("singlet population absent from SU(4) multiplet");
  std::size_t const total = shape[0] + shape[1] + shape[2] + shape[3] - singlets;
  auto const plus = std::min(shape[0], total - shape[2] - shape[3]);
  auto const remainder = total - plus;
  auto const minus = std::max(shape[3], remainder > shape[1] ? remainder - shape[1] : 0);
  Shape result{singlets, plus, remainder - minus, minus};
  if (reverse) std::swap(result[1], result[3]);
  return result;
}
inline std::int64_t magnetization(Shape const& populations)
{
  return std::int64_t(populations[1]) - std::int64_t(populations[3]);
}
// A four-color weight belongs to lambda iff every single entry is between
// lambda_3 and lambda_0 and every pair sums to at most lambda_0+lambda_1.
// With N_s fixed, these become identical lower/upper bounds on the three
// triplet populations. Intersect with (N_+,N_0,N_-)=(x+|M|,T-|M|-2x,x).
inline std::optional<Shape> fixed_magnetization_populations(Shape const& shape, std::size_t singlets,
                                                            std::int64_t projection)
{
  auto const l = shape[0] + shape[1] + shape[2] + shape[3];
  if (singlets < shape[3] || singlets > shape[0] || projection < -std::int64_t(l - singlets) ||
      projection > std::int64_t(l - singlets))
    return {};
  auto const m = std::size_t(projection < 0 ? -projection : projection), total = l - singlets;
  auto const pair = shape[0] + shape[1];
  auto const lower = std::max(shape[3], total > pair ? total - pair : 0);
  auto const upper = std::min(shape[0], pair - singlets);
  if (m > upper || total - m < lower) return {};
  auto const excess = total - m > upper ? total - m - upper : 0;
  auto const first = std::max(lower, excess / 2 + excess % 2);
  auto const last = std::min(upper - m, (total - m - lower) / 2);
  if (first > last) return {};
  Shape weight{singlets, first + m, total - m - 2 * first, first};
  if (projection < 0) std::swap(weight[1], weight[3]);
  return weight;
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
  s.magnetization = magnetization(s.populations);
  bethe::detail::CompensatedSum<Real> sum;
  sum.add(branch.energy);
  sum.add(-Real(s.rungs) / Real{4});
  sum.add(s.rung_coupling * (Real(s.rungs) / Real{4} - Real(s.singlets)));
  sum.add(-s.magnetic_field * Real(s.magnetization));
  s.energy = sum.value();
  if (!uni20::isfinite(*s.energy)) throw std::runtime_error("ladder energy is not representable at this precision");
}
template <uni20::Real Real> State<Real> product(std::size_t l, Real rung, Real field, std::size_t color)
{
  State<Real> s;
  s.rungs = l;
  s.singlets = color == 0 ? l : 0;
  s.rung_coupling = rung;
  s.magnetic_field = field;
  s.populations[color] = l;
  Branch<Real> b;
  b.shape = {l, 0, 0, 0};
  b.energy = Real(l);
  b.converged = true;
  b.status = BranchStatus::converged;
  set_energy(s, b);
  s.converged = s.analytic = true;
  s.status = SolveStatus::converged;
  return s;
}
template <uni20::Real Real>
SectorScan<Real> scan(std::size_t l, Real rung, Real field, std::optional<std::size_t> only,
                      SolverOptions<Real> const& options, std::optional<std::int64_t> projection = {})
{
  validate(l, rung, options);
  if (!uni20::isfinite(field)) throw std::invalid_argument("ladder magnetic field must be finite");
  if (only && *only > l) throw std::invalid_argument("singlet count must not exceed the number of rungs");
  if (projection && (*projection < -std::int64_t(l) || *projection > std::int64_t(l)))
    throw std::invalid_argument("ladder magnetization must be an integer in [-L,L]");
  auto const magnitude = projection ? std::size_t(*projection < 0 ? -*projection : *projection) : 0;
  if (only && *only > l - magnitude)
    throw std::invalid_argument("singlet count is incompatible with the requested magnetization");
  if (only && *only == l)
    return {.sectors = {product(l, rung, field, 0)}, .complete = true, .status = SolveStatus::converged};
  if (magnitude == l)
    return {.sectors = {product(l, rung, field, *projection > 0 ? 1 : 3)},
            .complete = true,
            .status = SolveStatus::converged};
  auto const order = l + l / 2, elements = std::size_t(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(Real);
  if (order > elements / order) throw std::length_error("ladder Newton matrix is too large");
  SectorScan<Real> out;
  auto const begin = only ? *only : 0, end = only ? *only : l - magnitude;
  for (std::size_t ns = begin; ns <= end; ++ns)
  {
    State<Real> s;
    s.rungs = l;
    s.singlets = ns;
    s.rung_coupling = rung;
    s.magnetic_field = field;
    s.populations = populations(l, ns);
    out.sectors.push_back(std::move(s));
  }
  Shape shape{};
  // Enumerate all highest weights, including those more dominant than a
  // requested population. Their descendants can be lower than that weight's
  // own filled-sea state on a ring; omitting them gives wrong sector minima.
  auto visit = [&]() {
    std::vector<std::optional<Shape>> weights;
    weights.reserve(out.sectors.size());
    bool relevant = false;
    for (auto const& s : out.sectors)
    {
      std::optional<Shape> weight;
      if (projection)
        weight = fixed_magnetization_populations(shape, s.singlets, *projection);
      else if (dominates(shape, populations(l, s.singlets)))
        weight =
            field == Real{0} ? populations(l, s.singlets) : polarized_populations(shape, s.singlets, field < Real{0});
      relevant = relevant || weight.has_value();
      weights.push_back(weight);
    }
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
      for (std::size_t index = 0; index < out.sectors.size(); ++index)
        if (weights[index])
        {
          auto& s = out.sectors[index];
          auto const& weight = *weights[index];
          // Cancel common shifts before comparing: a large field must not
          // erase differences between branches with identical magnetization.
          if (!s.energy ||
              branch.energy - s.highest_weight.energy < field * Real(magnetization(weight) - s.magnetization))
          {
            s.populations = weight;
            set_energy(s, branch);
          }
        }
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

/// Select the lowest retained candidate from one scan. Its convergence flag
/// still describes the whole scan; a partial scan does not certify a minimum.
template <uni20::Real Real> [[nodiscard]] State<Real> minimum_candidate(SectorScan<Real> const& scan)
{
  if (scan.sectors.empty()) throw std::invalid_argument("cannot select from an empty ladder scan");
  auto best = scan.sectors.front();
  for (auto const& state : scan.sectors)
    if (state.energy)
    {
      bethe::detail::CompensatedSum<Real> shift;
      shift.add(state.rung_coupling * Real(std::int64_t(state.singlets) - std::int64_t(best.singlets)));
      shift.add(state.magnetic_field * Real(state.magnetization - best.magnetization));
      if (!best.energy || state.highest_weight.energy - best.highest_weight.energy < shift.value()) best = state;
    }
  return best;
}

/// Minimize triplet populations too: balanced representatives at zero field,
/// extremal physical Sz in each compatible SU(4) multiplet at nonzero field.
template <uni20::Real Real = double>
[[nodiscard]] SectorScan<Real> sector_ground_states(std::size_t rungs, Real rung_coupling, Real magnetic_field,
                                                    SolverOptions<Real> const& options = {})
{
  return detail::scan(rungs, rung_coupling, magnetic_field, std::nullopt, options);
}

template <uni20::Real Real = double>
[[nodiscard]] State<Real> sector_ground_state(std::size_t rungs, std::size_t singlets, Real rung_coupling,
                                              Real magnetic_field, SolverOptions<Real> const& options = {})
{
  return detail::scan(rungs, rung_coupling, magnetic_field, singlets, options).sectors.front();
}

template <uni20::Real Real = double>
[[nodiscard]] State<Real> ground_state(std::size_t rungs, Real rung_coupling, Real magnetic_field,
                                       SolverOptions<Real> const& options = {})
{
  detail::validate(rungs, rung_coupling, options);
  if (!uni20::isfinite(magnetic_field)) throw std::invalid_argument("ladder magnetic field must be finite");
  // Operator bound sum(P-1)>=-4*N_triplet proves the rung product is a
  // ground state for J_r-|h|>=4, independently of parity and iteration budgets.
  auto const magnitude = std::abs(magnetic_field);
  if (rung_coupling >= Real{4} && rung_coupling - magnitude >= Real{4})
    return detail::product(rungs, rung_coupling, magnetic_field, 0);
  // The same bound about the favored triplet applies when all other rung
  // colors cost at least 4: |h|>=4+max(J_r,0). Sufficient, not necessary.
  if (magnitude >= Real{4} && magnitude - std::max(rung_coupling, Real{0}) >= Real{4})
    return detail::product(rungs, rung_coupling, magnetic_field, magnetic_field > Real{0} ? 1 : 3);
  auto const scan = sector_ground_states(rungs, rung_coupling, magnetic_field, options);
  return minimum_candidate(scan);
}
// Preserve the zero-field API and its balanced-triplet representative.
template <uni20::Real Real = double>
[[nodiscard]] SectorScan<Real> sector_ground_states(std::size_t rungs, Real rung,
                                                    SolverOptions<Real> const& options = {})
{
  return sector_ground_states(rungs, rung, Real{0}, options);
}
template <uni20::Real Real = double>
[[nodiscard]] State<Real> sector_ground_state(std::size_t rungs, std::size_t singlets, Real rung,
                                              SolverOptions<Real> const& options = {})
{
  return sector_ground_state(rungs, singlets, rung, Real{0}, options);
}
template <uni20::Real Real = double>
[[nodiscard]] State<Real> ground_state(std::size_t rungs, Real rung, SolverOptions<Real> const& options = {})
{
  return ground_state(rungs, rung, Real{0}, options);
}

/// All allowed singlet-count sectors at fixed total physical Sz. No parity
/// restriction on integer Sz: a rung triplet also has a zero-projection state.
template <uni20::Real Real = double>
[[nodiscard]] SectorScan<Real> magnetization_sector_ground_states(std::size_t rungs, std::int64_t projection, Real rung,
                                                                  Real field = Real{0},
                                                                  SolverOptions<Real> const& options = {})
{
  return detail::scan(rungs, rung, field, std::nullopt, options, projection);
}
/// Lowest energy at fixed total Sz, optionally also fixing the singlet count.
/// Field only adds -h*Sz here; it cannot change the selected multiplet.
template <uni20::Real Real = double>
[[nodiscard]] State<Real> magnetization_ground_state(std::size_t rungs, std::int64_t projection, Real rung,
                                                     Real field = Real{0}, SolverOptions<Real> const& options = {},
                                                     std::optional<std::size_t> singlets = {})
{
  return minimum_candidate(detail::scan(rungs, rung, field, singlets, options, projection));
}
} // namespace bethe::ladder
