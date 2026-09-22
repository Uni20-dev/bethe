// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/xxz_open_state.hpp>
#include <uni20/linalg/ops/linear_solve.hpp>

namespace bethe::xxz::open::massive
{
using SolveStatus = open::GroundSolveStatus;
template <uni20::Real Real> using BoundaryRoot = open::BoundaryRoot<Real>;
template <uni20::Real Real> using State = open::GroundState<Real>;

namespace detail
{
template <uni20::Real Real> struct Parameters
{
    explicit Parameters(Real delta)
        : inverse(Real{1} / delta), p(Real{1} + inverse), q(-(delta - Real{1}) / delta),
          r((delta - Real{1}) / (delta + Real{1})), one_minus_r2(Real{4} * inverse / (p * p))
    {}
    Real inverse, p, q, r, one_minus_r2;
    Real y(Real w) const { return r * r * std::expm1(-w); }
    Real one_plus_y(Real w) const { return one_minus_r2 + r * r * std::exp(-w); }
};

template <uni20::Real Real> class GroundSystem {
  public:
    GroundSystem(std::size_t sites, std::size_t roots)
        : sites(sites), order(roots), boundary(2 * roots == sites), bulk(roots - std::size_t(boundary)),
          pi(Real{4} * std::atan(Real{1}))
    {}

    std::vector<Real> seed(Real delta) const
    {
      std::vector<Real> x(order);
      for (std::size_t i = 0; i < order; ++i)
        x[i] = pi * Real(i + 1) / (Real{2} * Real(sites + 1));
      if (boundary)
      {
        Real const z = std::tan(x.back()), r = Parameters<Real>(delta).r;
        x.back() = -std::log1p(Real{1} / (z * z * r * r));
      }
      return x;
    }

    bool physical(std::span<Real const> x, Real delta) const
    {
      if (x.size() != order) return false;
      for (std::size_t i = 0; i < bulk; ++i)
        if (!uni20::isfinite(x[i]) || x[i] <= Real{0} || x[i] >= pi / Real{2} || (i && x[i] <= x[i - 1])) return false;
      if (boundary)
      {
        if (!uni20::isfinite(x.back())) return false;
        Real const y = Parameters<Real>(delta).y(x.back());
        if (!uni20::isfinite(y)) return false;
        // While real, the boundary root is the largest positive root. This
        // excludes an unphysical branch with interchanged/colliding roots.
        if (bulk && y > Real{0})
        {
          Real const z = std::tan(x[bulk - 1]);
          if (y * z * z >= Real{1}) return false;
        }
      }
      return true;
    }

    struct Evaluation
    {
        std::vector<Real> residual;
        Real norm = Real{0};
    };

    Evaluation evaluate(std::span<Real const> x, Real delta) const
    {
      using std::atan;
      using std::atan2;
      using std::log;
      using std::log1p;
      using std::sqrt;
      Parameters<Real> const a(delta);
      std::vector<Real> z(bulk);
      for (std::size_t i = 0; i < bulk; ++i)
        z[i] = std::tan(x[i]);
      Real const y = boundary ? a.y(x.back()) : Real{0};
      Evaluation result{.residual = std::vector<Real>(order)};
      for (std::size_t i = 0; i < bulk; ++i)
      {
        Real const u = z[i];
        bethe::detail::CompensatedSum<Real> f;
        f.add(Real{2} * Real(sites) * x[i] - pi * Real(i + 1));
        f.add(-Real{2} * atan(a.r * u));
        for (std::size_t j = 0; j < bulk; ++j)
          if (i != j)
          {
            f.add(-atan2(u - z[j], a.p - a.q * u * z[j]));
            f.add(-atan2(u + z[j], a.p + a.q * u * z[j]));
          }
        if (boundary)
        {
          // The combined phase's naive arguments both cancel near the pole
          // at strong coupling. Retain the 1/Delta pieces and the deviation.
          Real const deviation = a.r * a.r * std::exp(-x.back());
          Real const numerator = Real{2} * u * (Real{2} * a.r * a.inverse + a.p * deviation);
          Real const denominator = a.inverse * (Real{2} - a.inverse - a.r * a.r * (Real{2} + a.inverse) * u * u) +
                                   (a.p * a.p - u * u) * deviation;
          f.add(-atan2(numerator, denominator));
        }
        result.residual[i] = f.value() / Real(sites);
      }
      if (boundary)
      {
        Real const w = x.back();
        Real basic, boundary_phase;
        if (w < Real{0})
        {
          Real const v = sqrt(y), t = sqrt(std::expm1(-w));
          basic = atan(v) / v;
          boundary_phase = atan(t) / t;
        }
        else if (w > Real{0})
        {
          Real const t = sqrt(-std::expm1(-w)), v = a.r * t;
          // log(1-v^2) is evaluated without subtracting nearly equal values.
          basic = v < Real{1} / Real{2} ? std::atanh(v) / v : (log1p(v) - log(a.one_plus_y(w)) / Real{2}) / v;
          // atanh(t) = w/2+log1p(t), even when t rounds to one.
          boundary_phase = (w / Real{2} + log1p(t)) / t;
        }
        else
          basic = boundary_phase = Real{1};
        bethe::detail::CompensatedSum<Real> f;
        f.add(-Real(sites) * basic);
        f.add(boundary_phase / a.r);
        for (Real u : z)
        {
          Real const den = Real{1} + a.q * a.q * u * u;
          Real const re = (a.p - a.q * u * u) / den;
          Real const im = a.inverse * a.inverse * u / den;
          Real term;
          if (y > Real{0})
          {
            Real const v = sqrt(y), mod2 = (a.p * a.p + u * u) / den;
            term = atan2(Real{2} * v * re, Real{1} - y * mod2) / (Real{2} * v);
          }
          else if (y < Real{0})
          {
            Real const t = sqrt(-std::expm1(-w)), v = a.r * t;
            // 1-r*Re(a_j) has a positive rational form. Keep it, and
            // 1-t=exp(-w)/(1+t), instead of subtracting numbers near one.
            Real const minus =
                a.inverse * (Real{1} + a.r * a.r * a.p * u * u) / den + a.r * re * std::exp(-w) / (Real{1} + t);
            if (v < Real{1} / Real{2})
              term = log1p(Real{4} * v * re / (minus * minus + v * v * im * im)) / (Real{4} * v);
            else
              term = (log(std::hypot(Real{1} + v * re, v * im)) - log(std::hypot(minus, v * im))) / (Real{2} * v);
          }
          else
            term = re;
          f.add(term);
        }
        result.residual.back() = f.value() / Real(sites);
      }
      for (Real f : result.residual)
      {
        if (!uni20::isfinite(f)) throw std::runtime_error("nonfinite massive open XXZ residual");
        result.norm = std::max(result.norm, std::abs(f));
      }
      return result;
    }

    std::size_t sites, order;
    bool boundary;
    std::size_t bulk;
    Real pi;
};
} // namespace detail

/// Damped Newton in bulk atan(z) and a logarithmic boundary-root deviation.
/// Continuation follows the ground-state branch from near Delta=1; every
/// accepted Newton update counts against the single global iteration budget.
/// Dense storage is O(M^2); each finite-difference Newton update is O(M^3).
template <uni20::Real Real = double>
[[nodiscard]] State<Real> sector_ground_state(std::size_t sites, Real delta, uni20::half_int sz,
                                              SolverOptions<Real> const& options = {})
{
  auto const m = xxz::detail::sector_roots(sites, sz);
  if (!uni20::isfinite(delta) || delta <= Real{1})
    throw std::invalid_argument("massive open XXZ requires finite Delta > 1");
  if (!uni20::isfinite(options.residual_tolerance) || options.residual_tolerance <= Real{0})
    throw std::invalid_argument("residual tolerance must be finite and positive");
  State<Real> state;
  state.residual_convention = open::GroundResidualConvention::massive_regularized;
  state.delta = delta;
  state.sz = sz;
  state.spin_reversed = sz.twice() < 0;
  if (m == 0)
  {
    state.energy = Real(sites - 1) * (delta / Real{4});
    if (!uni20::isfinite(state.energy)) throw std::overflow_error("massive open XXZ energy overflow");
    state.root_delta = delta;
    state.converged = true;
    state.status = SolveStatus::converged;
    return state;
  }
  auto const elements = std::size_t(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(Real);
  if (m > elements / m) throw std::length_error("massive open XXZ Newton matrix is too large");
  detail::GroundSystem<Real> const system(sites, m);
  Real d = std::min(delta, Real{11} / Real{10});
  auto x = system.seed(d);
  uni20::DenseMatrix<Real> jacobian(m, m), step(m, 1);
  Real const difference_scale = std::cbrt(uni20::numeric_limits<Real>::epsilon());
  for (;;)
  {
    if (!system.physical(x, d)) throw std::runtime_error("massive open XXZ roots left the ground-state branch");
    auto const evaluation = system.evaluate(x, d);
    if (evaluation.norm <= options.residual_tolerance)
    {
      if (d == delta)
      {
        state.converged = true;
        state.status = SolveStatus::converged;
        break;
      }
      if (state.iterations == options.max_iterations) break;
      Real const next = std::min(delta, d + (d - Real{1}) / Real{4});
      // Before crossing infinity, keep the real z fixed. Afterwards, keep
      // the logarithmic distance to the moving boundary pole as the guess;
      // fixing y would discard exponentially small deviations at each step.
      if (system.boundary && x.back() < Real{0})
      {
        Real const y = detail::Parameters<Real>(d).y(x.back()), r = detail::Parameters<Real>(next).r;
        x.back() = -std::log1p(y / (r * r));
      }
      d = next;
      continue;
    }
    if (state.iterations == options.max_iterations) break;
    for (std::size_t j = 0; j < m; ++j)
    {
      auto plus = x, minus = x;
      Real const h = difference_scale * (Real{1} + std::abs(x[j]));
      plus[j] += h;
      minus[j] -= h;
      auto const fp = system.evaluate(plus, d), fm = system.evaluate(minus, d);
      for (std::size_t i = 0; i < m; ++i)
        jacobian[i, j] = (fp.residual[i] - fm.residual[i]) / (Real{2} * h);
      step[j, 0] = -evaluation.residual[j];
    }
    uni20::linalg::solve_inplace(jacobian, step);
    auto trial = x;
    Real damping = Real{1};
    bool accepted = false;
    for (int backtrack = 0; backtrack <= uni20::numeric_limits<Real>::digits; ++backtrack)
    {
      for (std::size_t j = 0; j < m; ++j)
        trial[j] = x[j] + damping * step[j, 0];
      if (system.physical(trial, d))
      {
        Real const norm = system.evaluate(trial, d).norm;
        if (norm <= options.residual_tolerance || norm < (Real{1} - damping / Real{10000}) * evaluation.norm)
        {
          accepted = true;
          break;
        }
      }
      damping /= Real{2};
    }
    if (!accepted)
    {
      state.status = SolveStatus::stalled;
      break;
    }
    x = std::move(trial);
    ++state.iterations;
  }
  state.root_delta = d;
  if (system.boundary && d != delta)
  {
    Real const y = detail::Parameters<Real>(d).y(x.back()), r = detail::Parameters<Real>(delta).r;
    x.back() = -std::log1p(y / (r * r));
  }
  // Even on budget exhaustion, roots, residuals and energy refer to the
  // requested Hamiltonian; root_delta says how far continuation progressed.
  state.residual_norm = system.evaluate(x, delta).norm;
  bethe::detail::CompensatedSum<Real> energy;
  energy.add((Real(sites - 1) / Real{4} - Real(m)) * delta);
  for (std::size_t i = 0; i < system.bulk; ++i)
  {
    Real const z = std::tan(x[i]);
    state.rapidities.push_back(z);
    state.quantum_numbers.emplace_back(std::int64_t(i + 1));
    energy.add(-std::cos(Real{2} * x[i]));
  }
  if (system.boundary)
  {
    detail::Parameters<Real> const parameters(delta);
    state.boundary_root = BoundaryRoot<Real>{x.back(), parameters.y(x.back()), uni20::half_int(std::int64_t(m))};
    energy.add(Real{2} / parameters.one_plus_y(x.back()) - Real{1});
  }
  state.energy = energy.value();
  if (!uni20::isfinite(state.energy)) throw std::overflow_error("massive open XXZ energy overflow");
  return state;
}

template <uni20::Real Real = double>
[[nodiscard]] State<Real> ground_state(std::size_t sites, Real delta, SolverOptions<Real> const& options = {})
{
  return massive::sector_ground_state<Real>(sites, delta, uni20::from_twice(std::int64_t(sites % 2)), options);
}
} // namespace bethe::xxz::open::massive
