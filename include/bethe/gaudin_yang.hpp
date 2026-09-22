// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <algorithm>
#include <array>
#include <bethe/solver.hpp>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <uni20/common/half_int.hpp>
#include <uni20/linalg/ops/linear_solve.hpp>
#include <vector>

namespace bethe::gaudin_yang
{
struct QuantumNumbers
{
    std::vector<uni20::half_int> charge, spin;
};
template <uni20::Real Real> using SolverOptions = bethe::SolverOptions<Real>;
enum class SolveStatus
{
  converged,
  iteration_limit,
  stalled
};

/// Periodic spin-1/2 fermions, H=-sum d_j^2+2c sum_(i<j) delta(x_i-x_j).
/// Interacting ground states currently require odd populations of both spins.
template <uni20::Real Real> struct State
{
    std::size_t up = 0, down = 0;
    Real length = Real{1}, interaction = Real{0};
    QuantumNumbers quantum_numbers;
    /// Physical k and lambda, not the dimensionless coordinates used internally.
    std::vector<Real> momenta, spin_rapidities;
    /// Populated only for exact free states; not interacting Bethe labels.
    std::array<std::vector<std::int64_t>, 2> free_modes;
    Real energy = Real{0}, momentum = Real{0};
    std::int64_t momentum_index = 0;
    Real charge_residual = Real{0}, spin_residual = Real{0}, residual_norm = Real{0};
    /// Physical coupling reached by continuation; diagnostics use requested c.
    Real root_interaction = Real{0};
    std::size_t iterations = 0;
    bool free = false, spin_reversed = false, converged = false;
    SolveStatus status = SolveStatus::iteration_limit;
};

namespace detail
{
inline std::size_t check_counts(std::size_t up, std::size_t down)
{
  auto const limit = std::size_t(std::numeric_limits<std::int64_t>::max() / 4);
  if (up > limit || down > limit - up)
    throw std::invalid_argument("Gaudin-Yang particle count exceeds the quantum-number range");
  return up + down;
}
inline void check_branch(std::size_t up, std::size_t down)
{
  check_counts(up, down);
  if (!up || !down || up % 2 == 0 || down % 2 == 0)
    throw std::invalid_argument("interacting Gaudin-Yang root branch requires odd N_up and N_down; other periodic "
                                "shell branches are not implemented");
}
inline QuantumNumbers centered_numbers(std::size_t n, std::size_t m)
{
  QuantumNumbers result;
  for (std::size_t j = 0; j < n; ++j)
    result.charge.push_back(uni20::from_twice(2 * std::int64_t(j) - std::int64_t(n - 1)));
  for (std::size_t j = 0; j < m; ++j)
    result.spin.push_back(uni20::from_twice(2 * std::int64_t(j) - std::int64_t(m - 1)));
  return result;
}

template <uni20::Real Real> class GroundSystem {
  public:
    GroundSystem(std::size_t n, std::size_t m)
        : n(n), m(m), nk(n / 2), ns(m / 2), order(nk + ns), pi(Real{4} * std::atan(Real{1})),
          labels(centered_numbers(n, m))
    {}
    static Real scale(Real g) { return std::max(Real{1}, g); }
    std::vector<Real> seed() const
    {
      std::vector<Real> x(order);
      for (std::size_t j = 0; j < nk; ++j)
        x[j] = pi * Real(labels.charge[nk + j].twice());
      for (std::size_t a = 0; a < ns; ++a)
        x[nk + a] = std::tan(pi * Real(a + 1) / Real(n)) / Real{2};
      return x;
    }
    std::array<std::vector<Real>, 2> expand(std::span<Real const> x, Real g) const
    {
      std::array<std::vector<Real>, 2> roots;
      for (std::size_t j = nk; j-- > 0;)
        roots[0].push_back(-x[j]);
      for (std::size_t j = 0; j < nk; ++j)
        roots[0].push_back(x[j]);
      for (std::size_t a = ns; a-- > 0;)
        roots[1].push_back(-scale(g) * x[nk + a]);
      roots[1].push_back(Real{0});
      for (std::size_t a = 0; a < ns; ++a)
        roots[1].push_back(scale(g) * x[nk + a]);
      return roots;
    }
    bool physical(std::span<Real const> x, Real g) const
    {
      if (x.size() != order) return false;
      Real const limit = uni20::numeric_limits<Real>::max() / Real{8};
      for (std::size_t i = 0; i < order; ++i)
      {
        if (!uni20::isfinite(x[i]) || x[i] <= Real{0} || x[i] > (i < nk ? limit : limit / scale(g)) ||
            (i && i != nk && x[i] <= x[i - 1]))
          return false;
      }
      return true;
    }
    // w/(w*w+d*d), without squaring very large or very small values.
    static Real kernel(Real d, Real w)
    {
      if (std::abs(d) > w)
      {
        Real const r = w / d;
        return (r / d) / (Real{1} + r * r);
      }
      Real const r = d / w;
      return (Real{1} / w) / (Real{1} + r * r);
    }
    struct Evaluation
    {
        std::vector<Real> residual;
        Real charge = Real{0}, spin = Real{0};
        Real norm() const { return std::max(charge, spin); }
    };
    Evaluation evaluate(std::span<Real const> x, Real g, uni20::DenseMatrix<Real>* jacobian = nullptr) const
    {
      auto const roots = expand(x, g);
      Evaluation result{.residual = std::vector<Real>(order)};
      if (jacobian)
        for (std::size_t i = 0; i < order; ++i)
          for (std::size_t j = 0; j < order; ++j)
            (*jacobian)[i, j] = Real{0};
      for (std::size_t row = 0; row < order; ++row)
      {
        bool const spin = row >= nk;
        std::size_t const a = spin ? 1 : 0, full = spin ? row - nk + ns + 1 : row + nk;
        Real const value = roots[a][full];
        auto const twice_label = spin ? labels.spin[full].twice() : labels.charge[full].twice();
        bethe::detail::CompensatedSum<Real> sum, diagonal;
        // At weak coupling, subtract exact pi-sized rank contributions before
        // summing the small complementary phases. This resolves the central
        // charge pair whose separation is O(sqrt(g)).
        std::int64_t rank = -twice_label;
        if (g >= Real{1}) sum.add(-pi * Real(twice_label));
        if (!spin)
        {
          sum.add(value);
          diagonal.add(Real{1});
        }
        for (std::size_t b = 0; b < 2; ++b)
        {
          if (!spin && b == 0) continue;
          Real const sign = spin && b == 1 ? -Real{1} : Real{1};
          Real const width = (spin && b == 1) ? g : g / Real{2};
          for (std::size_t j = 0; j < roots[b].size(); ++j)
          {
            if (a == b && j == full) continue;
            Real const d = value - roots[b][j];
            Real phase;
            if (g < Real{1} && std::abs(d) > width)
            {
              rank += (d > Real{0} ? 1 : -1) * (sign > Real{0} ? 1 : -1);
              phase = -Real{2} * std::atan(width / d);
            }
            else
              phase = Real{2} * std::atan(d / width);
            sum.add(sign * phase);
            if (jacobian)
            {
              Real const derivative = sign * Real{2} * kernel(d, width);
              diagonal.add(derivative * (spin ? scale(g) : Real{1}));
              auto const positives = b ? ns : nk;
              auto const first = positives + (b ? 1 : 0);
              if (j < positives)
                (*jacobian)[row, (b ? nk : 0) + positives - 1 - j] += derivative * (b ? scale(g) : Real{1});
              else if (j >= first)
                (*jacobian)[row, (b ? nk : 0) + j - first] -= derivative * (b ? scale(g) : Real{1});
            }
          }
        }
        if (g < Real{1}) sum.add(pi * Real(rank));
        if (jacobian) (*jacobian)[row, row] += diagonal.value();
        result.residual[row] = sum.value();
        if (!uni20::isfinite(sum.value())) throw std::runtime_error("nonfinite Gaudin-Yang residual");
        // Central weak-coupling charge roots need relative sqrt(g) resolution.
        // Other components use the physical root or unit equation scale.
        Real const normalization = g < Real{1} && !spin ? std::max(std::sqrt(g), std::abs(value)) : Real(n);
        auto& norm = spin ? result.spin : result.charge;
        norm = std::max(norm, std::abs(sum.value()) / normalization);
      }
      if (jacobian)
        for (std::size_t i = 0; i < order; ++i)
          for (std::size_t j = 0; j < order; ++j)
            if (!uni20::isfinite((*jacobian)[i, j]))
              throw std::overflow_error("Gaudin-Yang Jacobian exceeds this precision");
      return result;
    }
    std::size_t const n, m, nk, ns, order;
    Real const pi;
    QuantumNumbers const labels;
};

template <uni20::Real Real> void observables(State<Real>& state)
{
  bethe::detail::CompensatedSum<Real> energy;
  for (Real k : state.momenta)
  {
    if (k != Real{0} && k * k == Real{0})
      throw std::underflow_error("Gaudin-Yang kinetic energy underflows this precision; rescale the units");
    energy.add(k * k);
  }
  state.energy = energy.value();
  Real const pi = Real{4} * std::atan(Real{1});
  state.momentum = Real{2} * pi * Real(state.momentum_index) / state.length;
  if (!uni20::isfinite(state.energy) || !uni20::isfinite(state.momentum))
    throw std::overflow_error("Gaudin-Yang observables exceed this precision");
}

template <uni20::Real Real> Real physical_root(Real root, Real length)
{
  Real const result = root / length;
  if (!uni20::isfinite(result)) throw std::overflow_error("Gaudin-Yang physical root overflow; rescale the units");
  if (root != Real{0} && result == Real{0})
    throw std::underflow_error("Gaudin-Yang physical root underflow; rescale the units");
  return result;
}
} // namespace detail

/// Labels in the supported interacting branch; spin roots refer to the minority
/// species, with the majority as reference vacuum. Exact free states use modes.
inline QuantumNumbers ground_quantum_numbers(std::size_t up, std::size_t down)
{
  detail::check_branch(up, down);
  return detail::centered_numbers(up + down, std::min(up, down));
}

template <uni20::Real Real>
State<Real> ground_state(std::size_t up, std::size_t down, Real length, Real c, SolverOptions<Real> const& options = {})
{
  auto const n = detail::check_counts(up, down), m = std::min(up, down);
  if (!uni20::isfinite(length) || length <= Real{0} || !uni20::isfinite(c) || c < Real{0})
    throw std::invalid_argument("repulsive Gaudin-Yang requires finite length > 0 and c >= 0");
  if (!uni20::isfinite(options.residual_tolerance) || options.residual_tolerance <= Real{0})
    throw std::invalid_argument("residual tolerance must be finite and positive");
  State<Real> state;
  state.up = up;
  state.down = down;
  state.length = length;
  state.interaction = c;
  state.spin_reversed = down > up;
  if (c == Real{0} || m == 0)
  {
    state.free = state.converged = true;
    state.status = SolveStatus::converged;
    state.root_interaction = c;
    Real const pi = Real{4} * std::atan(Real{1});
    for (std::size_t species = 0; species < 2; ++species)
    {
      auto const count = species ? down : up;
      for (std::size_t j = 0; j < count; ++j)
      {
        auto const mode = std::int64_t(j) - std::int64_t((count - 1) / 2);
        state.free_modes[species].push_back(mode);
        state.momenta.push_back(detail::physical_root(Real{2} * pi * Real(mode), length));
      }
      if (count % 2 == 0) state.momentum_index += std::int64_t(count / 2);
    }
    std::sort(state.momenta.begin(), state.momenta.end());
    detail::observables(state);
    return state;
  }
  detail::check_branch(up, down);
  Real const target = c * length;
  if (!uni20::isfinite(target) || target / Real{2} <= Real{0})
    throw std::invalid_argument("c*length/2 must be finite and nonzero in the selected precision");
  auto const order = n / 2 + m / 2;
  auto const elements = std::size_t(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(Real);
  if (order > elements / order) throw std::length_error("Gaudin-Yang Newton matrix is too large");
  detail::GroundSystem<Real> const system(n, m);
  state.quantum_numbers = system.labels;
  Real g = std::max(target, Real{64} * Real(n));
  auto x = system.seed();
  uni20::DenseMatrix<Real> jacobian(order, order), step(order, 1);
  for (;;)
  {
    if (!system.physical(x, g)) throw std::runtime_error("Gaudin-Yang roots cannot be represented at this precision");
    auto const evaluation = system.evaluate(x, g, &jacobian);
    if (evaluation.norm() <= options.residual_tolerance)
    {
      if (g == target)
      {
        state.converged = true;
        state.status = SolveStatus::converged;
        break;
      }
      Real const next = std::max(target, g / Real{2});
      // Keep physical lambda*L fixed when its numerical scale changes. At
      // weak coupling the spin roots approach occupied free momenta, not zero.
      for (std::size_t a = 0; a < system.ns; ++a)
        x[system.nk + a] *= system.scale(g) / system.scale(next);
      g = next;
      continue;
    }
    if (state.iterations == options.max_iterations) break;
    for (std::size_t j = 0; j < order; ++j)
      step[j, 0] = -evaluation.residual[j];
    uni20::linalg::solve_inplace(jacobian, step);
    auto trial = x;
    Real damping = Real{1};
    bool accepted = false;
    for (int backtrack = 0; backtrack <= uni20::numeric_limits<Real>::digits; ++backtrack)
    {
      for (std::size_t j = 0; j < order; ++j)
        trial[j] = x[j] + damping * step[j, 0];
      if (system.physical(trial, g))
      {
        Real const norm = system.evaluate(trial, g).norm();
        if (norm <= options.residual_tolerance || norm < (Real{1} - damping / Real{10000}) * evaluation.norm())
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
  auto const roots = system.expand(x, g);
  for (Real q : roots[0])
    state.momenta.push_back(detail::physical_root(q, length));
  for (Real l : roots[1])
    state.spin_rapidities.push_back(detail::physical_root(l, length));
  // A stopped continuation must report residuals at the REQUESTED coupling.
  for (std::size_t a = 0; a < system.ns; ++a)
    x[system.nk + a] *= system.scale(g) / system.scale(target);
  auto const final = system.evaluate(x, target);
  state.charge_residual = final.charge;
  state.spin_residual = final.spin;
  state.residual_norm = final.norm();
  state.root_interaction = g / length;
  detail::observables(state);
  return state;
}
} // namespace bethe::gaudin_yang
