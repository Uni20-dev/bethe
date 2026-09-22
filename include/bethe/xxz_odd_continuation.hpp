// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/xxz_odd_bounds.hpp>
#include <bethe/xxz_polynomial.hpp>
#include <uni20/linalg/ops/linear_solve.hpp>

namespace bethe::xxz::detail
{
enum class PolynomialContinuationStatus
{
  equations_converged,
  iteration_limit,
  stalled,
  ill_conditioned,
  branch_rejected
};

/// Internal numerical continuation, NOT a physical-state/ground-state
/// certificate. The intended branch starts at the odd-ring free ground sea.
template <uni20::Real Real> struct OddPolynomialBranch
{
    std::vector<Real> coefficients; // Q(x), physical z=center+coordinate_scale*x
    Real delta = Real{0}, root_delta = Real{0};
    Real center = Real{0}, coordinate_scale = Real{1};
    Real energy = Real{0}, residual_norm = Real{0}, momentum_error = Real{0};
    Real reciprocal_condition = Real{1}, correction_norm = Real{0};
    Real variational_upper_bound = Real{0}; // At requested Delta, even on failure.
    std::size_t iterations = 0, rejected_stages = 0;
    std::size_t variational_rejections = 0;
    std::size_t momentum_index = 0;
    bool equations_converged = false;
    PolynomialContinuationStatus status = PolynomialContinuationStatus::iteration_limit;
};

template <uni20::Real Real> class OddPolynomialCoordinates {
  public:
    using Complex = std::complex<Real>;
    OddPolynomialCoordinates(std::size_t sites, std::size_t roots) : sites(sites), order(roots)
    {
      checked_sites(sites);
      if (sites % 2 == 0 || roots > sites / 2)
        throw std::invalid_argument("odd XXZ continuation requires odd N and M <= N/2");
      Real const pi = Real{4} * std::atan(Real{1});
      center = -std::tan(pi / (Real{2} * Real(sites)));
      momentum_index = roots % 2 ? (sites + roots) / 2 : roots / 2;
      Real const angle = Real{2} * pi * Real(momentum_index) / Real(sites);
      target_phase = {std::cos(angle), std::sin(angle)};
    }

    Real scale(Real delta) const { return std::sqrt((Real{1} + delta) / (Real{1} - delta) + center * center); }

    std::vector<Real> seed() const
    {
      std::vector<Real> c{Real{1}};
      Real const pi = Real{4} * std::atan(Real{1}), s = scale(Real{0});
      for (std::size_t j = 0; j < order; ++j)
      {
        Real const z = std::tan(pi * (Real(j) - Real(order) / Real{2}) / Real(sites));
        Real const x = (z - center) / s;
        std::vector<Real> next(c.size() + 1, Real{0});
        for (std::size_t k = 0; k < c.size(); ++k)
        {
          next[k] -= x * c[k];
          next[k + 1] += c[k];
        }
        c = std::move(next);
      }
      c.pop_back();
      impose_momentum(c, weights(s));
      return c;
    }

    // With this center, (i-center)^M already has the desired momentum phase.
    // Hence Im[1+sum c[k]*t^(M-k)]=0, t=scale/(i-center), is a LINEAR
    // constraint. Eliminate c[M-1], whose weight Im(t) never vanishes.
    std::vector<Real> weights(Real s) const
    {
      std::vector<Real> out(order ? order - 1 : 0);
      Complex const t = s / Complex{-center, Real{1}};
      Complex power = t;
      for (std::size_t degree = 2; degree <= order; ++degree)
      {
        power *= t;
        out[order - degree] = -power.imag() / t.imag();
      }
      return out;
    }

    static void impose_momentum(std::vector<Real>& c, std::span<Real const> weights)
    {
      if (c.empty()) return;
      bethe::detail::CompensatedSum<Real> sum;
      for (std::size_t j = 0; j < weights.size(); ++j)
        sum.add(c[j] * weights[j]);
      c.back() = sum.value();
    }

    static void rescale(std::vector<Real>& c, Real old_scale, Real new_scale)
    {
      Real const ratio = old_scale / new_scale;
      Real power = Real{1};
      for (std::size_t j = c.size(); j-- > 0;)
      {
        power *= ratio;
        c[j] *= power;
      }
    }

    std::size_t sites, order, momentum_index;
    Real center;
    Complex target_phase;
};

/// Follow the free odd-ring ground-sea polynomial towards -1 < Delta <= 0.
/// Newton uses M-1 coefficient equations and the exact momentum constraint;
/// acceptance checks ALL M equations, momentum, energy continuity/concavity,
/// and independent variational upper bounds on the sector minimum.
/// Every accepted Newton update, including rejected continuation stages,
/// consumes the one global budget. Failed outputs are NOT eigenstate energies.
template <uni20::Real Real>
OddPolynomialBranch<Real> continue_odd_polynomial(std::size_t sites, Real delta, uni20::half_int sz,
                                                  SolverOptions<Real> const& options = {})
{
  auto const m = sector_roots(sites, sz);
  OddPolynomialCoordinates<Real> const coordinates(sites, m);
  if (!uni20::isfinite(delta) || delta <= -Real{1} || delta > Real{0})
    throw std::invalid_argument("odd XXZ continuation requires -1 < Delta <= 0");
  if (!uni20::isfinite(options.residual_tolerance) || options.residual_tolerance <= Real{0})
    throw std::invalid_argument("residual tolerance must be finite and positive");
  auto const elements = std::size_t(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(Real);
  if (m && m > elements / m) throw std::length_error("odd XXZ Newton workspace is too large");
  Real const eps = uni20::numeric_limits<Real>::epsilon(), tolerance = options.residual_tolerance;
  OddPolynomialBranch<Real> state;
  OddSectorVariationalBounds<Real> const variational(sites, m);
  state.variational_upper_bound = variational.evaluate(delta).upper_bound();
  state.delta = delta;
  state.center = coordinates.center;
  state.momentum_index = coordinates.momentum_index;
  auto c = coordinates.seed(), saved = c;
  Real d = m <= 1 ? delta : Real{0}, saved_delta = d;
  Real s = coordinates.scale(d), saved_scale = s;
  Real saved_energy = PolynomialBetheSystem<Real>(sites, m, state.center, s).energy(c, d);
  Real preceding_delta = d, preceding_energy = saved_energy;
  bool have_secant = false;
  Real increment = Real{1} / Real{32};
  std::size_t stage_iterations = 0;
  for (;;)
  {
    auto failure = PolynomialContinuationStatus::stalled;
    PolynomialBetheSystem<Real> const system(sites, m, state.center, s);
    auto const weights = coordinates.weights(s);
    coordinates.impose_momentum(c, weights);
    bool numerical_failure = false;
    bool accepted_stage = false;
    std::vector<Real> correction(m, Real{0});
    try
    {
      uni20::DenseMatrix<Real> jac(m, m);
      auto const f = system.evaluate(c, d, &jac);
      state.momentum_error = std::abs(system.momentum_phase(c) - coordinates.target_phase);
      state.reciprocal_condition = Real{1};
      state.correction_norm = Real{0};
      if (m > 1)
      {
        // The last raw equation is not used to construct the Newton step,
        // but is retained in f.norm and may NEVER be omitted from acceptance.
        auto const r = m - 1;
        uni20::DenseMatrix<Real> matrix(r, r), rhs(r, r + 1);
        Real matrix_norm = Real{0};
        for (std::size_t i = 0; i < r; ++i)
        {
          Real row = Real{0};
          for (std::size_t j = 0; j < r; ++j)
          {
            matrix[i, j] = (jac[i, j] + jac[i, r] * weights[j]) / f.scale;
            if (!uni20::isfinite(matrix[i, j])) throw std::runtime_error("nonfinite reduced polynomial Jacobian");
            row += std::abs(matrix[i, j]);
            rhs[i, j + 1] = i == j ? Real{1} : Real{0};
          }
          matrix_norm = std::max(matrix_norm, row);
          rhs[i, 0] = -f.residual[i] / f.scale;
        }
        // The same factorization gives a condition estimate; no normal
        // equations or loss of long-double arithmetic is introduced.
        uni20::linalg::solve_inplace(matrix, rhs);
        Real inverse_norm = Real{0};
        for (std::size_t i = 0; i < r; ++i)
        {
          Real row = Real{0};
          correction[i] = rhs[i, 0];
          if (!uni20::isfinite(correction[i])) throw std::runtime_error("nonfinite polynomial Newton correction");
          for (std::size_t j = 0; j < r; ++j)
          {
            if (!uni20::isfinite(rhs[i, j + 1])) throw std::runtime_error("nonfinite inverse polynomial Jacobian");
            row += std::abs(rhs[i, j + 1]);
          }
          inverse_norm = std::max(inverse_norm, row);
        }
        coordinates.impose_momentum(correction, weights);
        state.reciprocal_condition = (Real{1} / matrix_norm) / inverse_norm;
        for (std::size_t i = 0; i < m; ++i)
        {
          if (!uni20::isfinite(correction[i])) throw std::runtime_error("nonfinite constrained Newton correction");
          state.correction_norm = std::max(state.correction_norm, std::abs(correction[i]) / (Real{1} + std::abs(c[i])));
        }
        if (!uni20::isfinite(state.reciprocal_condition) || state.reciprocal_condition <= Real{64} * eps ||
            !uni20::isfinite(state.correction_norm))
        {
          failure = PolynomialContinuationStatus::ill_conditioned;
          numerical_failure = true;
        }
      }
      if (!numerical_failure && f.norm <= tolerance && state.momentum_error <= tolerance &&
          state.correction_norm <= std::sqrt(tolerance))
      {
        Real const energy = system.energy(c, d);
        Real const upper_bound = variational.evaluate(d).upper_bound();
        Real const variational_roundoff =
            Real{256} * eps * (Real{1} + Real(sites) + std::abs(upper_bound) + std::abs(energy));
        bool const below_trial = uni20::isfinite(energy) && energy <= upper_bound + variational_roundoff;
        // ||dH/dDelta|| <= N/4 on a ring. Every continuously tracked
        // eigenvalue obeys this bound; it catches some, not all, branch jumps.
        Real const bound =
            Real(sites) * std::abs(d - saved_delta) / Real{4} + Real{256} * eps * (Real{1} + std::abs(saved_energy));
        // A sector's lowest energy is concave in Delta (minimum of affine
        // Rayleigh quotients). Its continuation to smaller Delta must lie
        // below the extrapolated preceding secant. The Lipschitz bound alone
        // can admit jumps to excited branches with exactly correct momentum.
        bool concave = true;
        if (have_secant && d != saved_delta)
        {
          Real const ratio = (d - saved_delta) / (saved_delta - preceding_delta);
          Real const upper = saved_energy + (saved_energy - preceding_energy) * ratio;
          Real const roundoff =
              Real{256} * eps * (Real{1} + std::abs(saved_energy) + std::abs(preceding_energy)) * (Real{1} + ratio);
          concave = energy <= upper + roundoff;
        }
        if (std::abs(energy - saved_energy) <= bound && concave && below_trial)
        {
          if (d != saved_delta)
          {
            preceding_delta = saved_delta;
            preceding_energy = saved_energy;
            have_secant = true;
          }
          saved = c;
          saved_delta = d;
          saved_scale = s;
          saved_energy = energy;
          accepted_stage = true;
        }
        else
        {
          if (!below_trial) ++state.variational_rejections;
          failure = PolynomialContinuationStatus::branch_rejected;
          numerical_failure = true;
        }
      }
      if (!numerical_failure && !accepted_stage && state.iterations < options.max_iterations && stage_iterations < 24)
      {
        bool accepted = false;
        Real damping = Real{1};
        for (int backtrack = 0; backtrack <= uni20::numeric_limits<Real>::digits; ++backtrack)
        {
          auto trial = c;
          bool finite = true;
          for (std::size_t i = 0; i < m; ++i)
          {
            trial[i] += damping * correction[i];
            finite = finite && uni20::isfinite(trial[i]);
          }
          if (finite)
          {
            coordinates.impose_momentum(trial, weights);
            finite = std::all_of(trial.begin(), trial.end(), [](Real x) { return uni20::isfinite(x); });
            try
            {
              if (finite)
              {
                auto const candidate = system.evaluate(trial, d);
                (void)system.momentum_phase(trial); // Reject observable poles during backtracking.
                if (candidate.norm < (Real{1} - damping / Real{10000}) * f.norm)
                {
                  c = std::move(trial);
                  ++state.iterations;
                  ++stage_iterations;
                  accepted = true;
                  break;
                }
              }
            }
            catch (std::runtime_error const&)
            {} // A shorter step can still be representable.
          }
          damping /= Real{2};
        }
        if (accepted) continue;
        numerical_failure = true;
        failure = PolynomialContinuationStatus::stalled;
      }
    }
    catch (std::runtime_error const&)
    {
      numerical_failure = true;
      failure = PolynomialContinuationStatus::ill_conditioned;
    }
    if (accepted_stage && d == delta)
    {
      state.equations_converged = true;
      state.status = PolynomialContinuationStatus::equations_converged;
      break;
    }
    if (state.iterations == options.max_iterations) break;
    if (!accepted_stage)
    {
      ++state.rejected_stages;
      increment = (saved_delta - d) / Real{2};
      if (d == saved_delta || increment <= Real{64} * eps)
      {
        state.status = failure;
        break;
      }
    }
    else if (stage_iterations < 8)
      increment = std::min(Real{1} / Real{16}, increment * Real{3} / Real{2});
    stage_iterations = 0;
    c = saved;
    d = std::max(delta, saved_delta - increment);
    s = coordinates.scale(d);
    coordinates.rescale(c, saved_scale, s);
  }
  state.root_delta = d;
  // Even on failure, reported observables and residual refer to the requested
  // Hamiltonian, in the returned coordinate system. root_delta records how
  // far the continuation attempt reached; it does not claim convergence there.
  state.coordinate_scale = s;
  state.coefficients = std::move(c);
  PolynomialBetheSystem<Real> const output(sites, m, state.center, s);
  state.energy = uni20::numeric_limits<Real>::quiet_NaN();
  state.residual_norm = state.momentum_error = uni20::numeric_limits<Real>::infinity();
  try
  {
    state.energy = output.energy(state.coefficients, delta);
    state.momentum_error = std::abs(output.momentum_phase(state.coefficients) - coordinates.target_phase);
    state.residual_norm = output.evaluate(state.coefficients, delta).norm;
  }
  catch (std::runtime_error const&)
  {
    state.equations_converged = false;
    state.status = PolynomialContinuationStatus::ill_conditioned;
  }
  return state;
}
} // namespace bethe::xxz::detail
