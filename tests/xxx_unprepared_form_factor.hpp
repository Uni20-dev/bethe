// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/xxx_form_factors.hpp>

// Unprepared public implementation at fe81a93, retained only as a regression
// and timing control for the prepared-context refactor. It deliberately
// recomputes all ground factors and retains the original summation order.
// The unchanged Gaudin/log-determinant kernels are shared; coordinate-wave
// and bit-basis ED tests provide the separate physical normalization oracles.
namespace test_support
{
using bethe::heisenberg::FormFactor;
using bethe::heisenberg::FormFactorStatus;
using bethe::heisenberg::RealState;
using bethe::heisenberg::sector_ground_quantum_numbers;

template <uni20::Real Real>
[[nodiscard]] FormFactor<Real> unprepared_raising_form_factor(std::size_t sites, RealState<Real> const& ground,
                                                              RealState<Real> const& excited)
{
  using Complex = std::complex<Real>;
  using std::abs;
  using std::exp;
  using std::log;
  using std::sqrt;
  bethe::heisenberg::detail::checked_sites(sites);
  auto const m = sites / 2;
  if (sites % 2 || ground.rapidities.size() != m || excited.rapidities.size() != m - 1 ||
      ground.quantum_numbers != sector_ground_quantum_numbers(sites, uni20::half_int{0}) ||
      excited.quantum_numbers.size() != m - 1 || ground.spin_reversed || excited.spin_reversed)
    throw std::invalid_argument("XXX form factors require the even-ring singlet and a real S=1 highest-weight state");
  bethe::heisenberg::detail::validate_quantum_numbers(sites, excited.quantum_numbers);
  FormFactor<Real> result;
  if (!ground.converged || !excited.converged)
  {
    result.status = FormFactorStatus::roots_unconverged;
    return result;
  }
  // Re-evaluate the equations: a loose caller tolerance (or stale state flag)
  // must not turn inaccurate roots into apparently converged form factors.
  Real const root_gate = Real{256} * Real(sites) * uni20::numeric_limits<Real>::epsilon();
  for (auto const* state : {&ground, &excited})
  {
    std::vector<Real> angles;
    for (Real z : state->rapidities)
      if (!uni20::isfinite(z)) return result;
    Real const residual =
        bethe::heisenberg::detail::residual<Real>(sites, state->quantum_numbers, state->rapidities, angles);
    if (!uni20::isfinite(residual) || residual > root_gate) return result;
  }
  std::vector<Real> mu(m), lambda(m - 1);
  for (std::size_t j = 0; j < m; ++j)
    mu[j] = ground.rapidities[j] / Real{2};
  for (std::size_t j = 0; j + 1 < m; ++j)
    lambda[j] = excited.rapidities[j] / Real{2};
  auto const norm_ground = bethe::heisenberg::detail::gaudin_log_determinant<Real>(sites, mu);
  auto const norm_excited = bethe::heisenberg::detail::gaudin_log_determinant<Real>(sites, lambda);
  result.pivot_margin = std::min(norm_ground.pivot_margin, norm_excited.pivot_margin);
  if (!norm_ground.log_absolute || !norm_excited.log_absolute) return result;
  std::vector<Complex> h(m * m);
  bethe::detail::CompensatedSum<Real> logarithm;
  logarithm.add(log(Real(sites)) - *norm_ground.log_absolute - *norm_excited.log_absolute);
  for (std::size_t b = 0; b + 1 < m; ++b)
  {
    Complex product{1};
    for (Real x : mu)
    {
      Real const d = x - lambda[b], scale = sqrt(Real{1} + d * d);
      product *= Complex(d, Real{-1}) / scale;
      // H column divided by prod sqrt(1+d^2); restore |det H|^2 in log space.
      logarithm.add(Real{2} * log(scale));
    }
    Complex const phase =
        bethe::heisenberg::detail::integer_power(Complex(lambda[b], Real{0.5}) / Complex(lambda[b], Real{-0.5}), sites);
    for (std::size_t a = 0; a < m; ++a)
    {
      Real const d = mu[a] - lambda[b];
      // The coincident-root limit needs a different determinant representation.
      // Do not perturb roots or silently divide through an unresolved pole.
      if (abs(d) <= root_gate * (Real{1} + abs(mu[a]) + abs(lambda[b]))) return result;
      h[a * m + b] = (product / Complex(d, Real{-1}) - phase * std::conj(product) / Complex(d, Real{1})) / d;
    }
  }
  for (std::size_t a = 0; a < m; ++a)
    h[a * m + m - 1] = Real{1} / (mu[a] * mu[a] + Real{0.25});
  auto const determinant = bethe::detail::log_absolute_determinant<Real>(std::move(h), m);
  result.pivot_margin = std::min(result.pivot_margin, determinant.pivot_margin);
  if (!determinant.log_absolute) return result;
  logarithm.add(Real{2} * *determinant.log_absolute);
  for (Real x : mu)
    logarithm.add(log(x * x + Real{0.25}));
  for (Real x : lambda)
    logarithm.add(-log(x * x + Real{0.25}));
  for (auto const* roots : {&mu, &lambda})
    for (std::size_t a = 0; a < roots->size(); ++a)
      for (std::size_t b = 0; b < a; ++b)
      {
        Real const d = (*roots)[a] - (*roots)[b];
        logarithm.add(-log(Real{1} + d * d));
      }
  Real const weight = exp(logarithm.value());
  // A local spin matrix element is bounded by one. Underflow is not an exact zero.
  if (!uni20::isfinite(weight) || weight <= Real{0} || weight > Real(sites)) return result;
  result.weight = weight;
  result.log_weight = logarithm.value();
  result.status = FormFactorStatus::converged;
  return result;
}
} // namespace test_support
