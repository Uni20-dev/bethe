// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/detail/log_determinant.hpp>
#include <bethe/heisenberg.hpp>

namespace bethe::heisenberg
{
enum class FormFactorStatus
{
  converged,
  roots_unconverged,
  precision_limit
};

template <uni20::Real Real> struct FormFactor
{
    /// N |<excited|S_0^+|ground>|^2, i.e. the normalized Fourier weight
    /// at the momentum transfer. Not multiplied by 2*pi.
    std::optional<Real> weight;
    std::optional<Real> log_weight;
    Real pivot_margin = Real{0};
    FormFactorStatus status = FormFactorStatus::precision_limit;
    bool converged() const { return status == FormFactorStatus::converged; }
};

namespace detail
{
/// Standard ABA lambda=z/2; Phi is the Jacobian of the unscaled logarithmic
/// equations with respect to lambda. Prefactors of the full ABA norm are
/// accounted for separately in raising_form_factor, not in this determinant.
template <uni20::Real Real>
bethe::detail::LogDeterminant<Real> gaudin_log_determinant(std::size_t sites, std::span<Real const> lambda)
{
  auto const m = lambda.size();
  std::vector<Real> matrix(m * m);
  for (std::size_t a = 0; a < m; ++a)
  {
    bethe::detail::CompensatedSum<Real> diagonal;
    diagonal.add(Real(sites) / (lambda[a] * lambda[a] + Real{0.25}));
    for (std::size_t b = 0; b < m; ++b)
      if (a != b)
      {
        Real const d = lambda[a] - lambda[b];
        matrix[a * m + b] = Real{2} / (Real{1} + d * d);
        diagonal.add(-matrix[a * m + b]);
      }
    matrix[a * m + a] = diagonal.value();
  }
  return bethe::detail::log_absolute_determinant<Real>(std::move(matrix), m);
}

template <typename Scalar> Scalar integer_power(Scalar base, std::size_t exponent)
{
  Scalar result{1};
  while (exponent)
  {
    if (exponent & 1U) result *= base;
    exponent >>= 1;
    if (exponent) base *= base;
  }
  return result;
}
} // namespace detail

/// Zero-field singlet -> finite-real-root S=Sz=1 state of an even periodic XXX
/// chain. Uses the rational limit of Caux-Hagemans-Maillet (2005), Eqs. 11-13,
/// with all norm/product prefactors. See docs/xxx-structure-factor.md.
/// Input states must come from this chain's solver; no strings or descendants.
template <uni20::Real Real>
[[nodiscard]] FormFactor<Real> raising_form_factor(std::size_t sites, RealState<Real> const& ground,
                                                   RealState<Real> const& excited)
{
  using Complex = std::complex<Real>;
  using std::abs;
  using std::exp;
  using std::log;
  using std::sqrt;
  detail::checked_sites(sites);
  auto const m = sites / 2;
  if (sites % 2 || ground.rapidities.size() != m || excited.rapidities.size() != m - 1 ||
      ground.quantum_numbers != sector_ground_quantum_numbers(sites, uni20::half_int{0}) ||
      excited.quantum_numbers.size() != m - 1 || ground.spin_reversed || excited.spin_reversed)
    throw std::invalid_argument("XXX form factors require the even-ring singlet and a real S=1 highest-weight state");
  detail::validate_quantum_numbers(sites, excited.quantum_numbers);
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
    Real const residual = detail::residual<Real>(sites, state->quantum_numbers, state->rapidities, angles);
    if (!uni20::isfinite(residual) || residual > root_gate) return result;
  }
  std::vector<Real> mu(m), lambda(m - 1);
  for (std::size_t j = 0; j < m; ++j)
    mu[j] = ground.rapidities[j] / Real{2};
  for (std::size_t j = 0; j + 1 < m; ++j)
    lambda[j] = excited.rapidities[j] / Real{2};
  auto const norm_ground = detail::gaudin_log_determinant<Real>(sites, mu);
  auto const norm_excited = detail::gaudin_log_determinant<Real>(sites, lambda);
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
    Complex const phase = detail::integer_power(Complex(lambda[b], Real{0.5}) / Complex(lambda[b], Real{-0.5}), sites);
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
} // namespace bethe::heisenberg
