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
/// accounted for separately in PreparedRaisingFormFactors, not in this determinant.
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

/// Prepared zero-field singlet -> finite-real-root S=Sz=1 raising form factors
/// on an even periodic XXX chain. Caux-Hagemans-Maillet (2005), Eqs. 11-13,
/// rational limit, with all norm/product factors. No strings or descendants.
///
/// Input RealState coordinates are z=2*lambda_ABA. This owning context validates
/// the supplied ground roots once and stores standard ABA lambda, its Gaudin
/// norm and ground-only factors. Later mutation/destruction of the input cannot
/// invalidate it. Each evaluation validates its excited roots independently.
/// Concurrent const calls share no scratch or mutable cache. See docs/xxx-structure-factor.md.
template <uni20::Real Real> class PreparedRaisingFormFactors {
  public:
    PreparedRaisingFormFactors(std::size_t sites, RealState<Real> const& ground) : sites_(sites)
    {
      using std::log;
      detail::checked_sites(sites);
      auto const m = sites / 2;
      if (sites % 2 || ground.rapidities.size() != m || ground.sz != uni20::half_int{0} || ground.spin_reversed ||
          ground.quantum_numbers != sector_ground_quantum_numbers(sites, uni20::half_int{0}))
        throw std::invalid_argument("XXX preparation requires the even-ring singlet ground state");
      root_tolerance_ = Real{256} * Real(sites) * uni20::numeric_limits<Real>::epsilon();
      status_ = validate_roots(ground, ground_residual_);
      if (status_ != FormFactorStatus::converged) return;
      status_ = FormFactorStatus::precision_limit;
      aba_rapidities_.reserve(m);
      for (Real z : ground.rapidities)
        aba_rapidities_.push_back(z / Real{2});
      auto const norm = detail::gaudin_log_determinant<Real>(sites, aba_rapidities_);
      ground_pivot_margin_ = norm.pivot_margin;
      if (!norm.log_absolute) return;
      bethe::detail::CompensatedSum<Real> logarithm;
      logarithm.add(log(Real(sites)) - *norm.log_absolute);
      last_column_.reserve(m);
      for (Real x : aba_rapidities_)
      {
        Real const denominator = x * x + Real{0.25};
        last_column_.push_back(Real{1} / denominator);
        logarithm.add(log(denominator));
      }
      for (std::size_t a = 0; a < m; ++a)
        for (std::size_t b = 0; b < a; ++b)
        {
          Real const d = aba_rapidities_[a] - aba_rapidities_[b];
          logarithm.add(-log(Real{1} + d * d));
        }
      ground_log_prefactor_ = logarithm.value();
      if (uni20::isfinite(ground_log_prefactor_)) status_ = FormFactorStatus::converged;
    }

    std::size_t sites() const { return sites_; }
    /// Owned standard ABA lambda, not the solver's scaled coordinate z.
    std::span<Real const> aba_rapidities() const { return aba_rapidities_; }
    /// Measured max|F|/N at the supplied ground roots, absent if not evaluated.
    std::optional<Real> ground_residual() const { return ground_residual_; }
    Real root_tolerance() const { return root_tolerance_; }
    FormFactorStatus status() const
    {
      // A moved-from context may have lost its snapshot and must not evaluate it.
      if (status_ == FormFactorStatus::converged &&
          (aba_rapidities_.size() != sites_ / 2 || last_column_.size() != sites_ / 2))
        return FormFactorStatus::precision_limit;
      return status_;
    }
    bool ready() const { return status() == FormFactorStatus::converged; }

    [[nodiscard]] FormFactor<Real> operator()(RealState<Real> const& excited) const
    {
      using Complex = std::complex<Real>;
      using std::abs;
      using std::exp;
      using std::log;
      using std::sqrt;
      auto const m = sites_ / 2;
      if (excited.rapidities.size() != m - 1 || excited.quantum_numbers.size() != m - 1 ||
          excited.sz != uni20::half_int{1} || excited.spin_reversed)
        throw std::invalid_argument("XXX raising form factors require a real S=Sz=1 highest-weight state");
      detail::validate_quantum_numbers(sites_, excited.quantum_numbers);
      FormFactor<Real> result;
      result.status = status();
      if (result.status != FormFactorStatus::converged) return result;
      std::optional<Real> residual;
      result.status = validate_roots(excited, residual);
      if (result.status != FormFactorStatus::converged) return result;
      result.status = FormFactorStatus::precision_limit;
      std::vector<Real> lambda;
      lambda.reserve(m - 1);
      for (Real z : excited.rapidities)
        lambda.push_back(z / Real{2});
      auto const norm = detail::gaudin_log_determinant<Real>(sites_, lambda);
      result.pivot_margin = std::min(ground_pivot_margin_, norm.pivot_margin);
      if (!norm.log_absolute) return result;
      auto const& mu = aba_rapidities_;
      std::vector<Complex> h(m * m);
      bethe::detail::CompensatedSum<Real> logarithm;
      logarithm.add(ground_log_prefactor_);
      logarithm.add(-*norm.log_absolute);
      for (std::size_t b = 0; b + 1 < m; ++b)
      {
        Complex product{1};
        for (Real x : mu)
        {
          Real const d = x - lambda[b], scale = sqrt(Real{1} + d * d);
          product *= Complex(d, Real{-1}) / scale;
          // Restore the H column's scaling in log space.
          logarithm.add(Real{2} * log(scale));
        }
        Complex const phase =
            detail::integer_power(Complex(lambda[b], Real{0.5}) / Complex(lambda[b], Real{-0.5}), sites_);
        for (std::size_t a = 0; a < m; ++a)
        {
          Real const d = mu[a] - lambda[b];
          // Coincident roots require a limiting representation, not perturbation.
          if (abs(d) <= root_tolerance_ * (Real{1} + abs(mu[a]) + abs(lambda[b]))) return result;
          h[a * m + b] = (product / Complex(d, Real{-1}) - phase * std::conj(product) / Complex(d, Real{1})) / d;
        }
      }
      for (std::size_t a = 0; a < m; ++a)
        h[a * m + m - 1] = last_column_[a];
      auto const determinant = bethe::detail::log_absolute_determinant<Real>(std::move(h), m);
      result.pivot_margin = std::min(result.pivot_margin, determinant.pivot_margin);
      if (!determinant.log_absolute) return result;
      logarithm.add(Real{2} * *determinant.log_absolute);
      for (Real x : lambda)
        logarithm.add(-log(x * x + Real{0.25}));
      for (std::size_t a = 0; a < lambda.size(); ++a)
        for (std::size_t b = 0; b < a; ++b)
        {
          Real const d = lambda[a] - lambda[b];
          logarithm.add(-log(Real{1} + d * d));
        }
      Real const weight = exp(logarithm.value());
      // A local spin matrix element is bounded by one. Underflow is not a zero.
      if (!uni20::isfinite(weight) || weight <= Real{0} || weight > Real(sites_)) return result;
      result.weight = weight;
      result.log_weight = logarithm.value();
      result.status = FormFactorStatus::converged;
      return result;
    }

  private:
    FormFactorStatus validate_roots(RealState<Real> const& state, std::optional<Real>& residual) const
    {
      if (!state.converged) return FormFactorStatus::roots_unconverged;
      for (Real z : state.rapidities)
        if (!uni20::isfinite(z)) return FormFactorStatus::precision_limit;
      std::vector<Real> angles;
      Real const measured = detail::residual<Real>(sites_, state.quantum_numbers, state.rapidities, angles);
      if (!uni20::isfinite(measured)) return FormFactorStatus::precision_limit;
      residual = measured;
      // Neither the caller's convergence flag nor its stored residual is a certificate.
      return measured <= root_tolerance_ ? FormFactorStatus::converged : FormFactorStatus::precision_limit;
    }

    std::size_t sites_;
    std::vector<Real> aba_rapidities_, last_column_;
    std::optional<Real> ground_residual_;
    Real root_tolerance_{}, ground_log_prefactor_{}, ground_pivot_margin_{};
    FormFactorStatus status_ = FormFactorStatus::precision_limit;
};
} // namespace bethe::heisenberg
