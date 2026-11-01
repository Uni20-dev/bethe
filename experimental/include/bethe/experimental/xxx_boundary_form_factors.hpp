// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/experimental/detail/xxx_open_gaudin.hpp>
#include <bethe/xxx_form_factors.hpp>
#include <memory>

namespace bethe::experimental::xxx::open
{
namespace model = bethe::heisenberg::open;
using BoundaryFormFactorStatus = bethe::heisenberg::FormFactorStatus;

template <uni20::Real R> struct BoundaryFormFactor
{
    /// Local |<S=Sz=1|S_1^+|singlet>|^2, without N or 2*pi. The zz weight is half.
    std::optional<R> weight, log_weight;
    /// Real scaled-determinant sign, not by itself the physical amplitude phase.
    int determinant_sign = 0;
    /// Pivot rejection diagnostic, not a condition number or error bound.
    std::optional<R> pivot_margin;
    BoundaryFormFactorStatus status = BoundaryFormFactorStatus::precision_limit;
    bool converged() const { return status == BoundaryFormFactorStatus::converged; }
};

namespace boundary_detail
{
template <uni20::Real R> R pair_log(R const& u, R const& v)
{
  R const one = uni20::scalar_like(u, 1), d = u - v, s = u + v;
  return uni20::math::log(one + d * d) + uni20::math::log(one + s * s);
}
} // namespace boundary_detail

/// Immutable even-chain free-end singlet normalization. Copies share its owned
/// snapshot; moves invalidate the source. Inputs use z=2*lambda_ABA. Only positive
/// finite real roots, zero field and singlet -> highest triplet are supported.
/// Uses the boundary scalar product at the exact off-shell root -i/2; no finite
/// approximation to an infinite descendant root or Fourier normalization.
template <uni20::Real R> class BoundaryFormFactors {
    struct Snapshot
    {
        explicit Snapshot(R const& value) : zero(value), prefactor(value), energy(value) {}
        R zero, prefactor, energy;
        std::vector<R> roots;
        std::optional<R> residual, pivot;
        BoundaryFormFactorStatus status = BoundaryFormFactorStatus::precision_limit;
    };
    std::size_t sites_;
    std::shared_ptr<Snapshot const> snapshot_;

  public:
    BoundaryFormFactors(std::size_t sites, model::RealState<R> const& ground) : sites_(sites)
    {
      boundary_detail::validate_singlet_shape(sites, ground);
      R const zero = bethe::detail::precision_zero(ground.rapidities.front());
      for (auto const& z : ground.rapidities)
        bethe::detail::require_same_precision(z, zero);
      auto saved = std::make_shared<Snapshot>(zero);
      snapshot_ = saved;
      if (!ground.converged)
      {
        saved->status = BoundaryFormFactorStatus::roots_unconverged;
        return;
      }
      saved->residual = boundary_detail::measured_residual(sites, ground, zero);
      R const one = uni20::scalar_like(zero, 1), n = uni20::scalar_like(zero, sites);
      if (!saved->residual || *saved->residual > 256 * n * uni20::numeric_limits<R>::epsilon(zero)) return;
      for (auto const& z : ground.rapidities)
        saved->roots.push_back(z / 2);
      auto norm = bethe::detail::log_determinant<R>(boundary_detail::gaudin_matrix<R>(sites, saved->roots, zero),
                                                    saved->roots.size(), zero);
      saved->pivot = norm.pivot_margin;
      if (!norm.log_absolute) return;
      bethe::detail::CompensatedSum<R> prefactor(zero), energy(zero);
      prefactor.add(-*norm.log_absolute);
      energy.add((n - one) / 4);
      for (std::size_t i = 0; i < saved->roots.size(); ++i)
      {
        R const u = saved->roots[i], denominator = u * u + one / 4;
        prefactor.add(uni20::math::log(4 * u * u) + 2 * uni20::math::log(denominator));
        energy.add(-one / (2 * denominator));
        for (std::size_t j = 0; j < i; ++j)
          prefactor.add(-boundary_detail::pair_log(u, saved->roots[j]));
      }
      saved->prefactor = prefactor.value();
      saved->energy = energy.value();
      if (uni20::isfinite(saved->prefactor) && uni20::isfinite(saved->energy))
        saved->status = BoundaryFormFactorStatus::converged;
    }
    std::size_t sites() const { return sites_; }
    BoundaryFormFactorStatus status() const
    {
      return snapshot_ ? snapshot_->status : BoundaryFormFactorStatus::precision_limit;
    }
    bool ready() const { return status() == BoundaryFormFactorStatus::converged; }
    std::optional<R> ground_residual() const { return snapshot_ ? snapshot_->residual : std::optional<R>{}; }
    std::optional<R> ground_energy() const
    {
      return ready() ? std::optional<R>(snapshot_->energy) : std::optional<R>{};
    }
    std::optional<R> working_zero() const { return ready() ? std::optional<R>(snapshot_->zero) : std::optional<R>{}; }
    std::span<R const> aba_rapidities() const
    {
      return ready() ? std::span<R const>(snapshot_->roots) : std::span<R const>{};
    }

    [[nodiscard]] BoundaryFormFactor<R> operator()(model::RealState<R> const& excited) const
    {
      using C = uni20::complex<R>;
      using uni20::math::abs;
      using uni20::math::log;
      auto const m = sites_ / 2;
      if (excited.rapidities.size() != m - 1 || excited.quantum_numbers.size() != m - 1 ||
          excited.sz != uni20::half_int{1} || excited.spin_reversed)
        throw std::invalid_argument("boundary form factors require a finite-real-root S=Sz=1 state");
      model::detail::validate_quantum_numbers(sites_, excited.quantum_numbers);
      BoundaryFormFactor<R> result;
      result.status = status();
      if (snapshot_) result.pivot_margin = snapshot_->pivot;
      if (!excited.converged) result.status = BoundaryFormFactorStatus::roots_unconverged;
      if (result.status != BoundaryFormFactorStatus::converged) return result;
      result.status = BoundaryFormFactorStatus::precision_limit;
      auto const& saved = *snapshot_;
      R const zero = saved.zero, one = uni20::scalar_like(zero, 1), n = uni20::scalar_like(zero, sites_);
      R const epsilon = uni20::numeric_limits<R>::epsilon(zero), pole_gate = 256 * n * epsilon;
      auto const residual = boundary_detail::measured_residual(sites_, excited, zero);
      if (!residual || *residual > pole_gate) return result;
      std::vector<R> roots;
      for (auto const& z : excited.rapidities)
        roots.push_back(z / 2);
      auto norm =
          bethe::detail::log_determinant<R>(boundary_detail::gaudin_matrix<R>(sites_, roots, zero), roots.size(), zero);
      result.pivot_margin = std::min(*result.pivot_margin, norm.pivot_margin);
      if (!norm.log_absolute) return result;
      bethe::detail::CompensatedSum<R> logarithm(zero);
      logarithm.add(saved.prefactor - *norm.log_absolute);
      std::vector<R> matrix(m * m, zero);
      for (std::size_t b = 0; b < roots.size(); ++b)
      {
        R const v = roots[b];
        C product(one, zero);
        for (auto const& u : saved.roots)
        {
          C const shifted = C(v, one) * C(v, one) - u * u;
          product *= shifted / abs(shifted);
          logarithm.add(boundary_detail::pair_log(v, u));
        }
        C const half_phase = bethe::heisenberg::detail::integer_power(C(v, one / 2) / C(v, -one / 2), sites_);
        product *= uni20::conj(half_phase);
        for (std::size_t a = 0; a < m; ++a)
        {
          R const u = saved.roots[a];
          // Coincident roots need another representation: do not perturb them.
          if (abs(v - u) <= pole_gate * (one + abs(v) + abs(u))) return result;
          C const shifted = C(v, one) * C(v, one) - u * u;
          matrix[a * m + b] = 2 * (product / shifted).imag() / ((v - u) * (v + u));
        }
        logarithm.add(-2 * log(v * v + one / 4));
        for (std::size_t j = 0; j < b; ++j)
          logarithm.add(-boundary_detail::pair_log(v, roots[j]));
      }
      for (std::size_t a = 0; a < m; ++a)
      {
        R const denominator = saved.roots[a] * saved.roots[a] + one / 4;
        matrix[a * m + m - 1] = one / (denominator * denominator);
      }
      auto determinant = bethe::detail::log_determinant<R>(std::move(matrix), m, zero);
      result.pivot_margin = std::min(*result.pivot_margin, determinant.pivot_margin);
      if (!determinant.log_absolute || !determinant.phase) return result;
      logarithm.add(2 * *determinant.log_absolute);
      R const weight = uni20::math::exp(logarithm.value());
      R const norm_gate = 4096 * uni20::scalar_like(zero, m) * epsilon;
      // Local singlet raising norm is 1/2. Rounding allowance is not an error bound.
      if (!uni20::isfinite(weight) || weight <= zero || weight > one / 2 + norm_gate) return result;
      result.weight = weight;
      result.log_weight = logarithm.value();
      result.determinant_sign = *determinant.phase > zero ? 1 : -1;
      result.status = BoundaryFormFactorStatus::converged;
      return result;
    }
};

template <uni20::Real R>
[[nodiscard]] BoundaryFormFactor<R> boundary_raising_form_factor(std::size_t sites, model::RealState<R> const& ground,
                                                                 model::RealState<R> const& excited)
{
  return BoundaryFormFactors<R>(sites, ground)(excited);
}
} // namespace bethe::experimental::xxx::open
