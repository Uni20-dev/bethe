// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/xxx_form_factors.hpp>
#include <optional>
#include <span>
#include <utility>

namespace bethe::experimental::xxx
{
/// Immutable prepared even-ring singlet for finite-root form factors, J=1.
/// Input RealState roots use z=2*lambda_ABA. The existing public preparation
/// validates those roots and owns the ground norm/products; no second LU is
/// performed. Energy is recomputed from those roots, not copied from metadata.
/// This is an in-memory snapshot, not a reader for historical research caches.
template <uni20::Real Real> class PreparedSinglet {
  public:
    PreparedSinglet(std::size_t sites, heisenberg::RealState<Real> const& ground) : factors_(sites, ground)
    {
      if (!factors_.ready()) return;
      Real const zero = bethe::detail::precision_zero(factors_.aba_rapidities().front());
      Real const one = uni20::scalar_like(zero, 1);
      bethe::detail::CompensatedSum<Real> sum(zero);
      sum.add(uni20::scalar_like(zero, sites) / (4 * one));
      for (Real const& root : factors_.aba_rapidities())
        sum.add(-one / (2 * (root * root + one / 4)));
      Real const value = sum.value();
      if (uni20::isfinite(value)) energy_ = value;
    }

    PreparedSinglet(PreparedSinglet const&) = default;
    PreparedSinglet& operator=(PreparedSinglet const& other)
    {
      if (this == &other) return *this;
      PreparedSinglet snapshot(other);
      return *this = std::move(snapshot);
    }
    PreparedSinglet(PreparedSinglet&&) = default;
    PreparedSinglet& operator=(PreparedSinglet&& other)
    {
      if (this == &other) return *this;
      // Do not expose the old energy if a subsequent scalar transfer throws.
      energy_.reset();
      factors_ = std::move(other.factors_);
      energy_ = std::move(other.energy_);
      return *this;
    }

    std::size_t sites() const { return factors_.sites(); }
    bool ready() const { return factors_.ready() && energy_.has_value(); }
    heisenberg::FormFactorStatus status() const
    {
      if (!factors_.ready()) return factors_.status();
      return energy_ ? heisenberg::FormFactorStatus::converged : heisenberg::FormFactorStatus::precision_limit;
    }

    /// Standard ABA lambda. Failed and moved-from snapshots expose no usable roots.
    std::span<Real const> aba_rapidities() const
    {
      if (!ready()) return {};
      return factors_.aba_rapidities();
    }
    std::optional<Real> energy() const
    {
      if (!ready()) return {};
      return energy_;
    }
    /// Finite-precision exemplar for independent real/complex scratch allocation.
    std::optional<Real> working_zero() const
    {
      if (!ready()) return {};
      return bethe::detail::precision_zero(factors_.aba_rapidities().front());
    }
    std::optional<Real> ground_log_prefactor() const
    {
      if (!ready()) return {};
      return factors_.ground_log_prefactor();
    }
    std::optional<Real> ground_pivot_margin() const
    {
      if (!ready()) return {};
      return factors_.ground_pivot_margin();
    }
    /// Measured residual diagnoses failed preparation; absent after moving the snapshot.
    std::optional<Real> ground_residual() const { return factors_.ground_residual(); }

    /// Existing supported real-root triplets, using the same prepared norm.
    /// Complex-root family evaluators must implement and validate their own matrices.
    heisenberg::FormFactor<Real> real_raising_form_factor(heisenberg::RealState<Real> const& excited) const
    {
      if (ready()) return factors_(excited);
      heisenberg::FormFactor<Real> result;
      result.status = status();
      return result;
    }

  private:
    heisenberg::PreparedRaisingFormFactors<Real> factors_;
    std::optional<Real> energy_;
};
} // namespace bethe::experimental::xxx
