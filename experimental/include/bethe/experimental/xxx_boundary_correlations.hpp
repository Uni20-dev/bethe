// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/experimental/detail/xxx_open_gaudin.hpp>
#include <uni20/linalg/ops/linear_solve.hpp>
#include <uni20/tensor/mdspec_tensor_view.hpp>

#include <array>
#include <optional>

namespace bethe::experimental::xxx::open
{
namespace model = bethe::heisenberg::open;
/// Reflection eigenvalue of a supported finite-real-root highest-weight state.
/// r=(-1)^(N*M + sum I_j + M*(M-1)/2). SU(2) descendants inherit r.
/// Evaluate modulo two to avoid overflowing products or sums of labels.
[[nodiscard]] inline int real_root_reflection_parity(std::size_t sites, std::span<uni20::half_int const> numbers)
{
  model::detail::validate_quantum_numbers(sites, numbers);
  auto const m = numbers.size();
  bool odd = ((sites % 2) && (m % 2)) != (m % 4 >= 2);
  for (auto number : numbers)
    odd = odd != ((number.twice() / 2) % 2 != 0);
  return odd ? -1 : 1;
}

enum class BoundaryCorrelationStatus
{
  converged,
  roots_unconverged,
  linear_solve_failed,
  precision_limit
};

template <uni20::Real Real> struct BoundaryTwoSiteDensity
{
    /// Full <S_1.S_2>, NOT a restricted spectral estimate.
    std::optional<Real> bond_correlation;
    std::optional<Real> ground_residual;
    /// The actual linear-solve result, absent if roots were rejected first.
    std::optional<uni20::linalg::SolveInfo> solve_info;
    /// Normwise infinity backward error of G*response=b; not a forward error bar.
    std::optional<Real> solve_backward_error;
    BoundaryCorrelationStatus status = BoundaryCorrelationStatus::precision_limit;
    bool converged() const { return status == BoundaryCorrelationStatus::converged; }

    [[nodiscard]] std::optional<Real> zz() const
    {
      if (!converged() || !bond_correlation) return {};
      return *bond_correlation / Real{3};
    }
    [[nodiscard]] std::optional<Real> singlet_probability() const
    {
      if (!converged() || !bond_correlation) return {};
      return (Real{1} / Real{4}) - *bond_correlation;
    }
    /// Probability per triplet component, not the combined triplet probability.
    [[nodiscard]] std::optional<Real> triplet_component_probability() const
    {
      if (auto z = zz()) return (Real{1} / Real{4}) + *z;
      return {};
    }
    /// Row-major rho in |up up>, |up down>, |down up>, |down down> order.
    [[nodiscard]] std::optional<std::array<Real, 16>> matrix() const
    {
      auto z = zz();
      if (!z) return {};
      std::array<Real, 16> rho;
      rho.fill(uni20::scalar_like(*z, 0));
      rho[0] = rho[15] = (Real{1} / Real{4}) + *z;
      rho[5] = rho[10] = (Real{1} / Real{4}) - *z;
      rho[6] = rho[9] = Real{2} * *z;
      return rho;
    }
};

/// Ground-root-only boundary bond and complete SU(2)-invariant two-site density
/// matrix. Differentiate the integrable boundary-impurity coupling 1/(1+theta^2)
/// at theta=0 (Wang, PRB 56, 14045; cond-mat/9805253).
/// With u=z/2, b_i=2*u_i/(u_i^2+1/4)^2, solve G*x=b and return 1/4-b^T*x/2.
/// One native-precision square solve through normal Uni20 dispatch; no finite
/// differences, excitation scan, string truncation or boundary-block sum.
template <uni20::Real Real>
[[nodiscard]] BoundaryTwoSiteDensity<Real> boundary_two_site_density(std::size_t sites,
                                                                     model::RealState<Real> const& ground)
{
  using uni20::math::abs;
  boundary_detail::validate_singlet_shape(sites, ground);
  Real const zero = bethe::detail::precision_zero(ground.rapidities.front());
  for (auto const& root : ground.rapidities)
    bethe::detail::require_same_precision(root, zero);
  BoundaryTwoSiteDensity<Real> result;
  if (!ground.converged)
  {
    result.status = BoundaryCorrelationStatus::roots_unconverged;
    return result;
  }
  result.ground_residual = boundary_detail::measured_residual(sites, ground, zero);
  if (!result.ground_residual ||
      *result.ground_residual > 256 * uni20::scalar_like(zero, sites) * uni20::numeric_limits<Real>::epsilon(zero))
    return result;
  std::vector<Real> u;
  for (Real z : ground.rapidities)
    u.push_back(z / Real{2});
  auto const n = u.size();
  auto const matrix = boundary_detail::gaudin_matrix<Real>(sites, u, zero);
  std::vector<Real> rhs(n, zero), response(n, zero), workspace = matrix;
  for (std::size_t i = 0; i < n; ++i)
  {
    Real const denominator = u[i] * u[i] + (Real{1} / Real{4});
    rhs[i] = response[i] = Real{2} * u[i] / (denominator * denominator);
  }
  // Column-major owned workspace permits LAPACK/MPLAPACK dispatch without a
  // backend override. Keep the original row-major matrix for residual checking.
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = i + 1; j < n; ++j)
      std::swap(workspace[i * n + j], workspace[j * n + i]);
  using Span = stdex::mdspan<Real, stdex::dextents<std::size_t, 2>, stdex::layout_left>;
  using ConstSpan = stdex::mdspan<Real const, stdex::dextents<std::size_t, 2>, stdex::layout_left>;
  using View = uni20::MdspecTensorView<Span, ConstSpan, uni20::HostStorage>;
  View a(Span(workspace.data(), n, n), ConstSpan(workspace.data(), n, n));
  View b(Span(response.data(), n, 1), ConstSpan(response.data(), n, 1));
  uni20::linalg::SolveOptions<Real> options{.relative_pivot_tolerance =
                                                64 * uni20::numeric_limits<Real>::epsilon(zero)};
  result.solve_info = [&] {
#if UNI20_HAS_MPREAL
    if constexpr (uni20::has_runtime_precision_v<Real>)
      return uni20::linalg::solve_inplace_with_info(a, b, zero.precision(), options);
    else
#endif
      return uni20::linalg::solve_inplace_with_info(a, b, options);
  }();
  if (!result.solve_info->succeeded())
  {
    result.status = BoundaryCorrelationStatus::linear_solve_failed;
    return result;
  }
  Real matrix_norm = zero, rhs_norm = zero, response_norm = zero, residual_norm = zero;
  bethe::detail::CompensatedSum<Real> contraction(zero);
  for (std::size_t i = 0; i < n; ++i)
  {
    if (!uni20::isfinite(response[i])) return result;
    contraction.add(rhs[i] * response[i]);
    rhs_norm = std::max(rhs_norm, abs(rhs[i]));
    response_norm = std::max(response_norm, abs(response[i]));
    bethe::detail::CompensatedSum<Real> row_norm(zero), residual(zero);
    residual.add(-rhs[i]);
    for (std::size_t j = 0; j < n; ++j)
    {
      row_norm.add(abs(matrix[i * n + j]));
      residual.add(matrix[i * n + j] * response[j]);
    }
    matrix_norm = std::max(matrix_norm, row_norm.value());
    residual_norm = std::max(residual_norm, abs(residual.value()));
  }
  Real const scale = matrix_norm * response_norm + rhs_norm;
  if (!uni20::isfinite(scale) || scale <= 0 || !uni20::isfinite(residual_norm)) return result;
  result.solve_backward_error = residual_norm / scale;
  Real const gate = Real{4096} * Real(n) * uni20::numeric_limits<Real>::epsilon(zero);
  Real const bond = (Real{1} / Real{4}) - contraction.value() / Real{2};
  if (*result.solve_backward_error > gate || !uni20::isfinite(bond) || bond < (-Real{3} / Real{4}) - gate ||
      bond > (Real{1} / Real{4}) + gate)
    return result;
  result.bond_correlation = bond; // Raw result, never clamp density eigenvalues.
  result.status = BoundaryCorrelationStatus::converged;
  return result;
}
} // namespace bethe::experimental::xxx::open
