// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/open_string_solver.hpp>
#include <complex>

namespace bethe::xxz::quantum_group::detail
{
/// Regular (non-internal) direct and reflected scattering from Eq. (5.12).
/// Log phases are principal branches; callers choose continuous fused phases.
template <uni20::Real Real> struct OpenScattering
{
    std::complex<Real> value, first, second;
};
template <uni20::Real Real> OpenScattering<Real> open_scattering(std::complex<Real> u, std::complex<Real> v, Real eta)
{
  using C = std::complex<Real>;
  auto term = [](C w) { return std::pair{std::log(std::sinh(w)), std::cosh(w) / std::sinh(w)}; };
  auto const [ap, dap] = term(u - v + eta);
  auto const [am, dam] = term(u - v - eta);
  auto const [bp, dbp] = term(u + v + eta);
  auto const [bm, dbm] = term(u + v - eta);
  return {ap - am + bp - bm, dap - dam + dbp - dbm, -dap + dam + dbp - dbm};
}
} // namespace bethe::xxz::quantum_group::detail
