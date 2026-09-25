// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/continuum_newton.hpp>
#include <bethe/detail/elliptic_theta.hpp>
#include <optional>
#include <span>
#include <uni20/common/half_int.hpp>

namespace bethe::xyz
{
template <uni20::Real Real> struct Couplings
{
    Real x{}, y{}, z{};
};

/// Rectangular Hermitian parametrization: theta_j(pi*eta | i*t)/theta_j(0 | i*t).
/// H=sum(Jx*Sx*Sx+Jy*Sy*Sy+Jz*Sz*Sz), S=sigma/2. This is only the
/// coupling parametrization; it does not select roots or compute an eigenstate.
template <uni20::Real Real> Couplings<Real> couplings(Real eta, Real t)
{
  if (!uni20::isfinite(eta)) throw std::invalid_argument("XYZ eta must be finite and real");
  auto ratio = [&](unsigned kind) {
    auto const numerator = bethe::detail::elliptic_theta(kind, std::complex<Real>(eta, 0), t);
    auto const denominator = bethe::detail::elliptic_theta(kind, std::complex<Real>{}, t);
    return bethe::detail::theta_ratio(numerator, 0, denominator, 0).real();
  };
  return {ratio(4), ratio(3), ratio(2)};
}
/// Evaluate the regular-root energy expression for even N and M=N/2.
/// This does NOT solve/validate Bethe equations or certify a physical state.
/// A complex result is retained; poles are rejected rather than regularized.
template <uni20::Real Real>
std::complex<Real> candidate_energy(std::size_t sites, std::span<std::complex<Real> const> roots, Real eta, Real t)
{
  using Complex = std::complex<Real>;
  if (sites < 2 || sites % 2 || roots.size() != sites / 2 || !uni20::isfinite(eta) || !(eta > Real{0} && eta < Real{1}))
    throw std::invalid_argument("XYZ regular energy requires even N>=2, N/2 roots, and real 0<eta<1");
  auto const at_eta = bethe::detail::elliptic_theta(1, Complex(eta, 0), t);
  auto const at_zero = bethe::detail::elliptic_theta(1, Complex{}, t);
  auto const factor = bethe::detail::theta_ratio(at_eta, 0, at_zero, 1);
  auto g = [&](Complex u) {
    auto const jet = bethe::detail::elliptic_theta(1, u, t);
    return factor * bethe::detail::theta_ratio(jet, 1, jet, 0);
  };
  bethe::detail::CompensatedSum<Complex> energy;
  energy.add(Real(sites) * g(Complex(eta, 0)) / Real{4});
  for (auto root : roots)
    energy.add((g(root - eta / Real{2}) - g(root + eta / Real{2})) / Real{2});
  if (!bethe::detail::theta_detail::finite(energy.value()))
    throw std::overflow_error("XYZ candidate energy exceeds range");
  return energy.value();
}
enum class Status
{
  converged,
  iteration_limit,
  stalled,
  precision_limit
};
template <uni20::Real Real> using SolverOptions = bethe::SolverOptions<Real>;
template <uni20::Real Real> struct GroundState
{
    std::size_t sites = 0, iterations = 0;
    Real eta{}, t{}, residual_norm = uni20::numeric_limits<Real>::infinity(), momentum{};
    Couplings<Real> exchange;
    std::vector<uni20::half_int> quantum_numbers;
    std::vector<std::complex<Real>> roots;
    std::optional<Real> energy;
    bool converged = false;
    Status status = Status::iteration_limit;
};
namespace detail
{
// Symmetric imaginary roots lambda=i*x, phase parameter xi=0 and sum(lambda)=0.
template <uni20::Real Real> struct GroundSystem
{
    std::size_t sites;
    Real eta, t;
    Real const pi = Real{4} * std::atan(Real{1});
    std::size_t count() const { return sites / 2; }
    std::size_t order() const { return count() / 2; }
    std::vector<Real> expand(std::span<Real const> x) const
    {
      std::vector<Real> roots;
      for (auto j = x.size(); j-- > 0;)
        roots.push_back(-x[j]);
      if (count() % 2) roots.push_back(Real{0});
      roots.insert(roots.end(), x.begin(), x.end());
      return roots;
    }
    bool physical(std::span<Real const> x) const
    {
      for (std::size_t j = 0; j < x.size(); ++j)
        if (!uni20::isfinite(x[j]) || !(x[j] > Real{0} && x[j] < t / Real{2}) || (j && !(x[j] > x[j - 1])))
          return false;
      return true;
    }
    struct Phase
    {
        Real value{}, derivative{};
    };
    Phase phase(Real a, Real x) const
    {
      if (a == Real{0.5}) return {};
      auto const jet = bethe::detail::elliptic_theta(1, std::complex<Real>(a, x), t);
      auto const value = jet.derivative[0];
      return {Real{2} * std::atan2(value.imag(), value.real()),
              Real{2} * bethe::detail::theta_ratio(jet, 1, jet, 0).real()};
    }
    struct Evaluation
    {
        std::vector<Real> residual;
        Real norm{};
    };
    Evaluation evaluate(std::span<Real const> x, uni20::DenseMatrix<Real>* jacobian = nullptr) const
    {
      auto const roots = expand(x);
      Evaluation out{std::vector<Real>(order())};
      if (jacobian)
        for (std::size_t i = 0; i < order(); ++i)
          for (std::size_t j = 0; j < order(); ++j)
            (*jacobian)[i, j] = Real{0};
      for (std::size_t row = 0; row < order(); ++row)
      {
        auto const full = order() + count() % 2 + row;
        auto const drive = phase(eta / Real{2}, x[row]);
        bethe::detail::CompensatedSum<Real> sum, diagonal;
        sum.add(Real(sites) * drive.value);
        sum.add(-pi * Real(2 * row + 1 + count() % 2));
        diagonal.add(Real(sites) * drive.derivative);
        for (std::size_t j = 0; j < roots.size(); ++j)
          if (j != full)
          {
            auto const scattering = phase(eta, x[row] - roots[j]);
            sum.add(-scattering.value);
            diagonal.add(-scattering.derivative);
            if (jacobian)
            {
              if (j < order())
                (*jacobian)[row, order() - 1 - j] -= scattering.derivative;
              else if (j >= order() + count() % 2)
                (*jacobian)[row, j - order() - count() % 2] += scattering.derivative;
            }
          }
        if (jacobian) (*jacobian)[row, row] += diagonal.value();
        if (jacobian)
          for (std::size_t col = 0; col < order(); ++col)
            if (!uni20::isfinite((*jacobian)[row, col])) throw std::overflow_error("nonfinite XYZ Jacobian");
        out.residual[row] = sum.value();
        if (!uni20::isfinite(sum.value())) throw std::overflow_error("nonfinite XYZ residual");
        // The equations flatten as eta approaches one. Do not accept an
        // unresolved seed merely because the unscaled phases become small.
        out.norm = std::max(out.norm, std::abs(sum.value()) / (Real(sites) * (Real{1} - eta)));
      }
      return out;
    }
};
} // namespace detail

/// Even periodic ground branch, real eta in (0,1), rectangular tau=i*t.
/// Failed iterations never publish an energy; retained roots are diagnostics.
template <uni20::Real Real>
GroundState<Real> ground_state(std::size_t sites, Real eta, Real t, SolverOptions<Real> options = {})
{
  if (sites < 2 || sites % 2 || sites > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max() / 4) ||
      !uni20::isfinite(eta) || !(eta > Real{0} && eta < Real{1}) || !uni20::isfinite(t) || !(t > Real{0}) ||
      !uni20::isfinite(options.residual_tolerance) || !(options.residual_tolerance > Real{0}))
    throw std::invalid_argument("XYZ ground branch requires even N>=2, 0<eta<1, t>0 and finite positive tolerance");
  auto const m = sites / 2, order = m / 2;
  auto const elements = static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(Real);
  if (order && order > elements / order) throw std::length_error("XYZ Newton matrix too large");
  GroundState<Real> out;
  out.sites = sites;
  out.eta = eta;
  out.t = t;
  Real const pi = Real{4} * std::atan(Real{1});
  out.momentum = m % 2 ? pi : Real{0};
  for (std::size_t j = 0; j < m; ++j)
    out.quantum_numbers.push_back(
        uni20::from_twice(2 * static_cast<std::int64_t>(j) - static_cast<std::int64_t>(m - 1)));
  detail::GroundSystem<Real> system{sites, eta, t};
  std::vector<Real> x(order);
  for (std::size_t j = 0; j < order; ++j)
    x[j] = std::min(t / Real{2}, eta) * Real(2 * j + 1 + m % 2) / Real(m);
  auto sync = [&] {
    out.roots.clear();
    for (Real v : system.expand(x))
      out.roots.emplace_back(0, v);
  };
  try
  {
    out.exchange = couplings(eta, t);
    if (!system.physical(x))
    {
      out.status = Status::precision_limit;
      sync();
      return out;
    }
    out.iterations = bethe::detail::continuum_newton(system, x, options);
    sync();
    out.residual_norm = system.evaluate(x).norm;
    if (out.residual_norm > options.residual_tolerance)
    {
      out.status = out.iterations == options.max_iterations ? Status::iteration_limit : Status::stalled;
      return out;
    }
    auto const energy = candidate_energy<Real>(sites, out.roots, eta, t);
    if (std::abs(energy.imag()) >
        Real{256} * uni20::numeric_limits<Real>::epsilon() * (Real{1} + std::abs(energy.real())))
    {
      out.status = Status::precision_limit;
      return out;
    }
    out.energy = energy.real();
    out.converged = true;
    out.status = Status::converged;
  }
  catch (std::overflow_error const&)
  {
    sync();
    out.status = Status::precision_limit;
  }
  catch (std::domain_error const&)
  {
    sync();
    out.status = Status::precision_limit;
  }
  return out;
}
} // namespace bethe::xyz
