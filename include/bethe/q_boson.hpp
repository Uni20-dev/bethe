// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/continuum_newton.hpp>
#include <optional>
#include <span>
#include <uni20/common/half_int.hpp>

namespace bethe::q_boson
{
enum class Status
{
  converged,
  iteration_limit,
  stalled,
  precision_limit
};
template <uni20::Real Real> using SolverOptions = bethe::SolverOptions<Real>;
template <uni20::Real Real> struct State
{
    std::size_t sites = 0, particles = 0;
    Real eta{}, residual_norm{};
    std::vector<uni20::half_int> quantum_numbers;
    std::vector<Real> momenta;
    std::optional<Real> energy;
    std::size_t iterations = 0;
    bool converged = false;
    Status status = Status::iteration_limit;
};
namespace detail
{
template <uni20::Real Real> struct GroundSystem
{
    std::size_t sites, particles;
    Real t; // tanh(eta); never construct exp(eta)
    Real const pi = Real{4} * std::atan(Real{1});
    bool physical(std::span<Real const> k) const
    {
      for (std::size_t j = 0; j < k.size(); ++j)
        if (!uni20::isfinite(k[j]) || !(std::abs(k[j]) < pi) || (j && !(k[j] > k[j - 1]))) return false;
      return true;
    }
    struct Evaluation
    {
        std::vector<Real> residual;
        Real norm{};
    };
    Evaluation evaluate(std::span<Real const> k, uni20::DenseMatrix<Real>* jacobian = nullptr) const
    {
      Evaluation out{std::vector<Real>(particles)};
      for (std::size_t j = 0; j < particles; ++j)
      {
        bethe::detail::CompensatedSum<Real> sum, diagonal;
        Real const target = pi * (Real{2} * Real(j) - Real(particles - 1));
        sum.add(Real(sites) * k[j] - (t < Real{0.5} ? Real{0} : target));
        diagonal.add(Real(sites));
        for (std::size_t l = 0; l < particles; ++l)
          if (l != j)
          {
            Real const half = (k[j] - k[l]) / Real{2}, s = std::sin(half), c = std::cos(half);
            // Cancel the exact rank phase against consecutive ground labels
            // before floating-point arithmetic in the weak-deformation limit.
            sum.add(t < Real{0.5} ? -Real{2} * std::atan(t * (c / s)) : Real{2} * std::atan2(s, t * c));
            if (jacobian)
            {
              Real const h = std::hypot(s, t * c), kernel = (t / h) / h;
              if (!uni20::isfinite(kernel)) throw std::overflow_error("q-boson Jacobian exceeds scalar precision");
              (*jacobian)[j, l] = -kernel;
              diagonal.add(kernel);
            }
          }
        if (jacobian) (*jacobian)[j, j] = diagonal.value();
        out.residual[j] = sum.value();
        Real const scale = Real(sites) * std::max(std::sqrt(t / Real(sites)), std::abs(k[j]));
        out.norm = std::max(out.norm, std::abs(out.residual[j]) / scale);
        if (!uni20::isfinite(out.residual[j]) || !uni20::isfinite(out.norm))
          throw std::overflow_error("q-boson residual exceeds scalar precision");
      }
      return out;
    }
};
} // namespace detail

/// H=-sum_j(B_j^dagger B_{j+1}+h.c.-2N_j), [n]_q=(1-exp(-2 eta n))/(1-exp(-2 eta)).
/// eta=0 is free hopping; +infinity is the phase model. L=2 includes both periodic bonds.
/// Only the fixed-N ground state is selected; failed solves have no published energy.
template <uni20::Real Real>
State<Real> ground_state(std::size_t sites, std::size_t particles, Real eta, SolverOptions<Real> options = {})
{
  if (sites < 2 || !(eta >= Real{0}) || (!uni20::isfinite(eta) && eta != uni20::numeric_limits<Real>::infinity()) ||
      !uni20::isfinite(options.residual_tolerance) || !(options.residual_tolerance > Real{0}))
    throw std::invalid_argument("q-boson requires L>=2, eta>=0 (or +infinity), and finite positive tolerance");
  if (particles > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max() / 4))
    throw std::length_error("q-boson particle number exceeds label range");
  auto const elements = static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(Real);
  if (particles && particles > elements / particles) throw std::length_error("q-boson Newton matrix too large");
  State<Real> out;
  out.sites = sites;
  out.particles = particles;
  out.eta = eta;
  out.momenta.resize(particles);
  out.quantum_numbers.resize(particles);
  Real const pi = Real{4} * std::atan(Real{1}), t = std::tanh(eta);
  for (std::size_t j = 0; j < particles; ++j)
  {
    auto const twice = 2 * static_cast<std::int64_t>(j) - static_cast<std::int64_t>(particles - 1);
    out.quantum_numbers[j] = uni20::from_twice(twice);
    Real const phase_root = pi * Real(twice) / (Real(sites) + Real(particles));
    Real const weak_root =
        Real(twice) * std::sqrt(t) / std::sqrt(Real(sites)) / std::sqrt(Real(std::max(std::size_t{1}, particles)));
    out.momenta[j] = eta == Real{0} ? Real{0}
                     : t == Real{1} ? phase_root
                                    : (phase_root < Real{0} ? -Real{1} : Real{1}) *
                                          std::min(std::abs(phase_root), std::abs(weak_root));
  }
  if (particles < 2 || eta == Real{0})
  {
    out.energy = Real{0};
    out.converged = true;
    out.status = Status::converged;
    return out;
  }
  detail::GroundSystem<Real> system{sites, particles, t};
  if (!system.physical(out.momenta) || !(t / Real(sites) > Real{0}))
  {
    out.status = Status::precision_limit;
    return out;
  }
  try
  {
    out.iterations = bethe::detail::continuum_newton(system, out.momenta, options);
    out.residual_norm = system.evaluate(out.momenta).norm;
  }
  catch (std::overflow_error const&)
  {
    out.status = Status::precision_limit;
    return out;
  }
  if (out.residual_norm > options.residual_tolerance)
  {
    out.status = out.iterations == options.max_iterations ? Status::iteration_limit : Status::stalled;
    return out;
  }
  bethe::detail::CompensatedSum<Real> energy;
  for (Real k : out.momenta)
  {
    Real const s = std::sin(k / Real{2});
    energy.add((Real{4} * s) * s);
  }
  Real const e = energy.value();
  Real const floor = (uni20::numeric_limits<Real>::min() * uni20::numeric_limits<Real>::epsilon()) * Real{8};
  if (!(e > Real{0}) || !uni20::isfinite(e) || !(e * options.residual_tolerance >= floor))
  {
    out.status = Status::precision_limit;
    return out;
  }
  out.energy = e;
  out.converged = true;
  out.status = Status::converged;
  return out;
}
} // namespace bethe::q_boson
