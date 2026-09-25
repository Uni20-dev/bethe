// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/continuum_newton.hpp>
#include <bethe/real_excitations.hpp>
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
using QuantumNumbers = std::vector<uni20::half_int>;
template <uni20::Real Real> struct State
{
    std::size_t sites = 0, particles = 0;
    Real eta{}, residual_norm{};
    std::vector<uni20::half_int> quantum_numbers;
    std::vector<Real> momenta;
    /// Canonical free-boson mode labels; interacting roots are 2*pi*m/L+deviation.
    std::vector<std::size_t> modes;
    std::vector<Real> deviations;
    std::size_t momentum_index = 0; // sum(m) modulo L
    Real momentum{};                // physical momentum in (-pi,pi], reduced as integers first
    std::optional<Real> energy;
    std::size_t iterations = 0;
    bool converged = false;
    Status status = Status::iteration_limit;
};
namespace detail
{
template <uni20::Real Real>
void validate(std::size_t sites, std::size_t particles, Real eta, SolverOptions<Real> const& options)
{
  if (sites < 2 || !(eta >= Real{0}) || (!uni20::isfinite(eta) && eta != uni20::numeric_limits<Real>::infinity()) ||
      !uni20::isfinite(options.residual_tolerance) || !(options.residual_tolerance > Real{0}))
    throw std::invalid_argument("q-boson requires L>=2, eta>=0 (or +infinity), and finite positive tolerance");
  auto const limit = static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max() / 4);
  if (particles > limit || sites > limit) throw std::length_error("q-boson counts exceed label range");
  auto const elements = static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(Real);
  if (particles && particles > elements / particles) throw std::length_error("q-boson Newton matrix too large");
}
template <uni20::Real Real> struct System
{
    std::size_t sites;
    std::span<std::size_t const> modes;
    Real t; // tanh(eta); never construct exp(eta)
    Real const pi = Real{4} * std::atan(Real{1});
    Real difference(std::span<Real const> u, std::size_t j, std::size_t l) const
    {
      return Real{2} * pi * Real(static_cast<__int128>(modes[j]) - modes[l]) / Real(sites) + (u[j] - u[l]);
    }
    bool physical(std::span<Real const> u) const
    {
      for (std::size_t j = 0; j < u.size(); ++j)
        if (!uni20::isfinite(u[j]) || (j && !(difference(u, j, j - 1) > Real{0}))) return false;
      if (u.size() > 1 && !(difference(u, u.size() - 1, 0) < Real{2} * pi)) return false;
      return true;
    }
    struct Evaluation
    {
        std::vector<Real> residual;
        Real norm{};
    };
    Evaluation evaluate(std::span<Real const> k, uni20::DenseMatrix<Real>* jacobian = nullptr) const
    {
      auto const particles = modes.size();
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
            Real const half = difference(k, j, l) / Real{2}, s = std::sin(half), c = std::cos(half);
            // Cancel the exact rank phase against the label's rank contribution
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
/// Sorted free mode labels 0<=m_j<L specify a canonical Bethe branch.
/// Failed solves have no published energy. Rounded momenta can coincide at tiny eta;
/// use the mode/deviation representation to retain their separation.
template <uni20::Real Real>
State<Real> solve_modes(std::size_t sites, std::span<std::size_t const> modes, Real eta,
                        SolverOptions<Real> options = {})
{
  auto const particles = modes.size();
  detail::validate(sites, particles, eta, options);
  State<Real> out;
  out.sites = sites;
  out.particles = particles;
  out.eta = eta;
  out.momenta.resize(particles);
  out.modes.assign(modes.begin(), modes.end());
  out.deviations.resize(particles);
  out.quantum_numbers.resize(particles);
  Real const pi = Real{4} * std::atan(Real{1}), t = std::tanh(eta);
  __int128 total = 0;
  for (std::size_t j = 0; j < particles; ++j)
  {
    if (modes[j] >= sites || (j && modes[j] < modes[j - 1])) throw std::invalid_argument("require sorted modes 0<=m<L");
    total += modes[j];
  }
  out.momentum_index = static_cast<std::size_t>(total % sites);
  auto const centered_momentum = static_cast<std::int64_t>(out.momentum_index) -
                                 (out.momentum_index > sites / 2 ? static_cast<std::int64_t>(sites) : 0);
  out.momentum = Real{2} * pi * Real(centered_momentum) / Real(sites);
  for (std::size_t j = 0; j < particles; ++j)
  {
    auto const twice = 2 * static_cast<std::int64_t>(j) - static_cast<std::int64_t>(particles - 1);
    out.quantum_numbers[j] = uni20::from_twice(twice + 2 * static_cast<std::int64_t>(modes[j]));
    out.deviations[j] =
        pi * (Real(twice) + Real{2} * Real(total - static_cast<__int128>(particles) * modes[j]) / Real(sites)) /
        (Real(sites) + Real(particles));
  }
  if (t < Real{0.5})
  {
    for (std::size_t begin = 0; begin < particles;)
    {
      auto end = begin + 1;
      while (end < particles && modes[end] == modes[begin])
        ++end;
      auto const count = end - begin;
      Real const factor =
          std::min(std::sqrt(t) / std::sqrt(Real(sites)) / std::sqrt(Real(count)), pi / (Real(sites) * Real(count)));
      for (auto j = begin; j < end; ++j)
        out.deviations[j] = (Real{2} * Real(j - begin) - Real(count - 1)) * factor;
      begin = end;
    }
  }
  auto sync_momenta = [&] {
    for (std::size_t j = 0; j < particles; ++j)
      out.momenta[j] = Real{2} * pi * Real(modes[j]) / Real(sites) + out.deviations[j];
  };
  sync_momenta();
  detail::System<Real> system{sites, modes, t};
  if (eta > Real{0} && particles > 1 && (!system.physical(out.deviations) || !(t / Real(sites) > Real{0})))
  {
    out.status = Status::precision_limit;
    return out;
  }
  try
  {
    if (eta > Real{0} && particles > 1)
    {
      out.iterations = bethe::detail::continuum_newton(system, out.deviations, options);
      out.residual_norm = system.evaluate(out.deviations).norm;
    }
  }
  catch (std::overflow_error const&)
  {
    sync_momenta();
    out.status = Status::precision_limit;
    return out;
  }
  sync_momenta();
  if (out.residual_norm > options.residual_tolerance)
  {
    out.status = out.iterations == options.max_iterations ? Status::iteration_limit : Status::stalled;
    return out;
  }
  bethe::detail::CompensatedSum<Real> energy;
  for (std::size_t j = 0; j < particles; ++j)
  {
    auto const centered =
        static_cast<std::int64_t>(modes[j]) - (modes[j] > sites / 2 ? static_cast<std::int64_t>(sites) : 0);
    Real const s = std::sin(pi * Real(centered) / Real(sites) + out.deviations[j] / Real{2});
    energy.add((Real{4} * s) * s);
  }
  Real const e = energy.value();
  Real const floor = (uni20::numeric_limits<Real>::min() * uni20::numeric_limits<Real>::epsilon()) * Real{8};
  bool const exact_zero = total == 0 && (particles < 2 || eta == Real{0});
  if (!exact_zero && (!(e > Real{0}) || !uni20::isfinite(e) || !(e * options.residual_tolerance >= floor)))
  {
    out.status = Status::precision_limit;
    return out;
  }
  out.energy = e;
  out.converged = true;
  out.status = Status::converged;
  return out;
}

template <uni20::Real Real>
State<Real> ground_state(std::size_t sites, std::size_t particles, Real eta, SolverOptions<Real> options = {})
{
  detail::validate(sites, particles, eta, options);
  std::vector<std::size_t> modes(particles, 0);
  return solve_modes<Real>(sites, modes, eta, options);
}

template <uni20::Real Real>
State<Real> solve_real(std::size_t sites, Real eta, QuantumNumbers const& numbers, SolverOptions<Real> options = {})
{
  detail::validate(sites, numbers.size(), eta, options);
  std::vector<std::size_t> modes(numbers.size());
  for (std::size_t j = 0; j < numbers.size(); ++j)
  {
    __int128 const twice_mode = static_cast<__int128>(numbers[j].twice()) - 2 * static_cast<__int128>(j) +
                                static_cast<__int128>(numbers.size()) - 1;
    if (twice_mode < 0 || twice_mode % 2 || twice_mode / 2 >= sites)
      throw std::invalid_argument("invalid canonical q-boson labels");
    modes[j] = static_cast<std::size_t>(twice_mode / 2);
  }
  return solve_modes<Real>(sites, modes, eta, options);
}

template <uni20::Real Real> using ExcitationScan = bethe::RealExcitationScan<State<Real>>;
inline std::size_t excitation_count(std::size_t sites, std::size_t particles, std::size_t max_candidates = 10000)
{
  auto const limit = static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max() / 4);
  if (sites < 2 || sites > limit || particles > limit)
    throw std::invalid_argument("q-boson counts outside label range");
  return bethe::detail::bounded_binomial(sites + particles - 1, particles, max_candidates);
}
template <uni20::Real Real>
ExcitationScan<Real> real_excitations(std::size_t sites, std::size_t particles, Real eta,
                                      bethe::RealExcitationOptions scan = {}, SolverOptions<Real> options = {})
{
  (void)excitation_count(sites, particles, scan.max_candidates);
  return bethe::detail::scan_real_combinations<ExcitationScan<Real>>(
      sites + particles - 1, particles, particles ? -static_cast<std::int64_t>(particles - 1) : 0, scan,
      [&](auto const& numbers) { return solve_real(sites, eta, numbers, options); },
      [&] { return ground_state(sites, particles, eta, options); }, bethe::detail::EnergyOrder::ascending,
      [](auto const& state) { return *state.energy; });
}
} // namespace bethe::q_boson
