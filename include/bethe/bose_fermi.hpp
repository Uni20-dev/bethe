// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <array>
#include <bethe/lieb_liniger.hpp>
#include <optional>

namespace bethe::bose_fermi
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
    std::size_t bosons = 0, fermions = 0, iterations = 0;
    Real length{}, interaction{}, residual_norm{}, momentum{};
    std::int64_t momentum_index = 0;
    std::vector<Real> momenta, auxiliary;
    std::vector<uni20::half_int> charge_numbers, auxiliary_numbers;
    std::vector<std::int64_t> free_fermion_modes;
    std::optional<Real> energy;
    bool converged = false;
    Status status = Status::iteration_limit;
};
namespace detail
{
// Reflection-reduced bipartite Bethe system. There is deliberately no
// auxiliary-auxiliary scattering: auxiliary particles have fermionic grading.
template <uni20::Real Real> struct System
{
    std::size_t n, m;
    Real g;
    Real const pi = Real{4} * std::atan(Real{1});
    Real scale() const { return std::max(Real{1}, g); }
    std::size_t order() const { return n / 2 + m / 2; }
    std::array<std::vector<Real>, 2> expand(std::span<Real const> x) const
    {
      std::array<std::vector<Real>, 2> result;
      for (std::size_t species = 0; species < 2; ++species)
      {
        auto const count = species ? m : n, half = count / 2, offset = species ? n / 2 : 0;
        Real const s = species ? scale() : Real{1};
        for (auto j = half; j-- > 0;)
          result[species].push_back(-s * x[offset + j]);
        if (count % 2) result[species].push_back(Real{0});
        for (std::size_t j = 0; j < half; ++j)
          result[species].push_back(s * x[offset + j]);
      }
      return result;
    }
    bool physical(std::span<Real const> x) const
    {
      for (std::size_t j = 0; j < x.size(); ++j)
        if (!uni20::isfinite(x[j]) || !(x[j] > Real{0}) || (j && j != n / 2 && !(x[j] > x[j - 1])) ||
            (j >= n / 2 && !uni20::isfinite(x[j] * scale())))
          return false;
      return true;
    }
    std::vector<Real> seed() const
    {
      std::vector<Real> x(order());
      for (std::size_t j = 0; j < n / 2; ++j)
      {
        auto const rank = 2 * j + 1 + n % 2;
        x[j] = pi * Real(rank);
        if (g < Real{1})
          x[j] = j < (m + 1) / 2 ? Real(rank) * std::sqrt(g) / std::sqrt(Real(m + 1))
                                 : Real{2} * pi * Real(j - (m + 1) / 2 + 1);
      }
      for (std::size_t j = 0; j < m / 2; ++j)
      {
        auto const rank = 2 * j + 1 + m % 2;
        x[n / 2 + j] = g < Real{1} ? Real(rank) * std::sqrt(g) / std::sqrt(Real(m + 1))
                                   : g / scale() / Real{2} * std::tan(pi * Real(rank) / (Real{2} * Real(n)));
      }
      return x;
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
        bool const auxiliary = row >= n / 2;
        auto const count = auxiliary ? m : n, other = auxiliary ? n : m;
        auto const local = auxiliary ? row - n / 2 : row;
        auto const full = count / 2 + count % 2 + local;
        Real const value = roots[auxiliary][full], width = g / Real{2};
        std::int64_t rank = -static_cast<std::int64_t>(2 * local + 1 + count % 2);
        bethe::detail::CompensatedSum<Real> sum, diagonal;
        if (!auxiliary)
        {
          sum.add(value);
          diagonal.add(Real{1});
        }
        if (g >= Real{1}) sum.add(pi * Real(rank));
        for (std::size_t j = 0; j < other; ++j)
        {
          Real const d = value - roots[!auxiliary][j];
          if (g < Real{1} && std::abs(d) > width)
          {
            rank += d > Real{0} ? 1 : -1;
            sum.add(-Real{2} * std::atan(width / d));
          }
          else
            sum.add(Real{2} * std::atan(d / width));
          if (jacobian)
          {
            Real const kernel = Real{2} * bethe::detail::rational_scattering_kernel(d, width);
            diagonal.add(kernel * (auxiliary ? scale() : Real{1}));
            auto const offset = auxiliary ? 0 : n / 2, half = other / 2;
            Real const derivative = kernel * (auxiliary ? Real{1} : scale());
            if (j < half)
              (*jacobian)[row, offset + half - 1 - j] += derivative;
            else if (j >= half + other % 2)
              (*jacobian)[row, offset + j - half - other % 2] -= derivative;
          }
        }
        if (g < Real{1}) sum.add(pi * Real(rank));
        if (jacobian) (*jacobian)[row, row] += diagonal.value();
        out.residual[row] = sum.value();
        Real const normalization = !auxiliary && g < Real{1} ? std::max(std::sqrt(g), std::abs(value)) : Real(n);
        Real const residual = std::abs(sum.value()) / normalization;
        if (!uni20::isfinite(residual)) throw std::overflow_error("Bose-Fermi residual exceeds precision");
        out.norm = std::max(out.norm, residual);
      }
      if (jacobian)
        for (std::size_t i = 0; i < order(); ++i)
          for (std::size_t j = 0; j < order(); ++j)
            if (!uni20::isfinite((*jacobian)[i, j])) throw std::overflow_error("Bose-Fermi Jacobian exceeds precision");
      return out;
    }
};
} // namespace detail

/// Equal masses, equal repulsive BB/BF contact strength 2c, spinless fermions.
/// Mixed interacting shells require odd fermion number; pure/free limits do not.
template <uni20::Real Real>
State<Real> ground_state(std::size_t bosons, std::size_t fermions, Real length, Real c,
                         SolverOptions<Real> options = {})
{
  auto const limit = static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max() / 4);
  if (bosons > limit || fermions > limit - bosons) throw std::length_error("Bose-Fermi particle count too large");
  auto const n = bosons + fermions;
  if (!uni20::isfinite(length) || !(length > Real{0}) || !uni20::isfinite(c) || c < Real{0} ||
      !uni20::isfinite(options.residual_tolerance) || !(options.residual_tolerance > Real{0}))
    throw std::invalid_argument("Bose-Fermi requires finite length>0, c>=0 and positive finite tolerance");
  if (bosons && fermions && c > Real{0} && fermions % 2 == 0)
    throw std::invalid_argument("interacting mixed periodic shells currently require odd fermion number");
  auto const elements = static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(Real);
  if (n && n > elements / n) throw std::length_error("Bose-Fermi matrix too large");
  State<Real> out;
  out.bosons = bosons;
  out.fermions = fermions;
  out.length = length;
  out.interaction = c;
  Real const pi = Real{4} * std::atan(Real{1});
  if (c == Real{0} || !bosons)
  {
    out.momenta.assign(bosons, Real{0});
    for (std::size_t j = 0; j < fermions; ++j)
    {
      auto const mode = static_cast<std::int64_t>(j) - static_cast<std::int64_t>((fermions - 1) / 2);
      out.free_fermion_modes.push_back(mode);
      out.momenta.push_back(Real{2} * pi * Real(mode) / length);
    }
    out.momentum_index = fermions % 2 ? 0 : static_cast<std::int64_t>(fermions / 2);
    out.momentum = Real{2} * pi * Real(out.momentum_index) / length;
  }
  else if (!fermions)
  {
    auto const state = bethe::lieb_liniger::ground_state(n, length, c, options);
    out.momenta = state.momenta;
    out.charge_numbers = state.quantum_numbers;
    out.iterations = state.iterations;
    out.residual_norm = state.residual_norm;
    if (!state.converged)
    {
      out.status =
          state.status == bethe::lieb_liniger::SolveStatus::iteration_limit ? Status::iteration_limit : Status::stalled;
      return out;
    }
  }
  else
  {
    Real const g = c * length;
    if (!uni20::isfinite(g) || !(g / Real{2} > Real{0}))
      throw std::invalid_argument("c*length and c*length/2 must be finite and positive");
    detail::System<Real> system{n, bosons, g};
    auto x = system.seed();
    for (std::size_t j = 0; j < n; ++j)
      out.charge_numbers.push_back(
          uni20::from_twice(2 * static_cast<std::int64_t>(j) - static_cast<std::int64_t>(n - 1)));
    for (std::size_t j = 0; j < bosons; ++j)
      out.auxiliary_numbers.push_back(
          uni20::from_twice(2 * static_cast<std::int64_t>(j) - static_cast<std::int64_t>(bosons - 1)));
    try
    {
      out.iterations = bethe::detail::continuum_newton(system, x, options);
      auto roots = system.expand(x);
      for (Real q : roots[0])
        out.momenta.push_back(q / length);
      for (Real q : roots[1])
        out.auxiliary.push_back(q / length);
      // Recheck the returned physical roots after dimensional conversion.
      for (std::size_t j = 0; j < n / 2; ++j)
        x[j] = out.momenta[n / 2 + n % 2 + j] * length;
      for (std::size_t j = 0; j < bosons / 2; ++j)
        x[n / 2 + j] = out.auxiliary[bosons / 2 + bosons % 2 + j] * length / system.scale();
      if (!system.physical(x))
      {
        out.status = Status::precision_limit;
        return out;
      }
      out.residual_norm = system.evaluate(x).norm;
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
  }
  bethe::detail::CompensatedSum<Real> sum;
  for (Real k : out.momenta)
  {
    if (!uni20::isfinite(k) || (k != Real{0} && k * k == Real{0}))
    {
      out.status = Status::precision_limit;
      return out;
    }
    sum.add(k * k);
  }
  Real const floor = (uni20::numeric_limits<Real>::min() * uni20::numeric_limits<Real>::epsilon()) * Real{8};
  if (!uni20::isfinite(sum.value()) || !uni20::isfinite(out.momentum) ||
      (sum.value() > Real{0} && !(sum.value() * options.residual_tolerance >= floor)))
  {
    out.status = Status::precision_limit;
    return out;
  }
  out.energy = sum.value();
  out.converged = true;
  out.status = Status::converged;
  return out;
}
} // namespace bethe::bose_fermi
