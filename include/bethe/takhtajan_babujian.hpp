// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/newton.hpp>
#include <bethe/solver.hpp>
#include <complex>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <uni20/common/half_int.hpp>

namespace bethe::takhtajan_babujian
{
template <uni20::Real Real> using SolverOptions = bethe::SolverOptions<Real>;
enum class SolveStatus
{
  converged,
  iteration_limit,
  stalled,
  ill_conditioned
};

/// Spin 1, H=sum(S.S-(S.S)^2), periodic even L>=4, zero-field singlet.
template <uni20::Real Real> struct State
{
    std::size_t sites = 0, iterations = 0;
    std::vector<uni20::half_int> string_quantum_numbers;
    std::vector<Real> centers, deviations;      // lambda = center +/- i*(1/2+deviation).
    std::vector<std::complex<Real>> rapidities; // Adjacent conjugate pairs.
    Real energy{}, momentum{}, residual_norm{}, phase_residual{}, modulus_residual{};
    std::size_t momentum_index = 0;
    bool converged = false;
    SolveStatus status = SolveStatus::iteration_limit;
};

namespace detail
{
inline void check_sites(std::size_t n)
{
  if (n < 4 || n % 2 || n > std::size_t(std::numeric_limits<std::int64_t>::max() / 4))
    throw std::invalid_argument("spin-1 TB ground state requires even L>=4 and L<=INT64_MAX/4");
}

// Vlijm-Caux (3.8),(3.10), for the filled two-string sea. In this branch
// J_++J_-=0 and every deviation is positive. Work with delta itself, avoiding
// cancellation in 1-y_j-y_k and in the within-pair scattering denominator.
template <uni20::Real Real> class GroundSystem {
  public:
    explicit GroundSystem(std::size_t n)
        : sites(n), strings(n / 2), positive(strings / 2), first(positive + strings % 2), order(strings),
          pi(Real{4} * std::atan(Real{1}))
    {}
    std::size_t center_index(std::size_t j) const { return j < positive ? positive - 1 - j : j - first; }
    std::size_t deviation_index(std::size_t j) const { return j < positive ? strings - 1 - j : j; }
    Real center_sign(std::size_t j) const { return j < positive ? -Real{1} : j >= first ? Real{1} : Real{0}; }
    // Variables: positive centers, followed by deviations at nonnegative centers.
    std::pair<std::vector<Real>, std::vector<Real>> expand(std::span<Real const> v) const
    {
      std::vector<Real> x(strings), d(strings);
      for (std::size_t j = 0; j < strings; ++j)
      {
        if (center_sign(j) != Real{0}) x[j] = center_sign(j) * v[center_index(j)];
        d[j] = v[deviation_index(j)];
      }
      return {std::move(x), std::move(d)};
    }
    std::vector<Real> seed() const
    {
      std::vector<Real> v(order);
      for (std::size_t j = positive; j < strings; ++j)
      {
        Real const angle = pi * Real(2 * std::int64_t(j) - std::int64_t(strings - 1)) / Real(sites);
        Real const x = std::asinh(std::tan(angle)) / pi;
        if (j >= first) v[center_index(j)] = x;
        v[deviation_index(j)] = std::log(Real{2}) / (Real{2} * pi * Real(sites) * std::cos(angle));
      }
      return v;
    }
    bool physical(std::span<Real const> v) const
    {
      if (v.size() != order) return false;
      Real const limit = std::sqrt(uni20::numeric_limits<Real>::max()) / Real{16};
      for (std::size_t j = 0; j < order; ++j)
      {
        if (!uni20::isfinite(v[j]) || v[j] <= Real{0}) return false;
        if (j < positive)
        {
          if (v[j] > limit || (j && v[j] <= v[j - 1])) return false;
        }
        else if (Real{1} / Real{2} + v[j] >= Real{1} || Real{1} / Real{2} + v[j] == Real{1} / Real{2})
          return false;
      }
      return true;
    }
    struct Evaluation
    {
        std::vector<Real> residual;
        Real phase{}, modulus{};
        Real norm() const { return std::max(phase, modulus); }
    };
    Evaluation evaluate(std::span<Real const> v, std::vector<Real>* jac = nullptr) const
    {
      auto const [x, d] = expand(v);
      Evaluation out{.residual = std::vector<Real>(order)};
      if (jac) jac->assign(order * order, Real{0});
      Real const invn = Real{1} / Real(sites);
      for (std::size_t j = positive; j < strings; ++j)
      {
        // The central string's phase vanishes by reflection; only its
        // modulus equation remains independent.
        bool const has_phase = j >= first;
        auto const pr = has_phase ? center_index(j) : 0, mr = deviation_index(j);
        bethe::detail::CompensatedSum<Real> phase, modulus;
        auto add = [&](Real a, Real b, Real pc, Real mc, std::size_t k, Real bj, Real bk, bool scatter) {
          Real const square = a * a + b * b;
          phase.add(pc * std::atan2(a, b));
          modulus.add(mc * std::log(square) / Real{2});
          if (!jac) return;
          Real const pa = pc * b / square, pb = -pc * a / square;
          Real const ma = mc * a / square, mb = mc * b / square;
          auto entry = [&](std::size_t row, Real da, Real db) {
            if (center_sign(j) != Real{0}) (*jac)[row * order + center_index(j)] += da * center_sign(j);
            (*jac)[row * order + deviation_index(j)] += db * bj;
            if (scatter)
            {
              if (center_sign(k) != Real{0}) (*jac)[row * order + center_index(k)] -= da * center_sign(k);
              (*jac)[row * order + deviation_index(k)] += db * bk;
            }
          };
          if (has_phase) entry(pr, pa, pb);
          entry(mr, ma, mb);
        };
        add(x[j], Real{3} / Real{2} + d[j], Real{1}, Real{1}, 0, Real{1}, Real{0}, false);
        add(x[j], Real{1} / Real{2} - d[j], Real{1}, -Real{1}, 0, -Real{1}, Real{0}, false);
        modulus.add(-(std::log1p(d[j]) - std::log(d[j])) * invn);
        if (jac) (*jac)[mr * order + deviation_index(j)] += invn / (d[j] * (Real{1} + d[j]));
        for (std::size_t k = 0; k < strings; ++k)
        {
          if (j == k) continue;
          Real const a = x[j] - x[k];
          add(a, Real{2} + d[j] + d[k], -invn, -invn, k, Real{1}, Real{1}, true);
          add(a, -d[j] - d[k], -invn, invn, k, -Real{1}, -Real{1}, true);
          add(a, Real{1} + d[j] - d[k], -invn, -invn, k, Real{1}, -Real{1}, true);
          add(a, Real{1} - d[j] + d[k], -invn, invn, k, -Real{1}, Real{1}, true);
        }
        if (has_phase)
        {
          out.residual[pr] = phase.value();
          out.phase = std::max(out.phase, std::abs(phase.value()));
        }
        out.residual[mr] = modulus.value();
        out.modulus = std::max(out.modulus, std::abs(modulus.value()));
      }
      for (Real r : out.residual)
        if (!uni20::isfinite(r)) throw std::runtime_error("nonfinite TB equation residual");
      return out;
    }
    std::size_t const sites, strings, positive, first, order;
    Real const pi;
};
} // namespace detail

template <uni20::Real Real = double>
[[nodiscard]] State<Real> ground_state(std::size_t sites, SolverOptions<Real> const& options = {})
{
  detail::check_sites(sites);
  if (!uni20::isfinite(options.residual_tolerance) || options.residual_tolerance <= Real{0})
    throw std::invalid_argument("TB residual tolerance must be finite and positive");
  auto const order = sites / 2, elements = std::size_t(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(Real);
  if (order > elements / order) throw std::length_error("TB Newton matrix is too large");
  detail::GroundSystem<Real> const system(sites);
  auto x = system.seed();
  if (!system.physical(x)) throw std::runtime_error("TB seed is not representable at this precision");
  State<Real> out;
  out.sites = sites;
  std::vector<Real> jac;
  for (;;)
  {
    auto const evaluation = system.evaluate(x, &jac);
    out.phase_residual = evaluation.phase;
    out.modulus_residual = evaluation.modulus;
    out.residual_norm = evaluation.norm();
    if (out.residual_norm <= options.residual_tolerance)
    {
      out.converged = true;
      out.status = SolveStatus::converged;
      break;
    }
    if (out.iterations == options.max_iterations) break;
    auto step = evaluation.residual;
    for (auto& v : step)
      v = -v;
    if (!bethe::detail::newton_step(std::move(jac), step))
    {
      out.status = SolveStatus::ill_conditioned;
      break;
    }
    bool accepted = false;
    Real damping{1};
    auto trial = x;
    for (int backtrack = 0; backtrack <= uni20::numeric_limits<Real>::digits; ++backtrack)
    {
      for (std::size_t j = 0; j < x.size(); ++j)
        trial[j] = x[j] + damping * step[j];
      if (system.physical(trial))
      {
        Real const norm = system.evaluate(trial).norm();
        if (norm <= options.residual_tolerance || norm < (Real{1} - damping / Real{10000}) * out.residual_norm)
        {
          accepted = true;
          break;
        }
      }
      damping /= Real{2};
    }
    if (!accepted)
    {
      out.status = SolveStatus::stalled;
      break;
    }
    x = std::move(trial);
    ++out.iterations;
  }
  auto [centers, deviations] = system.expand(x);
  out.centers = std::move(centers);
  out.deviations = std::move(deviations);
  bethe::detail::CompensatedSum<Real> energy;
  for (std::size_t j = 0; j < order; ++j)
  {
    Real const c = out.centers[j], y = Real{1} / Real{2} + out.deviations[j];
    out.string_quantum_numbers.push_back(uni20::from_twice(2 * std::int64_t(j) - std::int64_t(order - 1)));
    out.rapidities.emplace_back(c, y);
    out.rapidities.emplace_back(c, -y);
    // Real part of 1/(1+(c+i*y)^2), paired before summation.
    Real const a = Real{1} + c * c - y * y, b = Real{2} * c * y;
    Real const contribution = std::abs(a) >= std::abs(b) ? (Real{1} / a) / (Real{1} + (b / a) * (b / a))
                                                         : ((a / b) / b) / (Real{1} + (a / b) * (a / b));
    energy.add(-Real{8} * contribution);
  }
  out.energy = energy.value();
  // A reflection-paired two-string sea has P=0, including the central pair.
  return out;
}
} // namespace bethe::takhtajan_babujian
