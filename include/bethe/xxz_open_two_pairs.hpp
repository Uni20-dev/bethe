// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/xxz_open_two_string.hpp>

namespace bethe::xxz::quantum_group::two_pairs
{
template <uni20::Real Real> struct State
{
    std::size_t sites{};
    Real delta{}, energy{}, energy_shift{}, residual_norm{}, phase_residual{}, modulus_residual{};
    std::array<std::size_t, 2> string_labels{};
    std::array<Real, 2> centers{}, log_deviations{};
    std::array<int, 2> deviation_signs{};
    std::size_t iterations{};
    bool converged = false;
    SolveStatus status = SolveStatus::iteration_limit;
};

namespace detail
{
template <uni20::Real Real> class System {
  public:
    using Pair = two_string::detail::System<Real>;
    static constexpr std::size_t order = 4;
    // For ordered centers, the other pair contributes an extra 3pi phase
    // to the upper pair's unsquared equation. Its deviation parity flips.
    System(std::size_t n, Real delta, std::array<std::size_t, 2> labels)
        : sites(n), labels(labels), eta(std::acosh(delta)), pi(Real{4} * std::atan(Real{1})),
          pairs{Pair(n, delta, 2, labels[0], (n - labels[0]) % 2 ? 1 : -1),
                Pair(n, delta, 2, labels[1], (n - labels[1]) % 2 ? -1 : 1)}
    {}
    bool physical(std::span<Real const> x) const
    {
      return x.size() == order && pairs[0].physical(x.first(2)) && pairs[1].physical(x.subspan(2, 2)) && x[0] < x[2];
    }
    std::vector<Real> coordinate_scales(std::span<Real const> x) const
    {
      return {std::min(x[0], pi - x[0]), Real{2} * Real(sites), std::min(x[2], pi - x[2]), Real{2} * Real(sites)};
    }
    std::vector<Real> seed() const
    {
      std::vector<Real> x{pi * Real(labels[0]) / Real(sites - 5), Real{8} + std::abs(std::log(eta)),
                          pi * Real(labels[1]) / Real(sites - 5), Real{8} + std::abs(std::log(eta))};
      if (!physical(x)) throw std::overflow_error("two-pair seed is not resolvable at the selected precision");
      auto const f = evaluate(x, nullptr, true);
      for (std::size_t i = 0; i < 2; ++i)
        x[2 * i + 1] = std::max(x[2 * i + 1] + Real{2} * Real(sites) * f.residual[2 * i + 1], Real{2} - std::log(eta));
      return x;
    }
    struct Evaluation
    {
        std::vector<Real> residual;
        Real norm{}, phase_norm{}, modulus_norm{};
    };
    Evaluation evaluate(std::span<Real const> x, std::vector<Real>* jac = nullptr, bool ideal = false) const
    {
      Real const scale = Real{1} / (Real{2} * Real(sites));
      if (jac) jac->assign(16, Real{0});
      auto add = [&](std::size_t row, std::size_t col, Real value) {
        if (jac) (*jac)[4 * row + col] += value;
      };
      Evaluation out{.residual = std::vector<Real>(4)};
      for (std::size_t i = 0; i < 2; ++i)
      {
        auto const j = 1 - i, a = 2 * i, l = a + 1, b = 2 * j, k = b + 1;
        std::vector<Real> block;
        auto const base = pairs[i].evaluate(x.subspan(a, 2), jac ? &block : nullptr, ideal);
        out.residual[a] = base.residual[0];
        out.residual[l] = base.residual[1];
        if (jac)
          for (std::size_t r = 0; r < 2; ++r)
            for (std::size_t c = 0; c < 2; ++c)
              add(a + r, a + c, block[2 * r + c]);
        Real const d = ideal ? Real{0} : Real(pairs[i].deviation_sign) * std::exp(-x[l]);
        Real const e = ideal ? Real{0} : Real(pairs[j].deviation_sign) * std::exp(-x[k]);
        Real const t = (d - e) / Real{2}, s = (d + e) / Real{2};
        for (int sign : {-1, 1})
        {
          Real const beta = x[a] + Real(sign) * x[b];
          auto A = Pair::phase(beta, eta + t), B = Pair::phase(beta, eta - t), C = Pair::phase(beta, Real{2} * eta + s);
          auto D = Pair::phase(beta, s, true);
          // Continuous small-width correction on BOTH sides of beta=0.
          // atan2's complementary phase jumps by 2pi for negative beta.
          D.value = Real{2} * std::atan(std::tanh(s) / std::tan(beta / Real{2}));
          Real const phase = A.value + B.value + C.value + D.value;
          Real const da = A.angle + B.angle + C.angle + D.angle;
          out.residual[a] -= scale * phase;
          add(a, a, -scale * da);
          add(a, b, -scale * Real(sign) * da);
          add(a, l, scale * d / Real{2} * (A.width - B.width + C.width + D.width));
          add(a, k, -scale * e / Real{2} * (A.width - B.width - C.width - D.width));
          A = Pair::log_sinh(beta, eta + t);
          B = Pair::log_sinh(beta, eta - t);
          C = Pair::log_sinh(beta, Real{2} * eta + s);
          D = Pair::log_sinh(beta, s);
          Real const mod = A.value - B.value + C.value - D.value;
          Real const ma = A.angle - B.angle + C.angle - D.angle;
          out.residual[l] -= scale * mod;
          add(l, a, -scale * ma);
          add(l, b, -scale * Real(sign) * ma);
          add(l, l, scale * d / Real{2} * (A.width + B.width + C.width - D.width));
          add(l, k, -scale * e / Real{2} * (A.width + B.width - C.width + D.width));
        }
      }
      for (std::size_t i = 0; i < 4; ++i)
      {
        if (!uni20::isfinite(out.residual[i])) throw std::overflow_error("nonfinite two-pair residual");
        out.norm = std::max(out.norm, std::abs(out.residual[i]));
        auto& norm = i % 2 ? out.modulus_norm : out.phase_norm;
        norm = std::max(norm, std::abs(out.residual[i]));
      }
      return out;
    }
    Real energy_shift(std::span<Real const> x, Real delta) const
    {
      return pairs[0].energy_shift(x.first(2), delta) + pairs[1].energy_shift(x.subspan(2, 2), delta);
    }
    std::size_t sites;
    std::array<std::size_t, 2> labels;
    Real eta, pi;
    std::array<Pair, 2> pairs;
};
} // namespace detail

/// Two ordered two-strings, with no real roots: ell=N-8, N>=8, Delta>1.
/// Candidate labels 1<=J1<J2<=N-6; no completeness/convergence guarantee.
template <uni20::Real Real = double>
[[nodiscard]] State<Real> solve(std::size_t sites, Real delta, std::array<std::size_t, 2> labels,
                                SolverOptions<Real> const& options = {})
{
  (void)xxz::detail::checked_sites(sites);
  if (sites < 8) throw std::invalid_argument("two bound pairs require N>=8");
  if (!uni20::isfinite(delta) || delta <= Real{1})
    throw std::invalid_argument("two bound pairs require finite Delta>1");
  if (labels[0] == 0 || labels[0] >= labels[1] || labels[1] > sites - 6)
    throw std::invalid_argument("two-pair labels require 1<=J1<J2<=N-6");
  detail::System<Real> system(sites, delta, labels);
  auto const it = quantum_group::detail::solve_log_string<Real>(system, options);
  auto const f = system.evaluate(it.x);
  State<Real> state{.sites = sites, .delta = delta};
  state.string_labels = labels;
  state.energy_shift = system.energy_shift(it.x, delta);
  state.energy = Real(sites - 1) * delta / Real{4} + state.energy_shift;
  if (!uni20::isfinite(state.energy) || !uni20::isfinite(state.energy_shift))
    throw std::overflow_error("two-pair energy overflow");
  state.residual_norm = f.norm;
  state.phase_residual = f.phase_norm;
  state.modulus_residual = f.modulus_norm;
  state.iterations = it.iterations;
  state.converged = it.converged;
  state.status = it.status;
  for (std::size_t i = 0; i < 2; ++i)
  {
    state.centers[i] = it.x[2 * i];
    state.log_deviations[i] = it.x[2 * i + 1];
    state.deviation_signs[i] = system.pairs[i].deviation_sign;
  }
  return state;
}
} // namespace bethe::xxz::quantum_group::two_pairs
