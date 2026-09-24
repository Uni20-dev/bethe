// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/xxz_quantum_group.hpp>
#include <utility>

namespace bethe::xxz::quantum_group::two_string
{
/// One two-string, optionally above a real sea.
/// u_pair=(eta+sign*exp(-log_deviation))/2 +/- i*center/2; u_real=i*alpha/2.
/// log_deviation is authoritative: adding exp(-L) to eta may round it away.
/// NOT a spectrum enumeration or a guarantee of global excitation ordering.
template <uni20::Real Real> struct State
{
    std::size_t sites = 0;
    Real delta{}, center{}, log_deviation{};
    std::vector<Real> rapidities;
    Real energy{}, residual_norm{}, phase_residual{}, modulus_residual{};
    std::size_t iterations = 0;
    bool converged = false;
    SolveStatus status = SolveStatus::iteration_limit;
    std::size_t string_label = 1;
    int deviation_sign = 1;
    /// E_ref-(N-1)*Delta/4, evaluated without subtracting extensive energies.
    Real energy_shift{};
};

namespace detail
{
template <uni20::Real Real> class System {
  public:
    System(std::size_t n, Real delta) : System(n, delta, n / 2, 1, 1) {}
    System(std::size_t n, Real delta, std::size_t roots, std::size_t label, int sign)
        : sites(n), order(roots), sea(order - 2), string_label(label), deviation_sign(sign), eta(std::acosh(delta)),
          pi(Real{4} * std::atan(Real{1}))
    {}
    struct Function
    {
        Real value, angle, width;
    };
    // Theta(beta;w) and its two derivatives, without coth(w) near w=0.
    static Function phase(Real beta, Real w, bool complement = false)
    {
      Real const s = std::sin(beta / Real{2}), c = std::cos(beta / Real{2}), t = std::tanh(w);
      Real const den = s * s + t * t * c * c;
      Real const da = t / den, dw = -Real{2} * s * c * (Real{1} - t * t) / den;
      if (complement) return {Real{2} * std::atan2(t * c, s), -da, -dw};
      return {Real{2} * std::atan2(s, t * c), da, dw};
    }
    // log|sinh(w+i*beta/2)| and derivatives. Widths here are O(eta), not L.
    static Function log_sinh(Real beta, Real w)
    {
      Real const s = std::sin(beta / Real{2}), c = std::cos(beta / Real{2}), sh = std::sinh(w);
      Real const den = sh * sh + s * s;
      return {std::log(den) / Real{2}, s * c / (Real{2} * den), sh * std::cosh(w) / den};
    }
    // log(sinh(d)/d), and d*coth(d), including d=0 after underflow.
    static std::pair<Real, Real> sinh_correction(Real d)
    {
      if (std::abs(d) < std::sqrt(uni20::numeric_limits<Real>::epsilon()))
        return {d * d / Real{6}, Real{1} + d * d / Real{3}};
      return {std::log(std::sinh(d) / d), d / std::tanh(d)};
    }
    bool physical(std::span<Real const> x) const
    {
      if (x.size() != order) return false;
      for (Real v : x)
        if (!uni20::isfinite(v)) return false;
      for (std::size_t i = 0; i < sea; ++i)
        if (!(x[i] > Real{0} && x[i] < pi) || (i && x[i] <= x[i - 1])) return false;
      return x[sea] > Real{0} && x[sea] < pi && x[sea + 1] > -std::log(eta);
    }
    std::vector<Real> seed() const
    {
      std::vector<Real> x(order);
      for (std::size_t i = 0; i < sea; ++i)
        x[i] = Real{2} * std::atan(std::tanh(eta / Real{2}) * std::tan(pi * Real(i + 1) / (Real{2} * Real(sites))));
      x[sea] = sea ? pi / Real{2} : pi * (Real(string_label) / Real(sites - 2));
      x[sea + 1] = Real{8} + std::abs(std::log(eta));
      // Solve the ideal-string modulus equation for L at the bare sea seed.
      x[sea + 1] += Real{2} * Real(sites) * evaluate(x, nullptr, true).residual.back();
      x[sea + 1] = std::max(x[sea + 1], Real{2} - std::log(eta));
      return x;
    }
    struct Evaluation
    {
        std::vector<Real> residual;
        Real norm{}, phase_norm{}, modulus_norm{};
    };
    // All equations divided by 2N. Ideal mode is ONLY a Newton initializer;
    // convergence always uses the finite-deviation system (ideal=false).
    Evaluation evaluate(std::span<Real const> x, std::vector<Real>* jac = nullptr, bool ideal = false) const
    {
      Real const scale = Real{1} / (Real{2} * Real(sites));
      Real const a = x[sea], L = x[sea + 1], d = ideal ? Real{0} : Real(deviation_sign) * std::exp(-L);
      Real const wp = (Real{3} * eta + d) / Real{2}, wm = (eta - d) / Real{2};
      if (jac) jac->assign(order * order, Real{0});
      auto add = [&](std::size_t i, std::size_t j, Real v) {
        if (jac) (*jac)[i * order + j] += v * scale;
      };
      Evaluation out{.residual = std::vector<Real>(order)};
      for (std::size_t i = 0; i < sea; ++i)
      {
        auto drive = phase(x[i], eta / Real{2});
        bethe::detail::CompensatedSum<Real> f;
        f.add(drive.value - Real{2} * pi * Real(i + 1) * scale);
        add(i, i, Real{2} * Real(sites) * drive.angle);
        for (std::size_t j = 0; j < sea; ++j)
          if (i != j)
          {
            auto minus = phase(x[i] - x[j], eta), plus = phase(x[i] + x[j], eta);
            f.add(-scale * minus.value);
            f.add(-scale * plus.value);
            add(i, i, -minus.angle - plus.angle);
            add(i, j, minus.angle - plus.angle);
          }
        for (int sign : {-1, 1})
          for (int side : {-1, 1})
          {
            auto p = phase(x[i] + Real(sign) * a, side == 1 ? wp : wm);
            f.add(-scale * p.value);
            add(i, i, -p.angle);
            add(i, sea, -Real(sign) * p.angle);
            add(i, sea + 1, Real(side) * d * p.width / Real{2});
          }
        out.residual[i] = f.value();
      }
      auto p = phase(a, eta + d / Real{2}), c = phase(a, d / Real{2}, true), self = phase(Real{2} * a, eta);
      bethe::detail::CompensatedSum<Real> ph, mod;
      ph.add(p.value + c.value);
      ph.add(-scale * (Real{2} * self.value + Real{2} * pi * Real(string_label)));
      add(sea, sea, Real{2} * Real(sites) * (p.angle + c.angle) - Real{4} * self.angle);
      add(sea, sea + 1, -Real(sites) * d * (p.width + c.width));
      auto lp = log_sinh(a, eta + d / Real{2}), lm = log_sinh(a, d / Real{2});
      auto const [correction, dcoth] = sinh_correction(d);
      mod.add(lp.value - lm.value);
      mod.add(scale * (-std::log(std::sinh(Real{2} * eta + d)) - L + correction));
      add(sea + 1, sea, Real{2} * Real(sites) * (lp.angle - lm.angle));
      add(sea + 1, sea + 1, -Real(sites) * d * (lp.width - lm.width) + d / std::tanh(Real{2} * eta + d) - dcoth);
      for (std::size_t j = 0; j < sea; ++j)
        for (int sign : {-1, 1})
          for (int side : {-1, 1})
          {
            Real const beta = a + Real(sign) * x[j], w = side == 1 ? wp : wm;
            auto scatter = phase(beta, w), log = log_sinh(beta, w);
            ph.add(-scale * scatter.value);
            add(sea, sea, -scatter.angle);
            add(sea, j, -Real(sign) * scatter.angle);
            add(sea, sea + 1, Real(side) * d * scatter.width / Real{2});
            mod.add(-scale * Real(side) * log.value);
            add(sea + 1, sea, -Real(side) * log.angle);
            add(sea + 1, j, -Real(side * sign) * log.angle);
            add(sea + 1, sea + 1, d * log.width / Real{2});
          }
      out.residual[sea] = ph.value();
      out.residual[sea + 1] = mod.value();
      for (std::size_t i = 0; i < order; ++i)
      {
        if (!uni20::isfinite(out.residual[i])) throw std::overflow_error("nonfinite two-string residual");
        out.norm = std::max(out.norm, std::abs(out.residual[i]));
        if (i <= sea) out.phase_norm = std::max(out.phase_norm, std::abs(out.residual[i]));
      }
      out.modulus_norm = std::abs(mod.value());
      return out;
    }
    Real energy_sum(std::span<Real const> x, Real delta, Real offset) const
    {
      Real const sh = std::sinh(eta), s2 = sh * sh;
      bethe::detail::CompensatedSum<Real> sum;
      sum.add(offset);
      for (std::size_t i = 0; i < sea; ++i)
        sum.add(-s2 / (delta - std::cos(x[i])));
      Real const d = Real(deviation_sign) * std::exp(-x[sea + 1]), a = x[sea];
      Real const sa = std::sin(a / Real{2}), sd = std::sinh(d / Real{2});
      // delta-cosh(eta+d)*cos(a), without cancellation as a->0 or d->0.
      Real const re = Real{2} * delta * sa * sa - (Real{2} * delta * sd * sd + sh * std::sinh(d)) * std::cos(a);
      Real const im = std::sinh(eta + d) * std::sin(a);
      sum.add(-Real{2} * s2 * re / (re * re + im * im));
      return sum.value();
    }
    Real energy_shift(std::span<Real const> x, Real delta) const { return energy_sum(x, delta, Real{0}); }
    Real energy(std::span<Real const> x, Real delta) const
    {
      return energy_sum(x, delta, Real(sites - 1) * delta / Real{4});
    }
    std::size_t sites, order, sea, string_label;
    int deviation_sign;
    Real eta, pi;
};

template <uni20::Real Real>
State<Real> solve(System<Real> const& system, Real delta, SolverOptions<Real> const& options)
{
  if (!uni20::isfinite(delta) || delta <= Real{1}) throw std::invalid_argument("two-string requires finite Delta>1");
  if (!uni20::isfinite(options.residual_tolerance) || options.residual_tolerance <= Real{0})
    throw std::invalid_argument("residual tolerance must be finite and positive");
  quantum_group::detail::check_matrix_size<Real>(system.order);
  auto x = system.seed();
  if (!system.physical(x)) throw std::overflow_error("two-string seed is not resolvable at the selected precision");
  State<Real> state;
  state.sites = system.sites;
  state.delta = delta;
  state.string_label = system.string_label;
  state.deviation_sign = system.deviation_sign;
  std::vector<Real> jac, step(system.order);
  bool ideal = true;
  for (;;)
  {
    auto evaluation = system.evaluate(x, &jac, ideal);
    if (ideal && evaluation.norm <= std::sqrt(uni20::numeric_limits<Real>::epsilon()))
    {
      ideal = false;
      evaluation = system.evaluate(x, &jac);
    }
    if (!ideal && evaluation.norm <= options.residual_tolerance)
    {
      state.converged = true;
      state.status = SolveStatus::converged;
      break;
    }
    if (state.iterations == options.max_iterations) break;
    for (std::size_t i = 0; i < system.order; ++i)
      step[i] = -evaluation.residual[i];
    // Empty-sea coordinates have very different natural scales: the center
    // can be O(1/N), while L is O(N log N). Equilibrate this 2x2 correction
    // so the common pivot test does not mistake units for rank loss. Keep
    // acceptance and all published residuals in the original equation units.
    Real const center_scale = std::min(x[0], system.pi - x[0]);
    Real const log_scale = Real{2} * Real(system.sites);
    if (system.sea == 0)
      for (std::size_t i = 0; i < 2; ++i)
      {
        jac[2 * i] *= center_scale;
        jac[2 * i + 1] *= log_scale;
        Real const row_scale = std::max(std::abs(jac[2 * i]), std::abs(jac[2 * i + 1]));
        if (row_scale > Real{0})
        {
          jac[2 * i] /= row_scale;
          jac[2 * i + 1] /= row_scale;
          step[i] /= row_scale;
        }
      }
    if (!bethe::detail::newton_step(jac, step))
    {
      state.status = SolveStatus::singular_jacobian;
      break;
    }
    if (system.sea == 0)
    {
      step[0] *= center_scale;
      step[1] *= log_scale;
    }
    auto trial = x;
    Real damping{1};
    bool accepted = false;
    for (int backtrack = 0; backtrack <= uni20::numeric_limits<Real>::digits; ++backtrack)
    {
      for (std::size_t i = 0; i < system.order; ++i)
        trial[i] = x[i] + damping * step[i];
      if (system.physical(trial))
      {
        auto const norm = system.evaluate(trial, nullptr, ideal).norm;
        if (norm < (Real{1} - damping / Real{10000}) * evaluation.norm)
        {
          accepted = true;
          break;
        }
      }
      damping /= Real{2};
    }
    if (!accepted)
    {
      state.status = SolveStatus::stalled;
      break;
    }
    x = std::move(trial);
    ++state.iterations;
  }
  auto const final = system.evaluate(x); // Always finite-deviation diagnostics, including failed initialization.
  state.residual_norm = final.norm;
  state.phase_residual = final.phase_norm;
  state.modulus_residual = final.modulus_norm;
  state.center = x[system.sea];
  state.log_deviation = x.back();
  state.rapidities.assign(x.begin(), x.begin() + system.sea);
  state.energy = system.energy(x, delta);
  state.energy_shift = system.energy_shift(x, delta);
  if (!uni20::isfinite(state.energy)) throw std::overflow_error("two-string energy overflow");
  return state;
}
} // namespace detail

/// One positive-deviation two-string above a real sea, even N>=4, Delta>1.
/// Native precision, O(N^2) memory/O(N^3) per update; no enumeration or site cap.
/// Equations are the exact open-chain Bethe equations regularized in L=-log(d),
/// not an ideal-string approximation. See docs/xxz-open-two-string.md.
template <uni20::Real Real = double>
[[nodiscard]] State<Real> singlet(std::size_t sites, Real delta, SolverOptions<Real> const& options = {})
{
  (void)quantum_group::detail::sector_roots(sites, 0);
  if (sites < 4) throw std::invalid_argument("two-string singlet requires even N>=4");
  return detail::solve(detail::System<Real>(sites, delta), delta, options);
}

/// Target one empty-sea two-string in ell=N-4, odd/even N>=4.
/// mode=1,...,N-3 selects J=N-2-mode, with sign(d)=(-1)^(mode+1).
/// Two unknowns, constant storage/work per Newton step even for long chains.
/// This is NOT a full two-root spectrum or a guarantee of spectral ordering.
/// At other Delta, a chosen branch can collapse or fail to converge.
template <uni20::Real Real = double>
[[nodiscard]] State<Real> bound_pair(std::size_t sites, Real delta, std::size_t mode = 1,
                                     SolverOptions<Real> const& options = {})
{
  if (sites < 4) throw std::invalid_argument("bound pair requires N>=4");
  if (mode == 0 || mode > sites - 3) throw std::invalid_argument("bound-pair mode requires 1<=mode<=N-3");
  auto const label = sites - 2 - mode;
  auto const sign = mode % 2 ? 1 : -1;
  return detail::solve(detail::System<Real>(sites, delta, 2, label, sign), delta, options);
}
} // namespace bethe::xxz::quantum_group::two_string
