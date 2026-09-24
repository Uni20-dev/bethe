// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <array>
#include <bethe/detail/open_string_solver.hpp>
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
    /// Integer labels for the real roots; empty for an isolated bound pair.
    std::vector<std::size_t> real_labels;
};

namespace detail
{
template <uni20::Real Real> class System {
  public:
    System(std::size_t n, Real delta) : System(n, delta, n / 2, 1, 1) {}
    System(std::size_t n, Real delta, std::size_t roots, std::size_t label, int sign)
        : sites(n), order(roots), sea(order - 2), string_label(label), deviation_sign(sign), eta(std::acosh(delta)),
          pi(Real{4} * std::atan(Real{1}))
    {
      quantum_group::detail::check_matrix_size<Real>(order);
      real_labels.resize(sea);
      for (std::size_t i = 0; i < sea; ++i)
        real_labels[i] = i + 1;
    }
    System(std::size_t n, Real delta, std::size_t real_label, std::size_t pair_label)
        : System(n, delta, std::span<std::size_t const>{std::array<std::size_t, 1>{real_label}}, pair_label)
    {}
    System(std::size_t n, Real delta, std::span<std::size_t const> labels, std::size_t pair_label)
        : System(n, delta, labels.size() + 2, pair_label, (n - pair_label) % 2 ? 1 : -1)
    {
      // The product phase discards the sign of the singular reflected factor;
      // the unsquared pair equation requires sign(d)=(-1)^(N-J-1).
      real_labels.assign(labels.begin(), labels.end());
      selected_labels = true;
    }
    using Function = quantum_group::detail::StringPhase<Real>;
    // Theta(beta;w) and its two derivatives, without coth(w) near w=0.
    static Function phase(Real beta, Real w, bool complement = false)
    {
      return quantum_group::detail::string_phase(beta, w, complement);
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
        x[i] = Real{2} *
               std::atan(std::tanh(eta / Real{2}) * std::tan(pi * Real(real_labels[i]) / (Real{2} * Real(sites))));
      x[sea] = selected_labels ? pi * (Real(string_label) / Real(sites - 2 - 2 * sea))
               : sea           ? pi / Real{2}
                               : pi * (Real(string_label) / Real(sites - 2));
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
        f.add(drive.value - Real{2} * pi * Real(real_labels[i]) * scale);
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
    std::vector<Real> coordinate_scales(std::span<Real const> x) const
    {
      if (selected_labels)
      {
        std::vector<Real> scales(order);
        for (std::size_t i = 0; i <= sea; ++i)
          scales[i] = std::min(x[i], pi - x[i]);
        scales.back() = Real{2} * Real(sites);
        return scales;
      }
      if (sea) return {};
      return {std::min(x[0], pi - x[0]), Real{2} * Real(sites)};
    }
    std::size_t sites, order, sea, string_label;
    int deviation_sign;
    Real eta, pi;
    std::vector<std::size_t> real_labels;
    bool selected_labels = false;
};

template <uni20::Real Real>
State<Real> solve(System<Real> const& system, Real delta, SolverOptions<Real> const& options)
{
  if (!uni20::isfinite(delta) || delta <= Real{1}) throw std::invalid_argument("two-string requires finite Delta>1");
  auto const iteration = quantum_group::detail::solve_log_string<Real>(system, options);
  auto const& x = iteration.x;
  State<Real> state;
  state.sites = system.sites;
  state.delta = delta;
  state.string_label = system.string_label;
  state.deviation_sign = system.deviation_sign;
  state.real_labels = system.real_labels;
  state.iterations = iteration.iterations;
  state.converged = iteration.converged;
  state.status = iteration.status;
  auto const final = system.evaluate(x); // Always finite-deviation diagnostics, including failed initialization.
  state.residual_norm = final.norm;
  state.phase_residual = final.phase_norm;
  state.modulus_residual = final.modulus_norm;
  state.center = x[system.sea];
  state.log_deviation = x.back();
  state.rapidities.assign(x.begin(), x.begin() + system.sea);
  state.energy = system.energy(x, delta);
  state.energy_shift = system.energy_shift(x, delta);
  if (!uni20::isfinite(state.energy) || !uni20::isfinite(state.energy_shift))
    throw std::overflow_error("two-string energy overflow");
  return state;
}
} // namespace detail

/// One signed two-string plus selected real roots, M=labels.size()+2<=N/2.
/// Ordered distinct integer real labels in [1,N-M], pair label in [1,N-2M+1].
/// The string topology is an ansatz: convergence is not guaranteed at all Delta.
template <uni20::Real Real = double>
[[nodiscard]] State<Real> pair_with_real_roots(std::size_t sites, Real delta, std::span<std::size_t const> labels,
                                               std::size_t pair_label, SolverOptions<Real> const& options = {})
{
  (void)xxz::detail::checked_sites(sites);
  if (sites < 4 || labels.size() > sites / 2 - 2)
    throw std::invalid_argument("one pair plus real roots requires 2<=M<=N/2");
  auto const m = labels.size() + 2;
  quantum_group::detail::check_matrix_size<Real>(m);
  if (pair_label == 0 || pair_label > sites - 2 * m + 1)
    throw std::invalid_argument("pair label requires 1<=J<=N-2M+1");
  for (std::size_t i = 0; i < labels.size(); ++i)
    if (labels[i] == 0 || labels[i] > sites - m || (i && labels[i] <= labels[i - 1]))
      throw std::invalid_argument("real labels must be ordered distinct integers in [1,N-M]");
  return detail::solve(detail::System<Real>(sites, delta, labels, pair_label), delta, options);
}

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

/// Selected two-string plus one real root, ell=N-6, odd/even N>=6.
/// 1<=real_label<=N-3, 1<=pair_label<=N-5. These are logarithmic Bethe labels,
/// not energy ranks or momenta. Failure/collapse at other Delta is reported by
/// the solver; no full-spectrum completeness or convergence guarantee.
template <uni20::Real Real = double>
[[nodiscard]] State<Real> pair_defect(std::size_t sites, Real delta, std::size_t real_label, std::size_t pair_label,
                                      SolverOptions<Real> const& options = {})
{
  (void)xxz::detail::checked_sites(sites);
  if (sites < 6) throw std::invalid_argument("pair plus defect requires N>=6");
  if (real_label == 0 || real_label > sites - 3)
    throw std::invalid_argument("pair-defect real label requires 1<=I<=N-3");
  if (pair_label == 0 || pair_label > sites - 5)
    throw std::invalid_argument("pair-defect string label requires 1<=J<=N-5");
  return pair_with_real_roots<Real>(sites, delta, std::array<std::size_t, 1>{real_label}, pair_label, options);
}
} // namespace bethe::xxz::quantum_group::two_string
