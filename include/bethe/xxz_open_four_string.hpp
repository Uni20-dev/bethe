// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/xxz_open_three_string.hpp>
#include <bethe/xxz_open_two_string.hpp>

namespace bethe::xxz::quantum_group::four_string::detail
{
template <uni20::Real Real> class System {
  public:
    using C = uni20::complex<Real>;
    using Pair = two_string::detail::System<Real>;
    using Triple = three_string::detail::System<Real>;
    static constexpr std::size_t order = 4;
    System(std::size_t n, Real delta, std::size_t label)
        : sites(n), string_label(label), sign((n - label) % 2 ? 1 : -1), delta(delta),
          zero(bethe::detail::precision_zero(delta)), eta(uni20::math::acosh(delta)),
          pi(Real{4} * uni20::math::atan(uni20::scalar_like(delta, 1))), pair(n, delta, 2, label, sign)
    {}
    C deviation(std::span<Real const> x) const
    {
      return uni20::math::exp(-x[2]) * C{uni20::math::cos(x[3]), uni20::math::sin(x[3])};
    }
    void normalize(std::vector<Real>& x) const { x[3] = Triple::wrap(x[3]); }
    bool physical(std::span<Real const> x) const
    {
      if (x.size() != order) return false;
      for (Real v : x)
        if (!uni20::isfinite(v)) return false;
      if (!(x[0] > Real{0} && x[0] < pi && x[1] > uni20::math::log(Real{4} / eta) &&
            x[2] > uni20::math::log(Real{4} / eta)))
        return false;
      Real const b = x[0] + Real{2} * deviation(x).imag();
      return b > Real{0} && b < pi;
    }
    std::vector<Real> coordinate_scales(std::span<Real const> x) const
    {
      return {std::min(x[0], pi - x[0]), uni20::scalar_like(zero, 2) * Real(sites),
              uni20::scalar_like(zero, 2) * Real(sites), uni20::scalar_like(zero, 2) * Real(sites)};
    }
    std::vector<Real> seed() const
    {
      Real const L = Real{8} + uni20::math::abs(uni20::math::log(eta));
      std::vector<Real> x{pi * Real(string_label) / Real(sites - 6), L, L, zero};
      if (!physical(x)) throw std::overflow_error("four-string seed is not resolvable at the selected precision");
      auto const f = evaluate(x, nullptr, true);
      for (std::size_t j : {1, 2})
        x[j] = std::max(x[j] + uni20::scalar_like(zero, 2) * Real(sites) * f.residual[j],
                        Real{2} + uni20::math::log(Real{4} / eta));
      x[3] = Triple::wrap(-uni20::scalar_like(zero, 2) * Real(sites) * f.residual[3]);
      return x;
    }
    using Evaluation = typename Triple::Evaluation;
    Evaluation evaluate(std::span<Real const> x, std::vector<Real>* jac = nullptr, bool ideal = false) const
    {
      Real const a = x[0], d = ideal ? zero : Real(sign) * uni20::math::exp(-x[1]),
                 scale = uni20::scalar_like(zero, 1) / (uni20::scalar_like(zero, 2) * Real(sites));
      C const z = ideal ? C{zero, zero} : deviation(x), u{(eta + d) / Real{2}, a / Real{2}}, w = u + eta + z;
      auto log_sinh = [](C v) { return uni20::math::log(uni20::math::sinh(v)); };
      auto coth = [](C v) { return uni20::math::cosh(v) / uni20::math::sinh(v); };
      C const up = u + eta / Real{2}, um = u - eta / Real{2}, wp = w + eta / Real{2}, wm = w - eta / Real{2};
      C const av = Real{2} * eta + z;
      C const bp = Real{3} * eta + d + C{zero, a} + z, bm = eta + d + C{zero, a} + z;
      C const kp = Real{2} * eta + C{zero, a} + z, km = C{zero, a} + z;
      C const rp = Real{3} * eta + d + z, rm = eta + d + z;
      C const D0 = log_sinh(up) - log_sinh(um), D1 = log_sinh(wp) - log_sinh(wm), D = D0 + D1;
      C const A = log_sinh(av), B = log_sinh(bp) - log_sinh(bm), K = log_sinh(kp) - log_sinh(km),
              R = log_sinh(rp) - log_sinh(rm);
      Real const q = uni20::math::log(uni20::math::sinh(Real{2} * eta + d));
      Real const qp = Real{4} * eta + d + Real{2} * z.real(), qm = Real{2} * eta + d + Real{2} * z.real();
      Real const Q = uni20::math::log(uni20::math::sinh(qp)) - uni20::math::log(uni20::math::sinh(qm));
      auto const [cd, dcd] = Triple::correction(C{d});
      auto const [cz, dcz] = Triple::correction(z);
      auto const t0 = quantum_group::detail::string_phase(Real{2} * a, eta);
      auto const t1 = quantum_group::detail::string_phase(Real{2} * (a + Real{2} * z.imag()), eta);
      Evaluation out{
          .residual = {Real{2} * D.imag() + pi -
                           Real{2} * scale *
                               (t0.value + t1.value + Real{2} * B.imag() + Real{2} * K.imag() +
                                pi * (Real(string_label) + Real{1})),
                       D.real() - scale * (q + x[1] - cd.real() + Q + Real{2} * B.real() + Real{2} * R.real()),
                       D1.real() - scale * (A.real() + x[2] - cz.real() + B.real() + K.real() + R.real() + Q),
                       scale * Triple::wrap(D1.imag() / scale - A.imag() + x[3] + cz.imag() - B.imag() - K.imag() -
                                            R.imag() - t1.value + pi)},
          .norm = zero,
          .phase_norm = zero,
          .modulus_norm = zero};
      if (jac)
      {
        jac->assign(16, zero);
        for (std::size_t j = 0; j < 4; ++j)
        {
          Real const da = j == 0 ? Real{1} : Real{0}, dd = !ideal && j == 1 ? -d : Real{0};
          C const dz = ideal ? C{zero, zero} : j == 2 ? -z : j == 3 ? C{zero, Real{1}} * z : C{zero, zero};
          C const du{dd / Real{2}, da / Real{2}}, dw = du + dz;
          C const dD0 = (coth(up) - coth(um)) * du, dD1 = (coth(wp) - coth(wm)) * dw, dD = dD0 + dD1;
          C const dA = coth(av) * dz, dB = (coth(bp) - coth(bm)) * (C{dd, da} + dz);
          C const dK = (coth(kp) - coth(km)) * (C{zero, da} + dz), dR = (coth(rp) - coth(rm)) * (dd + dz);
          Real const dq = dd / uni20::math::tanh(Real{2} * eta + d),
                     dQ = (dd + Real{2} * dz.real()) *
                          (Real{1} / uni20::math::tanh(qp) - Real{1} / uni20::math::tanh(qm));
          Real const dt0 = Real{2} * t0.angle * da, dt1 = Real{2} * t1.angle * (da + Real{2} * dz.imag());
          C const dc = dcz * dz;
          (*jac)[j] = Real{2} * dD.imag() - Real{2} * scale * (dt0 + dt1 + Real{2} * dB.imag() + Real{2} * dK.imag());
          (*jac)[4 + j] = dD.real() - scale * (dq + (j == 1 ? Real{1} : Real{0}) - dcd.real() * dd + dQ +
                                               Real{2} * dB.real() + Real{2} * dR.real());
          (*jac)[8 + j] = dD1.real() - scale * (dA.real() + (j == 2 ? Real{1} : Real{0}) - dc.real() + dB.real() +
                                                dK.real() + dR.real() + dQ);
          (*jac)[12 + j] = dD1.imag() - scale * (dA.imag() - (j == 3 ? Real{1} : Real{0}) - dc.imag() + dB.imag() +
                                                 dK.imag() + dR.imag() + dt1);
        }
      }
      for (std::size_t i = 0; i < 4; ++i)
      {
        if (!uni20::isfinite(out.residual[i])) throw std::overflow_error("nonfinite four-string residual");
        out.norm = std::max(out.norm, uni20::math::abs(out.residual[i]));
        auto& norm = i == 0 || i == 3 ? out.phase_norm : out.modulus_norm;
        norm = std::max(norm, uni20::math::abs(out.residual[i]));
      }
      return out;
    }
    Real energy_shift(std::span<Real const> x) const
    {
      Real const d = Real(sign) * uni20::math::exp(-x[1]);
      C const w = C{(Real{3} * eta + d) / Real{2}, x[0] / Real{2}} + deviation(x);
      Real const sh = uni20::math::sinh(eta);
      return pair.energy_shift(x.first(2), delta) -
             Real{2} * sh * sh * (C{1} / (delta - uni20::math::cosh(Real{2} * w))).real();
    }
    std::size_t sites, string_label;
    int sign;
    Real delta, zero, eta, pi;
    Pair pair;
};
} // namespace bethe::xxz::quantum_group::four_string::detail

namespace bethe::xxz::quantum_group::four_string
{
template <uni20::Real Real> struct State
{
    std::size_t sites{}, string_label{};
    Real delta{}, center{}, inner_log_deviation{}, outer_log_deviation{}, outer_deviation_phase{};
    int inner_deviation_sign{};
    Real energy{}, energy_shift{}, residual_norm{}, phase_residual{}, modulus_residual{};
    std::size_t iterations{};
    bool converged = false;
    SolveStatus status = SolveStatus::iteration_limit;
};

/// Selected four-string droplet, ell=N-8, odd/even N>=8, Delta>1.
/// mode=1,...,N-7 labels a branch, not a complete-module excitation rank.
/// u=(eta+d+i*a)/2, w=u+eta+z, together with their conjugates;
/// d=sign*exp(-L_inner), z=exp(-L_outer+i*phi_outer).
template <uni20::Real Real = double>
[[nodiscard]] State<Real> bound_quartet(std::size_t sites, Real delta, std::size_t mode = 1,
                                        SolverOptions<Real> const& options = {})
{
  (void)xxz::detail::checked_sites(sites);
  if (sites < 8) throw std::invalid_argument("bound quartet requires N>=8");
  if (mode == 0 || mode > sites - 7) throw std::invalid_argument("bound-quartet mode requires 1<=mode<=N-7");
  (void)xxz::detail::solver_zero(delta, options);
  if (!uni20::isfinite(delta) || delta <= Real{1}) throw std::invalid_argument("four-string requires finite Delta>1");
  detail::System<Real> system(sites, delta, sites - 6 - mode);
  auto const it = quantum_group::detail::solve_log_string<Real>(system, options);
  auto const f = system.evaluate(it.x);
  Real const shift = system.energy_shift(it.x), energy = Real(sites - 1) * delta / Real{4} + shift;
  if (!uni20::isfinite(shift) || !uni20::isfinite(energy)) throw std::overflow_error("four-string energy overflow");
  return {.sites = sites,
          .string_label = system.string_label,
          .delta = delta,
          .center = it.x[0],
          .inner_log_deviation = it.x[1],
          .outer_log_deviation = it.x[2],
          .outer_deviation_phase = it.x[3],
          .inner_deviation_sign = system.sign,
          .energy = energy,
          .energy_shift = shift,
          .residual_norm = f.norm,
          .phase_residual = f.phase_norm,
          .modulus_residual = f.modulus_norm,
          .iterations = it.iterations,
          .converged = it.converged,
          .status = it.status};
}
} // namespace bethe::xxz::quantum_group::four_string
