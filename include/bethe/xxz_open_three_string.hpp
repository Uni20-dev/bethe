// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/open_string_solver.hpp>
#include <complex>

namespace bethe::xxz::quantum_group::three_string
{
template <uni20::Real Real> struct State
{
    std::size_t sites{}, string_label{};
    /// z=exp(-log_deviation+i*deviation_phase), even when exp(-L) underflows.
    Real delta{}, center{}, log_deviation{}, deviation_phase{};
    /// energy_shift is relative to the polarized reference, evaluated directly.
    Real energy{}, energy_shift{}, residual_norm{}, phase_residual{}, modulus_residual{};
    std::size_t iterations{};
    bool converged = false;
    SolveStatus status = SolveStatus::iteration_limit;
};

namespace detail
{
template <uni20::Real Real> class System {
  public:
    using C = std::complex<Real>;
    static constexpr std::size_t order = 3;
    System(std::size_t sites, Real delta, std::size_t label)
        : sites(sites), string_label(label), delta(delta), eta(std::acosh(delta)), pi(Real{4} * std::atan(Real{1}))
    {}
    static Real wrap(Real a) { return std::atan2(std::sin(a), std::cos(a)); }
    void normalize(std::vector<Real>& x) const { x[2] = wrap(x[2]); }
    C deviation(std::span<Real const> x) const { return std::exp(-x[1]) * C{std::cos(x[2]), std::sin(x[2])}; }
    bool physical(std::span<Real const> x) const
    {
      if (x.size() != order) return false;
      for (Real v : x)
        if (!uni20::isfinite(v)) return false;
      if (!(x[0] > Real{0} && x[0] < pi && x[1] > std::log(Real{4} / eta))) return false;
      Real const b = x[0] + Real{2} * deviation(x).imag();
      return b > Real{0} && b < pi;
    }
    std::vector<Real> coordinate_scales(std::span<Real const> x) const
    {
      // Near a branch edge, a or pi-a is O(1/N). The L/phi equations carry
      // a 1/(2N) coefficient; scale their Newton columns, not the tolerance.
      return {std::min(x[0], pi - x[0]), Real{2} * Real(sites), Real{2} * Real(sites)};
    }
    std::vector<Real> seed() const
    {
      std::vector<Real> x{pi * (Real(string_label) / Real(sites - 4)), Real{0}, Real{0}};
      auto const f = evaluate(x, nullptr, true);
      x[1] = std::max(Real{2} * Real(sites) * f.residual[1], Real{2} + std::log(Real{4} / eta));
      x[2] = wrap(-Real{2} * Real(sites) * f.residual[2]);
      return x;
    }
    // log(sinh(z)/z) and its complex derivative, including exact z=0 after
    // underflow. The omitted Taylor term is O(z^8), below native epsilon.
    static std::pair<C, C> correction(C z)
    {
      Real const cutoff = std::sqrt(std::sqrt(std::sqrt(uni20::numeric_limits<Real>::epsilon())));
      if (std::abs(z) < cutoff)
      {
        C const z2 = z * z;
        return {z2 * (Real{1} / Real{6} + z2 * (-Real{1} / Real{180} + z2 / Real{2835})),
                z * (Real{1} / Real{3} + z2 * (-Real{1} / Real{45} + Real{2} * z2 / Real{945}))};
      }
      return {std::log(std::sinh(z) / z), std::cosh(z) / std::sinh(z) - C{1} / z};
    }
    struct Evaluation
    {
        std::vector<Real> residual;
        Real norm{}, phase_norm{}, modulus_norm{};
    };
    Evaluation evaluate(std::span<Real const> x, std::vector<Real>* jac = nullptr, bool ideal = false) const
    {
      Real const a = x[0], L = x[1], phi = x[2], scale = Real{1} / (Real{2} * Real(sites));
      C const z = ideal ? C{} : deviation(x), u = C{eta, a / Real{2}} + z;
      Real const b = a + Real{2} * z.imag();
      auto log_sinh = [](C w) { return std::log(std::sinh(w)); };
      C const dp = u + eta / Real{2}, dm = u - eta / Real{2};
      C const ap = Real{2} * eta + z, bp = ap + C{0, a}, bm = z + C{0, a};
      C const cp{eta, b}, cm{-eta, b};
      Real const rp = Real{3} * eta + Real{2} * z.real(), rm = eta + Real{2} * z.real();
      C const drive = log_sinh(dp) - log_sinh(dm), A = log_sinh(ap);
      C const B = log_sinh(bp) - log_sinh(bm), Cphase = log_sinh(cp) - log_sinh(cm);
      Real const reflected = std::log(std::sinh(rp)) - std::log(std::sinh(rm));
      auto const [corr, dcorr] = correction(z);
      auto const theta = quantum_group::detail::string_phase(a, eta / Real{2});
      auto const self = quantum_group::detail::string_phase(Real{2} * b, eta);
      // Product phase (singular internal scattering cancelled), then the
      // original u+ modulus and phase. See docs/biquadratic-bound-triples.md.
      // log(sinh(z)) = -L+i*phi+corr is never formed by root subtraction.
      Evaluation result{
          .residual = {
              theta.value + Real{2} * drive.imag() -
                  Real{2} * scale * (Real{2} * B.imag() + self.value + pi * Real(string_label + 1)),
              drive.real() - scale * (A.real() + L - corr.real() + B.real() + reflected),
              wrap(Real{2} * Real(sites) * drive.imag() - A.imag() + phi + corr.imag() - B.imag() - Cphase.imag()) *
                  scale}};
      for (Real value : result.residual)
      {
        if (!uni20::isfinite(value)) throw std::overflow_error("nonfinite three-string residual");
        result.norm = std::max(result.norm, std::abs(value));
      }
      result.phase_norm = std::max(std::abs(result.residual[0]), std::abs(result.residual[2]));
      result.modulus_norm = std::abs(result.residual[1]);
      if (jac)
      {
        jac->assign(order * order, Real{0});
        auto coth = [](C w) { return std::cosh(w) / std::sinh(w); };
        for (std::size_t j = 0; j < order; ++j)
        {
          C const dz = ideal ? C{} : j == 1 ? -z : j == 2 ? C{0, 1} * z : C{};
          C const da{0, j == 0 ? Real{1} : Real{0}}, du = dz + da / Real{2};
          Real const db = (j == 0 ? Real{1} : Real{0}) + Real{2} * dz.imag();
          C const dD = (coth(dp) - coth(dm)) * du, dA = coth(ap) * dz;
          C const dB = (coth(bp) - coth(bm)) * (dz + da);
          C const dC = (coth(cp) - coth(cm)) * C{0, db};
          C const dc = dcorr * dz;
          Real const dR = Real{2} * dz.real() * (Real{1} / std::tanh(rp) - Real{1} / std::tanh(rm));
          (*jac)[j] =
              (j == 0 ? theta.angle : Real{0}) + Real{2} * dD.imag() - Real{4} * scale * (dB.imag() + self.angle * db);
          (*jac)[order + j] =
              dD.real() - scale * (dA.real() + (j == 1 ? Real{1} : Real{0}) - dc.real() + dB.real() + dR);
          (*jac)[2 * order + j] =
              dD.imag() - scale * (dA.imag() - (j == 2 ? Real{1} : Real{0}) - dc.imag() + dB.imag() + dC.imag());
        }
      }
      return result;
    }
    Real energy_shift(std::span<Real const> x) const
    {
      C const u = C{eta, x[0] / Real{2}} + deviation(x);
      Real const sh = std::sinh(eta);
      return -sh * sh *
             (Real{1} / (delta - std::cos(x[0])) + Real{2} * (C{1} / (delta - std::cosh(Real{2} * u))).real());
    }
    std::size_t sites, string_label;
    Real delta, eta, pi;
};
} // namespace detail

/// Selected empty-sea three-string in ell=N-6, odd/even N>=6, Delta>1.
/// z=exp(-L+i*phi), u0=i*a/2, u±=eta+Re(z) +/- i*(a/2+Im(z)).
/// The complex deviation, not rounded differences of roots, is authoritative.
/// mode=1,...,N-5 is a branch label, NOT a full-sector excitation rank.
template <uni20::Real Real = double>
[[nodiscard]] State<Real> bound_triple(std::size_t sites, Real delta, std::size_t mode = 1,
                                       SolverOptions<Real> const& options = {})
{
  if (sites < 6) throw std::invalid_argument("bound triple requires N>=6");
  if (mode == 0 || mode > sites - 5) throw std::invalid_argument("bound-triple mode requires 1<=mode<=N-5");
  if (!uni20::isfinite(delta) || delta <= Real{1}) throw std::invalid_argument("three-string requires finite Delta>1");
  detail::System<Real> system(sites, delta, sites - 4 - mode);
  auto const iteration = quantum_group::detail::solve_log_string<Real>(system, options);
  auto const& x = iteration.x;
  auto const final = system.evaluate(x);
  auto const shift = system.energy_shift(x);
  auto const energy = Real(sites - 1) * delta / Real{4} + shift;
  if (!uni20::isfinite(energy)) throw std::overflow_error("three-string energy overflow");
  return {.sites = sites,
          .string_label = system.string_label,
          .delta = delta,
          .center = x[0],
          .log_deviation = x[1],
          .deviation_phase = x[2],
          .energy = energy,
          .energy_shift = shift,
          .residual_norm = final.norm,
          .phase_residual = final.phase_norm,
          .modulus_residual = final.modulus_norm,
          .iterations = iteration.iterations,
          .converged = iteration.converged,
          .status = iteration.status};
}
} // namespace bethe::xxz::quantum_group::three_string
