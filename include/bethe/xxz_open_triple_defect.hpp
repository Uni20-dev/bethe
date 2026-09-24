// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/open_string_scattering.hpp>
#include <bethe/xxz_open_three_string.hpp>

namespace bethe::xxz::quantum_group::triple_defect::detail
{
template <uni20::Real Real> class System {
  public:
    using Triple = three_string::detail::System<Real>;
    using C = std::complex<Real>;
    static constexpr std::size_t order = 4;
    System(std::size_t n, Real delta, std::size_t real_label, std::size_t string_label)
        : triple(n, delta, string_label), real_label(real_label)
    {}
    bool physical(std::span<Real const> x) const
    {
      return x.size() == order && triple.physical(x.first(3)) && uni20::isfinite(x[3]) && x[3] > Real{0} &&
             x[3] < triple.pi && x[3] != x[0];
    }
    void normalize(std::vector<Real>& x, bool ideal) const
    {
      triple.normalize(x);
      if (!ideal && physical(x) && std::exp(-x[1]) == Real{0})
      {
        // Once z underflows, L and phi enter the exact residuals linearly.
        // Solve those two rows at the *rounded* trial angles. Near coincident
        // centers, a sub-ulp angular Newton step otherwise asks for a large
        // compensating L update although the represented angles never move.
        auto const f = evaluate(x);
        Real const L = x[1] + Real{2} * Real(triple.sites) * f.residual[1];
        if (uni20::isfinite(L) && std::exp(-L) == Real{0})
        {
          x[1] = L;
          x[2] = Triple::wrap(x[2] - Real{2} * Real(triple.sites) * f.residual[2]);
        }
      }
    }
    auto coordinate_scales(std::span<Real const> x) const
    {
      auto s = triple.coordinate_scales(x.first(3));
      s.push_back(std::min(x[3], triple.pi - x[3]));
      return s;
    }
    std::vector<Real> seed() const
    {
      // A linear string-center seed keeps nearby central and sea roots
      // separated better than independently inverting both bare phases.
      auto angle = [&](Real width, std::size_t label, std::size_t denominator) {
        return Real{2} *
               std::atan(std::tanh(width) * std::tan(triple.pi * Real(label) / (Real{2} * Real(denominator))));
      };
      std::vector<Real> x{triple.pi * Real(triple.string_label) / Real(triple.sites - 6), Real{0}, Real{0},
                          angle(triple.eta / Real{2}, real_label, triple.sites - 2)};
      auto const f = evaluate(x, nullptr, true);
      x[1] = std::max(Real{2} * Real(triple.sites) * f.residual[1], Real{2} + std::log(Real{4} / triple.eta));
      x[2] = Triple::wrap(-Real{2} * Real(triple.sites) * f.residual[2]);
      return x;
    }
    using Evaluation = typename Triple::Evaluation;
    Evaluation evaluate(std::span<Real const> x, std::vector<Real>* jac = nullptr, bool ideal = false) const
    {
      std::vector<Real> block;
      auto out = triple.evaluate(x.first(3), jac ? &block : nullptr, ideal);
      Real const scale = Real{1} / (Real{2} * Real(triple.sites)), eta = triple.eta;
      if (ideal)
      {
        // Initialize centers with the analytically fused phases, and L/phi
        // with the isolated triple. The individual external factors have
        // removable poles at alpha=a when z=0; they must not obstruct this
        // initializer. The final stage always uses the full finite equations.
        if (jac)
        {
          jac->assign(16, Real{0});
          for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j)
              (*jac)[4 * i + j] = block[3 * i + j];
        }
        auto const drive = quantum_group::detail::string_phase(x[3], eta / Real{2});
        out.residual.push_back(drive.value - Real{2} * scale * triple.pi * Real(real_label));
        if (jac) (*jac)[15] = drive.angle;
        for (int sign : {-1, 1})
          for (int width : {1, 2})
          {
            auto const t = quantum_group::detail::string_phase(x[0] + Real(sign) * x[3], Real(width) * eta);
            auto const r = quantum_group::detail::string_phase(x[3] + Real(sign) * x[0], Real(width) * eta);
            out.residual[0] -= scale * t.value;
            out.residual[3] -= scale * r.value;
            if (jac)
            {
              (*jac)[0] -= scale * t.angle;
              (*jac)[3] -= scale * Real(sign) * t.angle;
              (*jac)[12] -= scale * Real(sign) * r.angle;
              (*jac)[15] -= scale * r.angle;
            }
          }
        out.phase_norm = std::max({std::abs(out.residual[0]), std::abs(out.residual[2]), std::abs(out.residual[3])});
        out.norm = std::max(out.phase_norm, out.modulus_norm);
        return out;
      }
      C const z = ideal ? C{} : triple.deviation(x.first(3));
      C const u0{0, x[0] / Real{2}}, u = C{eta, x[0] / Real{2}} + z, v{0, x[3] / Real{2}};
      auto const t0 = quantum_group::detail::open_scattering(u0, v, eta);
      auto const tp = quantum_group::detail::open_scattering(u, v, eta);
      auto const r0 = quantum_group::detail::open_scattering(v, u0, eta);
      auto const rp = quantum_group::detail::open_scattering(v, u, eta);
      auto const rm = quantum_group::detail::open_scattering(v, std::conj(u), eta);
      auto fused = [&](Real a, Real b) {
        Real result{0};
        for (int sign : {-1, 1})
          for (int width : {1, 2})
            result += quantum_group::detail::string_phase(a + Real(sign) * b, Real(width) * eta).value;
        return result;
      };
      // The ideal fused phase fixes the 2pi branch, not the finite-string
      // value. Its derivative cancels inside/outside wrap locally, leaving
      // the derivatives of the original nonsingular scattering factors.
      auto continuous = [&](Real raw, Real target) { return target + Triple::wrap(raw - target); };
      out.residual[0] -= scale * continuous(t0.value.imag() + Real{2} * tp.value.imag(), fused(x[0], x[3]));
      out.residual[1] -= scale * tp.value.real();
      out.residual[2] = scale * Triple::wrap(out.residual[2] / scale - tp.value.imag());
      auto const theta = quantum_group::detail::string_phase(x[3], eta / Real{2});
      out.residual.push_back(theta.value -
                             scale * (Real{2} * triple.pi * Real(real_label) +
                                      continuous((r0.value + rp.value + rm.value).imag(), fused(x[3], x[0]))));
      if (jac)
      {
        jac->assign(16, Real{0});
        for (std::size_t i = 0; i < 3; ++i)
          for (std::size_t j = 0; j < 3; ++j)
            (*jac)[4 * i + j] = block[3 * i + j];
        for (std::size_t j = 0; j < 4; ++j)
        {
          C const dz = ideal ? C{} : j == 1 ? -z : j == 2 ? C{0, 1} * z : C{};
          C const du0{0, j == 0 ? Real{1} / Real{2} : Real{0}}, du = du0 + dz;
          C const dv{0, j == 3 ? Real{1} / Real{2} : Real{0}};
          C const dt0 = t0.first * du0 + t0.second * dv, dtp = tp.first * du + tp.second * dv;
          C const dr = r0.first * dv + r0.second * du0 + rp.first * dv + rp.second * du + rm.first * dv +
                       rm.second * std::conj(du);
          (*jac)[j] -= scale * (dt0 + Real{2} * dtp).imag();
          (*jac)[4 + j] -= scale * dtp.real();
          (*jac)[8 + j] -= scale * dtp.imag();
          (*jac)[12 + j] = (j == 3 ? theta.angle : Real{0}) - scale * dr.imag();
        }
      }
      out.norm = out.phase_norm = out.modulus_norm = Real{0};
      for (std::size_t i = 0; i < 4; ++i)
      {
        if (!uni20::isfinite(out.residual[i])) throw std::overflow_error("nonfinite triple-defect residual");
        out.norm = std::max(out.norm, std::abs(out.residual[i]));
        auto& norm = i == 1 ? out.modulus_norm : out.phase_norm;
        norm = std::max(norm, std::abs(out.residual[i]));
      }
      return out;
    }
    Real energy_shift(std::span<Real const> x) const
    {
      Real const sh = std::sinh(triple.eta);
      return triple.energy_shift(x.first(3)) - sh * sh / (triple.delta - std::cos(x[3]));
    }
    Triple triple;
    std::size_t real_label;
};
} // namespace bethe::xxz::quantum_group::triple_defect::detail

namespace bethe::xxz::quantum_group::triple_defect
{
template <uni20::Real Real> struct State : three_string::State<Real>
{
    std::size_t real_label{};
    Real rapidity{};
};

/// One three-string and one real root; selected branch, not a complete module.
/// 1<=I<=N-3, 1<=J<=N-7, N>=8, Delta>1. Labels are not energy ranks.
template <uni20::Real Real = double>
[[nodiscard]] State<Real> solve(std::size_t sites, Real delta, std::size_t real_label, std::size_t string_label,
                                SolverOptions<Real> const& options = {})
{
  (void)xxz::detail::checked_sites(sites);
  if (sites < 8) throw std::invalid_argument("triple plus defect requires N>=8");
  if (!uni20::isfinite(delta) || delta <= Real{1})
    throw std::invalid_argument("triple plus defect requires finite Delta>1");
  if (real_label == 0 || real_label > sites - 3 || string_label == 0 || string_label > sites - 7)
    throw std::invalid_argument("triple-defect labels require 1<=I<=N-3 and 1<=J<=N-7");
  detail::System<Real> system(sites, delta, real_label, string_label);
  auto const it = quantum_group::detail::solve_log_string<Real>(system, options);
  auto const f = system.evaluate(it.x);
  State<Real> s;
  s.sites = sites;
  s.delta = delta;
  s.real_label = real_label;
  s.string_label = string_label;
  s.center = it.x[0];
  s.log_deviation = it.x[1];
  s.deviation_phase = it.x[2];
  s.rapidity = it.x[3];
  s.energy_shift = system.energy_shift(it.x);
  s.energy = Real(sites - 1) * delta / Real{4} + s.energy_shift;
  if (!uni20::isfinite(s.energy) || !uni20::isfinite(s.energy_shift))
    throw std::overflow_error("triple-defect energy overflow");
  s.residual_norm = f.norm;
  s.phase_residual = f.phase_norm;
  s.modulus_residual = f.modulus_norm;
  s.iterations = it.iterations;
  s.converged = it.converged;
  s.status = it.status;
  return s;
}
} // namespace bethe::xxz::quantum_group::triple_defect
