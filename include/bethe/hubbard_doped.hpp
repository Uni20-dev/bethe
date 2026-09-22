// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/hubbard_kernel.hpp>
#include <bethe/detail/newton.hpp>
#include <bethe/hubbard_thermo.hpp>

namespace bethe::hubbard::thermo
{
enum class DopedBranch
{
  spinon,
  holon,
  charge_particle
};
enum class EnergyReference
{
  hamiltonian,
  fermi
};
enum class DopedStatus
{
  converged,
  mesh_limit,
  density_limit,
  linear_failure,
  momentum_limit,
  precision_limit
};
template <uni20::Real Real> struct DopedOptions
{
    Real tolerance = Real{4096} * uni20::numeric_limits<Real>::epsilon();
    std::size_t initial_nodes = 16, max_nodes = 256, max_background_iterations = 64, max_iterations = 160;
};
template <uni20::Real Real> struct DopedBackground
{
    Real interaction{}, density{};
    std::optional<Real> fermi_rapidity, mu_unshifted, mu_symmetric, energy_per_site_unshifted;
    Real mesh_error{}, density_error{};
    std::size_t nodes = 0, iterations = 0;
    bool converged = false;
    DopedStatus status = DopedStatus::mesh_limit;
};
template <uni20::Real Real> struct DopedPoint
{
    DopedBranch branch = DopedBranch::spinon;
    Convention convention = Convention::symmetric;
    EnergyReference reference = EnergyReference::hamiltonian;
    Real momentum{}, parameter{}, energy_error{}, momentum_error{};
    std::optional<Real> energy, symmetric_energy, fermi_energy;
    int delta_particles = 0;
    uni20::half_int spin{0};
    std::size_t iterations = 0;
    bool converged = false;
    DopedStatus status = DopedStatus::mesh_limit;
};
namespace detail
{
template <uni20::Real Real> struct DopedMesh
{
    Real u{}, q{}, density{}, mu{}, e0{}, density_derivative{};
    std::vector<Real> k, w, sine, cosine, rho, epsilon;

    // Returns (epsilon_c(k), p_c(k)); p_c grows from -pi to pi.
    std::array<Real, 2> charge(Real parameter) const
    {
      Real const sk = std::sin(parameter), p = pi<Real>();
      bethe::detail::CompensatedSum<Real> energy, momentum;
      for (std::size_t j = 0; j < k.size(); ++j)
      {
        Real const r = bethe::detail::hubbard_r(sk - sine[j], u) + bethe::detail::hubbard_r(sk + sine[j], u);
        energy.add(w[j] * cosine[j] * r * epsilon[j]);
        Real const phase =
            bethe::detail::hubbard_r_primitive(sk - sine[j], u) + bethe::detail::hubbard_r_primitive(sk + sine[j], u);
        momentum.add(w[j] * rho[j] * phase);
      }
      Real e = -Real{2} * std::cos(parameter) - mu + energy.value();
      Real pc = parameter + Real{2} * p * momentum.value();
      if (std::abs(parameter) == q)
      {
        e = Real{0};
        pc = std::copysign(p * density, parameter);
      }
      if (std::abs(parameter) == p) pc = parameter;
      return {e, pc};
    }
    // Hole momentum convention: p_s in [0,pi*n], matching the half-filled path.
    std::array<Real, 2> spinon(Real lambda) const
    {
      Real const l = std::abs(lambda), scale = Real{2} * u / pi<Real>();
      bethe::detail::CompensatedSum<Real> energy, momentum;
      for (std::size_t j = 0; j < k.size(); ++j)
      {
        Real const minus = (l - sine[j]) / scale, plus = (l + sine[j]) / scale;
        energy.add(-w[j] * cosine[j] * epsilon[j] * (sech(minus) + sech(plus)) / (Real{4} * u));
        momentum.add(Real{2} * w[j] * rho[j] * (atan_exp(-minus) + atan_exp(-plus)));
      }
      Real ps = momentum.value();
      if (l == Real{0})
        ps = pi<Real>() * density / Real{2};
      else if (lambda < Real{0})
        ps = pi<Real>() * density - ps;
      return {energy.value(), ps};
    }
};

template <uni20::Real Real>
std::optional<DopedMesh<Real>> doped_mesh(Real u, Real q, bethe::detail::GaussLegendreRule<Real> const& rule)
{
  std::size_t const n = rule.x.size();
  DopedMesh<Real> out;
  out.u = u;
  out.q = q;
  out.k.resize(n);
  out.w.resize(n);
  out.sine.resize(n);
  out.cosine.resize(n);
  out.rho.assign(n, Real{1} / (Real{2} * pi<Real>()));
  std::vector<Real> arho(n * n), ae(n * n), h(n), xi(n, Real{1});
  for (std::size_t j = 0; j < n; ++j)
  {
    out.k[j] = q * (rule.x[j] + Real{1}) / Real{2};
    out.w[j] = q * rule.w[j] / Real{2};
    out.sine[j] = std::sin(out.k[j]);
    out.cosine[j] = std::cos(out.k[j]);
    h[j] = -Real{2} * out.cosine[j];
  }
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
    {
      Real const kernel = bethe::detail::hubbard_r(out.sine[i] - out.sine[j], u) +
                          bethe::detail::hubbard_r(out.sine[i] + out.sine[j], u);
      arho[i * n + j] = Real(i == j) - out.cosine[i] * out.w[j] * kernel;
      ae[i * n + j] = Real(i == j) - out.cosine[j] * out.w[j] * kernel;
    }
  // Existing recoverable native-real pivoting solver: no process-global error
  // policy changes, including for native long double without BLAS support.
  if (!bethe::detail::newton_step(arho, out.rho) || !bethe::detail::newton_step(ae, h) ||
      !bethe::detail::newton_step(ae, xi))
    return std::nullopt;
  bethe::detail::CompensatedSum<Real> number, e0, hq, xiq, rhoq;
  hq.add(-Real{2} * std::cos(q));
  xiq.add(Real{1});
  rhoq.add(Real{1} / (Real{2} * pi<Real>()));
  for (std::size_t j = 0; j < n; ++j)
  {
    if (!(out.rho[j] > Real{0})) return std::nullopt;
    number.add(Real{2} * out.w[j] * out.rho[j]);
    e0.add(-Real{4} * out.w[j] * out.cosine[j] * out.rho[j]);
    Real const r =
        bethe::detail::hubbard_r(std::sin(q) - out.sine[j], u) + bethe::detail::hubbard_r(std::sin(q) + out.sine[j], u);
    hq.add(out.w[j] * out.cosine[j] * r * h[j]);
    xiq.add(out.w[j] * out.cosine[j] * r * xi[j]);
    rhoq.add(std::cos(q) * out.w[j] * r * out.rho[j]);
  }
  out.mu = hq.value() / xiq.value();
  out.density = number.value();
  out.e0 = e0.value();
  out.density_derivative = Real{2} * rhoq.value() * xiq.value();
  out.epsilon.resize(n);
  for (std::size_t j = 0; j < n; ++j)
    out.epsilon[j] = h[j] - out.mu * xi[j];
  if (!uni20::isfinite(out.mu) || !(out.density_derivative > Real{0})) return std::nullopt;
  return out;
}
} // namespace detail

/// Zero-field repulsive background below half filling. Mesh refinement is
/// native precision; mesh differences are estimates, not interval bounds.
template <uni20::Real Real> class DopedSolver {
  public:
    DopedSolver(Real interaction, Real density, DopedOptions<Real> options = {}) : options_(options)
    {
      info_.interaction = interaction;
      info_.density = density;
      if (!uni20::isfinite(interaction) || interaction <= Real{0} || !uni20::isfinite(density) || density <= Real{0} ||
          density >= Real{1})
        throw std::invalid_argument("doped Hubbard dispersions require finite U>0 and 0<density<1, at zero field");
      if (!uni20::isfinite(options.tolerance) || options.tolerance <= Real{0} || options.tolerance >= Real{1} ||
          options.initial_nodes < 4 || options.max_nodes > 512 || options.initial_nodes > options.max_nodes)
        throw std::invalid_argument(
            "invalid doped tolerance/mesh: require 0<tol<1 and 4<=initial_nodes<=max_nodes<=512");
      Real const u = interaction / Real{4}, p = detail::pi<Real>();
      if (!(u > Real{0}) || !uni20::isfinite(Real{1} / u))
      {
        info_.status = DopedStatus::precision_limit;
        return;
      }
      Real q = p * density * (Real{1} + interaction / (interaction + Real{4})) / Real{2};
      for (std::size_t n = options.initial_nodes;; n *= 2)
      {
        auto const rule = bethe::detail::gauss_legendre<Real>(n);
        Real lo{0}, hi = p;
        std::optional<detail::DopedMesh<Real>> candidate;
        for (;;)
        {
          if (info_.iterations == options.max_background_iterations)
          {
            info_.status = DopedStatus::density_limit;
            return;
          }
          ++info_.iterations;
          candidate = detail::doped_mesh(u, q, rule);
          if (!candidate)
          {
            info_.status = DopedStatus::linear_failure;
            return;
          }
          Real const residual = candidate->density - density;
          if (std::abs(residual) <= options.tolerance * density / Real{64}) break;
          if (residual > Real{0})
            hi = q;
          else
            lo = q;
          Real next = q - residual / candidate->density_derivative;
          if (!(next > lo && next < hi)) next = (lo + hi) / Real{2};
          if (next == q)
          {
            info_.status = DopedStatus::precision_limit;
            return;
          }
          q = next;
        }
        info_.nodes = n;
        info_.density_error = std::abs(candidate->density - density);
        if (fine_)
        {
          Real error = std::max({std::abs(candidate->q - fine_->q), std::abs(candidate->mu - fine_->mu),
                                 std::abs(candidate->e0 - fine_->e0)});
          for (int j = 0; j <= 8; ++j)
          {
            Real const k = p * Real(j) / Real{8};
            auto const a = candidate->charge(k), b = fine_->charge(k);
            for (int v = 0; v < 2; ++v)
              error = std::max(error, std::abs(a[v] - b[v]));
            Real const lambda = (Real{1} + u) * Real(j) / Real{4};
            auto const sa = candidate->spinon(lambda), sb = fine_->spinon(lambda);
            for (int v = 0; v < 2; ++v)
              error = std::max(error, std::abs(sa[v] - sb[v]));
          }
          info_.mesh_error = Real{8} * error;
          coarse_ = std::move(fine_);
        }
        fine_ = std::move(candidate);
        if (coarse_ && info_.mesh_error <= options.tolerance)
        {
          info_.converged = true;
          info_.status = DopedStatus::converged;
          info_.fermi_rapidity = fine_->q;
          info_.mu_unshifted = fine_->mu;
          info_.mu_symmetric = fine_->mu - interaction / Real{2};
          info_.energy_per_site_unshifted = fine_->e0;
          return;
        }
        if (n > options.max_nodes / 2) return;
      }
    }
    DopedBackground<Real> const& background() const { return info_; }

    /// Unwrapped momenta chosen to agree with the half-filled hole convention.
    std::pair<Real, Real> momentum_range(DopedBranch branch) const
    {
      Real const p = detail::pi<Real>(), pn = p * info_.density;
      switch (branch)
      {
        case DopedBranch::spinon:
          return {Real{0}, pn};
        case DopedBranch::holon:
          return {-pn / Real{2}, Real{3} * pn / Real{2}};
        case DopedBranch::charge_particle:
          return {pn / Real{2}, Real{2} * p - Real{3} * pn / Real{2}};
      }
      throw std::invalid_argument("invalid doped excitation branch");
    }

    /// Evaluate both converged meshes at the SAME physical momentum. A
    /// mesh-difference check is also made for this particular requested point.
    DopedPoint<Real> at_momentum(DopedBranch branch, Real momentum, Convention convention = Convention::symmetric,
                                 EnergyReference reference = EnergyReference::hamiltonian) const
    {
      auto const [low, high] = momentum_range(branch);
      if (!uni20::isfinite(momentum) || momentum < low || momentum > high)
        throw std::invalid_argument("momentum is outside this doped branch's unwrapped interval");
      if (convention != Convention::symmetric && convention != Convention::unshifted)
        throw std::invalid_argument("invalid Hubbard energy convention");
      if (reference != EnergyReference::hamiltonian && reference != EnergyReference::fermi)
        throw std::invalid_argument("invalid excitation energy reference");
      DopedPoint<Real> out;
      out.branch = branch;
      out.momentum = momentum;
      out.convention = convention;
      out.reference = reference;
      out.delta_particles = branch == DopedBranch::holon ? -1 : branch == DopedBranch::charge_particle ? 1 : 0;
      if (branch == DopedBranch::spinon) out.spin = uni20::from_twice(1);
      if (!info_.converged)
      {
        out.status = info_.status;
        return out;
      }
      Real ef{}, ec{};
      if (momentum == low || momentum == high)
      {
        out.parameter = branch == DopedBranch::spinon ? uni20::numeric_limits<Real>::infinity() : fine_->q;
        if (momentum == high) out.parameter = -out.parameter;
      }
      else
      {
        auto const fine = invert(*fine_, branch, momentum, out.iterations);
        if (fine.status != DopedStatus::converged)
        {
          out.status = fine.status;
          return out;
        }
        auto const coarse = invert(*coarse_, branch, momentum, out.iterations, fine.parameter);
        if (coarse.status != DopedStatus::converged)
        {
          out.status = coarse.status;
          return out;
        }
        ef = fine.energy;
        ec = coarse.energy;
        // An interior zero can be cancellation/underflow, not a gapless endpoint.
        if (!(ef > Real{0}) || !(ec > Real{0}))
        {
          out.status = DopedStatus::precision_limit;
          return out;
        }
        out.parameter = fine.parameter;
        out.momentum_error = std::max(fine.error, coarse.error);
        auto const mf = branch == DopedBranch::spinon ? fine_->spinon(out.parameter) : fine_->charge(out.parameter);
        auto const mc = branch == DopedBranch::spinon ? coarse_->spinon(out.parameter) : coarse_->charge(out.parameter);
        out.momentum_error += Real{8} * std::abs(mf[1] - mc[1]);
      }
      out.momentum_error +=
          Real{2} * detail::pi<Real>() * std::max(info_.density_error, std::abs(coarse_->density - info_.density));
      Real const dn = Real(out.delta_particles), mu = fine_->mu;
      Real const offset =
          reference == EnergyReference::fermi
              ? Real{0}
              : (mu - (convention == Convention::symmetric ? info_.interaction / Real{2} : Real{0})) * dn;
      Real const energy = ef + offset;
      out.energy_error =
          Real{8} * std::abs(ef - ec) +
          (reference == EnergyReference::fermi ? Real{0} : Real{8} * std::abs(fine_->mu - coarse_->mu) * std::abs(dn)) +
          Real{64} * uni20::numeric_limits<Real>::epsilon() * std::max(Real{1}, std::abs(energy));
      if (!uni20::isfinite(energy) || !uni20::isfinite(out.momentum_error) || out.momentum_error > options_.tolerance ||
          ef < Real{0} || out.energy_error > options_.tolerance * std::max(Real{1}, std::abs(energy)))
      {
        out.status = DopedStatus::mesh_limit;
        return out;
      }
      out.energy = energy;
      out.fermi_energy = ef;
      out.symmetric_energy = ef + (mu - info_.interaction / Real{2}) * dn;
      out.converged = true;
      out.status = DopedStatus::converged;
      return out;
    }

  private:
    struct Inversion
    {
        Real energy{}, parameter{}, error{};
        DopedStatus status = DopedStatus::precision_limit;
    };
    Inversion invert(detail::DopedMesh<Real> const& mesh, DopedBranch branch, Real momentum, std::size_t& iterations,
                     std::optional<Real> seed = {}) const
    {
      Inversion result;
      Real const p = detail::pi<Real>(), pn = p * info_.density;
      bool const spin = branch == DopedBranch::spinon;
      Real target;
      bool reflect;
      if (spin)
      {
        reflect = momentum > pn / Real{2};
        target = reflect ? pn - momentum : momentum;
      }
      else
      {
        target = branch == DopedBranch::holon ? pn / Real{2} - momentum : momentum + pn / Real{2};
        if (target > p) target -= Real{2} * p;
        reflect = target < Real{0};
        target = std::abs(target);
      }
      Real lo = spin || branch == DopedBranch::holon ? Real{0} : mesh.q;
      Real hi = spin ? Real{1} + Real{2} * mesh.u / p * std::log(Real{4} * info_.density / target)
                : branch == DopedBranch::holon ? mesh.q
                                               : p;
      Real flo = spin ? std::log((pn / Real{2}) / target) : target - (lo == Real{0} ? Real{0} : pn);
      Real fhi = spin ? -Real{1} : target - (hi == p ? p : pn);
      Real x = spin                           ? (target == pn / Real{2} ? Real{0} : (lo + hi) / Real{2})
               : branch == DopedBranch::holon ? mesh.q * target / pn
                                              : mesh.q + (p - mesh.q) * (target - pn) / (p - pn);
      if (seed && std::abs(*seed) >= lo && std::abs(*seed) <= hi) x = std::abs(*seed);
      int side = 0;
      for (;;)
      {
        auto const values = spin ? mesh.spinon(x) : mesh.charge(x);
        result.energy = spin || branch == DopedBranch::charge_particle ? values[0] : -values[0];
        result.error = std::abs(values[1] - target);
        result.parameter = reflect ? -x : x;
        if (!uni20::isfinite(result.energy) || !uni20::isfinite(result.error)) return result;
        if (result.error <= options_.tolerance * (spin ? target : Real{1}) / Real{16})
        {
          result.status = DopedStatus::converged;
          return result;
        }
        if (iterations == options_.max_iterations)
        {
          result.status = DopedStatus::momentum_limit;
          return result;
        }
        ++iterations;
        Real const f = spin ? std::log(values[1] / target) : target - values[1];
        if (f > Real{0})
        {
          lo = x;
          flo = f;
          if (side == 1) fhi /= Real{2};
          side = 1;
        }
        else
        {
          hi = x;
          fhi = f;
          if (side == -1) flo /= Real{2};
          side = -1;
        }
        Real fraction = flo / (flo - fhi);
        if (!(fraction > Real{0} && fraction < Real{1})) fraction = Real{.5};
        Real const next = lo + fraction * (hi - lo);
        if (next == x || next == lo || next == hi) return result;
        x = next;
      }
    }
    DopedOptions<Real> options_;
    DopedBackground<Real> info_;
    std::optional<detail::DopedMesh<Real>> fine_, coarse_;
};
} // namespace bethe::hubbard::thermo
