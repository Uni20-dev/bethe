// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/gauss_legendre.hpp>
#include <bethe/detail/hubbard_kernel.hpp>
#include <bethe/detail/newton.hpp>
#include <bethe/hubbard_thermo.hpp>
#include <memory>

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
    uni20::Precision precision() const { return uni20::native_precision<Real>(); }
};
#if UNI20_HAS_MPREAL
template <> struct DopedOptions<uni20::mpreal>
{
    explicit DopedOptions(uni20::Precision p) : tolerance(default_tolerance(p)), precision_(p) {}
    uni20::mpreal tolerance;
    std::size_t initial_nodes = 16, max_nodes = 256, max_background_iterations = 64, max_iterations = 160;
    uni20::Precision precision() const { return precision_; }

  private:
    static uni20::mpreal default_tolerance(uni20::Precision p)
    {
      if (p.is_exact()) throw std::invalid_argument("doped Hubbard dispersions require finite working precision");
      return 4096 * uni20::numeric_limits<uni20::mpreal>::epsilon(p);
    }
    uni20::Precision precision_;
};
#endif
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
    std::shared_ptr<bethe::detail::HubbardKernel<Real> const> kernel;

    // Returns (epsilon_c(k), p_c(k)); p_c grows from -pi to pi.
    std::array<Real, 2> charge(Real parameter) const
    {
      Real const sk = uni20::math::sin(parameter), p = pi(u), zero = uni20::scalar_like(u, 0);
      bethe::detail::CompensatedSum<Real> energy(zero), momentum(zero);
      for (std::size_t j = 0; j < k.size(); ++j)
      {
        Real const r = kernel->value(sk - sine[j]) + kernel->value(sk + sine[j]);
        energy.add(w[j] * cosine[j] * r * epsilon[j]);
        Real const phase = kernel->primitive(sk - sine[j]) + kernel->primitive(sk + sine[j]);
        momentum.add(w[j] * rho[j] * phase);
      }
      Real e = -Real{2} * uni20::math::cos(parameter) - mu + energy.value();
      Real pc = parameter + Real{2} * p * momentum.value();
      if (uni20::math::abs(parameter) == q)
      {
        e = zero;
        pc = uni20::math::copysign(p * density, parameter);
      }
      if (uni20::math::abs(parameter) == p) pc = parameter;
      return {e, pc};
    }
    // Hole momentum convention: p_s in [0,pi*n], matching the half-filled path.
    std::array<Real, 2> spinon(Real lambda) const
    {
      Real const l = uni20::math::abs(lambda), scale = Real{2} * u / pi(u), zero = uni20::scalar_like(u, 0);
      bethe::detail::CompensatedSum<Real> energy(zero), momentum(zero);
      for (std::size_t j = 0; j < k.size(); ++j)
      {
        Real const minus = (l - sine[j]) / scale, plus = (l + sine[j]) / scale;
        energy.add(-w[j] * cosine[j] * epsilon[j] * (sech(minus) + sech(plus)) / (Real{4} * u));
        momentum.add(Real{2} * w[j] * rho[j] * (atan_exp(-minus) + atan_exp(-plus)));
      }
      Real ps = momentum.value();
      if (l == Real{0})
        ps = pi(u) * density / Real{2};
      else if (lambda < Real{0})
        ps = pi(u) * density - ps;
      return {energy.value(), ps};
    }
};

template <uni20::Real Real>
std::optional<DopedMesh<Real>> doped_mesh(Real q, bethe::detail::GaussLegendreRule<Real> const& rule,
                                          std::shared_ptr<bethe::detail::HubbardKernel<Real> const> kernel)
{
  std::size_t const n = rule.x.size();
  Real const zero = uni20::scalar_like(q, 0), one = uni20::scalar_like(q, 1);
  auto const precision = uni20::precision_of(q);
  DopedMesh<Real> out;
  out.u = kernel->u();
  out.kernel = std::move(kernel);
  out.q = q;
  out.k.resize(n);
  out.w.resize(n);
  out.sine.resize(n);
  out.cosine.resize(n);
  out.rho.assign(n, Real{1} / (Real{2} * pi(q)));
  std::vector<Real> arho(n * n), ae(n * n), h(n), xi(n, one);
  for (std::size_t j = 0; j < n; ++j)
  {
    out.k[j] = q * (rule.x[j] + Real{1}) / Real{2};
    out.w[j] = q * rule.w[j] / Real{2};
    out.sine[j] = uni20::math::sin(out.k[j]);
    out.cosine[j] = uni20::math::cos(out.k[j]);
    h[j] = -Real{2} * out.cosine[j];
  }
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
    {
      Real const kernel = out.kernel->value(out.sine[i] - out.sine[j]) + out.kernel->value(out.sine[i] + out.sine[j]);
      arho[i * n + j] = Real(i == j) - out.cosine[i] * out.w[j] * kernel;
      ae[i * n + j] = Real(i == j) - out.cosine[j] * out.w[j] * kernel;
    }
  // Normal Uni20 dispatch selects the configured provider at working precision;
  // singular systems remain recoverable without a process-global error policy.
  if (!bethe::detail::newton_step(arho, out.rho, precision) || !bethe::detail::newton_step(ae, h, precision) ||
      !bethe::detail::newton_step(ae, xi, precision))
    return std::nullopt;
  bethe::detail::CompensatedSum<Real> number(zero), e0(zero), hq(zero), xiq(zero), rhoq(zero);
  hq.add(-Real{2} * uni20::math::cos(q));
  xiq.add(one);
  rhoq.add(Real{1} / (Real{2} * pi(q)));
  for (std::size_t j = 0; j < n; ++j)
  {
    if (!(out.rho[j] > Real{0})) return std::nullopt;
    number.add(Real{2} * out.w[j] * out.rho[j]);
    e0.add(-Real{4} * out.w[j] * out.cosine[j] * out.rho[j]);
    Real const r =
        out.kernel->value(uni20::math::sin(q) - out.sine[j]) + out.kernel->value(uni20::math::sin(q) + out.sine[j]);
    hq.add(out.w[j] * out.cosine[j] * r * h[j]);
    xiq.add(out.w[j] * out.cosine[j] * r * xi[j]);
    rhoq.add(uni20::math::cos(q) * out.w[j] * r * out.rho[j]);
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
/// at explicit working precision; mesh differences are estimates, not interval bounds.
template <uni20::Real Real> class DopedSolver {
  public:
    DopedSolver(Real interaction, Real density, DopedOptions<Real> options = {}) : options_(options)
    {
      auto const precision = options.precision();
      Real const zero = uni20::scalar_at_precision<Real>(0, precision);
      bethe::detail::require_same_precision(interaction, zero);
      bethe::detail::require_same_precision(density, zero);
      bethe::detail::require_same_precision(options.tolerance, zero);
      info_.mesh_error = info_.density_error = zero;
      info_.interaction = interaction;
      info_.density = density;
      if (!uni20::isfinite(interaction) || interaction <= Real{0} || !uni20::isfinite(density) || density <= Real{0} ||
          density >= Real{1})
        throw std::invalid_argument("doped Hubbard dispersions require finite U>0 and 0<density<1, at zero field");
      if (!uni20::isfinite(options.tolerance) || options.tolerance <= Real{0} || options.tolerance >= Real{1} ||
          options.initial_nodes < 4 || options.max_nodes > 512 || options.initial_nodes > options.max_nodes)
        throw std::invalid_argument(
            "invalid doped tolerance/mesh: require 0<tol<1 and 4<=initial_nodes<=max_nodes<=512");
      Real const u = interaction / Real{4}, p = detail::pi(info_.interaction);
      if (!(u > Real{0}) || !uni20::isfinite(Real{1} / u))
      {
        info_.status = DopedStatus::precision_limit;
        return;
      }
      auto const kernel = std::make_shared<bethe::detail::HubbardKernel<Real> const>(u, precision);
      Real q = p * density * (Real{1} + interaction / (interaction + Real{4})) / Real{2};
      for (std::size_t n = options.initial_nodes;; n *= 2)
      {
        auto const rule = bethe::detail::gauss_legendre<Real>(n, precision);
        Real lo = zero, hi = p;
        std::optional<detail::DopedMesh<Real>> candidate;
        for (;;)
        {
          if (info_.iterations == options.max_background_iterations)
          {
            info_.status = DopedStatus::density_limit;
            return;
          }
          ++info_.iterations;
          candidate = detail::doped_mesh(q, rule, kernel);
          if (!candidate)
          {
            info_.status = DopedStatus::linear_failure;
            return;
          }
          Real const residual = candidate->density - density;
          if (uni20::math::abs(residual) <= options.tolerance * density / Real{64}) break;
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
        info_.density_error = uni20::math::abs(candidate->density - density);
        if (fine_)
        {
          Real error = std::max({uni20::math::abs(candidate->q - fine_->q), uni20::math::abs(candidate->mu - fine_->mu),
                                 uni20::math::abs(candidate->e0 - fine_->e0)});
          for (int j = 0; j <= 8; ++j)
          {
            Real const k = p * Real(j) / Real{8};
            auto const a = candidate->charge(k), b = fine_->charge(k);
            for (int v = 0; v < 2; ++v)
              error = std::max(error, uni20::math::abs(a[v] - b[v]));
            Real const lambda = (Real{1} + u) * Real(j) / Real{4};
            auto const sa = candidate->spinon(lambda), sb = fine_->spinon(lambda);
            for (int v = 0; v < 2; ++v)
              error = std::max(error, uni20::math::abs(sa[v] - sb[v]));
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
      Real const p = detail::pi(info_.interaction), pn = p * info_.density;
      switch (branch)
      {
        case DopedBranch::spinon:
          return {uni20::scalar_like(p, 0), pn};
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
      bethe::detail::require_same_precision(momentum, info_.interaction);
      Real const zero = uni20::scalar_like(info_.interaction, 0);
      auto const [low, high] = momentum_range(branch);
      if (!uni20::isfinite(momentum) || momentum < low || momentum > high)
        throw std::invalid_argument("momentum is outside this doped branch's unwrapped interval");
      if (convention != Convention::symmetric && convention != Convention::unshifted)
        throw std::invalid_argument("invalid Hubbard energy convention");
      if (reference != EnergyReference::hamiltonian && reference != EnergyReference::fermi)
        throw std::invalid_argument("invalid excitation energy reference");
      DopedPoint<Real> out;
      out.parameter = out.energy_error = out.momentum_error = zero;
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
      Real ef = zero, ec = zero;
      if (momentum == low || momentum == high)
      {
        out.parameter = branch == DopedBranch::spinon
                            ? uni20::scalar_like(zero, std::numeric_limits<double>::infinity())
                            : fine_->q;
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
        out.momentum_error += Real{8} * uni20::math::abs(mf[1] - mc[1]);
      }
      out.momentum_error += Real{2} * detail::pi(info_.interaction) *
                            std::max(info_.density_error, uni20::math::abs(coarse_->density - info_.density));
      Real const dn = uni20::scalar_like(zero, out.delta_particles), mu = fine_->mu;
      Real const offset = reference == EnergyReference::fermi
                              ? zero
                              : (mu - (convention == Convention::symmetric ? info_.interaction / Real{2} : zero)) * dn;
      Real const energy = ef + offset;
      out.energy_error = Real{8} * uni20::math::abs(ef - ec) +
                         (reference == EnergyReference::fermi
                              ? zero
                              : Real{8} * uni20::math::abs(fine_->mu - coarse_->mu) * uni20::math::abs(dn)) +
                         Real{64} * uni20::numeric_limits<Real>::epsilon(options_.precision()) *
                             std::max(Real{1}, uni20::math::abs(energy));
      if (!uni20::isfinite(energy) || !uni20::isfinite(out.momentum_error) || out.momentum_error > options_.tolerance ||
          ef < Real{0} || out.energy_error > options_.tolerance * std::max(Real{1}, uni20::math::abs(energy)))
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
      Real const zero = uni20::scalar_like(info_.interaction, 0), one = uni20::scalar_like(info_.interaction, 1);
      Inversion result{zero, zero, zero};
      Real const p = detail::pi(info_.interaction), pn = p * info_.density;
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
        target = uni20::math::abs(target);
      }
      Real lo = spin || branch == DopedBranch::holon ? zero : mesh.q;
      Real hi = spin ? Real{1} + Real{2} * mesh.u / p * uni20::math::log(Real{4} * info_.density / target)
                : branch == DopedBranch::holon ? mesh.q
                                               : p;
      Real flo = spin ? uni20::math::log((pn / Real{2}) / target) : target - (lo == Real{0} ? zero : pn);
      Real fhi = spin ? -one : target - (hi == p ? p : pn);
      Real x = spin                           ? (target == pn / Real{2} ? zero : (lo + hi) / Real{2})
               : branch == DopedBranch::holon ? mesh.q * target / pn
                                              : mesh.q + (p - mesh.q) * (target - pn) / (p - pn);
      if (seed && uni20::math::abs(*seed) >= lo && uni20::math::abs(*seed) <= hi) x = uni20::math::abs(*seed);
      int side = 0;
      for (;;)
      {
        auto const values = spin ? mesh.spinon(x) : mesh.charge(x);
        result.energy = spin || branch == DopedBranch::charge_particle ? values[0] : -values[0];
        result.error = uni20::math::abs(values[1] - target);
        result.parameter = reflect ? -x : x;
        if (!uni20::isfinite(result.energy) || !uni20::isfinite(result.error)) return result;
        if (result.error <= options_.tolerance * (spin ? target : one) / Real{16})
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
        Real const f = spin ? uni20::math::log(values[1] / target) : target - values[1];
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
        if (!(fraction > Real{0} && fraction < Real{1})) fraction = one / Real{2};
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
