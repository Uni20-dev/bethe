// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/quadrature.hpp>
#include <optional>
#include <stdexcept>
#include <uni20/common/half_int.hpp>

namespace bethe::hubbard::thermo
{
enum class Branch
{
  spinon,
  holon,
  antiholon
};
enum class Convention
{
  symmetric,
  unshifted
};
enum class Status
{
  converged,
  quadrature_limit,
  momentum_limit,
  precision_limit
};
template <uni20::Real Real> struct Options
{
    // Convergence estimates, not rigorous error certificates. Budgets apply
    // to a whole point, including all evaluations during momentum inversion.
    Real relative_tolerance = Real{256} * uni20::numeric_limits<Real>::epsilon();
    std::size_t max_evaluations = 1000000, max_levels = 12, max_iterations = 160;
};
template <uni20::Real Real> struct Point
{
    Branch branch = Branch::spinon;
    Convention convention = Convention::symmetric;
    Real interaction{}, momentum{}, parameter{}, energy_offset{}, energy_error{}, momentum_error{};
    std::optional<Real> energy, symmetric_energy;
    int delta_particles = 0;
    uni20::half_int spin{0};
    std::size_t evaluations = 0, iterations = 0;
    bool converged = false;
    Status status = Status::precision_limit;
};
namespace detail
{
template <uni20::Real Real> Real pi() { return Real{4} * std::atan(Real{1}); }
template <uni20::Real Real> Real sech(Real x)
{
  Real const z = std::exp(-std::abs(x));
  return Real{2} * z / (Real{1} + z * z);
}
template <uni20::Real Real> Real atan_exp(Real x)
{
  return x <= Real{0} ? std::atan(std::exp(x)) : pi<Real>() / Real{2} - std::atan(std::exp(-x));
}
template <uni20::Real Real> Real atanh_exp_minus(Real x)
{
  if (x > Real{.5}) return std::atanh(std::exp(-x));
  return (std::log1p(std::exp(-x)) - std::log(-std::expm1(-x))) / Real{2};
}
template <uni20::Real Real> Real wrap(Real p)
{
  Real const half = pi<Real>(), period = Real{2} * half;
  if (p < -half) p += period;
  if (p >= half) p -= period;
  return p;
}
template <uni20::Real Real> void validate(Real interaction, Options<Real> const& o, Convention convention)
{
  if (!uni20::isfinite(interaction) || interaction <= Real{0})
    throw std::invalid_argument("half-filled Hubbard dispersions require finite U>0 (t=1, zero field)");
  if (!uni20::isfinite(o.relative_tolerance) || o.relative_tolerance <= Real{0} || o.relative_tolerance >= Real{1})
    throw std::invalid_argument("relative tolerance must be finite and in (0,1)");
  if (o.max_levels > 24) throw std::invalid_argument("quadrature levels must not exceed 24");
  if (convention != Convention::symmetric && convention != Convention::unshifted)
    throw std::invalid_argument("invalid Hubbard energy convention");
}
template <uni20::Real Real> struct IntegralPoint
{
    Real energy{}, momentum{}, energy_error{}, momentum_error{};
    bool converged = false;
};
template <uni20::Real Real> class Evaluator {
  public:
    Evaluator(Real interaction, Options<Real> const& options)
        : u(interaction / Real{4}), o(options), tolerance(options.relative_tolerance / Real{8}),
          s(Real{2} * u / pi<Real>())
    {}

    IntegralPoint<Real> spinon(Real lambda)
    {
      Real const l = std::abs(lambda), p = pi<Real>();
      IntegralPoint<Real> out;
      if (!(u > Real{0}) || !(s > Real{0})) return out;
      if (u >= Real{1} / Real{4})
      {
        auto f = [&](Real theta) {
          Real const sine = std::sin(theta), x = (l - std::cos(theta)) / s;
          return std::array<Real, 2>{sine * sine * sech(x) / u, Real{2} / p * atan_exp(-x)};
        };
        auto const integral = integrate(f, Real{0}, p);
        out = {integral.value[0], integral.value[1], integral.error[0], integral.error[1], integral.converged};
      }
      else
      {
        // Fourier-convolution form: nonoscillatory even as U tends to zero.
        // Split at x=+-1, where sqrt(1-x*x) and acos(x) change support.
        Real const cut = -std::log(tolerance) - std::log(u) / Real{2} + Real{12};
        std::vector<Real> breaks;
        bool const outside = l >= Real{1};
        Real const offset = outside ? (l - Real{1}) / s : Real{0};
        Real const prefactor = outside ? std::exp(-offset) : Real{1};
        if (prefactor == Real{0}) return out;
        breaks = {outside ? Real{0} : -cut, cut};
        for (Real edge : outside ? std::array<Real, 2>{Real{0}, Real{2} / s}
                                 : std::array<Real, 2>{(-Real{1} - l) / s, (Real{1} - l) / s})
          if (edge > breaks.front() && edge < cut) breaks.push_back(edge);
        std::sort(breaks.begin(), breaks.end());
        std::array<bethe::detail::CompensatedSum<Real>, 2> sums, errors;
        bool complete = true;
        auto f = [&](Real y) {
          Real const x = outside ? Real{1} - s * y : l + s * y;
          Real const weight =
              outside ? Real{2} * std::exp(-y) / (Real{1} + std::exp(-Real{2} * (offset + y))) : sech(y);
          Real const square = outside ? (s * y) * (Real{2} - s * y) : (Real{1} - x) * (Real{1} + x);
          Real const angle = outside         ? (s * y >= Real{2} ? p : Real{2} * std::asin(std::sqrt(s * y / Real{2})))
                             : x >= Real{1}  ? Real{0}
                             : x <= -Real{1} ? p
                                             : std::acos(x);
          return std::array<Real, 2>{Real{2} / p * weight * std::sqrt(std::max(Real{0}, square)), weight * angle / p};
        };
        for (std::size_t j = 1; j < breaks.size(); ++j)
        {
          auto const integral = integrate(f, breaks[j - 1], breaks[j]);
          complete = complete && integral.converged;
          for (std::size_t a = 0; a < 2; ++a)
          {
            sums[a].add(integral.value[a]);
            errors[a].add(integral.error[a]);
          }
          if (!complete) break;
        }
        // Both convolution integrands are bounded by constants times sech.
        Real const tail = Real{8} * std::exp(-cut);
        out = {prefactor * sums[0].value(), prefactor * sums[1].value(), prefactor * (errors[0].value() + tail),
               prefactor * (errors[1].value() + tail), complete};
      }
      if (l == Real{0})
      {
        out.momentum = p / Real{2};
        out.momentum_error = Real{0};
      }
      else if (lambda < Real{0})
        out.momentum = p - out.momentum;
      out.converged = out.converged && out.energy > Real{0} && uni20::isfinite(out.energy) &&
                      out.energy_error <= o.relative_tolerance * out.energy;
      return out;
    }

    IntegralPoint<Real> charge(Real k)
    {
      Real const p = pi<Real>();
      IntegralPoint<Real> out;
      if (!(s > Real{0})) return out;
      Real const a = Real{1} / s, cut = -std::log(tolerance) + Real{12};
      Real const ratio = cut * s;
      if (!uni20::isfinite(a) || !uni20::isfinite(ratio)) return out;
      Real const end = Real{2} * std::asinh(std::sqrt(ratio / Real{2}));
      Real dm = Real{2} * std::pow(std::sin((k - p / Real{2}) / Real{2}), 2);
      Real dp = Real{2} * std::pow(std::sin((k + p / Real{2}) / Real{2}), 2);
      if (k == Real{0} || std::abs(k) == p) dm = dp = Real{1};
      bool const outer = std::abs(k) < p / Real{2};
      auto f = [&](Real t) {
        Real const sh = std::sinh(t / Real{2}), h = Real{2} * sh * sh;
        Real const am = atanh_exp_minus(a * (dm + h)), ap = atanh_exp_minus(a * (dp + h));
        return std::array<Real, 2>{Real{4} / p * std::cosh(t) * (am + ap), Real{2} / p * (am - ap)};
      };
      // Resummation of Melzer's positive K_0/K_1 series; no Bessel library
      // or cancellation of order-one terms to obtain the Mott gap.
      auto const integral = integrate(f, Real{0}, end, {Real{0}, Real{1}});
      Real const sine = std::sqrt(ratio) * std::sqrt(Real{2} + ratio);
      Real const tail = std::exp(-a * std::min(dm, dp) - cut);
      out.energy = integral.value[0] + (outer ? Real{4} * std::cos(k) : Real{0});
      Real const base = outer ? p / Real{2} - Real{2} * k : k >= Real{0} ? -p / Real{2} : Real{3} * p / Real{2};
      out.momentum = base + integral.value[1]; // unwrapped: -pi/2 <= p_h <= 3*pi/2
      out.energy_error = integral.error[0] + Real{16} / p * s * (Real{1} + ratio) / sine * tail;
      out.momentum_error = integral.error[1] + Real{8} / p * s / sine * tail;
      out.converged = integral.converged && out.energy > Real{0} && uni20::isfinite(out.energy) &&
                      out.energy_error <= o.relative_tolerance * out.energy;
      return out;
    }
    std::size_t evaluations = 0;
    Status failure_status() const { return quadrature_failed ? Status::quadrature_limit : Status::precision_limit; }
    Real const u;
    Options<Real> const& o;
    Real const tolerance, s;

  private:
    bool quadrature_failed = false;
    template <typename Function>
    auto integrate(Function const& f, Real lo, Real hi, std::array<Real, 2> const& scales = {})
    {
      auto result =
          bethe::detail::tanh_sinh<Real, 2>(f, lo, hi, tolerance, evaluations, o.max_evaluations, o.max_levels, scales);
      quadrature_failed = quadrature_failed || !result.converged;
      return result;
    }
};

template <uni20::Real Real> Point<Real> initialize(Branch branch, Real interaction, Convention convention)
{
  Point<Real> out;
  out.branch = branch;
  out.interaction = interaction;
  out.convention = convention;
  if (branch == Branch::spinon)
    out.spin = uni20::from_twice(1);
  else if (branch == Branch::holon)
    out.delta_particles = -1;
  else if (branch == Branch::antiholon)
    out.delta_particles = 1;
  else
    throw std::invalid_argument("invalid Hubbard excitation branch");
  if (convention == Convention::unshifted) out.energy_offset = interaction / Real{2} * Real(out.delta_particles);
  return out;
}
template <uni20::Real Real> void finish(Point<Real>& out, IntegralPoint<Real> const& integral)
{
  out.energy_error = integral.energy_error;
  out.momentum_error = integral.momentum_error;
  if (!integral.converged) return;
  out.symmetric_energy = integral.energy;
  out.energy = integral.energy + out.energy_offset;
  out.energy_error +=
      Real{4} * uni20::numeric_limits<Real>::epsilon() * (std::abs(*out.energy) + std::abs(out.energy_offset));
  out.converged = true;
  out.status = Status::converged;
}
} // namespace detail

/// Parametric spinon line. Lambda is the conventional Lieb-Wu spin rapidity.
template <uni20::Real Real>
[[nodiscard]] Point<Real> spinon_at_rapidity(Real interaction, Real lambda,
                                             Convention convention = Convention::symmetric,
                                             Options<Real> const& options = {})
{
  detail::validate(interaction, options, convention);
  if (!uni20::isfinite(lambda)) throw std::invalid_argument("spin rapidity must be finite");
  auto out = detail::initialize<Real>(Branch::spinon, interaction, convention);
  out.parameter = lambda;
  detail::Evaluator<Real> eval(interaction, options);
  auto const integral = eval.spinon(lambda);
  out.status = eval.failure_status();
  out.momentum = integral.momentum;
  out.evaluations = eval.evaluations;
  detail::finish(out, integral);
  return out;
}
/// Parametric charge line, bare k in [-pi,pi]. The antiholon has the same
/// symmetric energy at k, and momentum p_h(k)-pi, not p_h(k).
template <uni20::Real Real>
[[nodiscard]] Point<Real> charge_at_bare_momentum(Branch branch, Real interaction, Real k,
                                                  Convention convention = Convention::symmetric,
                                                  Options<Real> const& options = {})
{
  detail::validate(interaction, options, convention);
  if (branch == Branch::spinon) throw std::invalid_argument("charge branch must be holon or antiholon");
  if (!uni20::isfinite(k) || std::abs(k) > detail::pi<Real>())
    throw std::invalid_argument("bare momentum must be in [-pi,pi]");
  auto out = detail::initialize<Real>(branch, interaction, convention);
  out.parameter = k;
  detail::Evaluator<Real> eval(interaction, options);
  auto const integral = eval.charge(k);
  out.status = eval.failure_status();
  out.momentum = detail::wrap(integral.momentum - (branch == Branch::antiholon ? detail::pi<Real>() : Real{0}));
  out.evaluations = eval.evaluations;
  detail::finish(out, integral);
  return out;
}

/// Evaluate at physical dressed momentum: spinon [0,pi], charge [-pi,pi].
/// Budgets cover ALL quadratures required by this point's momentum inversion.
template <uni20::Real Real>
[[nodiscard]] Point<Real> dispersion(Branch branch, Real interaction, Real momentum,
                                     Convention convention = Convention::symmetric, Options<Real> const& options = {})
{
  detail::validate(interaction, options, convention);
  auto out = detail::initialize<Real>(branch, interaction, convention);
  Real const p = detail::pi<Real>();
  out.momentum = momentum;
  if (!uni20::isfinite(momentum) || momentum > p || momentum < (branch == Branch::spinon ? Real{0} : -p))
    throw std::invalid_argument("dressed momentum must be in [0,pi] for spinons, [-pi,pi] for charge branches");
  if (branch == Branch::spinon && (momentum == Real{0} || momentum == p))
  {
    out.parameter =
        momentum == Real{0} ? uni20::numeric_limits<Real>::infinity() : -uni20::numeric_limits<Real>::infinity();
    detail::finish(out, {Real{0}, momentum, Real{0}, Real{0}, true});
    return out;
  }
  detail::Evaluator<Real> eval(interaction, options);
  bool const spin = branch == Branch::spinon, reflect = spin && momentum > p / Real{2};
  Real target = spin ? (reflect ? p - momentum : momentum) : momentum;
  if (!spin)
  {
    if (branch == Branch::antiholon) target += p;
    target = detail::wrap(target);
    if (target < -p / Real{2}) target += Real{2} * p;
  }
  // A rapidity bracket with p_s(upper)<target follows from cos(theta)<=1.
  Real lo = spin ? Real{0} : -p, hi = spin ? Real{1} + eval.s * std::log(Real{4} / target) : p;
  Real flo = spin ? std::log((p / Real{2}) / target) : Real{3} * p / Real{2} - target;
  Real fhi = spin ? -Real{1} : -p / Real{2} - target;
  Real x = spin ? (target == p / Real{2} ? Real{0} : (lo + hi) / Real{2}) : p / Real{2} - target;
  if (!spin && target == -p / Real{2}) x = p;
  int previous_side = 0;
  for (;;)
  {
    auto integral = spin ? eval.spinon(x) : eval.charge(x);
    out.evaluations = eval.evaluations;
    if (!integral.converged)
    {
      out.status = eval.failure_status();
      return out;
    }
    Real const residual = integral.momentum - target;
    Real const momentum_tolerance = options.relative_tolerance * (spin ? target : Real{1});
    out.momentum_error = std::abs(residual) + integral.momentum_error;
    if (out.momentum_error <= momentum_tolerance)
    {
      out.parameter = reflect ? -x : x;
      // Report horizontal and quadrature errors separately. The latter is not
      // an energy-error bound propagated through the momentum inversion.
      auto const error = out.momentum_error;
      detail::finish(out, integral);
      out.momentum_error = error;
      return out;
    }
    if (out.iterations == options.max_iterations)
    {
      out.status = Status::momentum_limit;
      return out;
    }
    ++out.iterations;
    Real const f = spin ? std::log(integral.momentum / target) : residual;
    if (f > Real{0})
    {
      lo = x;
      flo = f;
      if (previous_side == 1) fhi /= Real{2};
      previous_side = 1;
    }
    else
    {
      hi = x;
      fhi = f;
      if (previous_side == -1) flo /= Real{2};
      previous_side = -1;
    }
    Real fraction = flo / (flo - fhi);
    if (!(fraction > Real{0} && fraction < Real{1})) fraction = Real{.5};
    Real const next = lo + fraction * (hi - lo);
    if (next == x || next == lo || next == hi)
    {
      out.status = Status::precision_limit;
      return out;
    }
    x = next;
  }
}
} // namespace bethe::hubbard::thermo
