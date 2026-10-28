// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bethe/detail/quadrature.hpp>
#include <optional>
#include <stdexcept>
#include <uni20/common/half_int.hpp>
#include <uni20/core/finite_quadrature.hpp>

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
    uni20::Precision precision() const { return uni20::native_precision<Real>(); }
};
#if UNI20_HAS_MPREAL
template <> struct Options<uni20::mpreal>
{
    explicit Options(uni20::Precision p) : relative_tolerance(default_tolerance(p)), precision_(p) {}
    uni20::mpreal relative_tolerance;
    std::size_t max_evaluations = 1000000, max_levels = 12, max_iterations = 160;
    uni20::Precision precision() const { return precision_; }

  private:
    static uni20::mpreal default_tolerance(uni20::Precision p)
    {
      if (p.is_exact()) throw std::invalid_argument("Hubbard dispersions require finite working precision");
      return 256 * uni20::numeric_limits<uni20::mpreal>::epsilon(p);
    }
    uni20::Precision precision_;
};
#endif
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
template <uni20::Real Real> Real pi(Real const& value)
{
  return Real{4} * uni20::math::atan(uni20::scalar_like(value, 1));
}
template <uni20::Real Real>
  requires(!uni20::has_runtime_precision_v<Real>)
Real pi()
{
  return pi(Real{0});
}
template <uni20::Real Real> Real sech(Real x)
{
  Real const z = uni20::math::exp(-uni20::math::abs(x));
  return Real{2} * z / (Real{1} + z * z);
}
template <uni20::Real Real> Real atan_exp(Real x)
{
  return x <= Real{0} ? uni20::math::atan(uni20::math::exp(x))
                      : pi(x) / Real{2} - uni20::math::atan(uni20::math::exp(-x));
}
template <uni20::Real Real> Real atanh_exp_minus(Real x)
{
  if (x > Real{.5}) return uni20::math::atanh(uni20::math::exp(-x));
  return (uni20::math::log1p(uni20::math::exp(-x)) - uni20::math::log(-uni20::math::expm1(-x))) / Real{2};
}
template <uni20::Real Real> Real wrap(Real p)
{
  Real const half = pi(p), period = Real{2} * half;
  if (p < -half) p += period;
  if (p >= half) p -= period;
  return p;
}
template <uni20::Real Real> void validate(Real interaction, Options<Real> const& o, Convention convention)
{
  Real const zero = uni20::scalar_at_precision<Real>(0, o.precision());
  bethe::detail::require_same_precision(interaction, zero);
  bethe::detail::require_same_precision(o.relative_tolerance, zero);
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
          s(Real{2} * u / pi(interaction)), zero(uni20::scalar_like(interaction, 0))
    {}

    IntegralPoint<Real> spinon(Real lambda)
    {
      Real const l = uni20::math::abs(lambda), p = pi(u);
      IntegralPoint<Real> out{zero, zero, zero, zero};
      if (!(u > Real{0}) || !(s > Real{0})) return out;
      if (u >= Real{1} / Real{4})
      {
        auto f = [&](uni20::quadrature_point<Real> const& point) {
          Real const& theta = point.value;
          Real const sine = uni20::math::sin(theta), x = (l - uni20::math::cos(theta)) / s;
          return std::array<Real, 2>{sine * sine * sech(x) / u, Real{2} / p * atan_exp(-x)};
        };
        auto const integral = integrate(f, zero, p);
        out = {integral.value[0], integral.value[1], integral.error[0], integral.error[1], integral.converged};
      }
      else
      {
        // Fourier-convolution form: nonoscillatory even as U tends to zero.
        // Split at x=+-1, where sqrt(1-x*x) and acos(x) change support.
        Real const cut = -uni20::math::log(tolerance) - uni20::math::log(u) / Real{2} + Real{12};
        std::vector<Real> breaks;
        bool const outside = l >= Real{1};
        Real const offset = outside ? (l - Real{1}) / s : zero;
        Real const prefactor = outside ? uni20::math::exp(-offset) : uni20::scalar_like(zero, 1);
        if (prefactor == Real{0}) return out;
        breaks = {outside ? zero : -cut, cut};
        for (Real edge : outside ? std::array<Real, 2>{zero, Real{2} / s}
                                 : std::array<Real, 2>{(-Real{1} - l) / s, (Real{1} - l) / s})
          if (edge > breaks.front() && edge < cut) breaks.push_back(edge);
        std::sort(breaks.begin(), breaks.end());
        using Sum = bethe::detail::CompensatedSum<Real>;
        std::array<Sum, 2> sums{Sum(zero), Sum(zero)}, errors{Sum(zero), Sum(zero)};
        bool complete = true;
        for (std::size_t j = 1; j < breaks.size(); ++j)
        {
          Real const lo = breaks[j - 1], hi = breaks[j];
          auto f = [&](uni20::quadrature_point<Real> const& point) {
            Real const& y = point.value;
            Real const x = outside ? Real{1} - s * y : l + s * y;
            Real const weight =
                outside ? Real{2} * uni20::math::exp(-y) / (Real{1} + uni20::math::exp(-Real{2} * (offset + y)))
                        : sech(y);
            Real plus = Real{1} + x, minus = Real{1} - x;
            // The support edges are explicit interval boundaries. Use the
            // independent quadrature distances there: x may have rounded to
            // +/-1 even though the sampled point is strictly inside support.
            Real const left_edge = outside ? zero : (-Real{1} - l) / s;
            Real const right_edge = outside ? Real{2} / s : (Real{1} - l) / s;
            auto edge_distance = [&](Real const& edge) -> Real {
              if (edge == lo) return point.lower_distance;
              if (edge == hi) return -point.upper_distance;
              return y - edge;
            };
            if (left_edge == lo || left_edge == hi)
            {
              if (outside)
                minus = s * edge_distance(left_edge);
              else
                plus = s * edge_distance(left_edge);
            }
            if (right_edge == lo || right_edge == hi)
            {
              if (outside)
                plus = -s * edge_distance(right_edge);
              else
                minus = -s * edge_distance(right_edge);
            }
            Real const angle = minus <= zero  ? zero
                               : plus <= zero ? p
                               : minus < plus ? Real{2} * uni20::math::asin(uni20::math::sqrt(minus / Real{2}))
                                              : p - Real{2} * uni20::math::asin(uni20::math::sqrt(plus / Real{2}));
            return std::array<Real, 2>{Real{2} / p * weight * uni20::math::sqrt(std::max(zero, plus * minus)),
                                       weight * angle / p};
          };
          auto const integral = integrate(f, lo, hi);
          complete = complete && integral.converged;
          for (std::size_t a = 0; a < 2; ++a)
          {
            sums[a].add(integral.value[a]);
            errors[a].add(integral.error[a]);
          }
          if (!complete) break;
        }
        // Both convolution integrands are bounded by constants times sech.
        Real const tail = Real{8} * uni20::math::exp(-cut);
        out = {prefactor * sums[0].value(), prefactor * sums[1].value(), prefactor * (errors[0].value() + tail),
               prefactor * (errors[1].value() + tail), complete};
      }
      if (l == Real{0})
      {
        out.momentum = p / Real{2};
        out.momentum_error = zero;
      }
      else if (lambda < Real{0})
        out.momentum = p - out.momentum;
      out.converged = out.converged && out.energy > Real{0} && uni20::isfinite(out.energy) &&
                      out.energy_error <= o.relative_tolerance * out.energy;
      return out;
    }

    IntegralPoint<Real> charge(Real k)
    {
      Real const p = pi(u);
      IntegralPoint<Real> out{zero, zero, zero, zero};
      if (!(s > Real{0})) return out;
      Real const a = Real{1} / s, cut = -uni20::math::log(tolerance) + Real{12};
      Real const ratio = cut * s;
      if (!uni20::isfinite(a) || !uni20::isfinite(ratio)) return out;
      Real const end = Real{2} * uni20::math::asinh(uni20::math::sqrt(ratio / Real{2}));
      Real const sm = uni20::math::sin((k - p / Real{2}) / Real{2});
      Real const sp = uni20::math::sin((k + p / Real{2}) / Real{2});
      Real dm = Real{2} * (sm * sm);
      Real dp = Real{2} * (sp * sp);
      if (k == Real{0} || uni20::math::abs(k) == p) dm = dp = uni20::scalar_like(zero, 1);
      bool const outer = uni20::math::abs(k) < p / Real{2};
      auto f = [&](uni20::quadrature_point<Real> const& point) {
        Real const& t = point.value;
        Real const sh = uni20::math::sinh(t / Real{2}), h = Real{2} * sh * sh;
        Real const am = atanh_exp_minus(a * (dm + h)), ap = atanh_exp_minus(a * (dp + h));
        return std::array<Real, 2>{Real{4} / p * uni20::math::cosh(t) * (am + ap), Real{2} / p * (am - ap)};
      };
      // Resummation of Melzer's positive K_0/K_1 series; no Bessel library
      // or cancellation of order-one terms to obtain the Mott gap.
      auto const integral = integrate(f, zero, end, std::array<Real, 2>{zero, uni20::scalar_like(zero, 1)});
      Real const sine = uni20::math::sqrt(ratio) * uni20::math::sqrt(Real{2} + ratio);
      Real const tail = uni20::math::exp(-a * std::min(dm, dp) - cut);
      out.energy = integral.value[0] + (outer ? Real{4} * uni20::math::cos(k) : zero);
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
    Real const tolerance, s, zero;

  private:
    bool quadrature_failed = false;
    template <typename Function>
    auto integrate(Function const& f, Real lo, Real hi, std::optional<std::array<Real, 2>> scales = {})
    {
      uni20::tanh_sinh_options<Real, 2> controls(o.precision());
      controls.relative_tolerance = tolerance;
      controls.max_evaluations = o.max_evaluations - evaluations;
      controls.max_levels = o.max_levels;
      if (scales)
        for (std::size_t i = 0; i < 2; ++i)
          controls.absolute_tolerance[i] = tolerance * (*scales)[i];
      auto const integral = uni20::tanh_sinh_integrate<Real, 2>(f, lo, hi, controls);
      evaluations += integral.evaluations;
      bool const converged = integral.status == uni20::quadrature_status::converged;
      quadrature_failed = quadrature_failed || !converged;
      return bethe::detail::QuadratureResult<Real, 2>{integral.value, integral.error, converged};
    }
};

template <uni20::Real Real> Point<Real> initialize(Branch branch, Real interaction, Convention convention)
{
  Point<Real> out;
  Real const zero = uni20::scalar_like(interaction, 0);
  out.momentum = out.parameter = out.energy_offset = out.energy_error = out.momentum_error = zero;
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
  out.energy_error += Real{4} * uni20::numeric_limits<Real>::epsilon(out.interaction) *
                      (uni20::math::abs(*out.energy) + uni20::math::abs(out.energy_offset));
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
  bethe::detail::require_same_precision(interaction, lambda);
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
  bethe::detail::require_same_precision(interaction, k);
  if (branch == Branch::spinon) throw std::invalid_argument("charge branch must be holon or antiholon");
  if (!uni20::isfinite(k) || uni20::math::abs(k) > detail::pi(interaction))
    throw std::invalid_argument("bare momentum must be in [-pi,pi]");
  auto out = detail::initialize<Real>(branch, interaction, convention);
  out.parameter = k;
  detail::Evaluator<Real> eval(interaction, options);
  auto const integral = eval.charge(k);
  out.status = eval.failure_status();
  out.momentum = detail::wrap(
      integral.momentum - (branch == Branch::antiholon ? detail::pi(interaction) : uni20::scalar_like(interaction, 0)));
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
  bethe::detail::require_same_precision(interaction, momentum);
  auto out = detail::initialize<Real>(branch, interaction, convention);
  Real const zero = uni20::scalar_like(interaction, 0), one = uni20::scalar_like(interaction, 1);
  Real const p = detail::pi(interaction);
  out.momentum = momentum;
  if (!uni20::isfinite(momentum) || momentum > p || momentum < (branch == Branch::spinon ? Real{0} : -p))
    throw std::invalid_argument("dressed momentum must be in [0,pi] for spinons, [-pi,pi] for charge branches");
  if (branch == Branch::spinon && (momentum == Real{0} || momentum == p))
  {
    out.parameter = uni20::scalar_like(zero, momentum == zero ? std::numeric_limits<double>::infinity()
                                                              : -std::numeric_limits<double>::infinity());
    detail::finish(out, {zero, momentum, zero, zero, true});
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
  Real lo = spin ? zero : -p, hi = spin ? Real{1} + eval.s * uni20::math::log(Real{4} / target) : p;
  Real flo = spin ? uni20::math::log((p / Real{2}) / target) : Real{3} * p / Real{2} - target;
  Real fhi = spin ? -one : -p / Real{2} - target;
  Real x = spin ? (target == p / Real{2} ? zero : (lo + hi) / Real{2}) : p / Real{2} - target;
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
    Real const momentum_tolerance = options.relative_tolerance * (spin ? target : one);
    out.momentum_error = uni20::math::abs(residual) + integral.momentum_error;
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
    Real const f = spin ? uni20::math::log(integral.momentum / target) : residual;
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
    if (!(fraction > Real{0} && fraction < Real{1})) fraction = one / Real{2};
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
