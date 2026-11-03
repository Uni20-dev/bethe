// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <array>
#include <bethe/detail/newton.hpp>
#include <bethe/detail/newton_backtracking.hpp>
#include <bethe/experimental/xxx_complex_roots.hpp>
#include <optional>
#include <string>

namespace bethe::experimental::xxx::six_spinon
{
enum class Family
{
  two_pairs,
  triple
};

struct Configuration
{
    int sites;
    Family family;
    std::vector<uni20::half_int> holes;
    std::vector<uni20::half_int> branches;
};

struct Labels
{
    std::vector<int> lengths;
    std::vector<uni20::half_int> numbers;
    int momentum_index{};
};

inline Labels labels(Configuration const& c)
{
  if (c.sites < 12 || c.sites % 4 || (c.family != Family::two_pairs && c.family != Family::triple))
    throw std::invalid_argument("require N>=12 divisible by four and family 22 or 3");
  int const edge2 = c.sites / 2 + (c.family == Family::triple ? 1 : 0);
  auto holes = c.holes, branches = c.branches;
  std::sort(holes.begin(), holes.end());
  std::sort(branches.begin(), branches.end());
  if (holes.size() != 6 || std::adjacent_find(holes.begin(), holes.end()) != holes.end())
    throw std::invalid_argument("require six distinct holes");
  for (auto h : holes)
    if (h.twice() < -edge2 || h.twice() > edge2 || (h.twice() + edge2) % 2)
      throw std::invalid_argument("hole outside family counting window");
  bool const pairs = c.family == Family::two_pairs;
  if (branches.size() != (pairs ? 2U : 1U) || std::adjacent_find(branches.begin(), branches.end()) != branches.end())
    throw std::invalid_argument("invalid string quantum numbers");
  for (auto b : branches)
    if ((pairs && (b.twice() < -3 || b.twice() > 3 || b.is_integral())) ||
        (!pairs && (b.twice() < -2 || b.twice() > 2 || !b.is_integral())))
      throw std::invalid_argument("string quantum number outside family window");
  Labels result;
  for (int t = -edge2; t <= edge2; t += 2)
  {
    auto const i = uni20::half_int(t, uni20::half_int::twice_tag{});
    if (std::find(holes.begin(), holes.end(), i) == holes.end())
    {
      result.numbers.push_back(i);
      result.lengths.push_back(1);
    }
  }
  for (auto b : branches)
  {
    result.numbers.push_back(b);
    result.lengths.push_back(pairs ? 2 : 3);
  }
  // Both triplet families have N/2-3 strings: P=pi*Nstrings-2*pi*sum(I)/N.
  // N divisible by four has a zero-momentum ground state.
  long long twice_sum{};
  for (auto i : result.numbers)
    twice_sum += i.twice();
  result.momentum_index = int((c.sites / 2 - twice_sum / 2) % c.sites);
  if (result.momentum_index < 0) result.momentum_index += c.sites;
  return result;
}

// Deguchi-Giri (2015), Table I: these two N=12 three-string labels require
// self-conjugacy of the entire root set, not of each individual string.
// This is an explicit, bounded research mapping, not a general fallback.
inline int nonself_conjugate_sign(Configuration const& c, Labels const& l)
{
  if (c.sites != 12 || c.family != Family::triple || l.numbers.size() != 3) return 0;
  if (l.numbers[0].twice() == 3 && l.numbers[1].twice() == 5 && l.numbers[2].twice() == 2) return 1;
  if (l.numbers[0].twice() == -5 && l.numbers[1].twice() == -3 && l.numbers[2].twice() == -2) return -1;
  return 0;
}

enum class SolveStatus
{
  converged,
  seed_failed,
  unresolved_geometry,
  iteration_failed,
  residual_failed,
  momentum_mismatch,
  label_mismatch,
  energy_unresolved
};

template <uni20::Real R> struct State
{
    std::vector<uni20::complex<R>> roots; // Standard ABA lambda, with finite deviations.
    std::vector<R> variables;
    std::size_t real_count{};
    int momentum_index{};
    R energy, residual, pole_separation, label_error, momentum_error;
    bool nonself_conjugate{};
};

/// A failed solve preserves its request and attempted coordinates, with no usable state.
template <uni20::Real R> struct SolveResult
{
    Configuration request;
    SolveStatus status = SolveStatus::iteration_failed;
    std::optional<State<R>> state = {};
    std::vector<R> attempted_variables = {};
    std::size_t iterations = 0;
    std::string diagnostic = {};
    bool converged() const { return status == SolveStatus::converged && state.has_value(); }
};

namespace detail
{
struct SolveFailure : std::runtime_error
{
    SolveStatus status;
    SolveFailure(SolveStatus code, char const* message) : std::runtime_error(message), status(code) {}
};

template <uni20::Real R> R scattering(int a, int b, R x, bool derivative)
{
  R const zero = bethe::detail::precision_zero(x);
  auto theta = [&](int k) {
    return derivative
               ? uni20::scalar_like(zero, k) / (x * x + uni20::scalar_like(zero, k * k) / uni20::scalar_like(zero, 4))
               : uni20::scalar_like(zero, 2) *
                     uni20::math::atan(uni20::scalar_like(zero, 2) * x / uni20::scalar_like(zero, k));
  };
  R value = theta(a + b);
  if (a != b) value += theta(std::abs(a - b));
  for (int k = std::abs(a - b) + 2; k < a + b; k += 2)
    value += uni20::scalar_like(zero, 2) * theta(k);
  return value;
}

template <uni20::Real R>
std::vector<R> takahashi(int n, Labels const& labels, std::span<R const> centres, R const& zero,
                         bool derivative = false)
{
  R const one = uni20::scalar_like(zero, 1);
  auto const s = centres.size();
  std::vector<R> out(derivative ? s * s : s, zero);
  R const pi = uni20::scalar_like(zero, 4) * uni20::math::atan(one);
  for (std::size_t i = 0; i < s; ++i)
  {
    int const a = labels.lengths[i];
    complex_roots::Sum<R> value(zero);
    value.add(derivative
                  ? (uni20::scalar_like(zero, n) * a) /
                        (centres[i] * centres[i] + uni20::scalar_like(zero, a * a) / uni20::scalar_like(zero, 4))
                  : uni20::scalar_like(zero, 2) * uni20::scalar_like(zero, n) *
                            uni20::math::atan(uni20::scalar_like(zero, 2) * centres[i] / uni20::scalar_like(zero, a)) -
                        pi * uni20::scalar_like(zero, labels.numbers[i].twice()));
    for (std::size_t j = 0; j < s; ++j)
      if (i != j)
      {
        R const term = scattering(a, labels.lengths[j], centres[i] - centres[j], derivative);
        value.add(-term);
        if (derivative) out[i * s + j] = term;
      }
    out[derivative ? i * s + i : i] = value.value();
  }
  return out;
}

template <uni20::Real R> R norm(std::vector<R> const& values, R const& zero)
{
  R result = zero;
  for (R const& v : values)
  {
    if (!uni20::isfinite(v)) return uni20::scalar_like(zero, std::numeric_limits<double>::infinity());
    result = std::max(result, uni20::math::abs(v));
  }
  return result;
}

template <uni20::Real R> std::vector<R> ideal_centres(Configuration const& c, R const& zero)
{
  R const one = uni20::scalar_like(zero, 1);
  auto const l = labels(c);
  std::vector<R> x;
  R const pi = uni20::scalar_like(zero, 4) * uni20::math::atan(one),
          tolerance = uni20::scalar_like(zero, 128) * uni20::scalar_like(zero, c.sites) *
                      uni20::numeric_limits<R>::epsilon(zero);
  for (std::size_t i = 0; i < l.numbers.size(); ++i)
    x.push_back(uni20::scalar_like(zero, l.lengths[i]) *
                uni20::math::tan(pi * uni20::scalar_like(zero, l.numbers[i].twice()) / uni20::scalar_like(zero, 2) /
                                 uni20::scalar_like(zero, c.sites)) /
                uni20::scalar_like(zero, 2));
  for (int iteration = 0; iteration < 100; ++iteration)
  {
    auto rhs = takahashi<R>(c.sites, l, x, zero);
    R const error = norm(rhs, zero);
    if (error <= tolerance) return x;
    for (auto& v : rhs)
      v = -v;
    if (!bethe::detail::newton_step(takahashi<R>(c.sites, l, x, zero, true), rhs, uni20::precision_of(zero))) break;
    if (!bethe::detail::backtrack_newton(
            x, uni20::precision_of(zero), [&](std::size_t i) { return rhs[i]; },
            [&](auto const& trial, R damping) {
              return bethe::detail::newton_decreases(norm(takahashi<R>(c.sites, l, trial, zero), zero), error, damping,
                                                     tolerance);
            }))
      break;
  }
  throw SolveFailure(SolveStatus::seed_failed, "six-spinon Takahashi seed failed");
}

template <uni20::Real R> std::vector<uni20::complex<R>> unpack(std::span<R const> x, std::size_t reals)
{
  std::vector<uni20::complex<R>> roots(x.begin(), x.begin() + reals);
  for (std::size_t i = reals; i < x.size(); i += 2)
  {
    roots.emplace_back(x[i], x[i + 1]);
    roots.emplace_back(x[i], -x[i + 1]);
  }
  return roots;
}

// x: real roots, then (pair centre, positive pair height) for each conjugate
// pair. For a three-string, the central real root is the last real variable.
// No exact string is ever substituted for the full resolved equations.
template <uni20::Real R> std::vector<R> equations(int n, std::size_t reals, std::span<R const> x, R const& zero)
{
  auto const roots = unpack<R>(x, reals);
  std::vector<R> out(x.size(), zero);
  for (std::size_t i = 0; i < x.size(); ++i)
  {
    auto const l = complex_roots::log_equation<R>(n, roots, i, zero);
    if (i < reals)
      out[i] = complex_roots::phase(l.imag());
    else
    {
      out[i] = l.real();
      out[i + 1] = complex_roots::phase(l.imag());
      ++i;
    }
  }
  return out;
}

template <uni20::Real R> std::vector<R> jacobian(int n, std::size_t reals, std::span<R const> x, R const& zero)
{
  R const one = uni20::scalar_like(zero, 1);
  using C = uni20::complex<R>;
  auto const roots = unpack<R>(x, reals);
  auto const g = complex_roots::gaudin<R>(n, roots, zero);
  auto const s = x.size();
  std::vector<R> result(s * s, zero);
  for (std::size_t i = 0; i < s; ++i)
  {
    for (std::size_t j = 0; j < s; ++j)
    {
      C derivative(zero, zero);
      if (j < reals)
        derivative = C(zero, (-one)) * g[i * s + j];
      else if ((j - reals) % 2 == 0)
        derivative = C(zero, (-one)) * (g[i * s + j] + g[i * s + j + 1]);
      else
        derivative = g[i * s + j - 1] - g[i * s + j];
      if (i < reals)
        result[i * s + j] = derivative.imag();
      else
      {
        result[i * s + j] = derivative.real();
        result[(i + 1) * s + j] = derivative.imag();
      }
    }
    if (i >= reals) ++i;
  }
  return result;
}

// Recover Takahashi labels from the full Bethe quantum numbers, including
// branch cuts: Hagemans-Caux (2007), Eqs. (25) and (31). A momentum check
// alone cannot distinguish different configurations in the same sector.
template <uni20::Real R>
R quantum_number_error(Configuration const& c, std::span<uni20::complex<R> const> roots, R const& zero)
{
  R const one = uni20::scalar_like(zero, 1);
  auto const l = labels(c);
  if (roots.size() != static_cast<std::size_t>(c.sites / 2 - 1))
    throw std::invalid_argument("six-spinon root count differs from request");
  std::size_t const reals = std::count(l.lengths.begin(), l.lengths.end(), 1);
  R const pi = uni20::scalar_like(zero, 4) * uni20::math::atan(one);
  std::vector<R> j;
  for (std::size_t i = 0; i < roots.size(); ++i)
  {
    complex_roots::Sum<uni20::complex<R>> value(uni20::complex<R>(zero, zero));
    value.add(uni20::scalar_like(zero, 2) * uni20::scalar_like(zero, c.sites) *
              complex_roots::bethe_atan<R>(uni20::scalar_like(zero, 2) * roots[i]));
    for (std::size_t k = 0; k < roots.size(); ++k)
      if (k != i) value.add(-uni20::scalar_like(zero, 2) * complex_roots::bethe_atan<R>(roots[i] - roots[k]));
    R const number = value.value().real() / (uni20::scalar_like(zero, 2) * pi);
    if (!uni20::isfinite(number)) return uni20::scalar_like(zero, std::numeric_limits<double>::infinity());
    j.push_back(number);
  }
  auto sign = [&](R x) { return x > zero ? one : x < zero ? (-one) : zero; };
  R error = zero;
  if (int const reflection = nonself_conjugate_sign(c, l))
  {
    if (roots.size() != 5) throw std::invalid_argument("non-self-conjugate mapping requires five roots");
    // Canonical positive-momentum order: real, inner upper/lower, outer
    // upper/lower. Table I assigns inner lower to a one-string and inner
    // upper to the three-string. Its limiting corrections replace the
    // ambiguous sign of the equal real parts of the inner conjugate pair.
    std::array<std::size_t, 5> const order =
        reflection > 0 ? std::array<std::size_t, 5>{0, 1, 2, 3, 4} : std::array<std::size_t, 5>{0, 2, 1, 4, 3};
    std::array<R, 5> canonical;
    for (std::size_t i = 0; i < 5; ++i)
      canonical[i] = uni20::scalar_like(zero, reflection) * j[order[i]];
    std::array<R, 2> singles{uni20::scalar_like(zero, reflection) * (canonical[0] - (one / 2)),
                             uni20::scalar_like(zero, reflection) * (canonical[2] - (one / 2))};
    std::sort(singles.begin(), singles.end());
    for (std::size_t i = 0; i < 2; ++i)
      error = std::max(error, uni20::math::abs(singles[i] - uni20::scalar_like(zero, l.numbers[i].twice()) /
                                                                uni20::scalar_like(zero, 2)));
    R const triple = uni20::scalar_like(zero, reflection) *
                     (canonical[1] + canonical[3] + canonical[4] - uni20::scalar_like(zero, c.sites) + one);
    return std::max(error, uni20::math::abs(uni20::math::remainder(
                               triple - uni20::scalar_like(zero, l.numbers[2].twice()) / uni20::scalar_like(zero, 2),
                               uni20::scalar_like(zero, c.sites))));
  }
  for (std::size_t i = 0; i < reals; ++i)
  {
    R const correction =
        c.family == Family::triple ? sign(roots[i].real() - roots[reals].real()) / uni20::scalar_like(zero, 2) : zero;
    error =
        std::max(error, uni20::math::abs(j[i] + correction -
                                         uni20::scalar_like(zero, l.numbers[i].twice()) / uni20::scalar_like(zero, 2)));
  }
  if (c.family == Family::two_pairs)
  {
    for (std::size_t p = 0; p < 2; ++p)
    {
      std::size_t const i = reals + 2 * p, other = reals + 2 * (1 - p);
      R const number = (roots[i].imag() > (one / 2) ? one : zero) + uni20::scalar_like(zero, 2) * j[i] -
                       uni20::scalar_like(zero, c.sites) / uni20::scalar_like(zero, 2) +
                       sign(roots[i].real() - roots[other].real()) / uni20::scalar_like(zero, 2);
      error = std::max(
          error, uni20::math::abs(uni20::math::remainder(
                     number - uni20::scalar_like(zero, l.numbers[reals + p].twice()) / uni20::scalar_like(zero, 2),
                     uni20::scalar_like(zero, c.sites))));
    }
  }
  else
  {
    R number = j[reals] + j[reals + 1] + j[reals + 2];
    for (std::size_t i = 0; i < reals; ++i)
      number += sign(roots[reals].real() - roots[i].real()) / uni20::scalar_like(zero, 2);
    error =
        std::max(error, uni20::math::abs(uni20::math::remainder(
                            number - uni20::scalar_like(zero, l.numbers.back().twice()) / uni20::scalar_like(zero, 2),
                            uni20::scalar_like(zero, c.sites))));
  }
  return error;
}

template <uni20::Real R>
State<R> solve_state(Configuration const& c, R const& zero, std::vector<R>& x, std::size_t& iterations)
{
  R const one = uni20::scalar_like(zero, 1);
  using C = uni20::complex<R>;
  auto const l = labels(c);
  auto const seed = ideal_centres<R>(c, zero);
  int const nonself = nonself_conjugate_sign(c, l);
  std::size_t const real_strings = std::count(l.lengths.begin(), l.lengths.end(), 1);
  std::size_t const real_count = nonself ? 1 : real_strings + (c.family == Family::triple ? 1 : 0);
  x.assign(seed.begin(), seed.begin() + real_strings);
  if (nonself)
  {
    // Coarse root-geometry seed only. Full roots and determinants are solved
    // afresh in R; no published energy or weight enters the calculation.
    x = {uni20::scalar_like(zero, nonself) * uni20::scalar_like(zero, 18) / uni20::scalar_like(zero, 100),
         uni20::scalar_like(zero, nonself) * uni20::scalar_like(zero, 445) / uni20::scalar_like(zero, 1000),
         uni20::scalar_like(zero, 2) / uni20::scalar_like(zero, 100),
         uni20::scalar_like(zero, nonself) * uni20::scalar_like(zero, 492) / uni20::scalar_like(zero, 1000),
         uni20::scalar_like(zero, 96) / uni20::scalar_like(zero, 100)};
  }
  else if (c.family == Family::two_pairs)
  {
    for (std::size_t p = 0; p < 2; ++p)
    {
      R const v = seed[real_strings + p], other = seed[real_strings + 1 - p];
      std::vector<C> roots(seed.begin(), seed.begin() + real_strings);
      roots.emplace_back(other, (one / 2));
      roots.emplace_back(other, (-one / 2));
      roots.emplace_back(v, (one / 2));
      C const q = uni20::math::exp(complex_roots::log_equation<R>(c.sites, roots, roots.size() - 1, zero));
      x.push_back(v);
      x.push_back((one / 2) + (one / (q - one)).real());
    }
  }
  else
  {
    R const v = seed.back();
    std::vector<C> roots(seed.begin(), seed.begin() + real_strings);
    roots.emplace_back(v, (-one));
    roots.emplace_back(v, one);
    C const q = uni20::math::exp(complex_roots::log_equation<R>(c.sites, roots, roots.size() - 1, zero));
    C const displacement = C(zero, uni20::scalar_like(zero, 2)) / (q - one);
    x.push_back(v);
    x.push_back(v + displacement.real());
    x.push_back(one + displacement.imag());
  }
  auto domain = [&](std::vector<R> const& trial) {
    for (auto const& value : trial)
      if (!uni20::isfinite(value)) return false;
    if (nonself)
    {
      if (!(trial[2] > zero && trial[2] < (one / 2) && trial[4] > (one / 2) && trial[4] < (3 * one / 2))) return false;
      auto const roots = unpack<R>(trial, real_count);
      return complex_roots::separation<R>(roots, zero) >
             uni20::scalar_like(zero, 64) * uni20::math::sqrt(uni20::numeric_limits<R>::epsilon(zero));
    }
    for (std::size_t i = real_count + 1; i < trial.size(); i += 2)
      if (!(trial[i] > (c.family == Family::two_pairs ? zero : (one / 2))) ||
          trial[i] >= (c.family == Family::two_pairs ? one : (3 * one / 2)))
        return false;
    if (c.family == Family::two_pairs && !(trial[real_count] < trial[real_count + 2])) return false;
    for (std::size_t i = 1; i < real_strings; ++i)
      if (!(trial[i - 1] < trial[i])) return false;
    auto const roots = unpack<R>(trial, real_count);
    return complex_roots::separation<R>(roots, zero) >
           uni20::scalar_like(zero, 64) * uni20::math::sqrt(uni20::numeric_limits<R>::epsilon(zero));
  };
  if (!domain(x))
    throw SolveFailure(SolveStatus::unresolved_geometry,
                       "six-spinon seed unresolved; increase precision or reject state");
  for (iterations = 0; iterations < 100; ++iterations)
  {
    auto const roots = unpack<R>(x, real_count);
    R const tolerance = complex_roots::residual_tolerance<R>(c.sites, roots, zero);
    auto rhs = equations<R>(c.sites, real_count, x, zero);
    R const error = norm(rhs, zero);
    if (error <= tolerance)
    {
      R const residual = complex_roots::rational_residual<R>(c.sites, roots, zero);
      if (!uni20::isfinite(residual) || residual > uni20::scalar_like(zero, 4) * tolerance)
        throw SolveFailure(SolveStatus::residual_failed, "full six-spinon Bethe residual failed");
      complex_roots::Sum<C> total(C(zero, zero));
      for (auto z : roots)
        total.add(uni20::math::log(z + C(zero, (one / 2))) - uni20::math::log(z - C(zero, (one / 2))));
      R const pi = uni20::scalar_like(zero, 4) * uni20::math::atan(one);
      R const momentum_error = uni20::math::abs(complex_roots::phase(
          total.value().imag() - uni20::scalar_like(zero, 2) * pi * uni20::scalar_like(zero, l.momentum_index) /
                                     uni20::scalar_like(zero, c.sites)));
      if (!uni20::isfinite(momentum_error) || momentum_error > uni20::scalar_like(zero, 4) * tolerance)
        throw SolveFailure(SolveStatus::momentum_mismatch, "six-spinon state left expected momentum sector");
      R const label_error = quantum_number_error<R>(c, roots, zero);
      if (!uni20::isfinite(label_error) || label_error > uni20::scalar_like(zero, 4) * tolerance)
        throw SolveFailure(SolveStatus::label_mismatch, "six-spinon state left requested quantum numbers");
      R energy = zero;
      try
      {
        energy = complex_roots::energy<R>(c.sites, roots, zero);
      }
      catch (std::runtime_error const&)
      {
        throw SolveFailure(SolveStatus::energy_unresolved, "nonreal or unresolved energy");
      }
      return {roots,
              std::move(x),
              real_count,
              l.momentum_index,
              energy,
              residual,
              complex_roots::separation<R>(roots, zero),
              label_error,
              momentum_error,
              nonself != 0};
    }
    for (auto& v : rhs)
      v = -v;
    if (!bethe::detail::newton_step(jacobian<R>(c.sites, real_count, x, zero), rhs, uni20::precision_of(zero))) break;
    if (!bethe::detail::backtrack_newton(
            x, uni20::precision_of(zero), [&](std::size_t i) { return rhs[i]; },
            [&](auto const& trial, R damping) {
              return domain(trial) &&
                     bethe::detail::newton_decreases(norm(equations<R>(c.sites, real_count, trial, zero), zero), error,
                                                     damping, tolerance);
            }))
      break;
  }
  throw SolveFailure(SolveStatus::iteration_failed, "resolved six-spinon Newton solve failed");
}
} // namespace detail

/// Solve a selected triplet with resolved finite-size deviations, at the exemplar's
/// finite precision. Invalid requests/precision throw; numerical failure returns
/// an incomplete result. No claim of complete six-spinon-family coverage is made.
template <uni20::Real R> SolveResult<R> solve(Configuration const& request, R const& exemplar)
{
  (void)labels(request);
  R const zero = bethe::detail::precision_zero(exemplar);
  if (!uni20::isfinite(exemplar)) throw std::invalid_argument("finite working exemplar required");
  SolveResult<R> result{request};
  try
  {
    result.state = detail::solve_state<R>(request, zero, result.attempted_variables, result.iterations);
    result.status = SolveStatus::converged;
  }
  catch (detail::SolveFailure const& failure)
  {
    result.status = failure.status;
    result.diagnostic = failure.what();
  }
  return result;
}

template <uni20::Real R>
  requires(!uni20::has_runtime_precision_v<R>)
SolveResult<R> solve(Configuration const& request)
{
  return solve(request, R{0});
}

/// Revalidate original equations, family labels and momentum before evaluating
/// the physically normalized Fourier raising weight. Failed solves have no weight.
template <uni20::Real R>
heisenberg::FormFactor<R> raising_weight(PreparedSinglet<R> const& ground, SolveResult<R> const& solution)
{
  (void)labels(solution.request);
  if (ground.sites() != static_cast<std::size_t>(solution.request.sites))
    throw std::invalid_argument("six-spinon state and ground have different sizes");
  heisenberg::FormFactor<R> failure;
  failure.status = heisenberg::FormFactorStatus::roots_unconverged;
  if (!solution.converged()) return failure;
  if (!ground.ready())
  {
    failure.status = ground.status();
    return failure;
  }
  R const zero = *ground.working_zero(), one = uni20::scalar_like(zero, 1);
  auto const& state = *solution.state;
  auto const& request = solution.request;
  if (state.roots.size() != static_cast<std::size_t>(request.sites / 2 - 1))
    throw std::invalid_argument("six-spinon root count differs from request");
  // Shared helpers validate finite values and common precision before arithmetic.
  R const separation = complex_roots::separation<R>(state.roots, zero);
  R const eps = uni20::numeric_limits<R>::epsilon(zero);
  if (!(separation > 64 * uni20::math::sqrt(eps)))
  {
    failure.status = heisenberg::FormFactorStatus::precision_limit;
    return failure;
  }
  R const tolerance = 4 * complex_roots::residual_tolerance<R>(request.sites, state.roots, zero);
  R const label_error = detail::quantum_number_error<R>(request, state.roots, zero);
  if (!uni20::isfinite(label_error) || label_error > tolerance) return failure;
  complex_roots::Sum<uni20::complex<R>> momentum(uni20::complex<R>(zero, zero));
  uni20::complex<R> const imaginary(zero, one);
  for (auto const& root : state.roots)
    momentum.add(uni20::math::log(root + imaginary / (2 * one)) - uni20::math::log(root - imaginary / (2 * one)));
  R const error = uni20::math::abs(complex_roots::phase(
      momentum.value().imag() - 2 * bethe::detail::pi_like(zero) * labels(request).momentum_index / request.sites));
  if (!uni20::isfinite(error) || error > tolerance) return failure;
  return complex_roots::raising_weight<R>(ground, state.roots, labels(request).momentum_index);
}
} // namespace bethe::experimental::xxx::six_spinon
