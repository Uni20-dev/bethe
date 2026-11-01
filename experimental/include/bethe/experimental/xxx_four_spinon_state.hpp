// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <array>
#include <bethe/detail/newton.hpp>
#include <bethe/detail/newton_backtracking.hpp>
#include <bethe/experimental/xxx_complex_roots.hpp>
#include <memory>
#include <string>

namespace bethe::experimental::xxx
{
struct FourSpinonRequest
{
    int sites{};
    std::array<int, 4> holes{};
    int branch{};
    bool operator==(FourSpinonRequest const&) const = default;
};

enum class FourSpinonStatus
{
  converged,
  unsupported_geometry,
  unconverged,
  precision_limit
};

namespace four_detail
{
using complex_roots::C;
using complex_roots::Sum;

inline FourSpinonRequest canonical_request(FourSpinonRequest request)
{
  if (request.sites < 8 || request.sites % 4 || request.branch < -1 || request.branch > 1)
    throw std::invalid_argument("four-spinons require N>=8 divisible by four and branch -1,0,1");
  std::sort(request.holes.begin(), request.holes.end());
  if (request.holes.front() < -request.sites / 4 || request.holes.back() > request.sites / 4 ||
      std::adjacent_find(request.holes.begin(), request.holes.end()) != request.holes.end())
    throw std::invalid_argument("require four distinct integer holes in [-N/4,N/4]");
  return request;
}

inline std::vector<int> labels(FourSpinonRequest const& request)
{
  std::vector<int> result;
  for (int i = -request.sites / 4; i <= request.sites / 4; ++i)
    if (std::find(request.holes.begin(), request.holes.end(), i) == request.holes.end()) result.push_back(i);
  return result;
}

struct Failure : std::runtime_error
{
    FourSpinonStatus status;
    Failure(FourSpinonStatus value, char const* text) : std::runtime_error(text), status(value) {}
};

template <uni20::Real R> R kernel(int k, R const& x) { return k / (x * x + uni20::scalar_like(x, k * k) / 4); }
template <uni20::Real R> R theta(int k, R const& x) { return 2 * uni20::math::atan(2 * x / k); }
template <uni20::Real R> R norm(std::span<R const> values, R const& zero)
{
  R result = zero;
  for (auto const& value : values)
  {
    if (!uni20::isfinite(value)) return uni20::scalar_like(zero, std::numeric_limits<double>::infinity());
    result = std::max(result, uni20::math::abs(value));
  }
  return result;
}

// Ideal Takahashi equations are used only to seed the finite-deviation solve.
template <uni20::Real R>
std::vector<R> seed_equations(FourSpinonRequest const& request, std::span<int const> numbers, std::span<R const> x,
                              R const& zero)
{
  R const pi = bethe::detail::pi_like(zero);
  std::vector<R> result(x.size(), zero);
  for (std::size_t i = 0; i < numbers.size(); ++i)
  {
    Sum<R> sum(zero);
    sum.add(request.sites * theta(1, x[i]) - 2 * pi * numbers[i]);
    for (auto const& y : x.first(numbers.size()))
      sum.add(-theta(2, x[i] - y));
    sum.add(-theta(1, x[i] - x.back()) - theta(3, x[i] - x.back()));
    result[i] = sum.value() / request.sites;
  }
  Sum<R> sum(zero);
  sum.add(request.sites * theta(2, x.back()) - 2 * pi * request.branch);
  for (auto const& y : x.first(numbers.size()))
    sum.add(-theta(1, x.back() - y) - theta(3, x.back() - y));
  result.back() = sum.value() / request.sites;
  return result;
}

template <uni20::Real R> std::vector<R> seed_jacobian(int n, std::span<R const> x, R const& zero)
{
  auto const s = x.size(), last = s - 1;
  std::vector<R> matrix(s * s, zero);
  for (std::size_t i = 0; i < last; ++i)
  {
    Sum<R> diagonal(zero);
    diagonal.add(n * kernel(1, x[i]));
    for (std::size_t j = 0; j < last; ++j)
      if (i != j)
      {
        matrix[i * s + j] = kernel(2, x[i] - x[j]);
        diagonal.add(-matrix[i * s + j]);
      }
    matrix[i * s + last] = matrix[last * s + i] = kernel(1, x[i] - x[last]) + kernel(3, x[i] - x[last]);
    diagonal.add(-matrix[i * s + last]);
    matrix[i * s + i] = diagonal.value();
  }
  Sum<R> diagonal(zero);
  diagonal.add(n * kernel(2, x[last]));
  for (std::size_t j = 0; j < last; ++j)
    diagonal.add(-matrix[last * s + j]);
  matrix[last * s + last] = diagonal.value();
  for (auto& value : matrix)
    value /= n;
  return matrix;
}

template <uni20::Real R, typename Evaluate, typename Jacobian>
std::vector<R> newton(std::vector<R> x, R const& zero, Evaluate evaluate, Jacobian jacobian)
{
  R const tolerance = 32 * uni20::numeric_limits<R>::epsilon(zero);
  auto rhs = evaluate(x);
  for (int iteration = 0; iteration < 100; ++iteration)
  {
    R const error = norm<R>(rhs, zero);
    if (error <= tolerance) return x;
    if (!uni20::isfinite(error)) throw Failure(FourSpinonStatus::precision_limit, "nonfinite pair equations");
    for (auto& value : rhs)
      value = -value;
    if (!bethe::detail::newton_step(jacobian(x), rhs, uni20::precision_of(zero)))
      throw Failure(FourSpinonStatus::precision_limit, "finite-pair Newton factorization failed");
    std::vector<R> accepted;
    bool const success = bethe::detail::backtrack_newton(
        x, uni20::precision_of(zero), [&](std::size_t i) { return rhs[i]; },
        [&](auto const& trial, R const& damping) {
          auto values = evaluate(trial);
          if (!bethe::detail::newton_decreases(norm<R>(values, zero), error, damping, tolerance)) return false;
          accepted = std::move(values);
          return true;
        });
    if (!success) throw Failure(FourSpinonStatus::unconverged, "finite-pair backtracking failed");
    rhs = std::move(accepted);
  }
  throw Failure(FourSpinonStatus::unconverged, "finite-pair iteration limit");
}

template <uni20::Real R> C<R> drive(int n, C<R> const& z, R const& zero)
{
  R const one = uni20::scalar_like(zero, 1);
  C<R> const half_i(zero, one / 2);
  return (n * one) * (uni20::math::log(z + half_i) - uni20::math::log(z - half_i));
}

template <uni20::Real R> C<R> pair_log(int n, std::span<R const> reals, R const& centre, R const& delta, R const& zero)
{
  R const one = uni20::scalar_like(zero, 1);
  C<R> const u(centre, one / 2 + delta), imaginary(zero, one);
  Sum<C<R>> sum(C<R>(zero, zero));
  sum.add(drive(n, u, zero));
  for (auto const& x : reals)
    sum.add(uni20::math::log(u - x - imaginary) - uni20::math::log(u - x + imaginary));
  return sum.value();
}

template <uni20::Real R> bool in_domain(std::span<R const> x, int sign, R const& zero)
{
  for (auto const& a : x)
    if (!uni20::isfinite(a)) return false;
  R const delta = sign * uni20::math::exp(x.back()), bound = uni20::scalar_like(zero, 49) / 100;
  return uni20::isfinite(delta) && delta > -bound && delta < bound && delta != zero;
}

template <uni20::Real R>
std::vector<R> equations(int n, std::span<int const> numbers, int sign, std::span<R const> x, R const& zero)
{
  auto const r = numbers.size();
  R const one = uni20::scalar_like(zero, 1), pi = bethe::detail::pi_like(zero);
  std::vector<R> result(r + 2, uni20::scalar_like(zero, std::numeric_limits<double>::infinity()));
  if (x.size() != r + 2 || !in_domain<R>(x, sign, zero)) return result;
  R const v = x[r], delta = sign * uni20::math::exp(x.back()), a = one / 2 - delta, b = 3 * one / 2 + delta;
  for (std::size_t i = 0; i < r; ++i)
  {
    Sum<R> sum(zero);
    sum.add(n * theta(1, x[i]) - 2 * pi * numbers[i]);
    for (std::size_t j = 0; j < r; ++j)
      sum.add(-theta(2, x[i] - x[j]));
    sum.add(-2 * uni20::math::atan((x[i] - v) / a) - 2 * uni20::math::atan((x[i] - v) / b));
    result[i] = sum.value() / n;
  }
  auto const value = pair_log<R>(n, x.first(r), v, delta, zero);
  result[r] = (value.real() + x.back() - uni20::math::log1p(delta)) / n;
  result[r + 1] = complex_roots::phase(value.imag() + (sign < 0 ? pi : zero)) / n;
  return result;
}

template <uni20::Real R> std::vector<R> jacobian(int n, int sign, std::span<R const> x, R const& zero)
{
  auto const s = x.size(), r = s - 2;
  R const one = uni20::scalar_like(zero, 1), v = x[r], delta = sign * uni20::math::exp(x.back()), a = one / 2 - delta,
          b = 3 * one / 2 + delta;
  std::vector<R> matrix(s * s, zero);
  for (std::size_t i = 0; i < r; ++i)
  {
    Sum<R> diagonal(zero);
    diagonal.add(n * kernel(1, x[i]));
    for (std::size_t j = 0; j < r; ++j)
      if (i != j)
      {
        matrix[i * s + j] = kernel(2, x[i] - x[j]);
        diagonal.add(-matrix[i * s + j]);
      }
    R const d = x[i] - v, aa = d * d + a * a, bb = d * d + b * b;
    matrix[i * s + r] = 2 * a / aa + 2 * b / bb;
    matrix[i * s + r + 1] = 2 * delta * d * (one / bb - one / aa);
    diagonal.add(-matrix[i * s + r]);
    matrix[i * s + i] = diagonal.value();
  }
  C<R> const u(v, one / 2 + delta), minus_i(zero, -one);
  Sum<C<R>> diagonal(C<R>(zero, zero));
  diagonal.add(n * one / (u * u + one / 4));
  for (std::size_t j = 0; j < r; ++j)
  {
    C<R> const k = 2 * one / (one + (u - x[j]) * (u - x[j]));
    diagonal.add(-k);
    C<R> const derivative = minus_i * k;
    matrix[r * s + j] = derivative.real();
    matrix[(r + 1) * s + j] = derivative.imag();
  }
  C<R> const centre = minus_i * diagonal.value(), deviation = delta * diagonal.value() + one / (one + delta);
  matrix[r * s + r] = centre.real();
  matrix[(r + 1) * s + r] = centre.imag();
  matrix[r * s + r + 1] = deviation.real();
  matrix[(r + 1) * s + r + 1] = deviation.imag();
  for (auto& value : matrix)
    value /= n;
  return matrix;
}

template <uni20::Real R> R exponential_error(C<R> const& value, R const& zero)
{
  R const angle = complex_roots::phase(value.imag()), sine = uni20::math::sin(angle / 2),
          magnitude = uni20::math::exp(value.real());
  R const result = uni20::math::abs(
      C<R>(uni20::math::expm1(value.real()) - 2 * magnitude * sine * sine, magnitude * uni20::math::sin(angle)));
  return uni20::isfinite(result) ? result : uni20::scalar_like(zero, std::numeric_limits<double>::infinity());
}

// The pair self factor uses explicit delta, never the rounded imaginary parts.
template <uni20::Real R> R rational_error(int n, std::span<C<R> const> roots, R const& delta, R const& zero)
{
  R const one = uni20::scalar_like(zero, 1), pi = bethe::detail::pi_like(zero);
  C<R> const imaginary(zero, one);
  R error = zero;
  for (std::size_t i = 0; i < roots.size(); ++i)
  {
    Sum<C<R>> sum(C<R>(zero, zero));
    sum.add(drive(n, roots[i], zero));
    for (std::size_t j = 0; j < roots.size(); ++j)
      if (i != j)
      {
        if (i >= roots.size() - 2 && j >= roots.size() - 2)
        {
          R const self = uni20::math::log(uni20::math::abs(delta)) - uni20::math::log1p(delta);
          sum.add(C<R>(i < j ? self : -self, delta < zero ? pi : zero));
        }
        else
          sum.add(uni20::math::log(roots[i] - roots[j] - imaginary) -
                  uni20::math::log(roots[i] - roots[j] + imaginary));
      }
    error = std::max(error, exponential_error(sum.value(), zero));
  }
  return error;
}

template <uni20::Real R> struct StateData
{
    FourSpinonRequest request;
    std::vector<R> reals;
    R centre, deviation, log_abs_deviation, residual, quantum_residual, energy, zero;
    int momentum_index;
};
template <uni20::Real R> struct StateAccess;
} // namespace four_detail

// Immutable validated state: a shared snapshot also makes copy/move assignment
// atomic with respect to scalar allocation failures. Moved-from states are empty.
template <uni20::Real R> class FourSpinonState {
    std::shared_ptr<four_detail::StateData<R> const> data_;
    explicit FourSpinonState(four_detail::StateData<R> data)
        : data_(std::make_shared<four_detail::StateData<R> const>(std::move(data)))
    {}
    friend struct four_detail::StateAccess<R>;
    auto const& data() const
    {
      if (!data_) throw std::logic_error("moved-from four-spinon state");
      return *data_;
    }

  public:
    bool ready() const noexcept { return bool(data_); }
    FourSpinonRequest const& request() const { return data().request; }
    std::span<R const> real_rapidities() const { return data().reals; }
    R const& centre() const { return data().centre; }
    R const& deviation() const { return data().deviation; }
    R const& log_abs_deviation() const { return data().log_abs_deviation; }
    R const& rational_residual() const { return data().residual; }
    R const& quantum_residual() const { return data().quantum_residual; }
    R const& energy() const { return data().energy; }
    R const& working_zero() const { return data().zero; }
    int momentum_index() const { return data().momentum_index; }
    std::vector<complex_roots::C<R>> roots() const
    {
      auto const& s = data();
      std::vector<complex_roots::C<R>> result;
      result.reserve(s.reals.size() + 2);
      for (auto const& x : s.reals)
        result.emplace_back(x, s.zero);
      R const height = uni20::scalar_like(s.zero, 1) / 2 + s.deviation;
      result.emplace_back(s.centre, height);
      result.emplace_back(s.centre, -height);
      return result;
    }
};

namespace four_detail
{
template <uni20::Real R> struct StateAccess
{
    static FourSpinonState<R> make(StateData<R> data) { return FourSpinonState<R>(std::move(data)); }
};
} // namespace four_detail

template <uni20::Real R> struct FourSpinonSolve
{
    FourSpinonRequest request;
    FourSpinonStatus status = FourSpinonStatus::unconverged;
    std::optional<FourSpinonState<R>> state;
    std::string diagnostic;
    bool converged() const { return status == FourSpinonStatus::converged && state && state->ready(); }
};

template <uni20::Real R> FourSpinonSolve<R> solve_four_spinon(FourSpinonRequest request, R const& exemplar = R{0})
{
  using namespace four_detail;
  request = canonical_request(request);
  R const zero = bethe::detail::precision_zero(exemplar), one = uni20::scalar_like(zero, 1);
  R const pi = bethe::detail::pi_like(zero), eps = uni20::numeric_limits<R>::epsilon(zero);
  auto const numbers = labels(request);
  FourSpinonSolve<R> result{request};
  try
  {
    std::vector<R> x;
    for (int i : numbers)
      x.push_back(uni20::math::asinh(uni20::math::tan(2 * pi * i / (request.sites * one + 2))) / pi);
    x.push_back(request.branch * one / 2);
    x = newton<R>(
        std::move(x), zero, [&](auto const& v) { return seed_equations<R>(request, numbers, v, zero); },
        [&](auto const& v) { return seed_jacobian<R>(request.sites, v, zero); });
    if (uni20::math::abs(x.back()) <= 128 * eps)
      throw Failure(FourSpinonStatus::unsupported_geometry, "centred/singular pair is outside this family solver");
    auto const initial = pair_log<R>(request.sites, std::span<R const>(x).first(numbers.size()), x.back(), zero, zero);
    int const sign = uni20::math::cos(initial.imag()) < zero ? -1 : 1;
    x.push_back(-initial.real());
    R const initial_delta = sign * uni20::math::exp(x.back());
    if (!uni20::isfinite(initial_delta) || initial_delta == zero)
      throw Failure(FourSpinonStatus::precision_limit, "nonzero pair deviation is not representable");
    if (uni20::math::abs(initial_delta) >= 49 * one / 100)
      throw Failure(FourSpinonStatus::unsupported_geometry, "seed lies outside the narrow-pair domain");
    x = newton<R>(
        std::move(x), zero, [&](auto const& v) { return equations<R>(request.sites, numbers, sign, v, zero); },
        [&](auto const& v) { return jacobian<R>(request.sites, sign, v, zero); });
    auto const r = numbers.size();
    R const delta = sign * uni20::math::exp(x.back());
    StateData<R> data{request,
                      {x.begin(), x.begin() + r},
                      x[r],
                      delta,
                      x.back(),
                      zero,
                      norm<R>(equations<R>(request.sites, numbers, sign, x, zero), zero),
                      zero,
                      zero,
                      0};
    std::vector<C<R>> lambda;
    lambda.reserve(r + 2);
    for (auto const& value : data.reals)
      lambda.emplace_back(value, zero);
    lambda.emplace_back(data.centre, one / 2 + delta);
    lambda.emplace_back(data.centre, -one / 2 - delta);
    data.residual = rational_error<R>(request.sites, lambda, delta, zero);
    R const threshold = 512 * eps * request.sites;
    if (!uni20::isfinite(data.residual) || data.residual > threshold)
      throw Failure(FourSpinonStatus::unconverged, "original rational Bethe equations failed");
    int momentum = -request.branch;
    for (int hole : request.holes)
      momentum += hole;
    momentum %= request.sites;
    if (momentum < 0) momentum += request.sites;
    Sum<C<R>> phase_sum(C<R>(zero, zero));
    for (auto const& root : lambda)
      phase_sum.add(drive(1, root, zero));
    if (exponential_error(phase_sum.value() - C<R>(zero, 2 * pi * momentum / request.sites), zero) > threshold)
      throw Failure(FourSpinonStatus::unconverged, "resolved state left its momentum sector");
    data.momentum_index = momentum;
    try
    {
      data.energy = complex_roots::energy<R>(request.sites, lambda, zero);
    }
    catch (std::runtime_error const&)
    {
      throw Failure(FourSpinonStatus::precision_limit, "nonreal or unresolved pair energy");
    }
    result.state = StateAccess<R>::make(std::move(data));
    result.status = FourSpinonStatus::converged;
  }
  catch (Failure const& e)
  {
    result.status = e.status;
    result.diagnostic = e.what();
  }
  return result;
}
} // namespace bethe::experimental::xxx
