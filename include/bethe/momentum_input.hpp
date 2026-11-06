// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <algorithm>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <uni20/core/math.hpp>
#include <uni20/core/runtime_precision.hpp>
#include <uni20/core/scalar_io.hpp>
#include <utility>

namespace bethe
{
enum class MomentumInputKind
{
  radians,
  distance_from_pi,
  beta
};
enum class MomentumEndpoint
{
  interior,
  zero,
  pi,
  two_pi
};

/// Folded external momentum q_f=min(q,2*pi-q), b=log(tan(q_f/2)).
/// This is a coordinate for total momentum, not an individual particle rapidity.
/// Endpoints are explicit; a finite beta always denotes an interior point.
template <uni20::Real R> class MomentumCoordinate {
  public:
    static MomentumCoordinate from_beta(R beta, bool reflected = false)
    {
      if (!uni20::isfinite(beta)) throw std::invalid_argument("momentum beta must be finite");
      return {std::move(beta), MomentumEndpoint::interior, reflected};
    }
    static MomentumCoordinate at_endpoint(MomentumEndpoint endpoint, R const& exemplar)
    {
      if (endpoint == MomentumEndpoint::interior) throw std::invalid_argument("expected exact momentum endpoint");
      return {uni20::scalar_like(exemplar, 0), endpoint, endpoint == MomentumEndpoint::two_pi};
    }
    MomentumEndpoint endpoint() const { return endpoint_; }
    bool reflected() const { return reflected_; }
    std::optional<R> beta() const
    {
      return endpoint_ == MomentumEndpoint::interior ? std::optional(beta_) : std::nullopt;
    }

    /// Natural logarithm of distance from the nearest 0/pi/2pi endpoint.
    /// Exact endpoints have no finite logarithmic distance.
    std::optional<R> log_distance() const
    {
      if (!beta()) return std::nullopt;
      R const one = uni20::scalar_like(beta_, 1), t = uni20::math::exp(-uni20::math::abs(beta_));
      R const correction = t == R{0} ? uni20::scalar_like(beta_, 0) : uni20::math::log(uni20::math::atan(t) / t);
      return uni20::math::log(R{2} * one) - uni20::math::abs(beta_) + correction;
    }

    /// Explicit adapter for scalar-q backends. Never silently returns an endpoint
    /// for an interior coordinate. It does not certify downstream conditioning.
    R radians() const
    {
      R const one = uni20::scalar_like(beta_, 1), pi = R{4} * uni20::math::atan(one);
      if (endpoint_ == MomentumEndpoint::zero) return uni20::scalar_like(beta_, 0);
      if (endpoint_ == MomentumEndpoint::pi) return pi;
      if (endpoint_ == MomentumEndpoint::two_pi) return R{2} * pi;
      R const log_t = -uni20::math::abs(beta_);
      // For small t, atan(t)=t to working precision. Scale in the exponent
      // before rounding: t can underflow while the final distance 2*t is
      // still representable (including the smallest native subnormal).
      R const d = log_t < uni20::math::log(uni20::numeric_limits<R>::epsilon(one)) / R{2}
                      ? uni20::math::exp(log_t + uni20::math::log(R{2} * one))
                      : R{2} * uni20::math::atan(uni20::math::exp(log_t));
      R const folded = beta_ <= R{0} ? d : pi - d;
      R const q = reflected_ ? R{2} * pi - folded : folded;
      if (!(q > R{0} && q < R{2} * pi) || q == pi)
        throw std::range_error("interior momentum rounds to an endpoint in scalar-q geometry; increase geometry "
                               "precision or use a beta-coordinate backend");
      return q;
    }

  private:
    MomentumCoordinate(R beta, MomentumEndpoint endpoint, bool reflected)
        : beta_(std::move(beta)), endpoint_(endpoint), reflected_(reflected)
    {}
    R beta_;
    MomentumEndpoint endpoint_;
    bool reflected_;
};

template <uni20::Real R> struct ResolvedMomentum
{
    MomentumCoordinate<R> coordinate;
    std::string beta_absolute_error; ///< Estimated conversion error, not a rigorous bound.
    bool beta_accuracy_warning;
    std::optional<R> scalar_radians;
    std::string scalar_beta_absolute_error;
    bool scalar_accuracy_warning;
    std::int64_t input_bits;
    bool mpfr_input;

    R require_radians() const
    {
      if (!scalar_radians)
        throw std::range_error("interior momentum cannot be represented by this scalar-q backend; increase geometry "
                               "precision or use a beta-coordinate backend");
      return *scalar_radians;
    }
};

namespace detail
{
// Validate decimal syntax before parsing, and remember whether a zero result
// could represent exponent underflow rather than an exact endpoint.
inline bool nonzero_momentum_decimal(std::string_view text)
{
  std::size_t i = 0, digits = 0;
  bool nonzero = false;
  if (i < text.size() && (text[i] == '+' || text[i] == '-')) ++i;
  auto consume_digits = [&] {
    while (i < text.size() && text[i] >= '0' && text[i] <= '9')
    {
      nonzero = nonzero || text[i] != '0';
      ++digits;
      ++i;
    }
  };
  consume_digits();
  if (i < text.size() && text[i] == '.')
  {
    ++i;
    consume_digits();
  }
  if (!digits) throw std::invalid_argument("momentum input requires a decimal number");
  if (i < text.size() && (text[i] == 'e' || text[i] == 'E'))
  {
    ++i;
    if (i < text.size() && (text[i] == '+' || text[i] == '-')) ++i;
    auto const start = i;
    while (i < text.size() && text[i] >= '0' && text[i] <= '9')
      ++i;
    if (i == start) throw std::invalid_argument("momentum exponent requires digits");
  }
  if (i != text.size()) throw std::invalid_argument("invalid decimal momentum input");
  return nonzero;
}
template <uni20::Real R> R momentum_small_beta(R const& distance)
{
  R const one = uni20::scalar_like(distance, 1);
  // Avoid underflow of distance/2 and preserve its logarithm. The neglected
  // correction is below working precision on this branch.
  if (distance < uni20::math::sqrt(uni20::numeric_limits<R>::epsilon(one)))
    return uni20::math::log(distance) - uni20::math::log(R{2} * one);
  return uni20::math::log(uni20::math::tan(distance / R{2}));
}

// d(beta)/dq = 1/sin(q). Estimate amplification of angle/pi roundoff
// when forming an endpoint distance by subtraction. Tiny q or an explicitly
// supplied tiny distance has condition near one. Exact endpoints have no beta.
template <uni20::Real R> R momentum_angle_condition(R const& q, R const& one)
{
  R const pi = R{4} * uni20::math::atan(one);
  if (q > R{0} && q < R{2} * pi)
  {
    R const folded = q <= pi ? q : R{2} * pi - q;
    R const distance = std::min(folded, pi - folded);
    if (distance > R{0}) return std::max(one, q / distance);
  }
  return one;
}

template <uni20::Real R> MomentumCoordinate<R> momentum_from_radians(R const& q, R const& one)
{
  R const pi = R{4} * uni20::math::atan(one);
  if (!uni20::isfinite(q) || q < R{0} || q > R{2} * pi) throw std::invalid_argument("momentum must lie in [0,2*pi]");
  if (q == R{0}) return MomentumCoordinate<R>::at_endpoint(MomentumEndpoint::zero, one);
  if (q == pi) return MomentumCoordinate<R>::at_endpoint(MomentumEndpoint::pi, one);
  if (q == R{2} * pi) return MomentumCoordinate<R>::at_endpoint(MomentumEndpoint::two_pi, one);
  bool const reflected = q > pi;
  R const folded = reflected ? R{2} * pi - q : q;
  R const beta = folded <= pi / R{2} ? momentum_small_beta(folded) : -momentum_small_beta(pi - folded);
  return MomentumCoordinate<R>::from_beta(beta, reflected);
}

template <uni20::Real R> R momentum_parse(std::string_view text, R const& one)
{
  if constexpr (uni20::has_runtime_precision_v<R>)
    return uni20::parse_real<R>(text, one.precision());
  else
    return uni20::parse_real<R>(text);
}

template <uni20::Real R, uni20::Real From> R momentum_cast(R const& one, From const& value)
{
  if constexpr (uni20::has_runtime_precision_v<R>)
    return uni20::scalar_like(one, value);
  else
    // Native narrowing is intentional. GCC's C++23 extended-float rules
    // require an explicit cast for binary128 -> double/long double.
    return static_cast<R>(value);
}
} // namespace detail

/// Input text is retained until conversion. MPFR builds use a separate guarded
/// input precision, even for native working arithmetic; ambient defaults are
/// neither consulted nor changed. Without MPFR use the widest enabled native
/// format and report an estimated conversion error at that precision.
struct MomentumInput
{
    MomentumInputKind kind;
    std::string text;

    template <uni20::Real R>
    ResolvedMomentum<R> resolve(R const& exemplar, std::string_view beta_tolerance = "1e-12") const
    {
      R const one = uni20::scalar_like(exemplar, 1);
      if constexpr (uni20::has_runtime_precision_v<R>)
        if (one.precision().is_exact())
          throw std::invalid_argument("momentum conversion needs finite working precision");
#if UNI20_HAS_MPREAL
      if (text.size() > static_cast<std::size_t>((MPFR_PREC_MAX - 64) / 4) ||
          uni20::numeric_limits<R>::digits(one) > MPFR_PREC_MAX - 64)
        throw std::length_error("momentum input precision exceeds MPFR range");
      auto const bits = std::max<std::int64_t>(
          128, std::max<std::int64_t>(uni20::numeric_limits<R>::digits(one) + 64, 4 * text.size() + 64));
      uni20::mpreal const reference_one(1, uni20::Precision::bits(bits));
      return convert(one, reference_one, beta_tolerance, true);
#elif UNI20_HAS_FLOAT128
      return convert(one, uni20::float128{1}, beta_tolerance, false);
#else
      return convert(one, 1.L, beta_tolerance, false);
#endif
    }

  private:
    template <uni20::Real R, uni20::Real Ref>
    ResolvedMomentum<R> convert(R const& one, Ref const& reference_one, std::string_view tolerance_text,
                                bool mpfr) const
    {
      Ref const tolerance = detail::momentum_parse(tolerance_text, reference_one);
      if (!uni20::isfinite(tolerance) || tolerance <= Ref{0})
        throw std::invalid_argument("momentum beta tolerance must be positive and finite");
      Ref const pi = Ref{4} * uni20::math::atan(reference_one);
      Ref input_condition = reference_one;
      std::optional<Ref> input_angle;
      auto const reference = [&] {
        if (kind == MomentumInputKind::radians && (text == "pi" || text == "2*pi"))
          return MomentumCoordinate<Ref>::at_endpoint(text == "pi" ? MomentumEndpoint::pi : MomentumEndpoint::two_pi,
                                                      reference_one);
        if (kind == MomentumInputKind::distance_from_pi && text == "pi")
          return MomentumCoordinate<Ref>::at_endpoint(MomentumEndpoint::zero, reference_one);
        bool const nonzero = detail::nonzero_momentum_decimal(text);
        Ref const value = detail::momentum_parse(text, reference_one);
        if (nonzero && value == Ref{0})
          throw std::range_error("momentum input underflows the input format; specify beta instead");
        if (!uni20::isfinite(value)) throw std::invalid_argument("momentum input must be finite");
        if (kind == MomentumInputKind::beta) return MomentumCoordinate<Ref>::from_beta(value);
        input_condition = detail::momentum_angle_condition(value, reference_one);
        if (kind == MomentumInputKind::distance_from_pi)
        {
          if (value < Ref{0} || value > pi) throw std::invalid_argument("distance from pi must lie in [0,pi]");
          input_angle = pi - value;
          auto const distance = detail::momentum_from_radians(value, reference_one);
          if (distance.endpoint() == MomentumEndpoint::zero)
            return MomentumCoordinate<Ref>::at_endpoint(MomentumEndpoint::pi, reference_one);
          if (distance.endpoint() == MomentumEndpoint::pi)
            throw std::range_error(
                "decimal distance is indistinguishable from pi at input precision; specify pi or beta");
          return MomentumCoordinate<Ref>::from_beta(-*distance.beta());
        }
        input_angle = value;
        auto coordinate = detail::momentum_from_radians(value, reference_one);
        // Exact nonzero endpoints require symbolic input: a rounded decimal
        // equal to this precision's pi cannot identify its side or offset.
        if (coordinate.endpoint() == MomentumEndpoint::pi)
          throw std::range_error("decimal momentum is indistinguishable from pi at input precision; specify pi, an "
                                 "endpoint distance, or beta");
        if (coordinate.endpoint() == MomentumEndpoint::two_pi)
          throw std::range_error("decimal momentum is indistinguishable from 2*pi at input precision; specify 2*pi "
                                 "for the exact endpoint, or use wider input precision for an interior decimal");
        return coordinate;
      }();
      if (!reference.beta())
      {
        auto coordinate = MomentumCoordinate<R>::at_endpoint(reference.endpoint(), one);
        return {
            coordinate, "0", false, coordinate.radians(), "0", false, uni20::numeric_limits<Ref>::digits(reference_one),
            mpfr};
      }
      Ref const beta = *reference.beta();
      R const rounded = detail::momentum_cast(one, beta);
      if (!uni20::isfinite(rounded)) throw std::range_error("momentum beta exceeds the working format's range");
      auto coordinate = MomentumCoordinate<R>::from_beta(rounded, reference.reflected());
      Ref const reference_error = Ref{8} * uni20::numeric_limits<Ref>::epsilon(reference_one) *
                                  std::max(input_condition, uni20::math::abs(beta));
      Ref const error = uni20::math::abs(beta - uni20::scalar_like(reference_one, rounded)) + reference_error;
      std::optional<R> angle;
      std::string angle_error;
      bool angle_warning = false;
      try
      {
        // Preserve supplied angles, and convert beta/distance inputs before
        // rounding to R. Reconstructing from the rounded working beta can
        // shift ordinary q by several ULPs, or hide an endpoint collapse.
        angle = detail::momentum_cast(one, input_angle ? *input_angle : reference.radians());
        if (kind == MomentumInputKind::radians)
        {
          // Direct parsing avoids double rounding through a native reference.
          // Uni20's native parser rejects ERANGE even for representable
          // subnormals; the guarded conversion above supplies those values.
          try
          {
            angle = detail::momentum_parse(text, one);
          }
          catch (std::out_of_range const&)
          {}
        }
        R const working_pi = R{4} * uni20::math::atan(one);
        if (!(*angle > R{0} && *angle < R{2} * working_pi) || *angle == working_pi) angle.reset();
        if (angle)
        {
          Ref const reference_angle = uni20::scalar_like(reference_one, *angle);
          auto const roundtrip = detail::momentum_from_radians(reference_angle, reference_one);
          if (!roundtrip.beta() || roundtrip.reflected() != reference.reflected())
            angle.reset();
          else
          {
            // The roundtrip oracle has its own pi-subtraction uncertainty,
            // even when the original beta input needed no subtraction at all.
            Ref const roundtrip_error = Ref{8} * uni20::numeric_limits<Ref>::epsilon(reference_one) *
                                        std::max(detail::momentum_angle_condition(reference_angle, reference_one),
                                                 uni20::math::abs(*roundtrip.beta()));
            Ref const difference = uni20::math::abs(beta - *roundtrip.beta()) + reference_error + roundtrip_error;
            angle_error = uni20::format_scalar(difference);
            angle_warning = difference > tolerance;
          }
        }
      }
      catch (std::range_error const&)
      {
        angle.reset();
      }
      return {std::move(coordinate),
              uni20::format_scalar(error),
              error > tolerance,
              angle,
              angle_error,
              angle_warning,
              uni20::numeric_limits<Ref>::digits(reference_one),
              mpfr};
    }
};
} // namespace bethe
