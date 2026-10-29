// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <uni20/core/numeric_limits.hpp>
#include <uni20/core/scalar_io.hpp>
#if UNI20_HAS_MPREAL
#include <uni20/core/precision.hpp>
#endif
#include <charconv>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace bethe::cli
{
// Test the arithmetic format, not sizeof(T): extended80 commonly has padding.
template <typename T>
inline constexpr bool is_extended80 =
    uni20::numeric_limits<T>::radix == 2 && uni20::numeric_limits<T>::digits() == 64 &&
    uni20::numeric_limits<T>::min_exponent == -16381 && uni20::numeric_limits<T>::max_exponent == 16384;
inline constexpr bool has_fp80 = is_extended80<long double>;
inline constexpr bool has_arbitrary_structure_factor = UNI20_HAS_MPREAL;

// The frontend declares which operations it needs; arithmetic metadata itself
// stays separate from model support. Uni20 provides CPU QR and real/complex
// LU for enabled scalars, preferring an available LAPACK provider.
enum class ArbitraryPrecision
{
  none,
  real,
  real_lu,
  complex_linalg
};
constexpr bool has_arbitrary_precision(ArbitraryPrecision support)
{
  return support != ArbitraryPrecision::none && bool(UNI20_HAS_MPREAL);
}

enum class PrecisionKind
{
  fp64,
  long_double,
  fp80,
  fp128,
  bits,
  digits
};
struct PrecisionRequest
{
    PrecisionKind kind;
    std::int64_t count = 0;
    bool arbitrary() const { return kind == PrecisionKind::bits || kind == PrecisionKind::digits; }
#if UNI20_HAS_MPREAL
    uni20::Precision working_precision() const
    {
      if (kind == PrecisionKind::digits) return uni20::Precision::decimal_digits(count);
      if (kind != PrecisionKind::bits) throw std::logic_error("native precision has no runtime precision request");
      if (count > std::numeric_limits<mpfr_prec_t>::max())
        throw std::invalid_argument("requested bit precision is out of range");
      return uni20::Precision::bits(static_cast<mpfr_prec_t>(count));
    }
#endif
};

inline PrecisionRequest parse_precision(std::string_view text)
{
  if (text == "fp64") return {PrecisionKind::fp64};
  if (text == "long-double") return {PrecisionKind::long_double};
  if (text == "fp80") return {PrecisionKind::fp80};
  if (text == "fp128") return {PrecisionKind::fp128};
  for (auto const& [prefix, kind] :
       {std::pair{"bits-", PrecisionKind::bits}, std::pair{"digits-", PrecisionKind::digits}})
    if (text.starts_with(prefix))
    {
      auto const suffix = text.substr(std::string_view(prefix).size());
      std::int64_t count = 0;
      auto const [end, error] = std::from_chars(suffix.data(), suffix.data() + suffix.size(), count);
      if (suffix.empty() || suffix.front() < '0' || suffix.front() > '9' || error != std::errc{} ||
          end != suffix.data() + suffix.size() || count <= 0)
        throw std::invalid_argument("precision requires a positive integer in bits-N or digits-N: " +
                                    std::string(text));
      return {kind, count};
    }
  throw std::invalid_argument("unknown precision: " + std::string(text) +
                              "; choose fp64, long-double, fp80, fp128, bits-N or digits-N");
}

inline void validate_precision(PrecisionRequest request, ArbitraryPrecision support = ArbitraryPrecision::none)
{
  if (request.kind == PrecisionKind::fp80 && !has_fp80)
    throw std::invalid_argument("fp80 is unavailable: native long double is not the extended80 format; see bethe-info");
  if (request.kind == PrecisionKind::fp128 && !UNI20_HAS_FLOAT128)
    throw std::invalid_argument("fp128 is unavailable; configure with -DUNI20_ENABLE_FLOAT128=ON; see bethe-info");
  if (request.arbitrary())
  {
    if (support == ArbitraryPrecision::none)
      throw std::invalid_argument("this frontend does not support arbitrary precision");
    if (!has_arbitrary_precision(support))
      throw std::invalid_argument(
          support == ArbitraryPrecision::real ? "arbitrary precision requires -DUNI20_ENABLE_MPREAL=ON; see bethe-info"
          : support == ArbitraryPrecision::real_lu
              ? "arbitrary-precision linear solves require -DUNI20_ENABLE_MPREAL=ON; see bethe-info"
              : "arbitrary-precision form factors require -DUNI20_ENABLE_MPREAL=ON; see bethe-info");
#if UNI20_HAS_MPREAL
    (void)request.working_precision();
#endif
  }
}

struct ArithmeticInfo
{
    std::string format;
    std::int64_t bits;
    int decimal_digits, min_exponent, max_exponent;
    std::size_t storage_bytes;
};
template <typename T> ArithmeticInfo native_arithmetic()
{
  using Limits = uni20::numeric_limits<T>;
  std::string format = "native floating point";
  if constexpr (is_extended80<T>)
    format = "extended80";
  else if constexpr (Limits::radix == 2 && Limits::digits() == 53 && Limits::max_exponent == 1024)
    format = "binary64";
  else if constexpr (Limits::radix == 2 && Limits::digits() == 113 && Limits::max_exponent == 16384)
    format = "binary128";
  return {format, Limits::digits(), Limits::digits10(), Limits::min_exponent, Limits::max_exponent, sizeof(T)};
}
inline ArithmeticInfo arithmetic_info(PrecisionRequest request)
{
  validate_precision(request, ArbitraryPrecision::real);
  switch (request.kind)
  {
    case PrecisionKind::fp64:
      return native_arithmetic<double>();
    case PrecisionKind::long_double:
    case PrecisionKind::fp80:
      return native_arithmetic<long double>();
    case PrecisionKind::fp128:
#if UNI20_HAS_FLOAT128
      return native_arithmetic<uni20::float128>();
#else
      break;
#endif
    case PrecisionKind::bits:
    case PrecisionKind::digits:
#if UNI20_HAS_MPREAL
      return {"arbitrary binary", request.working_precision().bit_count(), 0, 0, 0, 0};
#else
      break;
#endif
  }
  throw std::logic_error("unavailable arithmetic");
}
} // namespace bethe::cli
