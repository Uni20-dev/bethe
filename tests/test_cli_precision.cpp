// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "cli-common.hpp"
#include "solver-precision.hpp"
#include "spinon-thermo.hpp"
#include <gtest/gtest.h>

namespace cli = bethe::cli;

TEST(CliPrecision, SpinonNativeComparisonUsesDistinctIntermediatePrecision)
{
  namespace spinon = cli::spinon;
  EXPECT_TRUE((std::is_same_v<spinon::DensityLow<double>, float>));
  EXPECT_TRUE((std::is_same_v<spinon::DensityLow<long double>, double>));
  auto const binary64 = spinon::recovery(spinon::Controls{}, 1.0);
  EXPECT_EQ(binary64.stages.front().bits, 24);
  EXPECT_EQ(binary64.stages.front().guard_bits, 53);
#if UNI20_HAS_FLOAT128
  auto const binary128 = spinon::recovery(spinon::Controls{}, uni20::float128{1});
  constexpr int native_bits = std::numeric_limits<long double>::digits;
  EXPECT_EQ(binary128.stages.front().bits, native_bits > 53 && native_bits < 113 ? native_bits : 53);
  EXPECT_EQ(binary128.stages.front().guard_bits, 113);
  EXPECT_LT(binary128.stages.front().bits, binary128.stages.front().guard_bits);
#endif
#if UNI20_HAS_MPREAL
  EXPECT_TRUE((std::is_same_v<spinon::DensityLow<uni20::mpreal>, uni20::mpreal>));
#endif
}

TEST(CliPrecision, ParseFormatsAndExplicitPrecision)
{
  EXPECT_EQ(cli::parse_precision("fp64").kind, cli::PrecisionKind::fp64);
  EXPECT_EQ(cli::parse_precision("fp80").kind, cli::PrecisionKind::fp80);
  EXPECT_EQ(cli::parse_precision("fp128").kind, cli::PrecisionKind::fp128);
  EXPECT_EQ(cli::parse_precision("long-double").kind, cli::PrecisionKind::long_double);
  EXPECT_EQ(cli::parse_precision("bits-512").count, 512);
  EXPECT_EQ(cli::parse_precision("digits-100").count, 100);
  EXPECT_TRUE(cli::parse_precision("bits-128").arbitrary());
  EXPECT_FALSE(cli::parse_precision("fp128").arbitrary());
  for (auto text : {"", "mpfr", "arbitrary", "FP80", "bits", "bits-", "bits-0", "bits--1", "bits-+512", "digits-1.5",
                    "digits-100x", "digits- 100", "bits-512 ", "digits-9223372036854775808"})
  {
    SCOPED_TRACE(text);
    EXPECT_THROW(cli::parse_precision(text), std::invalid_argument);
  }
}

TEST(CliPrecision, NativeFormatAndNoFallback)
{
  EXPECT_FALSE(cli::is_extended80<double>);
  auto bits = []<typename Real>() { return uni20::numeric_limits<Real>::digits(); };
  EXPECT_EQ(cli::dispatch_precision("fp64", bits), 53);
  EXPECT_EQ(cli::dispatch_precision("long-double", bits), std::numeric_limits<long double>::digits);
  if (cli::has_fp80)
  {
    EXPECT_EQ(cli::dispatch_precision("fp80", bits), 64);
    EXPECT_EQ(cli::arithmetic_info(cli::parse_precision("fp80")).format, "extended80");
  }
  else
    EXPECT_THROW(cli::dispatch_precision("fp80", bits), std::invalid_argument);
#if UNI20_HAS_FLOAT128
  EXPECT_FALSE(cli::is_extended80<uni20::float128>);
  EXPECT_EQ(cli::dispatch_precision("fp128", bits), 113);
#else
  EXPECT_THROW(cli::dispatch_precision("fp128", bits), std::invalid_argument);
#endif
  EXPECT_THROW(cli::dispatch_precision("bits-128", bits), std::invalid_argument);
  EXPECT_THROW(cli::validate_precision(cli::parse_precision("digits-100")), std::invalid_argument);
  if (!cli::has_arbitrary_structure_factor)
    EXPECT_THROW(cli::validate_precision(cli::parse_precision("bits-128"), cli::ArbitraryPrecision::complex_linalg),
                 std::invalid_argument);
  if (!UNI20_HAS_MPREAL)
    EXPECT_THROW(cli::validate_precision(cli::parse_precision("bits-128"), cli::ArbitraryPrecision::real),
                 std::invalid_argument);
  if (cli::has_arbitrary_precision(cli::ArbitraryPrecision::real_lu))
    EXPECT_NO_THROW(cli::validate_precision(cli::parse_precision("bits-128"), cli::ArbitraryPrecision::real_lu));
  else
    EXPECT_THROW(cli::validate_precision(cli::parse_precision("bits-128"), cli::ArbitraryPrecision::real_lu),
                 std::invalid_argument);
}

TEST(CliPrecision, ScalarDispatchAndParsing)
{
  auto native = []<typename Real>(uni20::Precision p) {
    static_assert(!uni20::has_runtime_precision_v<Real>);
    EXPECT_EQ(p, uni20::native_precision<Real>());
    EXPECT_EQ(cli::parse_scalar<Real>("1.25", p), Real{1.25});
    return int(p.bit_count());
  };
  EXPECT_EQ(cli::dispatch_scalar_precision<cli::ArbitraryPrecision::none>("fp64", native), 53);
  EXPECT_THROW(cli::dispatch_scalar_precision<cli::ArbitraryPrecision::none>("bits-128", native),
               std::invalid_argument);
  EXPECT_THROW(cli::parse_scalar<double>("1", uni20::Precision::bits(128)), std::invalid_argument);
  EXPECT_THROW(cli::parse_scalar<double>("1", uni20::Precision::exact()), std::invalid_argument);
  auto real = []<typename Real>(uni20::Precision p) {
    EXPECT_EQ(p.bit_count(), 333);
    auto value = cli::parse_scalar<Real>("1.00000000000000000000000000000000000000000000000001", p);
    EXPECT_GT(value, Real{1});
    if constexpr (uni20::has_runtime_precision_v<Real>) EXPECT_EQ(value.precision(), p);
    return 0;
  };
#if UNI20_HAS_MPREAL
  EXPECT_EQ(cli::dispatch_scalar_precision<cli::ArbitraryPrecision::real>("digits-100", real), 0);
  EXPECT_THROW(cli::parse_scalar<uni20::mpreal>("1", uni20::Precision::exact()), std::invalid_argument);
#else
  EXPECT_THROW(cli::dispatch_scalar_precision<cli::ArbitraryPrecision::real>("digits-100", real),
               std::invalid_argument);
#endif
}

#if UNI20_HAS_MPREAL
TEST(CliPrecision, DecimalAndBinaryWorkingPrecision)
{
  EXPECT_NO_THROW(cli::validate_precision(cli::parse_precision("bits-128"), cli::ArbitraryPrecision::real));
  EXPECT_EQ(cli::arithmetic_info(cli::parse_precision("digits-100")).bits, 333);
  EXPECT_EQ(cli::parse_precision("digits-100").working_precision().bit_count(), 333);
  EXPECT_EQ(cli::parse_precision("bits-512").working_precision().bit_count(), 512);
  EXPECT_EQ(cli::parse_precision("bits-128").working_precision().bit_count(), 128);
  EXPECT_THROW(cli::parse_precision("fp128").working_precision(), std::logic_error);
  EXPECT_THROW(cli::parse_precision("digits-9223372036854775807").working_precision(), std::invalid_argument);
}

TEST(CliPrecision, SolverControlsUseSelectedPrecision)
{
  auto check = []<typename Real>(bethe::SolverOptions<Real> const& options) {
    EXPECT_EQ(options.max_iterations, 17);
    if constexpr (uni20::has_runtime_precision_v<Real>)
    {
      EXPECT_EQ(options.precision(), uni20::Precision::bits(333));
      EXPECT_EQ(options.residual_tolerance.precision(), options.precision());
      EXPECT_EQ(options.residual_tolerance,
                Real("1.00000000000000000000000000000000000000001e-90", options.precision()));
    }
    return 0;
  };
  EXPECT_EQ(cli::dispatch_solver_precision<cli::ArbitraryPrecision::real>(
                "digits-100", "1.00000000000000000000000000000000000000001e-90", 17, check),
            0);
}
#endif
