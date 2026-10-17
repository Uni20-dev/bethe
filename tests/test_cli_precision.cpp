// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "cli-common.hpp"
#include <gtest/gtest.h>

namespace cli = bethe::cli;

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
  auto bits = []<typename Real>() { return uni20::numeric_limits<Real>::digits; };
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
  EXPECT_THROW(cli::validate_precision(cli::parse_precision("digits-100"), false), std::invalid_argument);
  if (!cli::has_arbitrary_structure_factor)
    EXPECT_THROW(cli::validate_precision(cli::parse_precision("bits-128"), true), std::invalid_argument);
}

#if UNI20_ENABLE_MPFR
TEST(CliPrecision, DecimalAndBinaryWorkingPrecision)
{
  EXPECT_EQ(cli::parse_precision("digits-100").working_precision().bit_count(), 333);
  EXPECT_EQ(cli::parse_precision("bits-512").working_precision().bit_count(), 512);
  EXPECT_EQ(cli::parse_precision("bits-128").working_precision().bit_count(), 128);
  EXPECT_THROW(cli::parse_precision("fp128").working_precision(), std::logic_error);
  EXPECT_THROW(cli::parse_precision("digits-9223372036854775807").working_precision(), std::invalid_argument);
}
#endif
