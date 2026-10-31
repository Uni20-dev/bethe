// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
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

