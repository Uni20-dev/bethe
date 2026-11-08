// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/detail/spinon_barnes.hpp>

template <class R> class SpinonBarnes : public ::testing::Test {};
TYPED_TEST_SUITE(SpinonBarnes, test_support::RealTypes, test_support::PrecisionNames);

template <class R> void check_barnes(uni20::Precision precision, uni20::log_gamma_options options = {})
{
  bethe::detail::SpinonBarnesContext<R> context(precision, options);
  auto n = [&](int x) { return uni20::scalar_at_precision<R>(x, precision); };
  // G(4) = Gamma(3) Gamma(2) Gamma(1) = 2, independent of the implementation.
  auto result = context.log_abs_barnes_g_right(uni20::complex<R>(n(4), n(0)));
  auto expected = uni20::math::log(n(2));
  EXPECT_LE(uni20::math::abs(result.value - expected),
            result.truncation_bound + result.roundoff_indicator);
  EXPECT_LT(result.truncation_bound + result.roundoff_indicator,
            n(10000000) * uni20::numeric_limits<R>::epsilon(precision));
}

TYPED_TEST(SpinonBarnes, DefaultAndExplicitCoefficientBudgets)
{
  using R = TypeParam;
  // New Uni20 leaves max_terms unset; the pinned version stores an integer.
  check_barnes<R>(uni20::native_precision<R>());
  uni20::log_gamma_options options;
  options.max_terms = 16;
  check_barnes<R>(uni20::native_precision<R>(), options);
  options.max_terms = 15;
  EXPECT_THROW((bethe::detail::SpinonBarnesContext<R>(options)), std::invalid_argument);
}

#if UNI20_HAS_MPREAL
TEST(SpinonBarnesMpfr, DefaultAndExplicitCoefficientBudgets)
{
  using R = uni20::mpreal;
  for (int bits : {128, 256})
  {
    auto precision = uni20::Precision::bits(bits);
    check_barnes<R>(precision);
    uni20::log_gamma_options options;
    options.max_terms = 256;
    check_barnes<R>(precision, options);
    options.max_terms = 16;
    EXPECT_THROW((bethe::detail::SpinonBarnesContext<R>(precision, options)), std::invalid_argument);
  }
}
#endif
