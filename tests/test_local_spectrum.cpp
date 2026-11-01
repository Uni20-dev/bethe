// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <array>
#include <bethe/dynamical_structure_factor.hpp>

namespace
{
template <uni20::Real R> void analytic_spectra(R const& zero)
{
  using C = uni20::complex<R>;
  auto equal = [](C const& a, C const& b) {
    EXPECT_EQ(a.real(), b.real());
    EXPECT_EQ(a.imag(), b.imag());
  };
  auto r = [&](int n) { return uni20::scalar_like(zero, n); };
  R const eps = uni20::numeric_limits<R>::epsilon(zero), tol = 64 * eps;
  R const time = bethe::detail::pi_like(zero) / 2;
  std::array<bethe::LocalSpectralLine<R>, 2> const local{{{0, r(0), r(2)}, {1, r(1), r(3)}}};
  auto const moments = bethe::local_spectral_moments<R>(local, zero);
  EXPECT_EQ(moments.weight, r(5));
  EXPECT_EQ(moments.first_moment, r(3));
  auto c = bethe::local_time_correlation<R>(local, time);
  EXPECT_REAL_NEAR(c.real(), r(2), tol);
  EXPECT_REAL_NEAR(c.imag(), r(-3), tol);
  auto backward = bethe::local_time_correlation<R>(local, -time);
  EXPECT_EQ(backward.real(), c.real());
  EXPECT_EQ(backward.imag(), -c.imag());
  auto initial = bethe::local_time_correlation<R>(local, zero);
  EXPECT_EQ(initial.real(), moments.weight);
  EXPECT_EQ(initial.imag(), zero);

  // A finite difference of the time signal independently checks its relation
  // to the moment operator. No common implementation of that operator is used.
  R const h = uni20::math::sqrt(uni20::math::sqrt(eps));
  auto const derivative =
      (bethe::local_time_correlation<R>(local, h) - bethe::local_time_correlation<R>(local, -h)) / (2 * h);
  EXPECT_REAL_NEAR(derivative.real(), zero, 8 * h * h);
  EXPECT_REAL_NEAR(derivative.imag(), -moments.first_moment, 8 * h * h);

  std::array<bethe::CrossSpectralLine<R>, 2> cross{{{0, r(1), C(r(1), r(2))}, {1, r(2), C(r(-3), r(1))}}};
  auto const cm = bethe::cross_spectral_moments<R>(cross, zero);
  equal(cm.value, C(r(-2), r(3)));
  equal(cm.first_moment, C(r(-5), r(4)));
  equal(bethe::cross_time_correlation<R>(cross, zero), cm.value);
  auto const ct = bethe::cross_time_correlation<R>(cross, time);
  EXPECT_REAL_NEAR(ct.real(), r(5), tol);
  EXPECT_REAL_NEAR(ct.imag(), r(-2), tol);
  auto const cd = (bethe::cross_time_correlation<R>(cross, h) - bethe::cross_time_correlation<R>(cross, -h)) / (2 * h);
  EXPECT_REAL_NEAR(cd.real(), cm.first_moment.imag(), 16 * h * h);
  EXPECT_REAL_NEAR(cd.imag(), -cm.first_moment.real(), 16 * h * h);
  auto const cb = bethe::cross_time_correlation<R>(cross, -time);
  EXPECT_GT(uni20::math::abs(cb - uni20::math::conj(ct)), r(1));
  for (auto& line : cross)
    line.coefficient = uni20::math::conj(line.coefficient);
  auto const reverse = bethe::cross_time_correlation<R>(cross, -time);
  EXPECT_REAL_NEAR(reverse.real(), ct.real(), tol);
  EXPECT_REAL_NEAR(reverse.imag(), -ct.imag(), tol);

  std::array<bethe::LehmannLine<R>, 3> const periodic{
      {{0, 0, zero, zero, r(2)}, {1, 0, zero, r(1), r(3)}, {2, 2, time, r(2), r(7)}}};
  auto const grouped = bethe::spectral_time_correlation<R>(4, periodic, time);
  ASSERT_EQ(grouped.size(), 4u);
  equal(grouped[0], c);
  equal(grouped[1], C(zero, zero));
  EXPECT_REAL_NEAR(grouped[2].real(), r(-7), tol);
  EXPECT_REAL_NEAR(grouped[2].imag(), zero, tol);
  equal(grouped[3], C(zero, zero));
  EXPECT_EQ(bethe::local_spectral_moments<R>({}, zero).weight, zero);
  equal(bethe::cross_spectral_moments<R>({}, zero).first_moment, C(zero, zero));
  equal(bethe::local_time_correlation<R>({}, time), C(zero, zero));
  equal(bethe::cross_time_correlation<R>({}, time), C(zero, zero));
  EXPECT_TRUE(bethe::spectral_time_correlation<R>(0, {}, time).empty());
  EXPECT_THROW((void)bethe::spectral_time_correlation<R>(2, periodic, time), std::invalid_argument);

  // Compensated signed cancellation. At this magnitude each individual +1
  // rounds away in a naive sum, but their accumulated contribution is eight.
  R const large = 8 / eps;
  std::vector<bethe::CrossSpectralLine<R>> cancel{{0, zero, C(large, -large)}};
  for (std::size_t i = 1; i <= 8; ++i)
    cancel.push_back({i, zero, C(r(1), r(-1))});
  cancel.push_back({9, zero, C(-large, large)});
  equal(bethe::cross_time_correlation<R>(cancel, zero), C(r(8), r(-8)));
}

template <typename T> class LocalSpectrum : public ::testing::Test {};
TYPED_TEST_SUITE(LocalSpectrum, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(LocalSpectrum, AnalyticMomentsTimeAndComplexCrossSpectra) { analytic_spectra(TypeParam{0}); }

TYPED_TEST(LocalSpectrum, InvalidInputsAndOverflow)
{
  using R = TypeParam;
  R const inf = uni20::numeric_limits<R>::infinity(), nan = uni20::numeric_limits<R>::quiet_NaN();
  R const max = uni20::numeric_limits<R>::max();
  std::array<bethe::LocalSpectralLine<R>, 1> local{{{0, R{1}, R{1}}}};
  for (R bad : {R{-1}, inf, nan})
  {
    local[0].gap = bad;
    EXPECT_THROW((void)bethe::local_spectral_moments<R>(local), std::invalid_argument);
    EXPECT_THROW((void)bethe::local_time_correlation<R>(local, R{0}), std::invalid_argument);
    local[0].gap = R{1};
    local[0].weight = bad;
    EXPECT_THROW((void)bethe::local_spectral_moments<R>(local), std::invalid_argument);
    EXPECT_THROW((void)bethe::local_time_correlation<R>(local, R{0}), std::invalid_argument);
    local[0].weight = R{1};
  }
  for (R bad : {inf, nan})
  {
    EXPECT_THROW((void)bethe::local_time_correlation<R>({}, bad), std::invalid_argument);
    EXPECT_THROW((void)bethe::cross_time_correlation<R>({}, bad), std::invalid_argument);
    EXPECT_THROW((void)bethe::spectral_time_correlation<R>(0, {}, bad), std::invalid_argument);
    std::array<bethe::CrossSpectralLine<R>, 1> cross{{{0, R{1}, {R{0}, bad}}}};
    EXPECT_THROW((void)bethe::cross_time_correlation<R>(cross, R{0}), std::invalid_argument);
    EXPECT_THROW((void)bethe::cross_spectral_moments<R>(cross), std::invalid_argument);
  }
  local[0] = {0, max, R{2}};
  EXPECT_THROW((void)bethe::local_spectral_moments<R>(local), std::overflow_error);
  EXPECT_THROW((void)bethe::local_time_correlation<R>(local, R{2}), std::overflow_error);
  std::array<bethe::LocalSpectralLine<R>, 2> overflow{{{0, R{0}, max}, {1, R{0}, max}}};
  EXPECT_THROW((void)bethe::local_time_correlation<R>(overflow, R{0}), std::overflow_error);
  EXPECT_THROW((void)bethe::local_spectral_moments<R>(overflow), std::overflow_error);
}

#if UNI20_HAS_MPREAL
TEST(LocalSpectrumMpfr, ExplicitPrecisionAndAnalyticChecks)
{
  using R = uni20::mpreal;
  auto const ambient = mpfr_get_default_prec();
  for (int bits : {80, 160, 256})
  {
    SCOPED_TRACE(bits);
    auto const p = uni20::Precision::bits(bits);
    R const zero(0, p), one(1, p);
    analytic_spectra(zero);
    auto const local = bethe::local_spectral_moments<R>({}, zero);
    auto const cross = bethe::cross_spectral_moments<R>({}, zero);
    for (auto const& value : {local.weight, local.first_moment, cross.value.real(), cross.first_moment.imag()})
      EXPECT_EQ(uni20::precision_of(value), p);
    EXPECT_EQ(uni20::precision_of(bethe::local_time_correlation<R>({}, one)), p);
    EXPECT_EQ(uni20::precision_of(bethe::cross_time_correlation<R>({}, one)), p);
    for (auto const& z : bethe::spectral_time_correlation<R>(2, {}, one))
      EXPECT_EQ(uni20::precision_of(z), p);
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

TEST(LocalSpectrumMpfr, RejectMixedExactAndUnsetScalars)
{
  using R = uni20::mpreal;
  using C = uni20::complex<R>;
  auto const p = uni20::Precision::bits(160);
  R const zero(0, p), one(1, p);
  std::array<bethe::LocalSpectralLine<R>, 1> local{{{0, one, one}}};
  std::array<bethe::CrossSpectralLine<R>, 1> cross{{{0, one, C(one, zero)}}};
  EXPECT_THROW((void)bethe::local_spectral_moments<R>({}), std::invalid_argument);
  EXPECT_THROW((void)bethe::cross_spectral_moments<R>({}), std::invalid_argument);
  for (R const& bad : {R{}, R{0}, R(0, uni20::Precision::bits(80))})
  {
    if (!bad.initialized() || bad.is_exact())
    {
      EXPECT_THROW((void)bethe::local_time_correlation<R>({}, bad), std::invalid_argument);
      EXPECT_THROW((void)bethe::cross_time_correlation<R>({}, bad), std::invalid_argument);
    }
    local[0].gap = bad;
    EXPECT_THROW((void)bethe::local_time_correlation<R>(local, one), std::invalid_argument);
    EXPECT_THROW((void)bethe::local_spectral_moments<R>(local, zero), std::invalid_argument);
    local[0].gap = one;
    local[0].weight = bad;
    EXPECT_THROW((void)bethe::local_time_correlation<R>(local, one), std::invalid_argument);
    local[0].weight = one;
  }
  for (C const& bad : {C{}, C(R{1}, R{0}), C(R(1, uni20::Precision::bits(80)), R(0, uni20::Precision::bits(80)))})
  {
    cross[0].coefficient = bad;
    EXPECT_THROW((void)bethe::cross_time_correlation<R>(cross, one), std::invalid_argument);
    EXPECT_THROW((void)bethe::cross_spectral_moments<R>(cross, zero), std::invalid_argument);
  }
}
#endif
} // namespace
