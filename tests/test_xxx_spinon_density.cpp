// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include "vertex_density_references.hpp"
#include <algorithm>
#include <array>
#include <bethe/prepared_observable_grid.hpp>
#include <bethe/xxx_spinon_density.hpp>
#include <bethe/xxx_thermodynamic_structure_factor.hpp>
#include <uni20/async/tbb_scheduler.hpp>

using bethe::heisenberg::SpinonDensityError;
using bethe::heisenberg::SpinonDensityOptions;
using bethe::heisenberg::XXXSpinonDensity;

template <class R> class SpinonDensity : public ::testing::Test {};
TYPED_TEST_SUITE(SpinonDensity, test_support::RealTypes, test_support::PrecisionNames);

template <class R> std::vector<R> bulk(XXXSpinonDensity<R> const& f)
{
  std::vector<R> beta;
  for (auto s : {"-2.1", "-1.3", "-.4", ".3", "1.2", "2"})
    beta.push_back(f.parse(s));
  return beta;
}

template <class R> void finite_approximation(XXXSpinonDensity<R> const& f)
{
  using uni20::math::abs;
  using C = uni20::complex<R>;
  auto beta = bulk(f);
  R eps = uni20::numeric_limits<R>::epsilon(f.precision());
  EXPECT_EQ(f.bits(), f.precision().bit_count());
  for (auto const& row : vertex_reference::densities)
  {
    SCOPED_TRACE(row.n);
    auto result = f.evaluate(std::span<R const>(beta).first(row.n));
    R expected = f.parse(row.value);
    EXPECT_LE(abs(result.value / expected - R{1}), R{1000000} * eps);
    EXPECT_LE(abs(result.value - expected), result.arithmetic_indicator);
    EXPECT_LE(abs(result.contour_norm / f.parse(row.norm) - R{1}), R{1000000} * eps);
    for (auto const& amplitude : vertex_reference::amplitudes)
      if (amplitude.n == row.n)
      {
        C reference(f.parse(amplitude.re), f.parse(amplitude.im));
        EXPECT_LE(abs(result.amplitudes[amplitude.index] - reference), R{1000000} * eps * abs(reference));
      }
  }
}

TYPED_TEST(SpinonDensity, IndependentFiniteContourFixtures) { finite_approximation(XXXSpinonDensity<TypeParam>{}); }

TYPED_TEST(SpinonDensity, ReconstructTwoSpinonSpectrum)
{
  using R = TypeParam;
  using uni20::math::abs;
  using uni20::math::atan;
  using uni20::math::cos;
  using uni20::math::exp;
  using uni20::math::sin;
  XXXSpinonDensity<R> f;
  R const pi = uni20::pi<R>;
  for (auto inputs : {std::array{"-2.1", "-1.3"}, std::array{"-.7", "1.2"}, std::array{"1.1", "2.3"}})
    for (R exchange : {R{1}, R{3}})
    {
      std::array<R, 2> b{f.parse(inputs[0]), f.parse(inputs[1])};
      R p = R{2} * atan(exp(b[0])), r = R{2} * atan(exp(b[1]));
      R q = p + r, energy = exchange * pi / R{2} * (sin(p) + sin(r));
      R jacobian = exchange * pi / R{2} * abs(cos(p) - cos(r));
      // Both labelled solutions of the conservation laws; W2 already has 1/2!.
      auto density = f.evaluate(b);
      R spectral = R{2} * density.value / jacobian;
      bethe::heisenberg::ThermodynamicTwoSpinonStructureFactor<R> direct(exchange);
      auto expected = direct(q, energy);
      ASSERT_TRUE(expected.value);
      ASSERT_TRUE(expected.error);
      EXPECT_REAL_NEAR(spectral, *expected.value,
                       R{2} * density.arithmetic_indicator / jacobian + *expected.error + f.parse("1e-24") * spectral);
    }
}

template <class R> void four_spinon_independent(XXXSpinonDensity<R> const& f)
{
  using uni20::math::abs;
  // Independent 8F7 sums and oscillatory pair integrals, multiplied by the
  // independently derived physical coefficient G(1/2)^16/3072. These values
  // contain no C++ prefactor, vertex-contour assembly or sum-rule fit.
  char const* cases[][5] = {{"-.6", "-.1", ".3", ".8", "0.0448573680600201527455311253220848271056983579181"},
                            {"-1", "-.4", ".2", "1.3", "1.45280512128837090622630719799604472657153619211"},
                            {"-2.204455496780781151440805630060717051913", "-2.029605920497355931277048454253291224376",
                             "-0.1789178381323832258356584261434473190251",
                             "0.1789157475173640394576444381598228210966",
                             "2159649.09516248006784124665512490013454103392657"}};
  R pi = bethe::detail::pi_like(f.parse("1"));
  for (auto const& c : cases)
  {
    std::array<R, 4> beta;
    for (int j = 0; j < 4; ++j)
      beta[j] = R{2} * pi * f.parse(c[j]);
    auto result = f.evaluate(beta, {64, 16});
    R expected = f.parse(c[4]);
    EXPECT_LE(abs(result.value - expected),
              result.contour_tail_error + result.arithmetic_indicator + abs(expected) * f.parse("1e-32"));
    EXPECT_LE(result.contour_tail_error / expected, f.parse("1e-8"));
    // Reported errors must be useful, not merely large enough to hide 2x errors.
    EXPECT_LE(result.arithmetic_indicator / expected, f.parse("1e-6"));
  }
}

TYPED_TEST(SpinonDensity, IndependentPhysicalFourSpinonNormalization)
{
  four_spinon_independent(XXXSpinonDensity<TypeParam>{});
}

TYPED_TEST(SpinonDensity, PermutationsAndUnresolvedInputs)
{
  using R = TypeParam;
  using uni20::math::abs;
  XXXSpinonDensity<R> f;
  auto beta = bulk(f);
  for (std::size_t n : {2u, 4u, 6u})
  {
    std::vector<R> b(beta.begin(), beta.begin() + n);
    auto first = f.evaluate(b);
    std::reverse(b.begin(), b.end());
    auto second = f.evaluate(b);
    EXPECT_LE(abs(first.value - second.value), first.contour_tail_error + second.contour_tail_error +
                                                   first.arithmetic_indicator + second.arithmetic_indicator);
  }
  EXPECT_THROW(f.evaluate(std::span<R const>(beta).first(3)), std::invalid_argument);
  EXPECT_THROW(f.evaluate(beta, {1, 12}), std::invalid_argument);
  EXPECT_THROW(f.evaluate(beta, {32, 0}), std::invalid_argument);
  beta[1] = beta[0];
  EXPECT_THROW(f.evaluate(beta), SpinonDensityError);
  beta[1] = uni20::numeric_limits<R>::infinity();
  EXPECT_THROW(f.evaluate(beta), std::invalid_argument);
}

TYPED_TEST(SpinonDensity, TailAgainstHurwitzZetaAndWorkspaceReuse)
{
  using R = TypeParam;
  using C = uni20::complex<R>;
  using uni20::math::abs;
  XXXSpinonDensity<R> f;
  bethe::detail::ControlledTailPlan<R> plan(f.parse("1"), 8);
  bethe::detail::ControlledTailWorkspace<R> workspace;
  for (int repeat = 0; repeat < 2; ++repeat)
    for (auto const& row : vertex_reference::tails)
    {
      C shift(f.parse(".25"), f.parse(".75"));
      std::vector<C> a(row.power, shift), b(row.power, shift + R{1});
      auto result = bethe::detail::controlled_hypergeometric_tail<R>(a, b, row.power, 16, plan, workspace);
      C expected(f.parse(row.re), f.parse(row.im));
      EXPECT_LE(abs(result.multiplier - expected), result.error_multiplier + result.arithmetic_indicator);
    }
}

template <class R> void compensated_components(R const& one)
{
  using C = uni20::complex<R>;
  R const zero = one - one;
  // Each component loses its unit in ordinary addition. The dominant
  // component changes between terms, so a complex magnitude comparison
  // cannot select the real Neumaier branch for both components.
  R const large = R{8} * one / uni20::numeric_limits<R>::epsilon(one);
  for (int sign : {-1, 1})
  {
    bethe::detail::SpinonSum<C> sum(C(zero, zero));
    sum.add(C(large, R(sign) * one));
    sum.add(C(one, R(sign) * large));
    sum.add(C(-large, -R(sign) * large));
    EXPECT_EQ(sum.value().real(), one);
    EXPECT_EQ(sum.value().imag(), R(sign) * one);
  }
}

TYPED_TEST(SpinonDensity, CompensatesComplexComponentsIndependently) { compensated_components(TypeParam{1}); }

TYPED_TEST(SpinonDensity, RejectsSubnormalWeightsWithRelativeErrorModel)
{
  using R = TypeParam;
  XXXSpinonDensity<R> f;
  std::array<R, 2> beta{R{0}, uni20::math::sqrt(uni20::numeric_limits<R>::min()) / R{100}};
  try
  {
    auto result = f.evaluate(beta);
    ADD_FAILURE() << "accepted subnormal density " << uni20::format_scalar(result.value)
                  << " with arithmetic indicator " << uni20::format_scalar(result.arithmetic_indicator);
  }
  catch (SpinonDensityError const& error)
  {
    EXPECT_EQ(error.failure(), bethe::heisenberg::SpinonDensityFailure::exponent_range);
  }
}

TEST(SpinonDensityFailures, FailuresHaveReasonsAndNeverBecomeZero)
{
  using bethe::heisenberg::SpinonDensityFailure;
  XXXSpinonDensity<double> f;
  auto beta = bulk(f);
  auto reason = [&](SpinonDensityOptions options) {
    try
    {
      (void)f.evaluate(beta, options);
      ADD_FAILURE() << "expected an explicit numerical failure";
    }
    catch (SpinonDensityError const& error)
    {
      return error.failure();
    }
    return SpinonDensityFailure::none;
  };
  beta[0] = -20.; // parameter modulus exceeds the tail start of two
  EXPECT_EQ(reason({2, 12}), SpinonDensityFailure::tail_domain);
  beta = bulk(f);
  beta[1] = beta[0];
  EXPECT_EQ(reason({32, 12}), SpinonDensityFailure::rounded_collision);
  beta = {-1000., -800., 1., 2.};
  EXPECT_NE(reason({32, 12}), SpinonDensityFailure::none);
  beta = {0., 1e-160}; // nonzero binary64 weight, but a zero relative-error indicator before the fix
  EXPECT_EQ(reason({32, 12}), SpinonDensityFailure::exponent_range);
  EXPECT_THROW(f.evaluate(beta, {std::numeric_limits<int>::max(), 12}), std::invalid_argument);
}

#if UNI20_HAS_MPREAL
TEST(SpinonDensityMpfr, BeyondBinary128WithoutAmbientPrecision)
{
  using R = uni20::mpreal;
  using uni20::math::abs;
  auto p160 = uni20::Precision::bits(160), p256 = uni20::Precision::bits(256);
  XXXSpinonDensity<R> low(p160), high(p256);
  compensated_components(low.parse("1"));
  compensated_components(high.parse("1"));
  std::array<R, 2> tiny_beta{high.parse("0"), high.parse("1e-160")};
  auto tiny = high.evaluate(tiny_beta);
  EXPECT_GT(tiny.value, high.parse("0"));
  EXPECT_GT(tiny.arithmetic_indicator, high.parse("0"));
  EXPECT_LT(tiny.arithmetic_indicator / tiny.value, high.parse("1e-65"));
  finite_approximation(low);
  finite_approximation(high);
  auto a = low.evaluate(bulk(low)), b = high.evaluate(bulk(high));
  R expected = high.parse(vertex_reference::densities[2].value);
  EXPECT_LT(abs(b.value - expected), abs(a.value.at(p256) - expected) / R{1000000});
  EXPECT_LT(abs(b.value / expected - R{1}), high.parse("1e-65"));
  auto mixed = bulk(high);
  mixed[0] = mixed[0].at(p160);
  EXPECT_THROW(high.evaluate(mixed), std::invalid_argument);
  mixed[0] = R{0};
  EXPECT_THROW(high.evaluate(mixed), std::invalid_argument);
  four_spinon_independent(high);
}

TEST(SpinonDensityMpfr, SharedContextAndIndependentPrecisionsOnOneScheduler)
{
  using R = uni20::mpreal;
  uni20::async::TbbScheduler scheduler(2);
  std::vector<std::vector<R>> expected(2);
  auto submit = [&](int lane, int bits) {
    auto p = uni20::Precision::bits(bits);
    XXXSpinonDensity<R> serial(p);
    for (int i = 0; i < 6; ++i)
    {
      auto beta = bulk(serial);
      beta.resize(2 * (1 + i % 3));
      beta.back() += serial.parse(".01") * R(i);
      expected[lane].push_back(serial.evaluate(beta).value);
    }
    auto context = bethe::prepare_background_async([p] { return XXXSpinonDensity<R>(p); }, {&scheduler});
    return bethe::evaluate_observable_grid_async(
        context.read(), 6,
        [](XXXSpinonDensity<R> const& f, std::size_t i) {
          auto beta = bulk(f);
          beta.resize(2 * (1 + i % 3));
          beta.back() += f.parse(".01") * R(i);
          return f.evaluate(beta).value;
        },
        [&, lane, p](std::size_t i, R value) {
          EXPECT_EQ(value.precision(), p);
          EXPECT_EQ(value, expected[lane][i]);
        },
        {.scheduler = &scheduler, .max_in_flight = 2});
  };
  auto a = submit(0, 160), b = submit(1, 256);
  EXPECT_EQ(a.get_wait(scheduler), 6u);
  EXPECT_EQ(b.get_wait(scheduler), 6u);
}
#endif
