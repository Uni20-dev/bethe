// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/detail/gauss_legendre.hpp>
#include <bethe/xxx_thermodynamic_structure_factor.hpp>
#include <bethe/xxz_thermodynamic_structure_factor.hpp>

namespace model = bethe::xxz;
template <typename T> class XXZThermodynamicDSF : public ::testing::Test {};
TYPED_TEST_SUITE(XXZThermodynamicDSF, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(XXZThermodynamicDSF, HighPrecisionReferences)
{
  using Real = TypeParam;
  char const* points[][4] = {{"1.5", "1.5", "2", "4.159749286581976836416246520552805055874094165803454070291"},
                             {"2", "1", "2.3", "0.5781035350458511934434787838145718029104220916152404592048"},
                             {"4", "2", "4", "3.470018047561235624535617914580890209465514431388532448227"}};
  for (auto const& p : points)
  {
    model::ThermodynamicTwoSpinonStructureFactor<Real> sf(uni20::parse_real<Real>(p[0]));
    auto const r = sf(uni20::parse_real<Real>(p[1]), uni20::parse_real<Real>(p[2]));
    ASSERT_TRUE(r.value) << p[0];
    Real const expected = uni20::parse_real<Real>(p[3]);
    EXPECT_REAL_NEAR(*r.value, expected, Real{50000} * uni20::numeric_limits<Real>::epsilon());
    EXPECT_LE(std::abs(*r.value - expected), *r.error);
  }
}

TYPED_TEST(XXZThermodynamicDSF, SymmetryScalingThresholdsAndBudget)
{
  using Real = TypeParam;
  Real const pi = bethe::detail::pi<Real>();
  model::ThermodynamicTwoSpinonStructureFactor<Real> sf(Real{1.5}), scaled(Real{1.5}, Real{3});
  auto const a = sf(Real{1.5}, Real{2}), b = sf(Real{2} * pi - Real{1.5}, Real{2}), c = scaled(Real{1.5}, Real{6});
  ASSERT_TRUE(a.value);
  ASSERT_TRUE(b.value);
  ASSERT_TRUE(c.value);
  Real const tol = Real{100000} * uni20::numeric_limits<Real>::epsilon();
  EXPECT_REAL_NEAR(*a.value, *b.value, tol);
  EXPECT_REAL_NEAR(*a.value, Real{3} * *c.value, tol);
  EXPECT_EQ(*sf(Real{1}, Real{-1}).value, Real{0});
  auto const band = sf.sheets(Real{1.5});
  EXPECT_EQ(sf(Real{1.5}, band[0].upper).status, model::SpectralDensityStatus::threshold);
  EXPECT_EQ(sf(Real{1.5}, band[0].lower).status, model::SpectralDensityStatus::threshold);
  auto const zero = sf(Real{0}, Real{2}), staggered = sf(pi, Real{2});
  ASSERT_TRUE(zero.value);
  ASSERT_TRUE(staggered.value);
  EXPECT_GT(*zero.value, Real{0});
  EXPECT_REAL_NEAR(*zero.value / *staggered.value, sf.gap() / sf.maximum_energy(), tol);
  model::ThermodynamicTwoSpinonStructureFactor<Real> failed(Real{1.5}, Real{1}, {.max_terms = 0});
  EXPECT_FALSE(failed(Real{1.5}, Real{2}).converged());
  EXPECT_FALSE(failed(Real{1.5}, Real{2}).value);
  EXPECT_THROW(sf(Real{-1}, Real{2}), std::invalid_argument);
  EXPECT_THROW(sf(Real{1}, uni20::numeric_limits<Real>::quiet_NaN()), std::invalid_argument);
  EXPECT_THROW((model::ThermodynamicTwoSpinonStructureFactor<Real>(Real{1})), std::invalid_argument);
  EXPECT_THROW((model::ThermodynamicTwoSpinonStructureFactor<Real>(Real{1.5}, Real{0})), std::invalid_argument);
}

TYPED_TEST(XXZThermodynamicDSF, TinyMassAtExactStaggeredMomentum)
{
  using Real = TypeParam;
  model::ThermodynamicTwoSpinonStructureFactor<Real> sf(uni20::parse_real<Real>("1.001"));
  Real const pi = bethe::detail::pi<Real>();
  auto const sheets = sf.sheets(pi);
  EXPECT_FALSE(sheets[0].present);
  ASSERT_TRUE(sheets[1].present);
  EXPECT_REAL_NEAR(sheets[1].lower, Real{2} * sf.gap(), Real{128} * uni20::numeric_limits<Real>::epsilon() * sf.gap());
}

TEST(XXZThermodynamicDSFLimit, ApproachesIndependentXXXEvaluator)
{
  bethe::heisenberg::ThermodynamicTwoSpinonStructureFactor<double> xxx;
  // Different evaluator and spectral convention: XXZ returns raising,
  // XXX returns zz. Keep J fixed; no fitting or intensity rescaling.
  for (double k : {.8, 1.5, 2.3})
  {
    auto [lo, hi] = xxx.boundaries(k);
    for (double fraction : {.25, .6})
    {
      double const w = lo + (hi - lo) * fraction;
      auto reference = xxx(k, w);
      ASSERT_TRUE(reference.value);
      double previous = 100;
      for (double delta : {1.01, 1.001, 1.0001})
      {
        auto const result = model::ThermodynamicTwoSpinonStructureFactor<double>(delta)(k, w);
        ASSERT_TRUE(result.value) << delta << " " << k;
        double const error = std::abs(*result.value / (2 * *reference.value) - 1);
        EXPECT_LT(error, previous);
        previous = error;
      }
      EXPECT_LT(previous, .02);
    }
  }
}

#if UNI20_HAS_MPREAL
TEST(XXZThermodynamicDSFMpfr, ExplicitPrecisionPointsBoundariesAndRejection)
{
  using R = uni20::mpreal;
  auto const ambient = mpfr_get_default_prec();
  char const* points[][4] = {
      {"1.5", "1.5", "2",
       "4."
       "1597492865819768364162465205528050558740941658034540702908692161256372658275335756250521366520670741337981378"},
      {"2", "1", "2.3",
       "0."
       "5781035350458511934434787838145718029104220916152404592047599422865903982329705536337680912851919544351941621"
       "9"},
      {"4", "2", "4",
       "3."
       "47001804756123562453561791458089020946551443138853244822721563090560035305566319127316704042331082732437839"}};
  for (int bits : {128, 256})
  {
    auto const p = uni20::Precision::bits(bits);
    auto parse = [p](char const* text) { return uni20::parse_real<R>(text, p); };
    R const zero = uni20::scalar_at_precision<R>(0, p), one = uni20::scalar_like(zero, 1);
    R const pi = bethe::detail::spinon_pi(zero);
    R const eps = uni20::numeric_limits<R>::epsilon(p);
    for (auto const& row : points)
    {
      model::ThermodynamicTwoSpinonStructureFactor<R> sf(parse(row[0]));
      R const k = parse(row[1]), w = parse(row[2]);
      auto value = sf(k, w);
      ASSERT_TRUE(value.value && value.error);
      EXPECT_EQ(uni20::precision_of(*value.value), p);
      EXPECT_EQ(uni20::precision_of(*value.error), p);
      R const reference = parse(row[3]);
      EXPECT_REAL_NEAR(*value.value, reference, 50000 * eps);
      EXPECT_LE(uni20::math::abs(*value.value - reference), *value.error);
      EXPECT_REAL_NEAR(*sf(2 * pi - k, w).value, *value.value, 100000 * eps);
      auto outside = sf(k, -one);
      ASSERT_TRUE(outside.value && outside.error);
      EXPECT_EQ(uni20::precision_of(*outside.value), p);
      EXPECT_EQ(uni20::precision_of(*outside.error), p);
      EXPECT_EQ(*outside.value, zero);
      for (auto const& sheet : sf.sheets(k))
      {
        EXPECT_EQ(uni20::precision_of(sheet.lower), p);
        EXPECT_EQ(uni20::precision_of(sheet.upper), p);
        if (sheet.present)
        {
          EXPECT_EQ(sf(k, sheet.lower).status, model::SpectralDensityStatus::threshold);
          EXPECT_EQ(sf(k, sheet.upper).status, model::SpectralDensityStatus::threshold);
        }
      }
      EXPECT_THROW((void)sf(k.at(uni20::Precision::bits(80)), w), std::invalid_argument);
      EXPECT_THROW((void)sf(k, w.at(uni20::Precision::bits(80))), std::invalid_argument);
      model::TwoSpinonStructureFactorOptions<R> options(p);
      options.max_terms = 0;
      auto fail = model::ThermodynamicTwoSpinonStructureFactor<R>(parse(row[0]), one, options)(k, w);
      EXPECT_FALSE(fail.value);
      EXPECT_EQ(fail.status, model::SpectralDensityStatus::numerical_failure);
      options.tolerance = parse("0");
      EXPECT_THROW((model::ThermodynamicTwoSpinonStructureFactor<R>(parse(row[0]), one, options)),
                   std::invalid_argument);
    }
    model::ThermodynamicTwoSpinonStructureFactor<R> sf(parse("1.5"));
    auto q0 = sf(zero, parse("2")), qpi = sf(pi, parse("2"));
    ASSERT_TRUE(q0.value && qpi.value);
    EXPECT_REAL_NEAR(*q0.value / *qpi.value, sf.gap() / sf.maximum_energy(), 100000 * eps);
    model::ThermodynamicTwoSpinonStructureFactor<R> tiny(parse("1.001"));
    auto edge = tiny.sheets(pi);
    EXPECT_FALSE(edge[0].present);
    EXPECT_EQ(uni20::precision_of(edge[0].lower), p);
    EXPECT_EQ(uni20::precision_of(edge[0].upper), p);
    EXPECT_REAL_NEAR(edge[1].lower, 2 * tiny.gap(), 128 * eps * tiny.gap());
    model::TwoSpinonStructureFactorOptions<R> other(uni20::Precision::bits(80));
    EXPECT_THROW((model::ThermodynamicTwoSpinonStructureFactor<R>(parse("1.5"), one, other)), std::invalid_argument);
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
#endif

namespace
{
std::pair<double, double> integrated_spectrum(double delta, int order)
{
  model::ThermodynamicTwoSpinonStructureFactor<double> sf(delta);
  auto const rule = bethe::detail::gauss_legendre<double>(order);
  double const pi = bethe::detail::pi<double>();
  double const critical = 2 * std::atan(std::sqrt(sf.gap() / sf.maximum_energy()));
  std::array<double, 4> const cuts{0, critical, pi - critical, pi};
  bethe::detail::CompensatedSum<double> intensity(0.), moment(0.);
  for (std::size_t panel = 0; panel + 1 < cuts.size(); ++panel)
    for (std::size_t i = 0; i < rule.x.size(); ++i)
    {
      double const width = cuts[panel + 1] - cuts[panel];
      double const k = cuts[panel] + width * (1 + rule.x[i]) / 2;
      double const kw = width * rule.w[i] / (2 * pi);
      std::vector<double> edges;
      for (auto const& band : sf.sheets(k))
        if (band.present)
        {
          edges.push_back(band.lower);
          edges.push_back(band.upper);
        }
      std::sort(edges.begin(), edges.end());
      edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
      for (std::size_t panel_w = 0; panel_w + 1 < edges.size(); ++panel_w)
        for (std::size_t j = 0; j < rule.x.size(); ++j)
        {
          double const t = pi * (1 + rule.x[j]) / 4;
          double const dw = edges[panel_w + 1] - edges[panel_w];
          double const w = edges[panel_w] + dw * std::sin(t) * std::sin(t);
          auto const value = sf(k, w);
          if (!value.value)
            throw std::runtime_error("unresolved integrated spectral node k=" + uni20::format_scalar(k) +
                                     " w=" + uni20::format_scalar(w));
          double const weight = kw * rule.w[j] * dw * std::sin(2 * t) / 8;
          intensity.add(weight * *value.value);
          moment.add(weight * w * *value.value);
        }
    }
  return {intensity.value(), moment.value()};
}
} // namespace

TEST(XXZThermodynamicDSFSumRules, OriginalVertexMeasure)
{
  // The reference values come from full-domain integration of original
  // two-vacuum q-product amplitudes, before imposing energy/momentum.
  // This test instead integrates the resolved DSF with both 2*pi measures.
  auto const coarse = integrated_spectrum(1.5, 32), fine = integrated_spectrum(1.5, 64);
  EXPECT_NEAR(coarse.first, fine.first, 1e-7);
  EXPECT_NEAR(fine.first, .4846429124574824, 1e-8);
  EXPECT_NEAR(fine.second, .7491062383225591, 1e-8);
  auto const ising_side = integrated_spectrum(4, 64);
  EXPECT_NEAR(ising_side.first, .4997481118452324, 1e-8);
  EXPECT_NEAR(ising_side.second, 1.9999994798394005, 1e-8);
}

TEST(XXZThermodynamicDSFSumRules, IsotropicIntensityAndFirstMoment)
{
  double const full_xxx_first = 4 * (std::log(2.) - .25) / 3;
  double previous_weight = 1, previous_moment = 1;
  for (double delta : {1.01, 1.001, 1.0001})
  {
    auto const [weight, moment] = integrated_spectrum(delta, 64);
    double const weight_error = std::abs(2 * weight - .7289);
    double const moment_error = std::abs(moment / full_xxx_first - .7130);
    std::cout << "Delta=" << delta << " intensity fraction=" << 2 * weight << " first/XXX=" << moment / full_xxx_first
              << '\n';
    EXPECT_LT(weight_error, previous_weight);
    EXPECT_LT(moment_error, previous_moment);
    previous_weight = weight_error;
    previous_moment = moment_error;
  }
  EXPECT_LT(previous_weight, .01);
  EXPECT_LT(previous_moment, .01);
}

TYPED_TEST(XXZThermodynamicDSF, InverseEllipticMapAndRange)
{
  using R = TypeParam;
  R const eps = uni20::numeric_limits<R>::epsilon(), pi = bethe::detail::pi<R>();
  EXPECT_REAL_NEAR(bethe::detail::carlson_rf_unit(R{0}, R{1}), pi / 2, 64 * eps);
  EXPECT_EQ(bethe::detail::carlson_rf_unit(R{1}, R{1}), R{1});
  for (R x : {R{1} / 4, eps * eps})
  {
    R const expected = uni20::math::acosh(1 / uni20::math::sqrt(x)) / uni20::math::sqrt(1 - x);
    EXPECT_REAL_NEAR(bethe::detail::carlson_rf_unit(x, x), expected, 256 * eps * expected);
  }
  EXPECT_THROW((void)bethe::detail::carlson_rf_unit(R{0}, R{0}), std::invalid_argument);
  EXPECT_THROW((void)bethe::detail::carlson_rf_unit(R{-1}, R{1}), std::invalid_argument);
  EXPECT_THROW((void)bethe::detail::carlson_rf_unit(R{2}, R{1}), std::invalid_argument);
}

namespace
{
template <uni20::Real R> void check_edge_points(R const& zero, bool wide)
{
  // 130-digit independent direct ellipf/theta evaluation, at fixed decimal
  // frequencies on either side of the internal sheet edge as well as the
  // exterior thresholds. These are finite offsets, never edge values.
  char const* rows[][3] = {
      {"1e-8",
       "1."
       "863643318437311673525067920463483091368239325314048593865926903987534241831156351600294764516976582383254007764"
       "207",
       "19895."
       "7778088541614359672895701435862680368524606170979936424697691982326756586659696718625021534778026504597864033"},
      {"1e-8",
       "2."
       "671358579744309470504687559157393093646582205378312937269028535979951556526904398932084180611707245546047785347"
       "117",
       "0."
       "000007809028696571039859826651022119614927952138463967208473668914421433881781550947185588912618246716090939664"
       "883778328"},
      {"1e-8",
       "2."
       "488827226915823926378113016393874899403988561369864543609081596074576836825929883702373072875094155250396678236"
       "337",
       "0."
       "049005979365788938900065157130817256600870055230121731359712523354229406995254636331078103342785012915341261826"
       "97986"},
      {"1e-8",
       "2."
       "488827239419502382528198127216297167373359963306814970328514595657064105513249491386403149213837701295599124885"
       "658",
       "0."
       "048807059744455220544353285516569975269985153108885343339190801123872674717861420600941113029573072961482534970"
       "26924"},
      {"1e-24",
       "1."
       "863643310360158898912035266121894924595172946688930914630932163599503178410705339335426038319022593224914924758"
       "903",
       "1989577799130."
       "844364920528533080515248309908220191305546890533457591971482694938692400792817844196776067583373229447"},
      {"1e-24",
       "2."
       "671358587821462245117720213498981260419648584003430616504023276367982619947355411196952906809661234704386868352"
       "42",
       "0."
       "000000000000078090286149719817076489588654348633916658588327206301114418378385744764655230793318654550309737064"
       "68900226852183689"},
      {"1e-24",
       "2."
       "488827233167663154453154946621163225884418721217226358500227999018299135197939708420024676678485544071494084623"
       "82",
       "0."
       "048807061114279030966457062489731484193239041002847580528819096384146255362750809072821665603488083995514331525"
       "95508"},
      {"1e-24",
       "2."
       "488827233167663154453156196989008840892929803459453155437368192713341807141239666668751545410446312474501718498"
       "175",
       "0."
       "048807061112289862119934245552592423074577967474516740452036479145391762772484397801439365930368284178164961910"
       "24441"}};
  auto parse = [&](char const* text) {
    if constexpr (uni20::has_runtime_precision_v<R>)
      return uni20::parse_real<R>(text, uni20::precision_of(zero));
    else
      return uni20::parse_real<R>(text);
  };
  R const eps = uni20::numeric_limits<R>::epsilon(zero);
  model::ThermodynamicTwoSpinonStructureFactor<R> sf(parse("1.5"));
  for (auto const& row : rows)
  {
    if (!wide && std::string_view(row[0]) == "1e-24") continue;
    SCOPED_TRACE(row[1]);
    auto const result = sf(parse("1.5"), parse(row[1]));
    ASSERT_EQ(result.status, model::SpectralDensityStatus::converged);
    ASSERT_TRUE(result.value && result.error);
    R const expected = parse(row[2]);
    R const difference = uni20::math::abs(*result.value - expected);
    EXPECT_LE(difference, *result.error + 4 * eps * expected);
    EXPECT_LE(difference / expected, 65536 * eps / parse(row[0]));
    EXPECT_EQ(uni20::precision_of(*result.value), uni20::precision_of(zero));
  }
}
} // namespace

TYPED_TEST(XXZThermodynamicDSF, IndependentNearEdgeReferences)
{
  check_edge_points(TypeParam{0}, uni20::numeric_limits<TypeParam>::digits() >= 100);
}
#if UNI20_HAS_MPREAL
TEST(XXZThermodynamicDSFMpfr, IndependentNearEdgeReferences)
{
  for (int bits : {128, 256})
    check_edge_points(uni20::mpreal(0, uni20::Precision::bits(bits)), true);
}
#endif

TEST(XXZThermodynamicDSFSumRules, IsingSaturation)
{
  double previous_weight = 1, previous_first = 1;
  for (double delta : {4., 8., 16.})
  {
    auto const [weight, moment] = integrated_spectrum(delta, 64);
    double const weight_error = std::abs(2 * weight - 1);
    double const first_error = std::abs(2 * moment / delta - 1);
    std::cout << "Delta=" << delta << " intensity fraction=" << 2 * weight << " first/Ising=" << 2 * moment / delta
              << '\n';
    EXPECT_LT(weight_error, previous_weight);
    EXPECT_LT(first_error, previous_first);
    previous_weight = weight_error;
    previous_first = first_error;
  }
  EXPECT_LT(previous_weight, 3e-6);
  EXPECT_LT(previous_first, 1e-9);
}
