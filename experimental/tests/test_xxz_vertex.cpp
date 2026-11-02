// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/experimental/detail/xxz_spectral_kinematics.hpp>
#include <bethe/experimental/xxz_vertex.hpp>
#include <bethe/experimental/xxz_vertex_spectrum.hpp>
#include <bethe/xxz_thermodynamic_structure_factor.hpp>

namespace
{
constexpr char const* amplitude_reference[8][2] = {
    {"0.23666140286824887438939348546561778634780377639490812694233051969606055220387128377315681876899",
     "1.1081814187441899962363340310547574734815524994460421592390996770497158080483117354970392105717"},
    {"0.64957156388327624178442602200584695464731891913534744379567627796796588175868326204365435469355",
     "-0.46200165256247321693293936983590122712312111460055753741769256172265301948057102186599218589656"},
    {"-0.39597187343727603574565947388127013313487860288453132553188433561980335090345410475969552845921",
     "0.4798479339005385173739509742072715134163924612693870920858018719025962802897607786153464430171"},
    {"-0.64176196456149270429025712983948618244241723532836997785362784222113927390985005148143437722765",
     "-0.94809406716071194969739632265230606656496082093857937442887649952757793554777509706094901248989"},
    {"0.30238041835285083606493963848854855056131374114879157843747524918411566882513404796029353286896",
     "0.17544199480758929735572259075967980155983482318726481871655642062813102565339594910121423726009"},
    {"0.44612381406382297502772420523639066873364369775278633080135601246404051348709368914645983149718",
     "-0.29575651465659417775598612140630175892109619512781439846663766103526811006154084481175194323426"},
    {"-0.66497769376873688546727560374532489667381348473604321810736684717510182301552932791364809102332",
     "0.13981804779106014248232286776582321807498154436765716586583059544928607367550472288006818734632"},
    {"-0.19364083111043985364970505414973310427079401886670403375646411914605992256238391476201915054496",
     "-0.31077357468231429858684972332368110580385247390609090126103971366066287080700273178915005573954"}};
}

template <typename T> class XXZVertex : public ::testing::Test {};
TYPED_TEST_SUITE(XXZVertex, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(XXZVertex, RejectsExcessiveAutomaticCounts)
{
  using R = TypeParam;
  // Check the work budget in Real before converting, rather than relying
  // on the host integer range as precision increases.
  R const delta = R{1} + uni20::numeric_limits<R>::epsilon();
  EXPECT_THROW((bethe::experimental::xxz::VertexEvaluator<R>(delta)), std::invalid_argument);
  EXPECT_THROW((bethe::experimental::xxz::VertexEvaluator<R>(R{1.5}, {20001, 8})), std::invalid_argument);
  EXPECT_THROW((bethe::experimental::xxz::VertexEvaluator<R>(R{1.5}, {24, 20001})), std::invalid_argument);
}

TYPED_TEST(XXZVertex, FixedRapidityNormalizationControls)
{
  using R = TypeParam;
  // Independent finite-q q-product values from the normalization audit,
  // before completeness/Fourier factors. Not a fitted spectral integral.
  std::array<R, 4> const rho{uni20::parse_real<R>("-.6"), uni20::parse_real<R>("-.1"), uni20::parse_real<R>(".3"),
                             uni20::parse_real<R>(".8")};
  R const pi = bethe::detail::pi<R>();
  char const* const controls[][2] = {{".1", "541.128035788733"}, {".05", "529.867275923635"}};
  for (auto const& control : controls)
  {
    R const eta = uni20::parse_real<R>(control[0]);
    std::array<R, 4> phi;
    R jacobian{1};
    for (std::size_t j = 0; j < rho.size(); ++j)
    {
      phi[j] = R{2} * eta * rho[j];
      jacobian *= eta * uni20::math::cosh(R{2} * pi * rho[j]);
    }
    auto const amplitudes = bethe::experimental::xxz::VertexEvaluator<R>(uni20::math::cosh(eta))(phi);
    R density{0};
    for (auto value : amplitudes.value[0])
      density += value.real() * value.real() + value.imag() * value.imag();
    EXPECT_REAL_NEAR(jacobian * density, uni20::parse_real<R>(control[1]), uni20::parse_real<R>("1e-8"));
  }
}

TYPED_TEST(XXZVertex, ScaledProductsCancelBeyondScalarRange)
{
  using R = TypeParam;
  using C = uni20::complex<R>;
  using S = bethe::experimental::detail::ScaledComplexProduct<R>;
  S large, small;
  R const factor = uni20::math::sqrt(uni20::numeric_limits<R>::max()) / R{2};
  for (int j = 0; j < 8; ++j)
  {
    large *= C{factor};
    small /= C{factor};
  }
  EXPECT_THROW(large.value(), std::overflow_error);
  large *= small;
  EXPECT_REAL_NEAR(large.value().real(), R{1}, R{64} * uni20::numeric_limits<R>::epsilon());
}
TYPED_TEST(XXZVertex, KinematicsRoundTripAndEndpoints)
{
  using R = TypeParam;
  R const pi = bethe::detail::pi<R>(), eps = uni20::numeric_limits<R>::epsilon();
  for (R delta : {R{1.5}, R{4}, R{40}})
  {
    bethe::experimental::xxz::SpectralKinematics<R> model(delta);
    for (int j = -7; j <= 7; ++j)
    {
      R const angle = R(j) / R{5};
      EXPECT_REAL_NEAR(model.angle(model.momentum(angle)), angle, R{200000} * eps);
    }
    EXPECT_EQ(model.momentum(pi / R{2}), R{0});
    EXPECT_EQ(model.angle(pi), -pi / R{2});
    EXPECT_EQ(model.energy(R{0}), model.gap());
  }
}
TYPED_TEST(XXZVertex, AmplitudeTruncationAndCollision)
{
  using R = TypeParam;
  std::array<R, 4> const phi{R{-1.2L}, R{-.2L}, R{.6L}, R{1.3L}};
  bethe::experimental::xxz::VertexEvaluator<R> eval(R{1.5});
  auto opts = eval.options();
  opts.products += 8;
  opts.residues += 4;
  auto const a = eval(phi), b = bethe::experimental::xxz::VertexEvaluator<R>(R{1.5}, opts)(phi);
  for (int v = 0; v < 2; ++v)
    for (int l = 0; l < 4; ++l)
      EXPECT_REAL_NEAR(uni20::math::abs(a.value[v][l] - b.value[v][l]), R{0},
                       R{50000} * uni20::numeric_limits<R>::epsilon());
  auto bad = phi;
  bad[1] = bad[0];
  EXPECT_THROW(eval(bad), std::domain_error);
}

TYPED_TEST(XXZVertex, SmallComplementSurvivesRoundedMainModulus)
{
  using R = TypeParam;
  bethe::experimental::xxz::SpectralKinematics<R> model(uni20::parse_real<R>("1.005"));
  R const eps = uni20::numeric_limits<R>::epsilon(), kp = model.complementary_modulus();
  ASSERT_GT(kp, R{0});
  EXPECT_EQ(R{1} - kp * kp, R{1}); // even fp128 cannot recover k' by subtraction
  R const pi = bethe::detail::pi<R>();
  EXPECT_EQ(model.jacobian(R{0}), model.jacobian(pi));
  EXPECT_REAL_NEAR(model.jacobian(R{0}) * kp, model.jacobian(pi / R{2}), R{64} * eps);
  for (R phi : {R{.5}, R{1}, R{1.4L}})
  {
    R const p = model.momentum(phi);
    ASSERT_GT(p, R{0});
    EXPECT_REAL_NEAR(model.angle(p), phi, R{200000} * eps);
  }
  if constexpr (std::is_same_v<R, double>) EXPECT_THROW(model.momentum(R{-1.4}), std::underflow_error);
}

TYPED_TEST(XXZVertex, IndependentHundredDigitComplexAmplitudes)
{
  using R = TypeParam;
  std::array<R, 4> const phi{uni20::parse_real<R>("-1.2"), uni20::parse_real<R>("-.2"), uni20::parse_real<R>(".6"),
                             uni20::parse_real<R>("1.3")};
  auto const& refs = amplitude_reference;
  auto const a = bethe::experimental::xxz::VertexEvaluator<R>(R{1.5})(phi);
  for (int v = 0; v < 2; ++v)
    for (int l = 0; l < 4; ++l)
    {
      EXPECT_REAL_NEAR(a.value[v][l].real(), uni20::parse_real<R>(refs[4 * v + l][0]),
                       R{50000} * uni20::numeric_limits<R>::epsilon());
      EXPECT_REAL_NEAR(a.value[v][l].imag(), uni20::parse_real<R>(refs[4 * v + l][1]),
                       R{50000} * uni20::numeric_limits<R>::epsilon());
    }
}

TYPED_TEST(XXZVertex, BroadenedExchangeScalingAndDeterministicReduction)
{
  using R = TypeParam;
  bethe::experimental::xxz::ThermodynamicVertexSpectrum<R> base(R{1.5}), scaled(R{1.5}, R{3});
  std::array<R, 2> const omega{R{2}, R{3}}, triple{R{6}, R{9}};
  auto const a = base.broadened_cut(4, R{1.5}, omega, R{.125}, 4);
  auto const b = scaled.broadened_cut(4, R{1.5}, triple, R{.375}, 4);
  auto const c = base.broadened_cut(4, R{1.5}, omega, R{.125}, 4, [](std::size_t n, auto body) {
    for (auto j = n; j > 0; --j)
      body(j - 1);
  });
  R const tolerance = R{1000000} * uni20::numeric_limits<R>::epsilon();
  for (unsigned i = 0; i < 2; ++i)
  {
    EXPECT_EQ(a.density[i], c.density[i]);
    EXPECT_REAL_NEAR(a.density[i], R{3} * b.density[i], tolerance);
  }
  EXPECT_REAL_NEAR(a.moments.weight, b.moments.weight, tolerance);
  EXPECT_REAL_NEAR(R{3} * a.moments.first_moment, b.moments.first_moment, tolerance);
  EXPECT_EQ(a.moments.nodes, std::size_t{256});
  EXPECT_THROW(base.broadened_cut(4, R{1.5}, omega, R{.1}, 3), std::invalid_argument);
  R const eps = uni20::numeric_limits<R>::epsilon();
  EXPECT_THROW(base.broadened_cut(2, eps * eps, omega, R{.1}, 4), std::underflow_error);
}

TYPED_TEST(XXZVertex, LargeExchangePreservesJacobianAndCut)
{
  using R = TypeParam;
  R const pi = bethe::detail::pi<R>(), eps = uni20::numeric_limits<R>::epsilon();
  R const maximum = uni20::numeric_limits<R>::max();
  bethe::experimental::xxz::SpectralKinematics<R> large(R{2}, maximum / R{4.5});
  ASSERT_TRUE(uni20::isfinite(R{2} * large.maximum_energy()));
  // Independent 100-digit theta-constant evaluation of 1/(2 K) at Delta=2.
  R const midpoint = uni20::parse_real<R>(
      "0.1331393188395339061872077240699224023339178280223141019043160248768945368247581602632319402406991944");
  EXPECT_REAL_NEAR(large.jacobian(pi / R{2}), midpoint, R{256} * eps);
  EXPECT_THROW(large.jacobian(-eps), std::invalid_argument);

  bethe::experimental::xxz::ThermodynamicVertexSpectrum<R> base(R{2});
  R const w = uni20::parse_real<R>("2.3"), width = uni20::parse_real<R>(".1");
  std::array<R, 1> const omega{w};
  auto const expected = base.broadened_cut(2, R{1.5}, omega, width, 8);
  ASSERT_GT(expected.moments.weight, R{0});
  for (int divisor : {8, 10, 12, 16})
  {
    R const exchange = maximum / R(divisor);
    bethe::experimental::xxz::ThermodynamicVertexSpectrum<R> scaled(R{2}, exchange);
    std::array<R, 1> const frequencies{exchange * w};
    auto const actual = scaled.broadened_cut(2, R{1.5}, frequencies, exchange * width, 8);
    // The existing band constructs its gap through exp(log(J)+constant).
    // At extreme J that contributes O(epsilon*abs(log(J))) relative error.
    R const tolerance = R{16} * eps * (R{1} + uni20::math::abs(uni20::math::log(exchange)));
    EXPECT_REAL_NEAR(actual.moments.weight / expected.moments.weight, R{1}, tolerance);
    EXPECT_REAL_NEAR((actual.moments.first_moment / exchange) / expected.moments.first_moment, R{1}, tolerance);
    EXPECT_REAL_NEAR((actual.density[0] * exchange) / expected.density[0], R{1}, tolerance);
    EXPECT_EQ(actual.moments.nodes, expected.moments.nodes);
  }
}

TEST(XXZVertexIntegration, IndependentPythonIntegratedWeight)
{
  bethe::experimental::xxz::ThermodynamicVertexSpectrum<double> sf(1.5, 1., {24, 8});
  auto const a = sf.integrated(4, 12);
  EXPECT_NEAR(a.weight, .015307326980575908, 1e-13);
  EXPECT_NEAR(a.first_moment, .03594894479635669, 1e-12);
  EXPECT_EQ(a.nodes, std::size_t{20736});
  std::array<double, 3> const omega{1.8, 2.2, 2.6};
  auto const two = sf.broadened_cut(2, bethe::detail::pi<double>() / 2, omega, .06, 128);
  std::array<double, 3> const expected{2.140377108167388, 2.197348192146067, .17242886353960657};
  for (unsigned j = 0; j < 3; ++j)
    EXPECT_NEAR(two.density[j], expected[j], 1e-9);
  std::array<double, 4> const grid{1.8, 2.2, 2.6, 3.2};
  auto const four = sf.broadened_cut(4, 1.5, grid, .06, 12);
  std::array<double, 4> const target{6.606952745684729e-06, .047679974335733656, .09266643564633772,
                                     .015852815526838926};
  for (unsigned j = 0; j < 4; ++j)
    EXPECT_NEAR(four.density[j], target[j], 1e-11);
  EXPECT_NEAR(four.moments.weight, .012093538495760878, 1e-12);
  EXPECT_NEAR(four.moments.first_moment, .031926812239310653, 1e-12);
}

TEST(XXZVertexIntegration, ScalarRangeFailuresAreNotSuccessfulNullRows)
{
  bethe::experimental::xxz::ThermodynamicVertexSpectrum<double> sf(1.5);
  std::array<double, 1> const omega{2.};
  EXPECT_THROW(sf.broadened_cut(4, 1.5, omega, uni20::numeric_limits<double>::min() / 4, 4), std::overflow_error);
  bethe::experimental::xxz::ThermodynamicVertexSpectrum<double> huge(1.5, uni20::numeric_limits<double>::max() / 4);
  EXPECT_THROW(huge.broadened_cut(4, 1.5, omega, .1, 4), std::overflow_error);
  EXPECT_THROW(huge.integrated(4, 4), std::overflow_error);
}

TEST(XXZVertexIntegration, NearIsotropicTwoSpinonsMatchIndependentReducedFormula)
{
  // Independent energy-sheet Jacobian versus original coherent amplitudes
  // and the momentum pair map, at the SAME Gaussian width and momentum.
  double const pi = bethe::detail::pi<double>(), k = pi / 2, width = .06;
  std::array<double, 3> const omega{1.8, 2., 2.3};
  auto const rule = bethe::detail::gauss_legendre<double>(256);
  for (double delta : {1.1, 1.025})
  {
    bethe::experimental::xxz::ThermodynamicVertexSpectrum<double> vertex(delta);
    bethe::xxz::ThermodynamicTwoSpinonStructureFactor<double> reduced(delta);
    auto const actual = vertex.broadened_cut(2, k, omega, width, 512);
    std::vector<double> edges;
    for (auto const& band : reduced.sheets(k))
      if (band.present)
      {
        edges.push_back(band.lower);
        edges.push_back(band.upper);
      }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    std::array<double, 3> expected{};
    for (std::size_t panel = 0; panel + 1 < edges.size(); ++panel)
      for (std::size_t i = 0; i < rule.x.size(); ++i)
      {
        double const t = pi * (1 + rule.x[i]) / 4, dw = edges[panel + 1] - edges[panel];
        double const energy = edges[panel] + dw * uni20::math::sin(t) * uni20::math::sin(t);
        auto const density = reduced(k, energy);
        ASSERT_TRUE(density.value);
        double const factor =
            *density.value * dw * uni20::math::sin(2 * t) * pi / 4 * rule.w[i] / (uni20::math::sqrt(2 * pi) * width);
        for (unsigned j = 0; j < omega.size(); ++j)
        {
          double const x = (omega[j] - energy) / width;
          expected[j] += factor * uni20::math::exp(-x * x / 2);
        }
      }
    for (unsigned j = 0; j < omega.size(); ++j)
      EXPECT_NEAR(actual.density[j], expected[j], 1e-7 * expected[j]) << delta << " " << omega[j];
  }
}

TEST(XXZVertexIntegration, IndependentLeadingIsingWeight)
{
  // Leading coefficient from the separate Ising/finite-chain analysis;
  // no local or frequency sum rule is used to normalize the evaluator.
  double const delta = 80;
  bethe::experimental::xxz::ThermodynamicVertexSpectrum<double> sf(delta);
  std::array<double, 1> const omega{2 * delta};
  // Match the independently refined Python pair rule: order 12 deviates
  // from the leading coefficient by 2.8e-4, outside the tolerance below.
  auto const cut = sf.broadened_cut(4, bethe::detail::pi<double>() / 2, omega, .1, 24);
  EXPECT_NEAR(cut.moments.weight * 16 * uni20::math::pow(delta, 4), 1., 1e-4);
  EXPECT_NEAR(cut.moments.first_moment / cut.moments.weight, 2 * delta, .02);
}

#if UNI20_HAS_MPREAL
TEST(XXZVertexMpfr, ExplicitPrecisionAmplitudesGeometryAndQuadrature)
{
  using R = uni20::mpreal;
  namespace model = bethe::experimental::xxz;
  for (auto bits : {128, 256})
  {
    auto const precision = uni20::Precision::bits(bits);
    auto parse = [&](char const* value) { return R(value, precision); };
    R const delta = parse("1.5"), zero = parse("0");
    R const tolerance = R{500000} * uni20::numeric_limits<R>::epsilon(precision);
    std::array<R, 4> phi{parse("-1.2"), parse("-.2"), parse(".6"), parse("1.3")};
    model::VertexEvaluator<R> evaluator(delta);
    auto const a = evaluator(phi);
    auto controls = evaluator.options();
    controls.products += 8;
    controls.residues += 4;
    auto const b = model::VertexEvaluator<R>(delta, controls)(phi);
    for (unsigned vacuum = 0; vacuum < 2; ++vacuum)
      for (unsigned polarization = 0; polarization < 4; ++polarization)
      {
        EXPECT_EQ(uni20::precision_of(a.value[vacuum][polarization]), precision);
        EXPECT_EQ(uni20::precision_of(a.arithmetic_indicator[vacuum][polarization]), precision);
        EXPECT_REAL_NEAR(uni20::math::abs(a.value[vacuum][polarization] - b.value[vacuum][polarization]), zero,
                         tolerance);
      }
    for (unsigned vacuum = 0; vacuum < 2; ++vacuum)
      for (unsigned polarization = 0; polarization < 4; ++polarization)
      {
        auto const& reference = amplitude_reference[4 * vacuum + polarization];
        EXPECT_REAL_NEAR(a.value[vacuum][polarization].real(), parse(reference[0]), tolerance);
        EXPECT_REAL_NEAR(a.value[vacuum][polarization].imag(), parse(reference[1]), tolerance);
      }
    model::SpectralKinematics<R> geometry(delta);
    EXPECT_REAL_NEAR(geometry.angle(geometry.momentum(phi[0])), phi[0], tolerance);
    EXPECT_EQ(uni20::precision_of(geometry.momentum(zero)), precision);
    EXPECT_EQ(uni20::precision_of(geometry.angle(bethe::detail::spinon_pi(zero) / R{2})), precision);
    model::ThermodynamicVertexSpectrum<R> spectrum(delta);
    auto const moments = spectrum.integrated(2, 4);
    EXPECT_EQ(moments.nodes, 16u);
    EXPECT_EQ(uni20::precision_of(moments.weight), precision);
    EXPECT_GT(moments.weight, zero);
    std::array<R, 1> omega{parse("2.2")};
    auto const cut = spectrum.broadened_cut(4, parse("1.5"), omega, parse(".125"), 2);
    EXPECT_EQ(uni20::precision_of(cut.density[0]), precision);
    EXPECT_GT(cut.moments.weight, zero);
    R const eps = uni20::numeric_limits<R>::epsilon(precision);
    EXPECT_THROW(spectrum.broadened_cut(2, eps * eps, omega, parse(".125"), 4), std::underflow_error);
    auto bad = phi;
    bad[0] = R("-1.2", uni20::Precision::bits(bits + 16));
    EXPECT_THROW(evaluator(bad), std::invalid_argument);
    EXPECT_THROW(geometry.momentum(bad[0]), std::invalid_argument);
    EXPECT_THROW(spectrum.broadened_cut(2, bad[0], omega, parse(".125"), 2), std::invalid_argument);
  }
}
#endif
