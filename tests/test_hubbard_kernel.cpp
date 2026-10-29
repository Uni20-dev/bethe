// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/detail/hubbard_kernel.hpp>
#include <bethe/observable_grid.hpp>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
struct Reference
{
    char const *u, *x, *kernel, *primitive;
};
// Independent mpmath 1.2.1 references, evaluated at 2048 bits (4096 for the last case), rounded
// to 140 significant decimal digits. With a=1+i*x/(2*u), b=x/(4*u):
// R = Re(digamma((a+1)/2)-digamma(a/2))/(4*pi*u),
// primitive = Im(loggamma(1+i*b)-loggamma(1/2+i*b))/pi.
// The real parts of the two digammas are subtracted only in this much higher
// precision oracle; the evaluated kernel uses the independent Euler series.
Reference const references[] = {
    {"1", "0",
     "0."
     "11031780007632579669822821605899884549134487436482773524463161685356088352334334124899655592698377794432415837089"
     "563016425375308391114052987",
     "0.0"},
    {"1", "0.5",
     "0."
     "10191806780841470105938227172551975601941710591006564317193565778518776249475860274048563745621920401530199456538"
     "039995954856078481021634187",
     "0."
     "05372207470990316480235901696017678489740151848575007470722230610109615666560836763033919778257636677332213435627"
     "2674707209288391846965760953"},
    {"1", "4",
     "0."
     "01165109450532912334773912305895696470276902368165092031445645226265826756482746592473302741669082899930966163799"
     "4438503342612805488961695088",
     "0."
     "20797330396366147321047668089876988468623211226602446531763713213649014208997147598278804534740724263812457053798"
     "552741577464609269784470814"},
    {"0.0625", "2",
     "0."
     "00249169160345381204657485098246890837709411859452175264296823996067639693389273737775145206032512699688978869316"
     "08290132370887853725564939521",
     "0."
     "24502315465044781119134400495527301838446630965624268827604883446410936914543483851067962588491399470754687407367"
     "751434627550127847435737936"},
    {"0.0001", "2",
     "0."
     "00000397887359719175167859638569537560585946105256568109702640453088016037112909670673729075299944789232506632495"
     "32032991904326847539354539878213",
     "0."
     "24999204225283214232112807625259391568122037055296221445784297693713374108640951951077026498848874892889546378646"
     "325555973759301107207121623"},
    {"1",
     "25822498780869085896559191720030118743297057928292235128306593565406476220168411946296453532801378314359031719727"
     "47493376",
     "2."
     "38684161070655894827501967680141115722561277698332825727080597390141514872257011096535163842162004118419359274077"
     "38775650900672073289380909e-242",
     "0."
     "24999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999"
     "999999993836578541740227624"},
    {"1",
     "3."
     "87259191484931827281803063328635184757021919204879086548776294134441634809768596486268223427701459690805754250755"
     "446753937083639899235031552231805065e-121",
     "0."
     "11031780007632579669822821605899884549134487436482773524463161685356088352334334124899655592698377794432415837089"
     "563016425375308391114052987",
     "4."
     "27215820639542786545076321256680315569007964200892070448899122470128147766389044430382042283166065381716172393134"
     "06453620542594705483749057e-122"},
    {"9."
     "33263618503218878990089544723817169617091446371708024621714339795966910975775634454440327097881102359594989930324"
     "24262421548752135403239484152081720393075623441066613832515027399507598590183151110049079626511311824051251479593"
     "37908051782711254151038106983788544264811194698142286609592220176629104427984561694488871474665280063283684526474"
     "2926182986e-302",
     "1",
     "1."
     "48533518092616119585428339731939273512461780437422358465723648189293647406936548218862708874588962106969104426112"
     "73298214329914824137756867e-302",
     "0.25"},
};

template <typename Real> class HubbardKernelNative : public ::testing::Test {};
TYPED_TEST_SUITE(HubbardKernelNative, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(HubbardKernelNative, IndependentCoreAndTailReferences)
{
  using Real = TypeParam;
  for (auto const& ref : references)
  {
    SCOPED_TRACE(std::string("u=") + ref.u + ", x=" + ref.x);
    Real const u = uni20::parse_real<Real>(ref.u), x = uni20::parse_real<Real>(ref.x);
    bethe::detail::HubbardKernel<Real> kernel(u);
    Real const expected = uni20::parse_real<Real>(ref.kernel), phi = uni20::parse_real<Real>(ref.primitive);
    Real const epsilon = uni20::numeric_limits<Real>::epsilon();
    EXPECT_REAL_NEAR(kernel.value(x), expected, Real{128} * epsilon * expected);
    if (phi == Real{0})
      EXPECT_EQ(kernel.primitive(x), phi);
    else
      EXPECT_REAL_NEAR(kernel.primitive(x), phi, Real{128} * epsilon * std::abs(phi));
    EXPECT_EQ(kernel.value(-x), kernel.value(x));
    EXPECT_EQ(kernel.primitive(-x), -kernel.primitive(x));
  }
}

TYPED_TEST(HubbardKernelNative, BudgetsAndDomainFailuresAreExplicit)
{
  using Real = TypeParam;
  using Kernel = bethe::detail::HubbardKernel<Real>;
  auto const p = uni20::native_precision<Real>();
  EXPECT_THROW(Kernel(Real{1}, uni20::Precision::exact()), std::invalid_argument);
  EXPECT_THROW(Kernel(Real{0}, p), std::invalid_argument);
  EXPECT_THROW(Kernel(Real{-1}, p), std::invalid_argument);
  Kernel const short_series(Real{1}, p, 1);
  EXPECT_THROW(short_series.value(Real{.5}), std::runtime_error);
  Kernel const kernel(Real{1}, p);
  EXPECT_THROW(kernel.value(uni20::numeric_limits<Real>::infinity()), std::invalid_argument);
  EXPECT_THROW(kernel.primitive(uni20::numeric_limits<Real>::quiet_NaN()), std::invalid_argument);
}

#if UNI20_HAS_MPREAL
using R = uni20::mpreal;
using uni20::Precision;
namespace math = uni20::math;

TEST(HubbardKernelMpfr, AccuracyGrowsBeyondBinary128IncludingRapidityTails)
{
  auto const ambient = mpfr_get_default_prec();
  auto const refp = Precision::bits(768);
  R previous(1, refp);
  for (int bits : {80, 160, 256, 333})
  {
    auto const p = Precision::bits(bits);
    R largest_error(0, refp);
    for (auto const& ref : references)
    {
      R const u(ref.u, p), x(ref.x, p), expected(ref.kernel, refp), phi(ref.primitive, refp);
      bethe::detail::HubbardKernel<R> const kernel(u, p);
      R const value = kernel.value(x), primitive = kernel.primitive(x);
      for (auto const* v : {&value, &primitive})
      {
        EXPECT_EQ(uni20::precision_of(*v), p);
        EXPECT_TRUE(uni20::isfinite(*v));
      }
      R const error = math::abs(uni20::scalar_at_precision<R>(value, refp) - expected) / expected;
      R const phase_error =
          math::abs(uni20::scalar_at_precision<R>(primitive, refp) - phi) / (phi == 0 ? R(1, refp) : math::abs(phi));
      largest_error = std::max({largest_error, error, phase_error});
      EXPECT_EQ(kernel.value(-x), value);
      EXPECT_EQ(kernel.primitive(-x), -primitive);
    }
    EXPECT_LT(largest_error, 128 * uni20::scalar_at_precision<R>(uni20::numeric_limits<R>::epsilon(p), refp));
    EXPECT_LT(largest_error, previous * R("1e-15", refp));
    if (bits >= 256) EXPECT_LT(largest_error, R("1e-70", refp));
    previous = largest_error;
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

TEST(HubbardKernelMpfr, RejectsMixedExactAndUnsetValues)
{
  auto const p = Precision::bits(160), other = Precision::bits(80);
  using Kernel = bethe::detail::HubbardKernel<R>;
  EXPECT_THROW(Kernel(R(1), p), std::invalid_argument);
  EXPECT_THROW(Kernel(R{}, p), std::invalid_argument);
  EXPECT_THROW(Kernel(R(1, other), p), std::invalid_argument);
  EXPECT_THROW(Kernel(R(1, p), Precision::exact()), std::invalid_argument);
  Kernel const kernel(R(1, p), p);
  for (R const& x : {R{}, R(1), R(1, other)})
  {
    EXPECT_THROW(kernel.value(x), std::invalid_argument);
    EXPECT_THROW(kernel.primitive(x), std::invalid_argument);
  }
}

TEST(HubbardKernelMpfr, SharedContextsRemainOwnedAcrossMixedPrecisionQueries)
{
  using Kernel = bethe::detail::HubbardKernel<R>;
  using Result = std::array<R, 2>;
  auto const ambient = mpfr_get_default_prec();
  for (int workers : {1, 2})
  {
    uni20::async::TbbScheduler scheduler(workers);
    std::array<std::shared_ptr<Kernel const>, 2> kernels;
    for (int j = 0; j < 2; ++j)
    {
      auto const p = Precision::bits(j == 0 ? 160 : 256);
      kernels[j] = std::make_shared<Kernel const>(R(1, p), p);
    }
    std::vector<Result> expected, actual;
    auto evaluate = [owned = kernels](std::size_t i) {
      auto const& kernel = *owned[i % 2];
      R const x = uni20::scalar_like(kernel.u(), i) / 8;
      return Result{kernel.value(x), kernel.primitive(x)};
    };
    for (std::size_t i = 0; i < 32; ++i)
      expected.push_back(evaluate(i));
    auto done = bethe::evaluate_observable_grid_async(
        expected.size(), std::move(evaluate),
        [&actual](std::size_t, Result value) { actual.push_back(std::move(value)); },
        {.scheduler = &scheduler, .max_in_flight = 4});
    kernels = {};
    (void)done.get_wait(scheduler);
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
      for (std::size_t j = 0; j < 2; ++j)
      {
        EXPECT_EQ(actual[i][j], expected[i][j]);
        EXPECT_EQ(uni20::precision_of(actual[i][j]), uni20::precision_of(expected[i][j]));
      }
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
#endif
} // namespace
