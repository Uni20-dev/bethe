// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/kondo.hpp>

namespace
{
template <typename Real> class Kondo : public ::testing::Test {};
TYPED_TEST_SUITE(Kondo, test_support::RealTypes, test_support::PrecisionNames);
TYPED_TEST(Kondo, IndependentReferences)
{
  using R = TypeParam;
  char const* refs[][3] = {
      {"0.1", "0.024120436714010866605544353146795416662520524221635417113550217572",
       "-0.0012079339556471688182054950998824299737601134567553014324079266758"},
      {"0.5", "0.11256604809692454792435668748526832773955993602502130129993816713",
       "-0.029147209883951833999521239464166941483327709815637216470321338254"},
      {"1", "0.19211164621093772781603846463983704837166745414312923905831121437",
       "-0.1067922423030353039230318678651420002441806867445409557659118577"},
      {"2", "0.27516363791067304687136724858153892615232627968484052219868642339",
       "-0.34629170027457867780906807893015178590470848543255021208916271878"},
      {"10", "0.38717932588109120947017766086536672158479856504446913615395010656",
       "-3.1798714372504321779665130965884524646973165211039948891430143624"},
      {"1.0001", "0.192124191899217789286549028908695952244968454", "-0.106811454094950626884725696525950107827116539"},
      {"1.01", "0.19336040352583326695234714279489446176410559", "-0.108719612319901181616841221803523583369859694"},
      {"1000000", "0.480548715930107962754613692984943018800281082",
       "-478804.124266610744890329181578815499893823977"}};
  for (auto const& ref : refs)
  {
    R const x = uni20::parse_real<R>(ref[0]);
    auto const s = bethe::kondo::ground_response(x, R{1});
    ASSERT_TRUE(s.converged) << ref[0] << " status " << int(s.status) << " lobes " << s.lobes << " evaluations "
                             << s.evaluations;
    R const tol = R{1048576} * uni20::numeric_limits<R>::epsilon();
    EXPECT_REAL_NEAR(*s.magnetization, uni20::parse_real<R>(ref[1]), tol);
    EXPECT_REAL_NEAR(*s.energy_change / x, uni20::parse_real<R>(ref[2]) / x, tol);
  }
}

TYPED_TEST(Kondo, ScaledGammaAgainstIndependentReferences)
{
  using R = TypeParam;
  char const* refs[][2] = {{"0.0001", "1.77391621018113358598956353900367318312051766294474770561283832938028189988"},
                           {"0.5", "2.33164398159712420336353606216840087638023629918758842300809644777601004941"},
                           {"1", "2.40901454734936102856076554562305940710651285559929926596630757129403823603"},
                           {"31.5", "2.50331501528966533704173696404450065640816112398761319492157427743204663599"},
                           {"32", "2.50336674533396979322532980037881789846154921353956920538722746446285382317"},
                           {"100", "2.50558406983204859993367118803367357618649823398872844965198108302382164123"},
                           {"1000", "2.5065238339681788333211283947765620246757417454316490213351903539708593187"}};
  for (auto const& ref : refs)
  {
    auto const a = bethe::kondo::detail::amplitude(uni20::parse_real<R>(ref[0]));
    EXPECT_REAL_NEAR(a.value, uni20::parse_real<R>(ref[1]), a.error);
  }
}

TYPED_TEST(Kondo, CrossoverAndScaleCovariance)
{
  using R = TypeParam;
  namespace model = bethe::kondo;
  model::Options<R> options;
  R const tol = options.tolerance;
  model::State<R> work;
  auto const anchor = model::detail::low(R{1}, options, work);
  ASSERT_TRUE(anchor);
  auto const high = model::detail::high(R{0}, *anchor, options, work);
  ASSERT_TRUE(high);
  EXPECT_REAL_NEAR(high->value[0], anchor->value[0], tol);
  EXPECT_REAL_NEAR(high->value[1], -anchor->value[1], tol);
  R const eps = uni20::numeric_limits<R>::epsilon();
  for (R x : {R{0.25}, R{1} - R{16} * eps, R{1} + R{16} * eps, R{3}})
  {
    auto const positive = model::ground_response(x, R{1});
    auto const negative = model::ground_response(-x * R{4}, R{4});
    ASSERT_TRUE(positive.converged && negative.converged) << uni20::format_scalar(x);
    EXPECT_REAL_NEAR(*positive.magnetization, -*negative.magnetization, tol);
    EXPECT_REAL_NEAR(*positive.energy_change, *negative.energy_change / R{4}, tol * std::max(x, R{1}));
    EXPECT_REAL_NEAR(*positive.zero_field_susceptibility, *negative.zero_field_susceptibility * R{4}, R{32} * eps);
  }
}

TYPED_TEST(Kondo, EnergyDerivativeAndPhysicalLimits)
{
  using R = TypeParam;
  namespace model = bethe::kondo;
  R const eps = uni20::numeric_limits<R>::epsilon(), h = std::pow(eps, R{1} / R{5});
  for (R x : {R{0.25}, R{1}, R{2}})
  {
    auto const a = model::ground_response(x - R{2} * h, R{1});
    auto const b = model::ground_response(x - h, R{1});
    auto const c = model::ground_response(x + h, R{1});
    auto const d = model::ground_response(x + R{2} * h, R{1});
    auto const center = model::ground_response(x, R{1});
    ASSERT_TRUE(a.converged && b.converged && c.converged && d.converged && center.converged);
    R const derivative =
        (*a.energy_change - R{8} * *b.energy_change + R{8} * *c.energy_change - *d.energy_change) / (R{12} * h);
    EXPECT_REAL_NEAR(-derivative, *center.magnetization, R{16777216} * eps / h);
  }
  R const small = std::sqrt(eps);
  auto const s = model::ground_response(small, R{1});
  ASSERT_TRUE(s.converged);
  EXPECT_REAL_NEAR(*s.magnetization / small, *s.zero_field_susceptibility, R{1048576} * eps);
  EXPECT_REAL_NEAR(*s.energy_change / (small * small), -*s.zero_field_susceptibility / R{2}, R{1048576} * eps);
  R previous{};
  for (R x : {R{10}, R{1000}, R{1000000}})
  {
    auto const state = model::ground_response(x, R{1});
    ASSERT_TRUE(state.converged);
    EXPECT_GT(*state.magnetization, previous);
    EXPECT_LT(*state.magnetization, R{0.5});
    EXPECT_GT(*state.energy_change / x, R{-0.5});
    previous = *state.magnetization;
  }
  // Ratio overflows, but all three requested physical observables fit.
  R const field = uni20::numeric_limits<R>::max() / R{16}, scale = uni20::numeric_limits<R>::min() * R{16};
  auto const huge = model::ground_response(field, scale);
  ASSERT_TRUE(huge.converged) << int(huge.status) << " evaluations " << huge.evaluations;
  EXPECT_GT(*huge.magnetization, R{0.49});
  EXPECT_LT(*huge.energy_change / field, R{-0.49});
}

TYPED_TEST(Kondo, BudgetsAndInvalidInputs)
{
  using R = TypeParam;
  namespace model = bethe::kondo;
  auto absent = [](auto const& s) {
    EXPECT_FALSE(s.converged);
    EXPECT_FALSE(s.magnetization);
    EXPECT_FALSE(s.energy_change);
    EXPECT_FALSE(s.zero_field_susceptibility);
  };
  model::Options<R> options;
  options.max_series_terms = 0;
  auto zero = model::ground_response(R{0}, R{1}, options);
  ASSERT_TRUE(zero.converged);
  EXPECT_EQ(*zero.magnetization, R{0});
  EXPECT_EQ(*zero.energy_change, R{0});
  auto series = model::ground_response(R{1}, R{1}, options);
  absent(series);
  EXPECT_EQ(series.status, model::Status::series_limit);
  options.max_series_terms = 16;
  series = model::ground_response(R{1}, R{1}, options);
  absent(series);
  EXPECT_EQ(series.status, model::Status::series_limit);
  EXPECT_EQ(series.series_terms, 16u);
  options = {};
  options.max_lobes = 0;
  auto tail = model::ground_response(R{2}, R{1}, options);
  absent(tail);
  EXPECT_EQ(tail.status, model::Status::tail_limit);
  options.max_lobes = 16;
  tail = model::ground_response(R{2}, R{1}, options);
  absent(tail);
  EXPECT_EQ(tail.status, model::Status::tail_limit);
  EXPECT_EQ(tail.lobes, 16u);
  options = {};
  options.max_evaluations = 0;
  auto quadrature = model::ground_response(R{2}, R{1}, options);
  absent(quadrature);
  EXPECT_EQ(quadrature.status, model::Status::quadrature_limit);
  options = {};
  options.max_quadrature_levels = 0;
  quadrature = model::ground_response(R{2}, R{1}, options);
  absent(quadrature);
  EXPECT_EQ(quadrature.status, model::Status::quadrature_limit);
  auto range = model::ground_response(uni20::numeric_limits<R>::min(), uni20::numeric_limits<R>::max());
  absent(range);
  EXPECT_EQ(range.status, model::Status::precision_limit);
  EXPECT_THROW((void)model::ground_response(R{0}, R{0}), std::invalid_argument);
  EXPECT_THROW((void)model::ground_response(uni20::numeric_limits<R>::quiet_NaN(), R{1}), std::invalid_argument);
  options = {};
  options.tolerance = R{0};
  EXPECT_THROW((void)model::ground_response(R{0}, R{1}, options), std::invalid_argument);
  options = {};
  options.max_lobes = 4097;
  EXPECT_THROW((void)model::ground_response(R{0}, R{1}, options), std::length_error);
}
} // namespace
