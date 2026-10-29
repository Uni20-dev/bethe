// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/kondo.hpp>
#include <bethe/observable_grid.hpp>
#include <tuple>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using R = uni20::mpreal;
using uni20::Precision;
namespace model = bethe::kondo;
namespace math = uni20::math;

void check_precision(model::State<R> const& state, Precision p)
{
  for (auto const* value : {&state.field, &state.scale, &state.magnetization_error, &state.scaled_energy_error})
    EXPECT_EQ(uni20::precision_of(*value), p);
  for (auto const* value : {&state.magnetization, &state.energy_change, &state.zero_field_susceptibility})
    if (*value)
    {
      EXPECT_EQ(uni20::precision_of(**value), p);
      EXPECT_TRUE(uni20::isfinite(**value));
    }
}
auto fields(model::State<R> const& s)
{
  return std::tie(s.field, s.scale, s.magnetization, s.energy_change, s.zero_field_susceptibility,
                  s.magnetization_error, s.scaled_energy_error, s.series_terms, s.lobes, s.evaluations, s.converged,
                  s.status);
}

TEST(KondoMpfr, IndependentSeriesAndIntegralGainAccuracy)
{
  // Independent mpmath nsum and quadosc at 105/130 decimal digits, with
  // direct Gamma rather than production log-Gamma/lobe tanh-sinh quadrature.
  char const* refs[][3] = {
      {"1",
       "0."
       "192111646210937727816038464639837048371667454143129239058311214369588527019054275778408325279344171282415032934"
       "79403595707395",
       "-0."
       "106792242303035303923031867865142000244180686744540955765911857704475886963720209710877374655759466399969423317"
       "7803816046091"},
      {"2",
       "0."
       "275163637910673046871367248581538926152326279684840522198686423385892399422730234552560359629415304183013842214"
       "5228546957692",
       "-0."
       "346291700274578677809068078930151785904708485432550212089162718779733803188613168667705536779580074463737734032"
       "95083008090387"},
  };
  auto oracle = Precision::bits(768);
  auto ambient = mpfr_get_default_prec();
  for (auto const& ref : refs)
  {
    R previous(1, oracle);
    for (int bits : {80, 160, 256, 333})
    {
      SCOPED_TRACE(ref[0]);
      SCOPED_TRACE(bits);
      auto p = Precision::bits(bits);
      model::Options<R> options(p);
      auto state = model::ground_response(R(ref[0], p), R(1, p), options);
      ASSERT_TRUE(state.converged) << int(state.status) << " terms " << state.series_terms << " lobes " << state.lobes
                                   << " evaluations " << state.evaluations;
      check_precision(state, p);
      R em = math::abs(state.magnetization->at(oracle) - R(ref[1], oracle));
      R ee = math::abs(state.energy_change->at(oracle) - R(ref[2], oracle)) / R(ref[0], oracle);
      EXPECT_LE(em, state.magnetization_error.at(oracle));
      EXPECT_LE(ee, state.scaled_energy_error.at(oracle));
      R error = std::max(em, ee);
      EXPECT_LT(error, options.tolerance.at(oracle));
      EXPECT_LT(error, previous * R("1e-14", oracle));
      if (bits == 333) EXPECT_LT(error, R("1e-92", oracle));
      previous = error;
    }
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

TEST(KondoMpfr, ScaledGammaHasNoFixedCoefficientAccuracyCeiling)
{
  char const* refs[][2] = {
      {"0", "1."
            "7724538509055160272981674833411451827975494561223871282138077898529112845910321813749506567385446654162268"
            "236242825706662362"},
      {"0.0001", "1."
                 "77391621018113358598956353900367318312051766294474770561283832938028189988417891593844653536617195048"
                 "36933354052427816916044"},
      {"0.5", "2."
              "33164398159712420336353606216840087638023629918758842300809644777601004941265734950262999179547775424833"
              "50263120733624483664"},
      {"1", "2."
            "4090145473493610285607655456230594071065128555992992659663075712940382360271737691653362282019443302033909"
            "073783245358704805"},
      {"32", "2."
             "503366745333969793225329800378817898461549213539569205387227464462853823166911336317954009298529361904564"
             "1715044971676113317"},
      {"4096", "2."
               "5066027760193815629020357359078366715086839782657902203585588647635684596574577156419379603397726599996"
               "687484212543221781969"},
  };
  auto p = Precision::bits(333), oracle = Precision::bits(768);
  for (auto const& ref : refs)
  {
    auto a = model::detail::amplitude(R(ref[0], p));
    EXPECT_EQ(uni20::precision_of(a.value), p);
    EXPECT_EQ(uni20::precision_of(a.error), p);
    auto error = math::abs(a.value.at(oracle) - R(ref[1], oracle));
    EXPECT_LE(error, a.error.at(oracle));
    EXPECT_LT(error, R("1e-93", oracle));
  }
}

TEST(KondoMpfr, ZerosScalingAndExtremeExponentsRetainPrecision)
{
  auto p = Precision::bits(160);
  model::Options<R> options(p);
  auto zero = model::ground_response(R(0, p), R(1, p), options);
  ASSERT_TRUE(zero.converged);
  check_precision(zero, p);
  EXPECT_EQ(*zero.magnetization, 0);
  EXPECT_EQ(*zero.energy_change, 0);
  auto state = model::ground_response(R(1, p), R(1, p), options);
  ASSERT_TRUE(state.converged);
  for (char const* scale : {"2", "1e1000", "1e-1000"})
  {
    R factor(scale, p);
    auto scaled = model::ground_response(-factor, factor, options);
    ASSERT_TRUE(scaled.converged);
    check_precision(scaled, p);
    EXPECT_EQ(*scaled.magnetization, -*state.magnetization);
    EXPECT_LT(math::abs(*scaled.energy_change / factor - *state.energy_change),
              16 * uni20::numeric_limits<R>::epsilon(p));
    EXPECT_LT(math::abs(*scaled.zero_field_susceptibility * factor - *state.zero_field_susceptibility),
              16 * uni20::numeric_limits<R>::epsilon(p));
  }
  auto tiny = model::ground_response(R("1e-1000", p), R(1, p), options);
  ASSERT_TRUE(tiny.converged);
  check_precision(tiny, p);
  EXPECT_GT(*tiny.magnetization, 0);
  EXPECT_LT(*tiny.energy_change, 0);
  EXPECT_LT(math::abs(*tiny.magnetization / R("1e-1000", p) - *tiny.zero_field_susceptibility), options.tolerance);
}

TEST(KondoMpfr, InvalidInputsAndBudgetFailuresStayExplicit)
{
  auto p = Precision::bits(160);
  model::Options<R> options(p);
  EXPECT_THROW(model::Options<R>(Precision::exact()), std::invalid_argument);
  EXPECT_THROW(model::ground_response(R(1), R(1, p), options), std::invalid_argument);
  EXPECT_THROW(model::ground_response(R(1, Precision::bits(80)), R(1, p), options), std::invalid_argument);
  EXPECT_THROW(model::ground_response(R(1, p), R(0, p), options), std::invalid_argument);
  options.tolerance = R(1, Precision::bits(80));
  EXPECT_THROW(model::ground_response(R(1, p), R(1, p), options), std::invalid_argument);
  for (int failure = 0; failure < 3; ++failure)
  {
    options = model::Options<R>(p);
    if (failure == 0) options.max_series_terms = 0;
    if (failure == 1) options.max_lobes = 0;
    if (failure == 2) options.max_evaluations = 0;
    auto state = model::ground_response(R(2, p), R(1, p), options);
    EXPECT_FALSE(state.converged);
    EXPECT_FALSE(state.magnetization);
    EXPECT_FALSE(state.energy_change);
    EXPECT_FALSE(state.zero_field_susceptibility);
    check_precision(state, p);
    EXPECT_EQ(state.status, failure == 0   ? model::Status::series_limit
                            : failure == 1 ? model::Status::tail_limit
                                           : model::Status::quadrature_limit);
  }
}

TEST(KondoMpfr, OwningMixedPrecisionAsyncMatchesCompleteStates)
{
  using namespace uni20::async;
  struct Request
  {
      R field, scale;
      model::Options<R> options;
  };
  std::vector<Request> requests;
  for (int bits : {80, 160, 256})
  {
    auto p = Precision::bits(bits);
    model::Options<R> options(p);
    options.tolerance = R("1e-12", p);
    for (int field : {-1, 0, 2})
      requests.push_back({R(field, p), R(1, p), options});
    options.max_lobes = 0;
    requests.push_back({R(2, p), R(1, p), options});
  }
  std::vector<model::State<R>> reference;
  for (auto const& r : requests)
    reference.push_back(model::ground_response(r.field, r.scale, r.options));
  auto check = [&](IAsyncScheduler& scheduler) {
    std::vector<model::State<R>> actual;
    auto completion = [&] {
      auto owned = requests;
      return bethe::evaluate_observable_grid_async(
          owned.size(),
          [owned](std::size_t i) {
            auto const& r = owned[i];
            return model::ground_response(r.field, r.scale, r.options);
          },
          [&](std::size_t i, model::State<R> state) {
            EXPECT_EQ(i, actual.size());
            actual.push_back(std::move(state));
          },
          {.scheduler = &scheduler, .max_in_flight = 2});
    }();
    EXPECT_EQ(completion.get_wait(scheduler), reference.size());
    ASSERT_EQ(actual.size(), reference.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
      EXPECT_TRUE(fields(actual[i]) == fields(reference[i]));
      check_precision(actual[i], requests[i].options.precision());
    }
  };
  DebugScheduler deferred({.order = DebugSchedulerOrder::reverse});
  check(deferred);
  for (int workers : {1, 2})
  {
    TbbScheduler scheduler(workers);
    check(scheduler);
  }
}
} // namespace
