// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "spinon-thermo.hpp"
#include <gtest/gtest.h>

namespace cli = bethe::cli;

TEST(CliPrecision, SpinonNativeComparisonUsesDistinctIntermediatePrecision)
{
  namespace spinon = cli::spinon;
  EXPECT_TRUE((std::is_same_v<spinon::DensityLow<double>, float>));
  constexpr int native_bits = std::numeric_limits<long double>::digits;
  using NativeLow = std::conditional_t<native_bits == 53, float, double>;
  EXPECT_TRUE((std::is_same_v<spinon::DensityLow<long double>, NativeLow>));
  auto const native = spinon::recovery(spinon::Controls{}, 1.L);
  EXPECT_EQ(native.stages.front().bits, native_bits == 53 ? 24 : 53);
  EXPECT_EQ(native.stages.front().guard_bits, native_bits);
  auto const binary64 = spinon::recovery(spinon::Controls{}, 1.0);
  EXPECT_EQ(binary64.stages.front().bits, 24);
  EXPECT_EQ(binary64.stages.front().guard_bits, 53);
#if UNI20_HAS_FLOAT128
  auto const binary128 = spinon::recovery(spinon::Controls{}, uni20::float128{1});
  EXPECT_EQ(binary128.stages.front().bits, native_bits > 53 && native_bits < 113 ? native_bits : 53);
  EXPECT_EQ(binary128.stages.front().guard_bits, 113);
  EXPECT_LT(binary128.stages.front().bits, binary128.stages.front().guard_bits);
#endif
#if UNI20_HAS_MPREAL
  EXPECT_TRUE((std::is_same_v<spinon::DensityLow<uni20::mpreal>, uni20::mpreal>));
#endif
}

TEST(CliPrecision, SampleTaskSizeIsIndependentOfTheSamplingBudget)
{
  cli::spinon::Controls controls;
  CLI::App app;
  cli::spinon::add_options(app, controls);
  app.parse("--samples-per-task 3 --samples 29 --replicates 5");
  EXPECT_EQ(controls.samples_per_task, 3u);
  EXPECT_EQ(controls.sampling.samples_per_replicate, 29u);
  EXPECT_EQ(controls.sampling.replicates, 5u);
  EXPECT_THROW(app.parse("--samples-per-task 0"), CLI::ValidationError);
}

TEST(CliPrecision, SixStageListIsExplicitAndOrdered)
{
  using Stage = cli::spinon::PrecisionStage;
  auto parse = [](std::string text) { return cli::spinon::precision_stages(std::move(text)); };
  EXPECT_TRUE(cli::spinon::precision_stages(std::nullopt).empty());
  EXPECT_EQ(parse("fp80,fp128,reference"), (std::vector{Stage::fp80, Stage::fp128, Stage::reference}));
  EXPECT_EQ(parse("fp128,reference"), (std::vector{Stage::fp128, Stage::reference}));
  EXPECT_EQ(parse("fp80,reference"), (std::vector{Stage::fp80, Stage::reference}));
  EXPECT_EQ(parse("reference"), (std::vector{Stage::reference}));
  for (auto text : {"", "fp64,reference", "long-double,reference", "fp80", "fp128", "fp80,fp128", "reference,fp128",
                    "fp128,fp80,reference", "fp80,fp80,reference", "reference,reference", ",reference",
                    "fp80,,reference", "fp128,reference,"})
  {
    SCOPED_TRACE(text);
    EXPECT_THROW(parse(text), std::invalid_argument);
  }
}

TEST(CliPrecision, SixStageValidationKeepsFixedPrecisionAndRequiresCapableReference)
{
  namespace spinon = cli::spinon;
  auto const reference = spinon::precision_stages("reference");
  auto const native = spinon::precision_stages("fp80,reference");
  EXPECT_NO_THROW(spinon::validate_precision_stages({}, cli::parse_precision("fp64"), false, std::nullopt));
  EXPECT_THROW(spinon::validate_precision_stages(reference, cli::parse_precision("fp64"), true, std::nullopt),
               std::invalid_argument);
  EXPECT_THROW(spinon::validate_precision_stages(reference, cli::parse_precision("bits-192"), false, std::nullopt),
               std::invalid_argument);
#if UNI20_HAS_MPREAL
  EXPECT_NO_THROW(spinon::validate_precision_stages(reference, cli::parse_precision("digits-100"), true, 192));
  for (auto guard : {24u, 193u})
    EXPECT_THROW(spinon::validate_precision_stages(reference, cli::parse_precision("bits-192"), true, guard),
                 std::invalid_argument);
  if (cli::has_fp80)
  {
    EXPECT_NO_THROW(spinon::validate_precision_stages(native, cli::parse_precision("bits-192"), true, 128));
    EXPECT_THROW(spinon::validate_precision_stages(native, cli::parse_precision("bits-192"), true, 63),
                 std::invalid_argument);
  }
  else
    EXPECT_THROW(spinon::validate_precision_stages(native, cli::parse_precision("bits-192"), true, 128),
                 std::invalid_argument);
  auto const wide = spinon::precision_stages("fp128,reference");
#if UNI20_HAS_FLOAT128
  EXPECT_NO_THROW(spinon::validate_precision_stages(wide, cli::parse_precision("bits-192"), true, 128));
  EXPECT_THROW(spinon::validate_precision_stages(wide, cli::parse_precision("bits-192"), true, 112),
               std::invalid_argument);
#else
  EXPECT_THROW(spinon::validate_precision_stages(wide, cli::parse_precision("bits-192"), true, 128),
               std::invalid_argument);
#endif
#else
  EXPECT_THROW(spinon::validate_precision_stages(reference, cli::parse_precision("bits-192"), true, std::nullopt),
               std::invalid_argument);
#endif
}

TEST(CliPrecision, SixStageOptionBelongsToTheSpectrumFrontend)
{
  cli::spinon::Controls controls;
  CLI::App spectrum;
  cli::spinon::add_options(spectrum, controls);
  spectrum.parse("--six-precision-stages fp80,fp128,reference");
  EXPECT_EQ(controls.precision_stages, "fp80,fp128,reference");
  CLI::App moments;
  cli::spinon::Controls moment_controls;
  cli::spinon::add_options(moments, moment_controls, true);
  EXPECT_THROW(moments.parse("--six-precision-stages reference"), CLI::ExtrasError);
}

TEST(CliPrecision, NativeStagePreservesPhysicalControlsAndUsesItsOwnDensityPair)
{
  if (!cli::has_fp80) GTEST_SKIP();
  namespace spinon = cli::spinon;
  spinon::Controls controls;
  controls.sectors = "6";
  controls.density_guard_bits = 320;
  controls.six_cutoff = "1040";
  controls.six_soft = "0.375";
  controls.relative_tolerance = "0.00003";
  controls.absolute_tolerance = "0";
  controls.terms = 96;
  controls.tail_order = 20;
  controls.recovery_seconds = 17;
  auto const native = spinon::native_six_stage<long double>(controls, 1.25L);
  // A physical boundary exposes the prepared controls without doing density work.
  auto const record = native.at(0.L, 2.L).evaluate(controls.sampling);
  EXPECT_EQ(record.exchange, 1.25L);
  EXPECT_EQ(record.rapidity_cutoff, 1040.L);
  EXPECT_EQ(record.proposal_soft_fraction, 0.375L);
  EXPECT_EQ(record.recovery.relative_tolerance, controls.relative_tolerance);
  EXPECT_EQ(record.recovery.absolute_tolerance, controls.absolute_tolerance);
  EXPECT_EQ(record.recovery.max_wall_seconds, controls.recovery_seconds);
  ASSERT_EQ(record.recovery.stages.size(), 1u);
  EXPECT_EQ(record.recovery.stages[0].bits, 53);
  EXPECT_EQ(record.recovery.stages[0].guard_bits, 64);
  EXPECT_EQ(record.recovery.stages[0].contours[0].residue_terms, 96);
  EXPECT_EQ(record.recovery.stages[0].contours[1].residue_terms, 112);
  EXPECT_EQ(record.recovery.stages[0].contours[0].tail_order, 20);
  EXPECT_EQ(record.recovery.stages[0].contours[1].tail_order, 24);
  EXPECT_EQ(controls.density_guard_bits, 320u);
  controls.sampling.uniform_bits = 64;
  EXPECT_THROW(spinon::native_six_stage<long double>(controls, 1.L), std::invalid_argument);
  controls.sampling.uniform_bits = 52;
  EXPECT_THROW(spinon::native_six_stage<long double>(controls, 0.L), std::invalid_argument);
}

#if UNI20_HAS_MPREAL
TEST(CliPrecision, SixPrecisionDispatchKeepsPhysicalResultsAndPerReplicateStageWork)
{
  if (!cli::has_fp80) GTEST_SKIP();
  namespace spinon = cli::spinon;
  using R = uni20::mpreal;
  auto const precision = uni20::Precision::bits(192);
  R const one(1, precision), q("1.5", precision), omega(2, precision);
  bethe::heisenberg::StructureFactorOptions<R> options(precision);
  spinon::Controls controls;
  controls.sectors = "6";
  controls.sampling = {2, 4, 42, 0, 52};
  controls.samples_per_task = 3;
  controls.recovery_seconds = 60;
  std::shared_ptr<spinon::Plan<R> const> fixed = std::make_shared<spinon::Plan<R>>(controls, one, options);
  controls.precision_stages = "fp80,reference";
  std::shared_ptr<spinon::Plan<R> const> staged = std::make_shared<spinon::Plan<R>>(controls, one, options);
  uni20::async::TbbScheduler scheduler(2);
  bethe::ExecutionOptions execution{.scheduler = &scheduler, .batch_size = 3, .max_in_flight = 2};
  auto reference = spinon::point_async(fixed, q, omega, 0, execution).get_wait(scheduler);
  auto mixed = spinon::point_async(staged, q, omega, 0, execution).get_wait(scheduler);
  ASSERT_TRUE(reference.six && mixed.six);
  ASSERT_TRUE(reference.six->value && mixed.six->value);
  EXPECT_LT(uni20::math::abs(*mixed.six->value - *reference.six->value),
            R("1e-12", precision) * uni20::math::abs(*reference.six->value));
  EXPECT_EQ(mixed.six->supported_roots, reference.six->supported_roots);
  EXPECT_EQ(mixed.six->completed_samples, 8u);
  EXPECT_TRUE(reference.six_stages.empty());
  ASSERT_EQ(mixed.six_stages.size(), 2u);
  EXPECT_EQ(mixed.six_stages[0].attempts, 8u);
  EXPECT_EQ(mixed.six_stages[0].accepted + mixed.six_stages[1].accepted, 8u);
  EXPECT_EQ(mixed.six_stages[0].retried, mixed.six_stages[1].attempts);
  ASSERT_EQ(mixed.six_replicate_stages.size(), 2u);
  for (auto const& row : mixed.six_replicate_stages)
  {
    ASSERT_TRUE(row);
    ASSERT_EQ(row->size(), 2u);
    EXPECT_EQ((*row)[0].attempts, 4u);
    EXPECT_EQ((*row)[0].accepted + (*row)[1].accepted, 4u);
  }
  R const seam = R{8} * uni20::math::atan(one);
  auto boundary = spinon::point_async(staged, seam, omega, 0, execution).get_wait(scheduler);
  ASSERT_TRUE(boundary.six);
  EXPECT_EQ(boundary.six->status, bethe::heisenberg::SixSpinonSpectrumStatus::outside_continuum);
  EXPECT_EQ(boundary.six->value, one - one);
  EXPECT_TRUE(boundary.six_replicate_stages.empty());
  ASSERT_EQ(boundary.six_stages.size(), 2u);
  EXPECT_EQ(boundary.six_stages[0].attempts + boundary.six_stages[1].attempts, 0u);
}
#endif
