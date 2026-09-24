// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "run-metadata.hpp"
#include "test_support.hpp"

namespace
{
namespace cli = bethe::cli;

template <typename Real> class RunMetadataPrecision : public ::testing::Test {};
TYPED_TEST_SUITE(RunMetadataPrecision, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(RunMetadataPrecision, NativeSnapshotAndLegacyProjection)
{
  using Real = TypeParam;
  Real const x = Real{1} + uni20::numeric_limits<Real>::epsilon();
  uni20::metadata_document values;
  values.group("model", "Model");
  values.add("model", "interaction", x, {.label = "Human label", .display = {.numeric = {.precision = 2}}});
  values.add("model", "spin", uni20::from_twice(1));
  values.add("model", "missing", std::optional<Real>{});
  values.add("model", "empty", "");
  values.add("model", "infinite", uni20::numeric_limits<Real>::infinity());
  auto snapshot = values;
  values.replace("interaction", Real{2});
  auto exported = cli::run_metadata(snapshot, {{"interaction", "U (t=1)"}, {"spin", "Spin"}});
  EXPECT_EQ(snapshot.find("interaction")->value.type(), typeid(Real));
  test_support::expect_exact(snapshot.find("interaction")->value.template get<Real>(), x, "native value");
  test_support::expect_exact(uni20::parse_real<Real>(exported.at("U (t=1)")), x, "round-trip metadata");
  EXPECT_EQ(exported.at("Spin"), "0.5");
  EXPECT_EQ(exported.at("missing"), "unavailable");
  EXPECT_TRUE(snapshot.find("missing")->value.missing());
  EXPECT_EQ(exported.at("empty"), "");
  EXPECT_EQ(exported.at("infinite"), "inf");
  EXPECT_THROW(cli::run_metadata(snapshot, {{"interaction", "Spin"}, {"spin", "Spin"}}), std::invalid_argument);
}

uni20::metadata_document result_summary(std::string status, std::size_t rows)
{
  uni20::metadata_document result;
  result.group("results");
  result.add("results", "status", std::move(status));
  result.add("results", "rows", rows);
  return result;
}

TEST(RunMetadata, ComputationRunAndElapsedTimingStayDistinct)
{
  uni20::run_clock_sample now{10, 1};
  uni20::run_context run({.name = "test"}, {.clock = [&] { return now; }, .utc = [] { return "fixed"; }});
  now = {11, 2};
  run.measure([&] { now = {13, 3.125L}; });
  now = {15, 4}; // Work outside computation intervals, such as rendering/export.
  run.measure([&] { now = {17, 6}; });
  now = {17.25L, 6.5L};
  auto const& summary = run.finish(uni20::run_outcome::partial, result_summary("incomplete", 3));
  auto exported = cli::run_summary(summary);
  EXPECT_EQ(summary.find("compute_cpu_seconds")->value.get<std::optional<long double>>(), 3.125L);
  EXPECT_EQ(summary.find("rows")->value.get<std::size_t>(), 3u);
  EXPECT_EQ(exported.at("CPU time"), "3.125000 s");
  EXPECT_EQ(exported.at("Run CPU seconds"), "5.5");
  EXPECT_EQ(exported.at("Elapsed seconds"), "7.25");
  EXPECT_EQ(exported.at("Outcome"), "partial");
  EXPECT_EQ(exported.at("Status"), "incomplete");
  EXPECT_EQ(exported.at("Rows"), "3");
  now = {100, 90};
  EXPECT_EQ(cli::run_summary(summary), exported); // Formatting never resamples clocks.
}

TEST(RunMetadata, UnavailableCpuAndFailedComputationRemainExplicit)
{
  uni20::run_clock_sample now{10, 1};
  uni20::run_context run({.name = "test"}, {.clock = [&] { return now; }});
  EXPECT_THROW(run.measure([&] {
    now = {12, std::nullopt};
    throw std::runtime_error("computation failed");
  }),
               std::runtime_error);
  auto const& summary = run.finish(uni20::run_outcome::failed, result_summary("aborted", 0));
  auto exported = cli::run_summary(summary);
  EXPECT_EQ(exported.at("CPU time"), "unavailable");
  EXPECT_EQ(exported.at("Run CPU seconds"), "unavailable");
  EXPECT_EQ(exported.at("Elapsed seconds"), "2");
  EXPECT_EQ(exported.at("Outcome"), "failed");
  EXPECT_EQ(exported.at("Status"), "aborted");
}
} // namespace
