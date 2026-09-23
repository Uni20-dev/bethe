// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "data-output.hpp"
#include "test_support.hpp"
#include <sstream>

namespace
{
namespace cli = bethe::cli;
namespace data = uni20::presentation;

TEST(DataOutput, ShellQuotingAndCommentEscaping)
{
  EXPECT_EQ(cli::quote_argument(""), "''");
  EXPECT_EQ(cli::quote_argument("a'b $HOME"), "'a'\\''b $HOME'");
  EXPECT_EQ(cli::comment_text("a\nb\rc\td\\e\x1b"), "a\\nb\\rc\\td\\\\e\\x1b");
  std::ostringstream out;
  cli::write_comments(out, {{"Command", "with\nnewline"}});
  EXPECT_EQ(out.str(), "# Command: with\\nnewline\n");
}

template <typename Real> class DataOutputPrecision : public ::testing::Test {};
TYPED_TEST_SUITE(DataOutputPrecision, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(DataOutputPrecision, MachineOutputRetainsNativePrecision)
{
  using Real = TypeParam;
  Real const x = Real{1} + uni20::numeric_limits<Real>::epsilon();
  std::ostringstream csv, tsv, json;
  auto table = data::make_data_table("test", {.metadata = {{"U", "4"}}}, data::data_column<Real>("x"),
                                     data::data_column<uni20::half_int>("spin"),
                                     data::data_column<std::optional<Real>>("missing"));
  table.attach(cli::CommentedSink(csv, data::csv_sink(csv), true));
  table.attach(cli::CommentedSink(tsv, data::tsv_sink(tsv), false));
  table.attach(data::json_sink(json));
  table.append(x, uni20::from_twice(1), std::nullopt);
  table.finish({{"Status", "converged"}});
  auto const start = csv.str().find("x,spin,missing\n") + std::string("x,spin,missing\n").size();
  auto const token = csv.str().substr(start, csv.str().find(',', start) - start);
  test_support::expect_exact(uni20::parse_real<Real>(token), x, "CSV full precision");
  EXPECT_NE(csv.str().find(token + ",0.5,\n# Status: converged\n"), std::string::npos);
  EXPECT_EQ(tsv.str(), "x\tspin\tmissing\n" + token + "\t0.5\t\n");
  EXPECT_NE(json.str().find("[\"" + token + "\",0.5,null]"), std::string::npos);
}

struct FailingSink
{
    template <typename... Args> void begin(Args const&...) {}
    template <typename... Args> void row(Args const&...) { throw std::runtime_error("disk full"); }
    void finish(data::table_metadata const&) {}
};

TEST(DataOutput, AcceptedRowIsNotRetriedAndHealthyJsonFinishesOnAbort)
{
  cli::DataOutputOptions options;
  options.quiet = true;
  cli::DataOutput output(options);
  std::ostringstream json;
  auto table = data::make_data_table("test", {}, data::data_column<int>("i"));
  // Put the failing sink first: the healthy sink must still receive the row.
  table.attach(cli::NamedSink("failed.csv", FailingSink{}));
  table.attach(data::json_sink(json));
  try
  {
    table.append(42);
    FAIL() << "expected required delivery error";
  }
  catch (data::data_delivery_error const& error)
  {
    std::ostringstream diagnostic;
    cli::print_output_error(diagnostic, error);
    EXPECT_NE(diagnostic.str().find("failed.csv: disk full"), std::string::npos);
  }
  EXPECT_EQ(table.size(), 1u);
  output.abort(table, {{"Status", "aborted"}});
  EXPECT_TRUE(table.finished());
  EXPECT_NE(json.str().find("\"rows\":[[42]]"), std::string::npos);
  EXPECT_NE(json.str().find("\"Status\":\"aborted\""), std::string::npos);
}

struct FlushFailure : std::stringbuf
{
    int sync() override { return -1; }
};

TEST(DataOutput, FinalFlushFailureIsReportedWithoutLosingHealthyOutput)
{
  cli::DataOutputOptions options;
  options.quiet = true;
  cli::DataOutput output(options);
  FlushFailure buffer;
  std::ostream broken(&buffer);
  std::ostringstream good;
  auto table = data::make_data_table("test", {}, data::data_column<int>("i"));
  table.attach(cli::NamedSink("failed.csv", cli::CommentedSink(broken, data::csv_sink(broken), true)));
  table.attach(data::json_sink(good));
  table.append(1);
  EXPECT_THROW(output.finish(table, {{"Status", "converged"}}), data::data_delivery_error);
  // An I/O failure does not retroactively turn computed rows into failed solves.
  output.abort(table, {{"Status", "aborted"}});
  EXPECT_EQ(table.summary()->at("Status"), "converged");
  EXPECT_NE(good.str().find("\"Status\":\"converged\""), std::string::npos);
}

TEST(DataOutput, FinalReportNeedsRetainedRows)
{
  cli::DataOutputOptions options;
  options.retain = false;
  EXPECT_THROW(options.validate(), std::invalid_argument);
  options.stream = true;
  EXPECT_NO_THROW(options.validate());
  options.stream = false;
  options.quiet = true;
  EXPECT_NO_THROW(options.validate());
  options.quiet = false;
  options.format = "json";
  EXPECT_NO_THROW(options.validate());
}

TEST(DataOutput, LateCommentedSinkReplaysRetainedRowsAndSummary)
{
  std::ostringstream early, late;
  auto table = data::make_data_table("test", {.metadata = {{"U", "4"}}}, data::data_column<int>("i"));
  table.attach(cli::CommentedSink(early, data::csv_sink(early), true));
  table.append(1);
  table.append(2);
  table.finish({{"Status", "converged"}});
  table.attach(cli::CommentedSink(late, data::csv_sink(late), true));
  EXPECT_EQ(late.str(), early.str());
  EXPECT_NE(late.str().find("# First row: 0\ni\n1\n2\n# Status: converged\n"), std::string::npos);
}

TEST(DataOutput, NonretainedFutureOnlySinkReportsItsOffset)
{
  std::ostringstream out;
  auto table = data::make_data_table("test", {.retain = data::retention::none}, data::data_column<int>("i"));
  table.append(1);
  table.attach(cli::CommentedSink(out, data::csv_sink(out), true), {.replay = data::sink_replay::future_only});
  table.append(2);
  table.finish({{"Rows", std::to_string(table.size())}});
  EXPECT_EQ(table.size(), 2u);
  EXPECT_EQ(table.retained_size(), 0u);
  EXPECT_EQ(out.str(), "# First row: 1\ni\n2\n# Rows: 2\n");
}
} // namespace
