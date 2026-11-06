// SPDX-License-Identifier: GPL-3.0-or-later
#include "momentum-input.hpp"
#include <gtest/gtest.h>

namespace cli = bethe::cli;
TEST(CliMomentumInput, FormsAreMutuallyExclusive)
{
  for (auto command : {"--momentum 1 --momentum-beta 2", "--momentum-from-pi 1 --momentum 1",
                       "--momentum-beta 1 --momentum-from-pi 1"})
  {
    CLI::App app;
    cli::MomentumArguments args;
    cli::add_momentum_options(app, args);
    EXPECT_THROW(app.parse(command), CLI::ExcludesError);
  }
}
TEST(CliMomentumInput, PreservesInputAndDetectsScalarCollapse)
{
  CLI::App app;
  cli::MomentumArguments args;
  cli::add_momentum_options(app, args);
  app.parse("--momentum-from-pi 1e-100");
  EXPECT_TRUE(bool(args));
  EXPECT_EQ(args.input().text, "1e-100");
  auto r = args.resolve(1.);
  ASSERT_TRUE(r.coordinate.beta());
  EXPECT_GT(*r.coordinate.beta(), 230.);
  EXPECT_THROW(cli::momentum_radians(r), std::range_error);
}
TEST(CliMomentumInput, WarnsWhenOnlyScalarGeometryLosesAccuracy)
{
  cli::MomentumArguments args;
  args.from_pi = "1e-10";
  auto r = args.resolve(1.);
  testing::internal::CaptureStderr();
  (void)cli::momentum_radians(r);
  auto text = testing::internal::GetCapturedStderr();
  EXPECT_NE(text.find("scalar momentum geometry"), std::string::npos);
  EXPECT_EQ(text.find("input to working beta"), std::string::npos);
}
