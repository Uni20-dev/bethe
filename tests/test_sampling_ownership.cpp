// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/xxx_six_spinon_moments.hpp>
#include <bethe/xxx_six_spinon_structure_factor.hpp>
#include <uni20/async/debug_scheduler.hpp>

namespace
{
using namespace uni20::async;
class RejectSubmission : public DebugScheduler {
  public:
    explicit RejectSubmission(std::size_t n) : reject(n) {}
    std::size_t calls = 0, reject;
    void schedule(AsyncTask&& task) override
    {
      if (++calls == reject) throw std::bad_alloc{};
      DebugScheduler::schedule(std::move(task));
    }
};

template <class Submit> void check_owner(Submit submit)
{
  RejectSubmission first(1);
  EXPECT_THROW((void)submit(first), std::bad_alloc);
  EXPECT_TRUE(first.done());

  // Only one orchestration task may be admitted by the public call. Later
  // observation admission failures belong to its returned result and drain
  // any blocks accepted before the failing admission.
  for (std::size_t reject : {2u, 3u})
  {
    SCOPED_TRACE(reject);
    RejectSubmission scheduler(reject);
    try
    {
      auto pending = submit(scheduler);
      EXPECT_EQ(scheduler.calls, 1u);
      EXPECT_THROW((void)pending.get_wait(scheduler), std::bad_alloc);
    }
    catch (std::bad_alloc const&)
    {
      ADD_FAILURE() << "public call admitted a dependency before failing";
    }
    EXPECT_TRUE(scheduler.done());
    if (!scheduler.done()) scheduler.run_all(); // Drain only a failing old-code test.
  }
}
} // namespace

TEST(SamplingOwnership, MonteCarloAdmissionFailuresRetainOneOwner)
{
  check_owner([](auto& scheduler) {
    return bethe::monte_carlo_async<1>(
        {2, 2, 42}, 1, 1.,
        [](auto const&, std::stop_token) { return bethe::MonteCarloObservation<double>{std::vector{1.}}; },
        {&scheduler, 1, 2});
  });
}

TEST(SamplingOwnership, ThrowingObserverMoveCannotAdmitWork)
{
  struct ThrowMove
  {
      ThrowMove() = default;
      ThrowMove(ThrowMove&&) { throw std::bad_alloc{}; }
      bethe::MonteCarloObservation<double> operator()(bethe::UniformSample<1> const&, std::stop_token) const
      {
        return {std::vector{1.}};
      }
  };
  RejectSubmission scheduler(99);
  EXPECT_THROW((void)bethe::monte_carlo_async<1>({2, 2}, 1, 1., ThrowMove{}, {&scheduler}), std::bad_alloc);
  EXPECT_EQ(scheduler.calls, 0u);
  EXPECT_TRUE(scheduler.done());
}

TEST(SamplingOwnership, SpectrumAndMomentsOwnTheirSamplingDependency)
{
  using namespace bethe::heisenberg;
  if (uni20::numeric_limits<long double>::digits() <= uni20::numeric_limits<double>::digits()) GTEST_SKIP();
  bethe::heisenberg::detail::SpinonJointProposal<long double> proposal(
      8.L, .5L, {{bethe::heisenberg::detail::SpinonProposalKind::background, 1.L, 1.L, 1.L}});
  SpinonRecoveryPolicy policy{
      "1e-5", "0", {{53, uni20::numeric_limits<long double>::digits(), {{{32, 12}, {64, 20}}}}}, 60};
  ThermodynamicSixSpinonStructureFactor<long double, double, long double> spectrum(proposal, policy, 1.L);
  ThermodynamicSixSpinonMoments<long double, double, long double> moments(proposal, policy, 1.L);
  auto point = spectrum.at(4 * std::atan(1.L) / 2, 2.L);
  check_owner([&](auto& scheduler) { return point.evaluate_async({2, 2, 42}, {&scheduler, 1, 2}); });
  check_owner([&](auto& scheduler) {
    return point.evaluate_observations_async({2, 2, 42},
                                             [](auto const&, std::stop_token) {
                                               return bethe::MonteCarloObservation<long double, SpinonSampleMetadata>{
                                                   std::vector{1.L, 0.L, 0.L}};
                                             },
                                             {&scheduler, 1, 2});
  });
  check_owner(
      [&](auto& scheduler) { return moments.fixed_moments(1.L).evaluate_async({2, 2, 42}, {&scheduler, 1, 2}); });
  check_owner([&](auto& scheduler) { return moments.global_async({2, 2, 42}, {&scheduler, 1, 2}); });
}
