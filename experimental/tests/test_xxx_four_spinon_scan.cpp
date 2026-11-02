// SPDX-License-Identifier: GPL-3.0-or-later
#include "test_support.hpp"
#include "xxx_complex_fixture.hpp"
#include <bethe/experimental/xxx_four_spinon_scan.hpp>
#include <condition_variable>
#include <set>
#include <thread>
#include <uni20/async/tbb_scheduler.hpp>

namespace f = bethe::experimental::xxx;
using namespace uni20::async;
namespace
{
template <uni20::Real R> R parse(std::string_view text, R const& zero)
{
  if constexpr (uni20::has_runtime_precision_v<R>)
    return uni20::parse_real<R>(text, uni20::precision_of(zero));
  else
    return uni20::parse_real<R>(text);
}
template <uni20::Real R> f::PreparedSinglet<R> background(R const& zero)
{
  bethe::heisenberg::RealState<R> ground;
  ground.converged = true;
  ground.sz = uni20::half_int{0};
  ground.quantum_numbers = bethe::heisenberg::sector_ground_quantum_numbers(12, ground.sz);
  for (auto text : complex_fixture::ground)
    ground.rapidities.push_back(2 * parse(text, zero));
  return {12, ground};
}
std::vector<f::FourSpinonRequest> selected()
{
  // Deliberately not canonical order; includes an unsupported centred pair.
  return {{12, {-2, 0, 1, 3}, -1}, {12, {2, 0, -3, -1}, 1}, {12, {-2, -1, 1, 2}, 0}};
}

template <uni20::Real R> void numerical_scan(R const& zero)
{
  auto ground = background(zero);
  ASSERT_TRUE(ground.ready());
  auto requests = selected();
  requests.push_back({12, {-2, -1, 0, 2}, -1}); // Exact q=0 selection-rule zero.
  DebugScheduler serial;
  auto reference = f::scan_four_spinons(ground, requests, {&serial});
  TbbScheduler concurrent(2);
  auto pending = f::scan_four_spinons_async(std::move(ground), requests, {&concurrent, 1, 2});
  EXPECT_FALSE(ground.ready());
  auto const& scan = pending.get_wait(concurrent);
  ASSERT_EQ(scan.records.size(), 4u);
  ASSERT_EQ(scan.lines.size(), 3u);
  EXPECT_FALSE(scan.converged());
  ASSERT_EQ(scan.moments.size(), 12u);
  R const one = uni20::scalar_like(zero, 1), eps = uni20::numeric_limits<R>::epsilon(zero);
  auto const expected = parse(complex_fixture::weight, zero);
  std::set<std::size_t> transfers;
  for (std::size_t i = 0; i < scan.records.size(); ++i)
  {
    auto const& row = scan.records[i];
    EXPECT_EQ(row.candidate_id, reference.records[i].candidate_id);
    EXPECT_EQ(row.solved.status, reference.records[i].solved.status);
    if (!row.line)
    {
      EXPECT_EQ(row.solved.status, f::FourSpinonStatus::unsupported_geometry);
      EXPECT_FALSE(row.factor);
      EXPECT_FALSE(row.diagnostic.empty());
      continue;
    }
    ASSERT_TRUE(reference.records[i].line);
    EXPECT_EQ(row.line->weight, reference.records[i].line->weight);
    EXPECT_EQ(row.line->gap, reference.records[i].line->gap);
    EXPECT_EQ(row.line->state_id, i + 1);
    EXPECT_EQ(uni20::precision_of(row.line->weight), uni20::precision_of(zero));
    EXPECT_EQ(row.line->momentum_index, f::four_spinon_transfer_index(row.solved.request));
    EXPECT_EQ(row.line->momentum_index, (12 - row.solved.state->momentum_index()) % 12);
    if (row.line->momentum_index == 0)
    {
      EXPECT_EQ(row.line->weight, zero);
      ASSERT_TRUE(row.factor);
      EXPECT_TRUE(row.factor->converged());
      EXPECT_FALSE(row.factor->log_weight);
      EXPECT_TRUE(row.diagnostic.empty());
    }
    else
      EXPECT_REAL_NEAR(row.line->weight / expected, one, 1000000 * eps);
    EXPECT_EQ(scan.moments[row.line->momentum_index].weight, row.line->weight);
    EXPECT_EQ(scan.moments[row.line->momentum_index].first_moment, row.line->gap * row.line->weight);
    transfers.insert(row.line->momentum_index);
  }
  EXPECT_EQ(transfers, (std::set<std::size_t>{0, 3, 9}));
}

AsyncTask count_lines(ReadBuffer<f::FourSpinonScan<double>> input, WriteBuffer<std::size_t> output)
{
  auto const& scan = co_await input;
  co_await output = scan.lines.size();
}
} // namespace

TEST(ExperimentalFourCatalogue, IndependentExhaustiveLabelsAndMomentumFilters)
{
  for (int n : {8, 12, 16, 20})
  {
    std::vector<f::FourSpinonRequest> brute;
    int const h = n / 4;
    for (int a = -h; a <= h; ++a)
      for (int b = -h; b <= h; ++b)
        for (int c = -h; c <= h; ++c)
          for (int d = -h; d <= h; ++d)
            if (a < b && b < c && c < d)
              for (int branch : {-1, 0, 1})
                brute.push_back({n, {a, b, c, d}, branch});
    auto const all = f::four_spinon_catalogue(n);
    EXPECT_EQ(all, brute);
    std::set<std::string> identities;
    std::size_t count = 0;
    for (int q = 0; q < n; ++q)
    {
      std::vector<f::FourSpinonRequest> expected;
      for (auto const& r : brute)
      {
        int transfer = r.branch;
        for (int hole : r.holes)
          transfer -= hole;
        while (transfer < 0)
          transfer += n;
        if (transfer % n == q) expected.push_back(r);
      }
      auto const filtered = f::four_spinon_catalogue(n, q);
      EXPECT_EQ(filtered, expected);
      count += filtered.size();
      for (auto const& r : filtered)
        EXPECT_TRUE(identities.insert(f::four_spinon_id(r)).second);
    }
    EXPECT_EQ(count, all.size());
  }
  EXPECT_EQ(f::four_spinon_catalogue(12).size(), 105u);
  EXPECT_EQ(f::four_spinon_id({12, {2, -3, 0, -1}, 1}), "N12-s4-h-3,-1,0,2-b1");
  EXPECT_EQ(f::four_spinon_transfer_index({12, {-3, -1, 0, 2}, 1}), 3);
  EXPECT_THROW(f::four_spinon_catalogue(10), std::invalid_argument);
  EXPECT_THROW(f::four_spinon_catalogue(12, 12), std::invalid_argument);
  EXPECT_THROW(f::four_spinon_catalogue(12, {}, 104), std::length_error);
  auto const size = f::four_spinon_catalogue(12, 3).size();
  EXPECT_EQ(f::four_spinon_catalogue(12, 3, size).size(), size);
  EXPECT_THROW(f::four_spinon_catalogue(12, 3, size - 1), std::length_error);
}

template <typename R> class ExperimentalFourScan : public ::testing::Test {};
TYPED_TEST_SUITE(ExperimentalFourScan, test_support::RealTypes, test_support::PrecisionNames);
TYPED_TEST(ExperimentalFourScan, IndependentWeightsFailuresAndSerialConcurrentAgreement)
{
  numerical_scan(TypeParam{0});
}

TEST(ExperimentalFourScanExecution, OwnershipCompositionPreflightAndEmpty)
{
  for (auto order : {DebugSchedulerOrder::fifo, DebugSchedulerOrder::reverse, DebugSchedulerOrder::random})
  {
    DebugScheduler scheduler({.order = order, .random_seed = 17});
    Async<std::size_t> final;
    {
      auto pending = f::scan_four_spinons_async(background(0.), selected(), {&scheduler});
      scheduler.schedule(count_lines(pending.read(), final.write()));
    }
    EXPECT_EQ(final.get_wait(scheduler), 2u);
  }
  DebugScheduler scheduler;
  auto const ground = background(0.);
  auto requests = selected();
  requests.push_back(requests.front());
  EXPECT_THROW(f::scan_four_spinons_async(ground, requests, {&scheduler}), std::invalid_argument);
  requests = selected();
  requests.front().sites = 16;
  EXPECT_THROW(f::scan_four_spinons_async(ground, requests, {&scheduler}), std::invalid_argument);
  EXPECT_THROW(f::scan_four_spinons_async(ground, selected(), {&scheduler}, 2), std::length_error);
  EXPECT_THROW(f::scan_four_spinons_async(ground, selected(), {&scheduler, 1, 0}), std::invalid_argument);
  auto empty = f::scan_four_spinons(ground, {}, {&scheduler});
  EXPECT_TRUE(empty.converged());
  EXPECT_TRUE(empty.records.empty());
  ASSERT_EQ(empty.moments.size(), 12u);
  for (auto const& m : empty.moments)
    EXPECT_EQ(m.weight, 0.);
  struct Reject : DebugScheduler
  {
      void schedule(AsyncTask&&) override { throw std::bad_alloc{}; }
  } reject;
  EXPECT_THROW(f::scan_four_spinons_async(ground, selected(), {&reject}), std::bad_alloc);
  std::stop_source stop;
  stop.request_stop();
  auto cancelled = f::scan_four_spinons_async(ground, selected(), {&scheduler, 1, 2, stop.get_token()});
  EXPECT_THROW(cancelled.get_wait(scheduler), bethe::ExecutionCancelled);
}

TEST(ExperimentalFourScanExecution, FailureAndCancellationWaitForHeldSibling)
{
  TbbScheduler scheduler(2);
  for (bool cancel : {false, true})
  {
    std::mutex mutex;
    std::condition_variable cv;
    bool held = false, failed = false, release = false, observed = false, finished = false;
    std::stop_source stop;
    auto requests = selected();
    std::sort(requests.begin(), requests.end(), f::four_detail::request_less);
    auto pending = f::four_detail::submit_scan(
        background(0.), requests, {&scheduler, 1, 2, stop.get_token()}, 10,
        [&](auto const& request, auto const&, std::stop_token) {
          std::unique_lock lock(mutex);
          if (request == requests.front())
          {
            held = true;
            cv.notify_all();
            EXPECT_TRUE(cv.wait_for(lock, std::chrono::seconds(5), [&] { return release; }));
            finished = true;
          }
          else
          {
            EXPECT_TRUE(cv.wait_for(lock, std::chrono::seconds(5), [&] { return held; }));
            failed = true;
            cv.notify_all();
            if (cancel)
              stop.request_stop();
            else
              throw std::bad_alloc{};
          }
          return f::FourSpinonRecord<double>{};
        });
    // The host joins get_wait below. A separate observer checks completion
    // while one admitted kernel is still deliberately held, then releases it.
    std::jthread observer([&] {
      std::unique_lock lock(mutex);
      observed = cv.wait_for(lock, std::chrono::seconds(5), [&] { return held && failed; });
      EXPECT_TRUE(observed);
      auto reader = pending.read();
      EXPECT_FALSE(reader.await_ready());
      release = true;
      cv.notify_all();
    });
    if (cancel)
      EXPECT_THROW(pending.get_wait(scheduler), bethe::ExecutionCancelled);
    else
      EXPECT_THROW(pending.get_wait(scheduler), std::bad_alloc);
    observer.join();
    EXPECT_TRUE(finished);
  }
}

#if UNI20_HAS_MPREAL
TEST(ExperimentalFourScanMpfr, ExplicitPrecisionInWorkers)
{
  auto const ambient = mpfr_get_default_prec();
  for (int bits : {128, 256})
    numerical_scan(uni20::mpreal(0, uni20::Precision::bits(bits)));
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
#endif
