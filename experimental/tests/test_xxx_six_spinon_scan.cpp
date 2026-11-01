// SPDX-License-Identifier: GPL-3.0-or-later
#include "test_support.hpp"
#include "xxx_complex_fixture.hpp"
#include "xxx_six_spinon_fixture.hpp"
#include <bethe/experimental/xxx_six_spinon_scan.hpp>
#include <condition_variable>
#include <functional>
#include <set>
#include <thread>
#include <uni20/async/tbb_scheduler.hpp>

namespace six = bethe::experimental::xxx::six_spinon;
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
six::Configuration request(six_fixture::Entry const& entry)
{
  six::Configuration c{entry.sites, entry.triple ? six::Family::triple : six::Family::two_pairs, {}, {}};
  for (int twice : entry.holes_twice)
    c.holes.push_back(uni20::from_twice(twice));
  for (int i = 0; i < (entry.triple ? 1 : 2); ++i)
    c.branches.push_back(uni20::from_twice(entry.branches_twice[i]));
  return c;
}
std::vector<six::Configuration> selected()
{
  return {request(six_fixture::entries[7]), request(six_fixture::entries[0]), request(six_fixture::entries[6]),
          request(six_fixture::entries[1])};
}

template <uni20::Real R> void numerical_scan(R const& zero)
{
  auto ground = background(zero);
  ASSERT_TRUE(ground.ready());
  DebugScheduler serial;
  auto reference = six::scan(ground, selected(), {&serial});
  TbbScheduler concurrent(2);
  auto pending = six::scan_async(std::move(ground), selected(), {&concurrent, 1, 2});
  EXPECT_FALSE(ground.ready());
  auto const& result = pending.get_wait(concurrent);
  ASSERT_TRUE(result.converged());
  ASSERT_EQ(result.lines.size(), 4u);
  ASSERT_EQ(result.records.size(), 4u);
  ASSERT_EQ(result.moments.size(), 12u);
  R const one = uni20::scalar_like(zero, 1), eps = uni20::numeric_limits<R>::epsilon(zero);
  for (std::size_t i = 0; i < result.records.size(); ++i)
  {
    auto const& row = result.records[i];
    ASSERT_TRUE(row.solved.converged());
    ASSERT_TRUE(row.factor);
    ASSERT_TRUE(row.line);
    EXPECT_EQ(row.candidate_id, reference.records[i].candidate_id);
    EXPECT_EQ(row.line->weight, reference.records[i].line->weight);
    EXPECT_EQ(row.line->gap, reference.records[i].line->gap);
    EXPECT_EQ(row.line->state_id, i + 1);
    EXPECT_EQ(row.line->momentum_index, (12 - row.solved.state->momentum_index) % 12);
    EXPECT_EQ(uni20::precision_of(row.line->weight), uni20::precision_of(zero));
    bool matched = false;
    for (auto index : {0, 1, 6, 7})
      if (row.candidate_id == six::candidate_id(request(six_fixture::entries[index])))
      {
        EXPECT_REAL_NEAR(row.line->weight / parse(six_fixture::entries[index].weight, zero), one, 1000000000 * eps);
        matched = true;
      }
    EXPECT_TRUE(matched);
  }
  for (std::size_t q = 0; q < 12; ++q)
  {
    R m0 = zero, m1 = zero;
    for (auto const& line : result.lines)
      if (line.momentum_index == q)
      {
        m0 += line.weight;
        m1 += line.weight * line.gap;
      }
    EXPECT_REAL_NEAR(result.moments[q].weight, m0, 8 * eps);
    EXPECT_REAL_NEAR(result.moments[q].first_moment, m1, 8 * eps);
  }
}

AsyncTask count_lines(ReadBuffer<six::ScanResult<double>> input, WriteBuffer<std::size_t> output)
{
  auto const& result = co_await input;
  co_await output = result.lines.size();
}
} // namespace

TEST(ExperimentalSixCatalogue, ExhaustiveLabelsAndMomentumIndependentOfSuffixCounts)
{
  for (int n : {12, 16, 20})
  {
    std::vector<std::vector<std::string>> expected(n);
    std::array<std::size_t, 2> family_counts{};
    for (bool triple : {false, true})
    {
      int const edge = n / 2 + triple;
      std::vector<std::vector<int>> branches;
      if (triple)
        branches = {{-2}, {0}, {2}};
      else
        for (int a : {-3, -1, 1, 3})
          for (int b : {-3, -1, 1, 3})
            if (a < b) branches.push_back({a, b});
      for (auto const& branch : branches)
      {
        six::Configuration c{n, triple ? six::Family::triple : six::Family::two_pairs, {}, {}};
        for (int b : branch)
          c.branches.push_back(uni20::from_twice(b));
        std::function<void(int)> choose = [&](int start) {
          if (c.holes.size() == 6)
          {
            // Independent momentum from occupied quantum numbers and string count.
            int sum_twice = 0;
            for (auto b : c.branches)
              sum_twice += b.twice();
            for (int t = -edge; t <= edge; t += 2)
              if (std::find(c.holes.begin(), c.holes.end(), uni20::from_twice(t)) == c.holes.end()) sum_twice += t;
            int p = ((n / 2) * (n / 2 - 3) - sum_twice / 2) % n;
            while (p < 0)
              p += n;
            int const q = (n - p) % n;
            expected[q].push_back(six::candidate_id(c));
            EXPECT_EQ(six::transfer_index(c), q);
            ++family_counts[triple];
            return;
          }
          for (int t = start; t <= edge; t += 2)
          {
            c.holes.push_back(uni20::from_twice(t));
            choose(t + 2);
            c.holes.pop_back();
          }
        };
        choose(-edge);
      }
    }
    std::array<std::size_t, 2> totals{};
    std::set<std::string> identities;
    for (int q = 0; q < n; ++q)
    {
      six::MomentumCatalogue const catalogue(n, q);
      ASSERT_EQ(catalogue.size(), expected[q].size());
      totals[0] += catalogue.family_size(six::Family::two_pairs);
      totals[1] += catalogue.family_size(six::Family::triple);
      for (std::size_t i = 0; i < catalogue.size(); ++i)
      {
        auto const c = catalogue.candidate(i);
        EXPECT_EQ(six::candidate_id(c), expected[q][i]);
        EXPECT_EQ(six::transfer_index(c), q);
        EXPECT_TRUE(identities.insert(six::candidate_id(c)).second);
      }
      EXPECT_THROW(catalogue.candidate(catalogue.size()), std::out_of_range);
    }
    EXPECT_EQ(totals, family_counts);
    if (n == 12)
    {
      EXPECT_EQ(totals, (std::array<std::size_t, 2>{42, 84}));
    }
  }
}

TEST(ExperimentalSixCatalogue, IndexedReferenceIdentityAndBudgets)
{
  six::MomentumCatalogue const c(64, 16);
  EXPECT_EQ(c.size(), 166849u);
  EXPECT_EQ(c.family_size(six::Family::two_pairs), 103817u);
  EXPECT_EQ(c.family_size(six::Family::triple), 63032u);
  EXPECT_EQ(six::candidate_id(c.candidate(0)), "N64-22-h-32,-30,-28,-26,-16,32-b-3,-1");
  EXPECT_EQ(six::candidate_id(c.candidate(83424)), "N64-22-h-16,-12,2,14,18,28-b-1,3");
  EXPECT_EQ(six::candidate_id(c.candidate(166848)), "N64-3-h21,23,25,29,31,33-b2");
  EXPECT_THROW(c.requests(), std::length_error);
  EXPECT_THROW(six::MomentumCatalogue(10, 0), std::invalid_argument);
  EXPECT_THROW(six::MomentumCatalogue(12, 12), std::invalid_argument);
  EXPECT_THROW(six::MomentumCatalogue(12, 0, 0), std::length_error);
  EXPECT_THROW(six::MomentumCatalogue(std::numeric_limits<int>::max() - 3, 0), std::length_error);
  auto small = six::MomentumCatalogue(12, 3);
  EXPECT_EQ(small.requests(small.size()).size(), small.size());
  EXPECT_THROW(small.requests(small.size() - 1), std::length_error);
  EXPECT_THROW(small.requests(0), std::invalid_argument);
  EXPECT_THROW(small.family_size(static_cast<six::Family>(4)), std::invalid_argument);
  EXPECT_THROW(six::detail::checked_count_sum(std::numeric_limits<std::size_t>::max(), 1), std::overflow_error);
  auto canonical = request(six_fixture::entries[0]);
  auto shuffled = canonical;
  std::reverse(shuffled.holes.begin(), shuffled.holes.end());
  std::reverse(shuffled.branches.begin(), shuffled.branches.end());
  EXPECT_EQ(six::candidate_id(canonical), six::candidate_id(shuffled));
}

template <typename R> class ExperimentalSixScan : public ::testing::Test {};
TYPED_TEST_SUITE(ExperimentalSixScan, test_support::RealTypes, test_support::PrecisionNames);
TYPED_TEST(ExperimentalSixScan, IndependentWeightsAndSerialConcurrentAgreement) { numerical_scan(TypeParam{0}); }

TEST(ExperimentalSixScanExecution, OwnershipCompositionPreflightAndEmpty)
{
  for (auto order : {DebugSchedulerOrder::fifo, DebugSchedulerOrder::reverse, DebugSchedulerOrder::random})
  {
    DebugScheduler scheduler({.order = order, .random_seed = 17});
    Async<std::size_t> final;
    {
      auto pending = six::scan_async(background(0.), selected(), {&scheduler});
      scheduler.schedule(count_lines(pending.read(), final.write()));
    }
    EXPECT_EQ(final.get_wait(scheduler), 4u);
  }
  DebugScheduler scheduler;
  auto ground = background(0.);
  auto requests = selected();
  requests.push_back(requests.front());
  std::reverse(requests.back().holes.begin(), requests.back().holes.end());
  EXPECT_THROW(six::scan_async(ground, requests, {&scheduler}), std::invalid_argument);
  requests = selected();
  requests.front().sites = 16;
  EXPECT_THROW(six::scan_async(ground, requests, {&scheduler}), std::invalid_argument);
  EXPECT_THROW(six::scan_async(ground, selected(), {&scheduler}, 3), std::length_error);
  EXPECT_THROW(six::scan_async(ground, selected(), {&scheduler, 1, 0}), std::invalid_argument);
  auto empty = six::scan(ground, {}, {&scheduler});
  EXPECT_TRUE(empty.converged());
  EXPECT_TRUE(empty.records.empty());
  std::stop_source stop;
  stop.request_stop();
  auto cancelled = six::scan_async(ground, selected(), {&scheduler, 1, 2, stop.get_token()});
  EXPECT_THROW((void)cancelled.get_wait(scheduler), bethe::ExecutionCancelled);
}

TEST(ExperimentalSixScanExecution, FailureAndCancellationDrainAdmittedKernels)
{
  TbbScheduler scheduler(2);
  for (bool cancel : {false, true})
  {
    std::mutex mutex;
    std::condition_variable cv;
    bool held = false, failed = false, release = false, finished = false;
    std::stop_source stop;
    auto requests = selected();
    std::sort(requests.begin(), requests.end(), six::detail::request_less);
    auto const first = six::candidate_id(requests.front());
    auto pending =
        six::detail::submit_scan(background(0.), requests, {&scheduler, 1, 2, stop.get_token()}, 10,
                                 [&](auto const& request, auto const&, std::stop_token) {
                                   std::unique_lock lock(mutex);
                                   if (six::candidate_id(request) == first)
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
                                   return six::ScanRecord<double>{};
                                 });
    std::jthread observer([&] {
      std::unique_lock lock(mutex);
      EXPECT_TRUE(cv.wait_for(lock, std::chrono::seconds(5), [&] { return held && failed; }));
      auto reader = pending.read();
      EXPECT_FALSE(reader.await_ready());
      release = true;
      cv.notify_all();
    });
    if (cancel)
      EXPECT_THROW((void)pending.get_wait(scheduler), bethe::ExecutionCancelled);
    else
      EXPECT_THROW((void)pending.get_wait(scheduler), std::bad_alloc);
    observer.join();
    EXPECT_TRUE(finished);
  }
}

TEST(ExperimentalSixScanExecution, UnresolvedCandidateRetainsIdentityAndNoLine)
{
  six::Configuration const request{
      20, six::Family::two_pairs, {-4, -3, -1, 0, 1, 4}, {uni20::from_twice(1), uni20::from_twice(3)}};
  auto ground_state = bethe::heisenberg::sector_ground_state<double>(20, uni20::half_int{0});
  f::PreparedSinglet<double> ground(20, ground_state);
  DebugScheduler scheduler;
  auto const result = six::scan(std::move(ground), {request}, {&scheduler});
  ASSERT_EQ(result.records.size(), 1u);
  EXPECT_FALSE(result.converged());
  EXPECT_EQ(result.records[0].candidate_id, six::candidate_id(request));
  EXPECT_EQ(result.records[0].solved.status, six::SolveStatus::unresolved_geometry);
  EXPECT_FALSE(result.records[0].factor);
  EXPECT_FALSE(result.records[0].line);
  EXPECT_FALSE(result.records[0].diagnostic.empty());
  EXPECT_TRUE(result.lines.empty());
}

TEST(ExperimentalSixScanExecution, SmallCatalogueRetainsMixedSuccessAndFailure)
{
  six::MomentumCatalogue const catalogue(12, 3);
  DebugScheduler scheduler;
  auto const result = six::scan(background(0.), catalogue.requests(), {&scheduler});
  ASSERT_EQ(result.records.size(), catalogue.size());
  ASSERT_FALSE(result.lines.empty());
  EXPECT_FALSE(result.converged());
  std::size_t accepted = 0, rejected = 0;
  double m0 = 0, m1 = 0;
  for (std::size_t i = 0; i < result.records.size(); ++i)
  {
    auto const& row = result.records[i];
    EXPECT_EQ(row.candidate_id, six::candidate_id(catalogue.candidate(i)));
    if (row.line)
    {
      ++accepted;
      EXPECT_EQ(row.line->state_id, i + 1);
      EXPECT_EQ(row.line->momentum_index, 3u);
      EXPECT_GT(row.line->weight, 0.);
      m0 += row.line->weight;
      m1 += row.line->weight * row.line->gap;
    }
    else
    {
      ++rejected;
      EXPECT_FALSE(row.diagnostic.empty());
      EXPECT_TRUE(!row.solved.converged() || !row.factor || !row.factor->converged());
    }
  }
  EXPECT_GT(rejected, 0u);
  EXPECT_EQ(accepted, result.lines.size());
  EXPECT_EQ(accepted + rejected, catalogue.size());
  EXPECT_NEAR(result.moments[3].weight, m0, 1e-20);
  EXPECT_NEAR(result.moments[3].first_moment, m1, 1e-20);
}

#if UNI20_HAS_MPREAL
TEST(ExperimentalSixScanMpfr, ExplicitPrecisionInWorkers)
{
  auto const ambient = mpfr_get_default_prec();
  for (int bits : {128, 256})
    numerical_scan(uni20::mpreal(0, uni20::Precision::bits(bits)));
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
#endif
