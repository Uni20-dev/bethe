// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/experimental/spinon_observation_store.hpp>
#include <condition_variable>
#include <fstream>
#include <oneapi/tbb/global_control.h>
#include <sys/wait.h>
#include <thread>
#include <uni20/async/tbb_scheduler.hpp>

namespace x = bethe::experimental::xxx;
namespace h = bethe::heisenberg;
namespace cd = x::checkpoint_detail;
namespace
{
struct Scratch
{
    std::filesystem::path path;
    Scratch()
    {
      auto name = (std::filesystem::temp_directory_path() / "bethe-observations-XXXXXX").string();
      auto result = ::mkdtemp(name.data());
      if (!result) throw std::runtime_error("test scratch creation failed");
      path = result;
    }
    ~Scratch()
    {
      std::error_code ignored;
      std::filesystem::remove_all(path, ignored);
    }
};
using Store = x::SpinonObservationStore<long double, double, long double>;
bool native_pair() { return uni20::numeric_limits<long double>::digits() > 53; }
template <uni20::Real G> auto group(G const& one, bool zero_only = false)
{
  using Group = x::SixSpinonObservationGroup<G>;
  G const zero = uni20::scalar_like(one, 0), pi = 4 * uni20::math::atan(one);
  typename Group::Proposal proposal(8 * one, one / 2,
                                    {{h::detail::SpinonProposalKind::background, one / 2, one, one},
                                     {h::detail::SpinonProposalKind::one_pair, one / 4, 4 * one, one / 2},
                                     {h::detail::SpinonProposalKind::two_pairs, one / 4, 8 * one, one}});
  std::vector<x::SpinonSectionRequest<G>> requests;
  for (int s = 0; s < 3; ++s)
  {
    G const omega = zero_only ? 20 * one : 2 * one;
    requests.push_back({"whole-" + std::to_string(s), pi / 2, omega, zero, 8 * one, s});
    requests.push_back({"compact-" + std::to_string(s), pi / 2, omega, zero, 4 * one, s});
    requests.push_back({"remainder-" + std::to_string(s), pi / 2, omega, 4 * one, 8 * one, s});
  }
  return Group(std::move(proposal), one, std::move(requests));
}
auto policy()
{
  return h::SpinonRecoveryPolicy{
      "1e-5", "0", {{53, uni20::numeric_limits<long double>::digits(), {{{32, 12}, {64, 20}}}}}, 60};
}
std::vector<x::SpinonObservationIdentity> samples(std::size_t count = 12)
{
  bethe::MonteCarloStream<5> stream(42, 0, 24);
  std::vector<x::SpinonObservationIdentity> result;
  for (std::size_t i = 0; i < count; ++i)
    result.push_back({"sample-" + std::to_string(i), stream.next()});
  return result;
}
void same(Store::Observation const& a, Store::Observation const& b)
{
  EXPECT_EQ(a.values, b.values);
  EXPECT_EQ(a.cancelled, b.cancelled);
  EXPECT_EQ(a.reason, b.reason);
  EXPECT_EQ(a.metadata.failure, b.metadata.failure);
  EXPECT_EQ(a.metadata.density_attempts, b.metadata.density_attempts);
  EXPECT_EQ(a.metadata.supported_configurations, b.metadata.supported_configurations);
}
std::size_t selected(Store const& store)
{
  for (std::size_t i = 0; i < store.size(); ++i)
    for (auto const& row : store.prepared(i).rows())
      if (row.support == x::SpinonSectionSupport::selected) return i;
  throw std::runtime_error("fixture has no selected roots");
}
void replace(std::filesystem::path const& p, std::string const& data)
{
  std::ofstream out(p, std::ios::binary | std::ios::trunc);
  out.write(data.data(), static_cast<std::streamsize>(data.size()));
  if (!out) throw std::runtime_error("test replacement failed");
}

// Exercise the production chunk executor with a held recovery call, without
// depending on the duration of a numerical kernel or filesystem operation.
struct HeldChunk
{
    std::mutex mutex;
    std::condition_variable changed;
    std::size_t held_index = 0, controls = 0;
    bool started = false, stopped = false, release = false, exited = false;
    std::size_t size() const { return 2; }
    int run(std::size_t index, std::stop_token stop)
    {
      if (index != held_index)
      {
        std::unique_lock lock(mutex);
        if (!changed.wait_for(lock, std::chrono::seconds(5), [&] { return started; }))
          throw std::runtime_error("sibling recovery did not start");
        return 0;
      }
      std::stop_callback callback(stop, [&] {
        std::lock_guard lock(mutex);
        stopped = true;
        changed.notify_all();
      });
      std::unique_lock lock(mutex);
      started = true;
      ++controls;
      changed.notify_all();
      if (!changed.wait_for(lock, std::chrono::seconds(10), [&] { return stopped || release; }))
        throw std::runtime_error("held recovery did not receive stop or test release");
      if (!stopped) ++controls; // A second control must not start after grid failure.
      exited = true;
      return 0;
    }
};

class RejectChunkSibling : public uni20::async::TbbScheduler {
    std::shared_ptr<HeldChunk> chunk_;
    bool reject_;
    std::atomic<std::size_t> calls_ = 0;

  public:
    RejectChunkSibling(std::shared_ptr<HeldChunk> chunk, bool reject)
        : TbbScheduler(3), chunk_(std::move(chunk)), reject_(reject)
    {}
    void schedule(uni20::async::AsyncTask&& task) override
    {
      if (++calls_ == 3 && reject_)
      {
        std::unique_lock lock(chunk_->mutex);
        if (!chunk_->changed.wait_for(lock, std::chrono::seconds(5), [&] { return chunk_->started; }))
          throw std::runtime_error("first recovery did not start");
        throw std::bad_alloc{};
      }
      TbbScheduler::schedule(std::move(task));
    }
};

void check_chunk_failure_stops_recovery(bool reject)
{
  oneapi::tbb::global_control workers(oneapi::tbb::global_control::max_allowed_parallelism, 3);
  auto chunk = std::make_shared<HeldChunk>();
  chunk->held_index = reject ? 0 : 1;
  RejectChunkSibling scheduler(chunk, reject);
  std::stop_source external;
  auto pending = x::observation_store_detail::evaluate_chunk_async(
      chunk, [](auto, auto) { throw std::runtime_error("collector failed"); },
      {.scheduler = &scheduler, .max_in_flight = 2, .stop = external.get_token()});
  {
    std::unique_lock lock(chunk->mutex);
    EXPECT_TRUE(chunk->changed.wait_for(lock, std::chrono::seconds(5), [&] { return chunk->stopped; }));
    chunk->release = true; // Bound an old-code failure without orphaning held work.
    chunk->changed.notify_all();
  }
  if (reject)
    EXPECT_THROW((void)pending.get_wait(scheduler), std::bad_alloc);
  else
    EXPECT_THROW((void)pending.get_wait(scheduler), std::runtime_error);
  EXPECT_TRUE(chunk->exited);
  EXPECT_EQ(chunk->controls, 1u);
  EXPECT_FALSE(external.stop_requested());
}
} // namespace

TEST(ObservationStore, CollectorFailureStopsRunningRecovery) { check_chunk_failure_stops_recovery(false); }

TEST(ObservationStore, AdmissionFailureStopsRunningRecovery) { check_chunk_failure_stops_recovery(true); }

TEST(ObservationStore, ExactDecimalsPreserveCrossPrecisionTies)
{
  if (!native_pair()) GTEST_SKIP();
  bool shortest_lost_tie = false;
  for (int j = 1; j < 32; j += 2)
    for (long double sign : {-1.L, 1.L})
    {
      long double const value = sign * (1.L + j * std::ldexp(1.L, -53));
      double const expected = static_cast<double>(value);
      EXPECT_EQ(uni20::parse_real<double>(cd::exact_decimal(value)), expected);
      shortest_lost_tie |= uni20::parse_real<double>(uni20::format_real(value)) != expected;
    }
  EXPECT_TRUE(shortest_lost_tie); // A round-trip-at-G implementation is distinguishable.
}

TEST(ObservationStore, FreshReplayBothRootsMasksAndOriginalDenominator)
{
  if (!native_pair()) GTEST_SKIP();
  Scratch scratch;
  auto model = group(1.L);
  h::XXXSpinonDensityRecovery<double, long double> recovery(policy());
  std::vector<Store::Observation> expected;
  std::size_t controls = 0, accepted = 0, failures = 0;
  {
    Store store(scratch.path, model, samples(), policy(), "test-v1", "fresh uniforms seed42/r0/24bits");
    EXPECT_THROW((Store(scratch.path, model, samples(), policy(), "test-v1", "fresh uniforms seed42/r0/24bits")),
                 std::system_error);
    for (std::size_t i = 0; i < store.size(); ++i)
    {
      EXPECT_FALSE(store.load(i));
      auto actual = store.run(i, [&](auto, auto, auto const&) { ++controls; });
      auto direct = model.observation(store.prepared(i), recovery);
      same(actual, direct);
      if (actual.values)
      {
        ++accepted;
        for (std::size_t field = 0; field < 3; ++field)
          for (std::size_t section = 0; section < 3; ++section)
          {
            auto j = 9 * field + 3 * section;
            EXPECT_EQ((*actual.values)[j], (*actual.values)[j + 1] + (*actual.values)[j + 2]);
          }
      }
      else
        ++failures;
      same(*store.load(i), actual);
      expected.push_back(std::move(actual));
    }
  }
  EXPECT_GT(controls, 0u);
  EXPECT_GT(accepted, 0u);
  EXPECT_EQ(accepted + failures, samples().size());
  {
    Store store(scratch.path, model, samples(), policy(), "test-v1", "fresh uniforms seed42/r0/24bits");
    std::size_t replay_controls = 0;
    for (std::size_t i = 0; i < store.size(); ++i)
      same(store.run(i, [&](auto, auto, auto const&) { ++replay_controls; }), expected[i]);
    EXPECT_EQ(replay_controls, 0u);
  }
  std::size_t files = 0;
  for (auto const& entry : std::filesystem::recursive_directory_iterator(scratch.path))
    files += entry.is_regular_file();
  EXPECT_EQ(files, 4u); // Two stable locks and two chunk snapshots, not files per sample/root.
}

TEST(ObservationStore, ProcessDeathResumesOnlyMissingControls)
{
  if (!native_pair()) GTEST_SKIP();
  Scratch scratch, fresh;
  std::size_t index;
  {
    Store store(scratch.path, group(1.L), samples(), policy(), "v1", "fixture");
    index = selected(store);
  }
  auto pid = ::fork();
  ASSERT_GE(pid, 0);
  if (!pid)
  {
    Store store(scratch.path, group(1.L), samples(), policy(), "v1", "fixture");
    store.run(index, [](auto, auto, auto const&) { ::_exit(17); });
    ::_exit(99);
  }
  int status = 0;
  ASSERT_EQ(::waitpid(pid, &status, 0), pid);
  ASSERT_TRUE(WIFEXITED(status));
  ASSERT_EQ(WEXITSTATUS(status), 17);
  Store resumed(scratch.path, group(1.L), samples(), policy(), "v1", "fixture");
  Store direct(fresh.path, group(1.L), samples(), policy(), "v1", "fixture");
  EXPECT_FALSE(resumed.load(index));
  std::size_t remaining = 0, all = 0;
  auto actual = resumed.run(index, [&](auto, auto, auto const&) { ++remaining; });
  auto expected = direct.run(index, [&](auto, auto, auto const&) { ++all; });
  same(actual, expected);
  EXPECT_EQ(remaining + 1, all);
}

TEST(ObservationStore, CancelledPartialAndFailedNumericsRemainDistinct)
{
  if (!native_pair()) GTEST_SKIP();
  Scratch scratch;
  std::size_t index;
  {
    Store store(scratch.path, group(1.L), samples(), policy(), "v1", "fixture");
    index = selected(store);
    std::stop_source stop;
    auto value = store.run(index, [&](auto, auto, auto const&) { stop.request_stop(); }, stop.get_token());
    EXPECT_TRUE(value.cancelled);
    EXPECT_FALSE(value.values);
    EXPECT_FALSE(store.load(index));
  }
  Store store(scratch.path, group(1.L), samples(), policy(), "v1", "fixture");
  std::size_t controls = 0;
  auto result = store.run(index, [&](auto, auto, auto const&) { ++controls; });
  EXPECT_FALSE(result.cancelled);
  EXPECT_EQ(result.metadata.density_attempts, controls + 1);
  Scratch failed;
  auto strict = policy();
  strict.relative_tolerance = "1e-100";
  Store rejected(failed.path, group(1.L), samples(), strict, "v1", "fixture");
  auto no_value = rejected.run(selected(rejected));
  EXPECT_FALSE(no_value.values);
  EXPECT_FALSE(no_value.cancelled);
  EXPECT_EQ(no_value.metadata.failure, h::SpinonObservationFailure::density);
  ASSERT_TRUE(rejected.load(selected(rejected)));
  same(*rejected.load(selected(rejected)), no_value);
}

TEST(ObservationStore, BindingCorruptionAndPublicationFailure)
{
  if (!native_pair()) GTEST_SKIP();
  Scratch scratch;
  std::size_t index;
  {
    int publications = 0;
    Store store(scratch.path, group(1.L), samples(), policy(), "v1", "fixture", [&](cd::PublicationPhase phase) {
      if (phase == cd::PublicationPhase::before_replace && ++publications == 2)
        throw std::runtime_error("injected publication failure");
    });
    index = selected(store);
    EXPECT_THROW(store.run(index), std::runtime_error);
    EXPECT_THROW(store.load(index), std::runtime_error);
  }
  // The observation commit failed; density controls already persisted safely.
  {
    Store restored(scratch.path, group(1.L), samples(), policy(), "v1", "fixture");
    EXPECT_FALSE(restored.load(index));
    std::size_t calls = 0;
    auto result = restored.run(index, [&](auto, auto, auto const&) { ++calls; });
    EXPECT_EQ(calls, 0u);
    ASSERT_TRUE(restored.load(index));
  }
  auto changed = samples();
  std::swap(changed[0], changed[1]);
  EXPECT_THROW((Store(scratch.path, group(1.L), changed, policy(), "v1", "fixture")), std::invalid_argument);
  EXPECT_THROW((Store(scratch.path, group(1.L), samples(), policy(), "v2", "fixture")), std::invalid_argument);
  EXPECT_THROW((Store(scratch.path, group(1.L, true), samples(), policy(), "v1", "fixture")), std::invalid_argument);
  auto bytes = *cd::read_snapshot(scratch.path / "snapshot.bin");
  replace(scratch.path / "snapshot.bin", bytes.substr(0, bytes.size() - 3));
  EXPECT_THROW((Store(scratch.path, group(1.L), samples(), policy(), "v1", "fixture")), std::invalid_argument);
}

TEST(ObservationStore, ZeroOnlyChunksAndOrderedAsyncReplayOwnTheirStorage)
{
  if (!native_pair()) GTEST_SKIP();
  Scratch scratch;
  using namespace uni20::async;
  std::vector<Store::Observation> original;
  for (std::size_t window : {1u, 4u})
  {
    DebugScheduler scheduler({.order = DebugSchedulerOrder::reverse});
    auto store = std::make_shared<Store>(scratch.path, group(1.L, true), samples(), policy(), "v1", "fixture");
    std::weak_ptr<Store> weak = store;
    std::vector<Store::Observation> results;
    auto pending = x::evaluate_spinon_observation_chunk_async(store,
                                                              [&](std::size_t i, Store::Observation value) {
                                                                EXPECT_EQ(i, results.size());
                                                                ASSERT_TRUE(value.values);
                                                                for (auto v : *value.values)
                                                                  EXPECT_EQ(v, 0);
                                                                results.push_back(std::move(value));
                                                              },
                                                              {.scheduler = &scheduler, .max_in_flight = window});
    store.reset();
    EXPECT_FALSE(weak.expired());
    EXPECT_EQ(pending.get_wait(scheduler), samples().size());
    EXPECT_TRUE(weak.expired());
    EXPECT_TRUE(scheduler.done());
    if (original.empty())
      original = results;
    else
      for (std::size_t i = 0; i < results.size(); ++i)
        same(results[i], original[i]);
  }
  EXPECT_FALSE(std::filesystem::exists(scratch.path / "controls"));
}

TEST(ObservationStore, SavedValuesAreRecomputedEvenWhenChecksumIsValid)
{
  if (!native_pair()) GTEST_SKIP();
  Scratch scratch;
  {
    Store store(scratch.path, group(1.L), samples(), policy(), "v1", "fixture");
    (void)store.run(selected(store));
  }
  auto bytes = *cd::read_snapshot(scratch.path / "snapshot.bin");
  cd::Decoder d{std::string_view(bytes).substr(0, bytes.size() - 8)};
  (void)d.text();
  (void)d.bytes(d.count(cd::max_bytes));
  auto count = d.count(cd::max_points);
  bool changed = false;
  for (std::size_t i = 0; i < count; ++i)
    if (d.flag())
    {
      auto size = d.count(cd::max_bytes);
      bytes[d.at] ^= 1; // Change the presence of the saved value; retain the rest.
      (void)d.bytes(size);
      changed = true;
      break;
    }
  ASSERT_TRUE(changed);
  cd::Encoder e;
  e.bytes(std::string_view(bytes).substr(0, bytes.size() - 8));
  replace(scratch.path / "snapshot.bin", e.finish());
  EXPECT_THROW((Store(scratch.path, group(1.L), samples(), policy(), "v1", "fixture")), std::invalid_argument);
}

TEST(ObservationStore, SectionUpperEdgeReplaysAsZeroWithoutDensityControls)
{
  if (!native_pair()) GTEST_SKIP();
  Scratch scratch;
  using Group = Store::Group;
  Group::Proposal proposal(8.L, .5L, {{h::detail::SpinonProposalKind::background, 1.L, 1.L, 1.L}});
  bethe::UniformSample<5>::Words words;
  words[0] = {0};
  for (int i = 1; i < 5; ++i)
    words[i] = {std::uint64_t{1} << 22};
  bethe::UniformSample<5> identity(24, words);
  auto draw = proposal.draw(identity.at_precision(1.L));
  long double const k = 12 * uni20::math::atan(uni20::math::exp(draw.beta[0]));
  long double const pi = 4 * uni20::math::atan(1.L), upper = 3 * pi * uni20::math::sin(k / 6);
  ASSERT_LT(k, 2 * pi);
  Group model(proposal, 1.L, {{"section-zero-edge", k, upper, 0.L, 8.L, 0}});
  for (int replay = 0; replay < 2; ++replay)
  {
    Store store(scratch.path, model, {{"edge", identity}}, policy(), "v1", "equal-momentum section edge");
    std::size_t calls = 0;
    auto value = store.run(0, [&](auto, auto, auto const&) { ++calls; });
    ASSERT_TRUE(value.values);
    EXPECT_EQ(*value.values, (std::vector<long double>{0, 0, 0}));
    EXPECT_EQ(calls, 0u);
    EXPECT_FALSE(std::filesystem::exists(scratch.path / "controls"));
  }
}

TEST(ObservationStore, PhysicalAsyncResultsSurviveReversedReplayOrder)
{
  if (!native_pair()) GTEST_SKIP();
  Scratch scratch;
  auto model = group(1.L);
  auto identities = samples();
  auto found = std::find_if(identities.begin(), identities.end(), [&](auto const& sample) {
    auto p = model.prepare(sample.id, sample.uniforms);
    return std::any_of(p.rows().begin(), p.rows().end(),
                       [](auto const& row) { return row.support == x::SpinonSectionSupport::selected; });
  });
  ASSERT_NE(found, identities.end());
  // Deliberately duplicate the coordinates under distinct observation IDs to
  // exercise concurrent distinct control histories. This is a replay fixture,
  // not a claim of independent Monte Carlo samples or a covariance estimate.
  std::vector<x::SpinonObservationIdentity> inputs{{"first", found->uniforms}, {"second", found->uniforms}};
  h::XXXSpinonDensityRecovery<double, long double> recovery(policy());
  auto expected = model.observation(model.prepare(inputs[0].id, inputs[0].uniforms), recovery);
  ASSERT_TRUE(expected.values);
  ASSERT_GT(expected.metadata.density_attempts, 0u);
  {
    uni20::async::TbbScheduler scheduler(2);
    auto store = std::make_shared<Store>(scratch.path, model, inputs, policy(), "v1", "correlated replay fixture");
    std::size_t collected = 0;
    auto pending = x::evaluate_spinon_observation_chunk_async(store,
                                                              [&](auto i, auto actual) {
                                                                EXPECT_EQ(i, collected++);
                                                                same(actual, expected);
                                                              },
                                                              {.scheduler = &scheduler, .max_in_flight = 2});
    EXPECT_EQ(pending.get_wait(scheduler), 2u);
    EXPECT_EQ(collected, 2u);
  }
  {
    uni20::async::DebugScheduler scheduler({.order = uni20::async::DebugSchedulerOrder::reverse});
    auto store = std::make_shared<Store>(scratch.path, model, inputs, policy(), "v1", "correlated replay fixture");
    std::size_t collected = 0;
    auto pending = x::evaluate_spinon_observation_chunk_async(store,
                                                              [&](auto i, auto actual) {
                                                                EXPECT_EQ(i, collected++);
                                                                same(actual, expected);
                                                              },
                                                              {.scheduler = &scheduler, .max_in_flight = 1});
    EXPECT_EQ(pending.get_wait(scheduler), 2u);
    EXPECT_EQ(collected, 2u);
  }
}

TEST(ObservationStore, AsyncSubmissionAndCollectorFailuresDrainOwnership)
{
  if (!native_pair()) GTEST_SKIP();
  using namespace uni20::async;
  struct Reject : DebugScheduler
  {
      int calls = 0, reject = 1;
      void schedule(AsyncTask&& task) override
      {
        if (++calls == reject) throw std::bad_alloc{};
        DebugScheduler::schedule(std::move(task));
      }
  };
  for (int reject : {1, 2, 4, 100})
  {
    Scratch scratch;
    Reject scheduler;
    scheduler.reject = reject;
    auto store = std::make_shared<Store>(scratch.path, group(1.L, true), samples(), policy(), "v1", "fixture");
    std::weak_ptr<Store> weak = store;
    auto start = [&] {
      return x::evaluate_spinon_observation_chunk_async(
          std::move(store), [](auto, auto) { throw std::runtime_error("collector failed"); },
          {.scheduler = &scheduler, .max_in_flight = 4});
    };
    if (reject == 1)
      EXPECT_THROW(start(), std::bad_alloc);
    else
    {
      auto pending = start();
      EXPECT_THROW((void)pending.get_wait(scheduler), std::exception);
    }
    EXPECT_TRUE(scheduler.done());
    EXPECT_TRUE(weak.expired());
    Store reopened(scratch.path, group(1.L, true), samples(), policy(), "v1", "fixture");
  }
}

TEST(ObservationStore, AsyncCancellationAndExceptionWaitForHeldPublication)
{
  if (!native_pair()) GTEST_SKIP();
  for (bool fail : {false, true})
  {
    Scratch scratch;
    std::mutex mutex;
    std::condition_variable changed;
    bool entered = false, release = false, armed = false;
    std::atomic<bool> exited = false;
    auto store = std::make_shared<Store>(
        scratch.path, group(1.L), samples(), policy(), "v1", "fixture", [&](cd::PublicationPhase phase) {
          std::unique_lock lock(mutex);
          if (!armed || phase != cd::PublicationPhase::after_replace || entered) return;
          entered = true;
          changed.notify_all();
          if (!changed.wait_for(lock, std::chrono::seconds(5), [&] { return release; }))
            throw std::runtime_error("held publication test timed out");
          exited = true;
          if (fail) throw std::runtime_error("publication completed but observer failed");
        });
    armed = true;
    uni20::async::TbbScheduler scheduler(2);
    std::stop_source stop;
    auto pending = x::evaluate_spinon_observation_chunk_async(
        store, [](auto, auto) {}, {.scheduler = &scheduler, .max_in_flight = 4, .stop = stop.get_token()});
    store.reset();
    {
      std::unique_lock lock(mutex);
      EXPECT_TRUE(changed.wait_for(lock, std::chrono::seconds(5), [&] { return entered; }));
      EXPECT_FALSE(pending.read().await_ready());
      if (!fail) stop.request_stop();
      release = true;
      changed.notify_all();
    }
    if (fail)
      EXPECT_THROW((void)pending.get_wait(scheduler), std::runtime_error);
    else
      EXPECT_THROW((void)pending.get_wait(scheduler), bethe::ExecutionCancelled);
    EXPECT_TRUE(exited);
    // Completed snapshots remain usable after cancellation or observer failure.
    scheduler.run_all(); // Test teardown only: final producer reference release.
    Store restored(scratch.path, group(1.L), samples(), policy(), "v1", "fixture");
  }
}

#if UNI20_HAS_MPREAL
TEST(ObservationStoreMpfr, ExactBinaryInputsAndPhysicalReplay)
{
  using R = uni20::mpreal;
  auto p = uni20::Precision::bits(160);
  R const one(1, p);
  for (R const value : {one / 3, uni20::math::ldexp(one, -400), one + uni20::math::ldexp(one, -79)})
  {
    auto exact = cd::exact_decimal(value);
    EXPECT_EQ(R(exact, p), value);
    for (int bits : {53, 80, 128})
    {
      h::XXXSpinonDensity<R> context(uni20::Precision::bits(bits));
      EXPECT_EQ(context.parse(exact), context.import(value));
    }
  }
  EXPECT_THROW(cd::exact_decimal(uni20::math::ldexp(one, -1000000)), std::invalid_argument);
  Scratch scratch;
  using MpStore = x::SpinonObservationStore<R, R>;
  auto model = group(one);
  h::SpinonRecoveryPolicy configuration{"1e-5", "0", {{80, 128, {{{32, 12}, {64, 20}}}}}, 60};
  std::vector<MpStore::Observation> expected;
  std::size_t attempted = 0;
  h::XXXSpinonDensityRecovery<R> recovery(configuration);
  {
    MpStore store(scratch.path, model, samples(12), configuration, "mpfr-v1", "fresh fixture");
    for (std::size_t i = 0; i < store.size(); ++i)
    {
      auto actual = store.run(i, [&](auto, auto, auto const&) { ++attempted; }),
           direct = model.observation(store.prepared(i), recovery);
      EXPECT_EQ(actual.values, direct.values);
      EXPECT_EQ(actual.metadata.failure, direct.metadata.failure);
      expected.push_back(std::move(actual));
    }
  }
  EXPECT_GT(attempted, 0u);
  MpStore store(scratch.path, model, samples(12), configuration, "mpfr-v1", "fresh fixture");
  std::size_t calls = 0;
  for (std::size_t i = 0; i < store.size(); ++i)
    EXPECT_EQ(store.run(i, [&](auto, auto, auto const&) { ++calls; }).values, expected[i].values);
  EXPECT_EQ(calls, 0u);
}
#endif
