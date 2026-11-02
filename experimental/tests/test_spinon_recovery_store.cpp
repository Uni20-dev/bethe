// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/detail/continuous_work.hpp>
#include <bethe/execution.hpp>
#include <bethe/experimental/spinon_recovery_store.hpp>
#include <csignal>
#include <fstream>
#include <sys/wait.h>
#include <thread>
#include <uni20/async/tbb_scheduler.hpp>

namespace x = bethe::experimental::xxx;
namespace cd = x::checkpoint_detail;
using namespace bethe::heisenberg;

namespace
{
struct Scratch
{
    std::filesystem::path path;
    Scratch()
    {
      auto name = (std::filesystem::temp_directory_path() / "bethe-recovery-XXXXXX").string();
      auto result = ::mkdtemp(name.data());
      if (!result) throw std::runtime_error("cannot create test scratch directory");
      path = result;
    }
    ~Scratch()
    {
      std::error_code error;
      std::filesystem::remove_all(path, error);
    }
};
struct StopAfterControl
{};
using Store = x::SpinonRecoveryStore<double, long double>;
bool native_pair() { return uni20::numeric_limits<long double>::digits() > 53; }
x::RecoveryChunkBinding binding()
{
  return {
      "test-producer-v1",
      "fresh distinct points; longitudinal labelled dp; no sample redraw",
      {"1e-6", "0", {{53, uni20::numeric_limits<long double>::digits(), {{{32, 12}, {64, 20}}}}}, 60},
      {{"two", {"-2", "1"}}, {"four", {"-2.1", "-.4", ".3", "2"}}, {"six", {"-2.1", "-1.3", "-.4", ".3", "1.2", "2"}}}};
}
void replace(std::filesystem::path const& p, std::string const& data)
{
  std::ofstream out(p, std::ios::binary | std::ios::trunc);
  out.write(data.data(), static_cast<std::streamsize>(data.size()));
  if (!out) throw std::runtime_error("test replacement failed");
}
template <uni20::Real R> void scalars(R const& zero)
{
  int const bits = uni20::precision_of(zero).bit_count();
  for (R const& value : {zero, -zero, uni20::scalar_like(zero, 1) / 3, uni20::scalar_like(zero, -2)})
  {
    cd::Encoder e;
    e.real(value);
    cd::Decoder d{e.data};
    auto restored = d.real<R>(bits);
    EXPECT_EQ(uni20::format_real(restored, {.normalize_negative_zero = false}),
              uni20::format_real(value, {.normalize_negative_zero = false}));
    EXPECT_EQ(uni20::precision_of(restored), uni20::precision_of(value));
    EXPECT_EQ(d.at, d.data.size());
  }
}
} // namespace

template <typename R> class RecoveryScalarCodec : public ::testing::Test {};
TYPED_TEST_SUITE(RecoveryScalarCodec, test_support::RealTypes, test_support::PrecisionNames);
TYPED_TEST(RecoveryScalarCodec, NativeExponentRangeAndSubnormalRoundTrip)
{
  using R = TypeParam;
  scalars(R{0});
  for (R const value : {uni20::numeric_limits<R>::max(), uni20::numeric_limits<R>::min(),
                        uni20::numeric_limits<R>::denorm_min(), -uni20::numeric_limits<R>::denorm_min()})
  {
    cd::Encoder e;
    e.real(value);
    cd::Decoder d{e.data};
    EXPECT_EQ(d.real<R>(uni20::numeric_limits<R>::digits()), value);
  }
}

TEST(RecoveryStore, PartialCompletedAndAllPolarizationRoundTrips)
{
  if (!native_pair()) GTEST_SKIP();
  Scratch scratch;
  auto b = binding();
  {
    Store store(scratch.path, b);
    EXPECT_EQ(store.size(), 3u);
    EXPECT_FALSE(store.load(0));
    EXPECT_THROW(store.run(0,
                           [&](auto const& progress) {
                             EXPECT_EQ(progress.attempts.size(), 1u);
                             EXPECT_EQ(store.load(0)->attempts.size(), 1u);
                             EXPECT_THROW(store.run(0), std::invalid_argument);
                             throw StopAfterControl{};
                           }),
                 StopAfterControl);
  }
  {
    Store store(scratch.path, b);
    ASSERT_EQ(store.load(0)->attempts.size(), 1u);
    int controls = 0;
    auto done = store.run(0, [&](auto const&) { ++controls; });
    ASSERT_TRUE(done.accepted());
    EXPECT_EQ(controls, 1);
    auto replay = store.run(0, [](auto const&) { ADD_FAILURE() << "completed point evaluated"; });
    EXPECT_EQ(replay.estimate->value, done.estimate->value);
    for (std::size_t i = 1; i < store.size(); ++i)
      ASSERT_TRUE(store.run(i).accepted());
  }
  Store reopened(scratch.path, b);
  Store::Driver direct(b.policy);
  for (std::size_t i = 0; i < b.points.size(); ++i)
  {
    auto result = *reopened.load(i);
    auto fresh = direct(b.points[i].beta);
    ASSERT_TRUE(result.accepted() && fresh.accepted());
    EXPECT_EQ(result.estimate->value, fresh.estimate->value);
    EXPECT_EQ(result.estimate->amplitudes, fresh.estimate->amplitudes);
  }
  EXPECT_LE(std::distance(std::filesystem::directory_iterator(scratch.path), std::filesystem::directory_iterator{}), 3);
}

TEST(RecoveryStore, BindingOwnershipAndCorruption)
{
  if (!native_pair()) GTEST_SKIP();
  Scratch scratch;
  auto b = binding();
  {
    Store first(scratch.path, b);
    EXPECT_THROW(Store(scratch.path, b), std::system_error);
    ASSERT_TRUE(first.run(0).accepted());
  }
  auto reject = [&](auto mutate) {
    auto changed = b;
    mutate(changed);
    EXPECT_THROW(Store(scratch.path, changed), std::invalid_argument);
  };
  reject([](auto& a) { a.producer += "new"; });
  reject([](auto& a) { a.provenance += "changed measure"; });
  reject([](auto& a) { a.points[0].beta[0] = "-2.0"; });
  reject([](auto& a) { std::swap(a.points[0], a.points[1]); });
  reject([](auto& a) { a.points[1].id = a.points[0].id; });
  reject([](auto& a) { a.policy.max_wall_seconds *= 2; });
  reject([](auto& a) { a.policy.relative_tolerance = "1e-3"; });
  auto saved = *cd::read_snapshot(scratch.path / "snapshot.bin");
  for (auto damaged : {saved.substr(0, saved.size() - 1), saved + "trailing", std::string("bad")})
  {
    replace(scratch.path / "snapshot.bin", damaged);
    EXPECT_THROW(Store(scratch.path, b), std::invalid_argument);
  }
  auto damaged = saved;
  damaged[damaged.size() / 2] ^= 1;
  replace(scratch.path / "snapshot.bin", damaged);
  EXPECT_THROW(Store(scratch.path, b), std::invalid_argument);
  replace(scratch.path / "snapshot.bin", saved);
  // A well-framed but internally inconsistent record must fail semantic
  // validation too. Reset consumed time without changing completion evidence.
  auto reset = saved.substr(0, saved.size() - 8);
  cd::Decoder fields{reset};
  (void)fields.text();
  (void)fields.bytes(fields.count(cd::max_bytes));
  (void)fields.clock();
  (void)fields.integer(); // number of points
  (void)fields.integer(); // first point index
  (void)fields.clock();   // original start
  (void)fields.clock();   // original deadline
  (void)fields.flag();
  reset.replace(fields.at, 8, 8, '\0');
  cd::Encoder reframed;
  reframed.data = reset;
  replace(scratch.path / "snapshot.bin", reframed.finish());
  EXPECT_THROW(Store(scratch.path, b), std::invalid_argument);
  replace(scratch.path / "snapshot.bin", saved);
  Store restored(scratch.path, b);
  EXPECT_TRUE(restored.load(0)->accepted());
}

TEST(RecoveryStore, CodecRetainsEarlierContradictionAndFailedControls)
{
  if (!native_pair()) GTEST_SKIP();
  auto p = binding().policy;
  auto next = p.stages[0];
  next.contours = {{{96, 22}, {128, 24}}};
  p.stages.push_back(next);
  Store::Driver driver(p);
  int calls = 0;
  // Deliberately inconsistent factory only for testing retained evidence.
  auto history = driver.recover([&](auto const& c) {
    ++calls;
    if (calls == 3) throw SpinonDensityError(SpinonDensityFailure::tail_domain, "injected failure");
    return std::vector{c.parse("-2"), c.parse(calls > 3 ? "2" : "1")};
  });
  ASSERT_EQ(history.status, SpinonRecoveryStatus::inconsistent_controls);
  cd::Encoder e;
  for (auto const& attempt : history.attempts)
    cd::encode_attempt(e, attempt);
  cd::Decoder d{e.data};
  auto decoded = history;
  decoded.attempts.clear();
  while (d.at < d.data.size())
    decoded.attempts.push_back(cd::decode_attempt<long double>(d, p));
  auto result = driver.restore(decoded);
  EXPECT_EQ(result.status, SpinonRecoveryStatus::inconsistent_controls);
  EXPECT_EQ(result.conflicting_attempts, history.conflicting_attempts);
  EXPECT_EQ(result.attempts[1].failure, SpinonDensityFailure::tail_domain);
  EXPECT_FALSE(result.attempts[1].checked);
  EXPECT_FALSE(result.estimate);
}

TEST(RecoveryStore, CancellationAndOriginalDeadlineSurviveDowntime)
{
  if (!native_pair()) GTEST_SKIP();
  Scratch scratch;
  auto b = binding();
  b.policy.max_wall_seconds = 1;
  {
    Store store(scratch.path, b);
    std::stop_source stop;
    auto r = store.run(0, [&](auto const&) { stop.request_stop(); }, stop.get_token());
    ASSERT_EQ(r.status, SpinonRecoveryStatus::cancelled);
    ASSERT_EQ(r.attempts.size(), 1u);
    ASSERT_TRUE(store.run(1).accepted());
  }
  std::this_thread::sleep_for(std::chrono::milliseconds(1100));
  Store store(scratch.path, b);
  auto expired = store.run(0, [](auto const&) { ADD_FAILURE() << "expired point ran a control"; });
  EXPECT_EQ(expired.status, SpinonRecoveryStatus::time_limit);
  EXPECT_EQ(expired.attempts.size(), 1u);
  EXPECT_FALSE(expired.estimate);
  EXPECT_GE(expired.wall_seconds, 1);
  EXPECT_TRUE(store.run(1, [](auto const&) { ADD_FAILURE() << "completed result reran"; }).accepted());
}

TEST(RecoveryStore, PublicationFailurePoisonsOwnerAndRetainsPreviousPrefix)
{
  if (!native_pair()) GTEST_SKIP();
  Scratch scratch;
  auto b = binding();
  {
    Store store(scratch.path, b);
    EXPECT_THROW(store.run(0,
                           [&](auto const& progress) {
                             if (progress.attempts.size() == 1)
                               std::filesystem::create_directory(scratch.path / "snapshot.pending");
                           }),
                 std::system_error);
    EXPECT_THROW(store.load(0), std::runtime_error);
    EXPECT_THROW(store.run(1), std::runtime_error);
  }
  std::filesystem::remove(scratch.path / "snapshot.pending");
  Store reopened(scratch.path, b);
  ASSERT_EQ(reopened.load(0)->attempts.size(), 1u);
  int controls = 0;
  EXPECT_TRUE(reopened.run(0, [&](auto const&) { ++controls; }).accepted());
  EXPECT_EQ(controls, 1);
}

TEST(RecoveryStore, InitialPublicationTimeCannotBuyExtraKernelBudget)
{
  if (!native_pair()) GTEST_SKIP();
  Scratch scratch;
  auto b = binding();
  b.policy.max_wall_seconds = 1;
  int publications = 0;
  Store store(scratch.path, b, [&](auto phase) {
    if (phase == cd::PublicationPhase::before_replace && ++publications == 2)
      std::this_thread::sleep_for(std::chrono::milliseconds(1100));
  });
  auto result = store.run(0, [](auto const&) { ADD_FAILURE() << "expired admission ran a kernel"; });
  EXPECT_EQ(result.status, SpinonRecoveryStatus::time_limit);
  EXPECT_TRUE(result.attempts.empty());
}

TEST(RecoveryStore, ProcessDeathBeforeAndAfterAtomicReplacement)
{
  if (!native_pair()) GTEST_SKIP();
  for (auto phase : {cd::PublicationPhase::before_replace, cd::PublicationPhase::after_replace})
  {
    Scratch scratch;
    auto b = binding();
    {
      Store initial(scratch.path, b);
    }
    pid_t pid = ::fork();
    ASSERT_GE(pid, 0);
    if (pid == 0)
    {
      int publications = 0;
      Store child(scratch.path, b, [&](auto event) {
        if (event == cd::PublicationPhase::before_replace) ++publications;
        if (publications == 2 && event == phase) ::_exit(42);
      });
      (void)child.run(0);
      ::_exit(1);
    }
    int status = 0;
    auto const deadline = std::chrono::steady_clock::now() + std::chrono::seconds(10);
    pid_t observed = 0;
    while ((observed = ::waitpid(pid, &status, WNOHANG)) == 0 && std::chrono::steady_clock::now() < deadline)
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
    if (observed == 0)
    {
      ::kill(pid, SIGKILL);
      (void)::waitpid(pid, &status, 0);
      FAIL() << "checkpoint child exceeded bounded wait";
    }
    ASSERT_EQ(observed, pid);
    ASSERT_TRUE(WIFEXITED(status));
    ASSERT_EQ(WEXITSTATUS(status), 42);
    Store recovered(scratch.path, b);
    auto prefix = recovered.load(0);
    ASSERT_TRUE(prefix);
    std::size_t retained = phase == cd::PublicationPhase::before_replace ? 0 : 1;
    EXPECT_EQ(prefix->attempts.size(), retained);
    int controls = 0;
    EXPECT_TRUE(recovered.run(0, [&](auto const&) { ++controls; }).accepted());
    EXPECT_EQ(controls, 2 - retained);
  }
}

TEST(RecoveryStore, ConcurrentPointsRetainIdentityAndDeterministicResults)
{
  if (!native_pair()) GTEST_SKIP();
  Scratch scratch;
  auto b = binding();
  {
    Store store(scratch.path, b);
    uni20::async::TbbScheduler scheduler(2);
    EXPECT_EQ(bethe::detail::continuous_work(
                  scheduler, 2, store.size(),
                  [&](std::size_t i, std::stop_token stop) { EXPECT_TRUE(store.run(i, stop).accepted()); }),
              store.size());
  }
  Store store(scratch.path, b);
  Store::Driver direct(b.policy);
  for (std::size_t i = 0; i < store.size(); ++i)
  {
    EXPECT_EQ(store.point(i).id, b.points[i].id);
    EXPECT_EQ(store.run(i).estimate->value, direct(b.points[i].beta).estimate->value);
  }
}

#if UNI20_HAS_MPREAL
TEST(RecoveryStoreMpfr, MixedPrecisionControlEvidenceRoundTrips)
{
  using R = uni20::mpreal;
  Scratch scratch;
  auto b = binding();
  b.points.resize(1);
  b.policy = {"1e-24", "0", {{32, 64, {{{32, 12}, {64, 20}}}}, {160, 256, {{{64, 20}, {128, 24}}}}}, 60};
  using S = x::SpinonRecoveryStore<R>;
  auto ambient = mpfr_get_default_prec();
  {
    S store(scratch.path, b);
    EXPECT_THROW(store.run(0,
                           [](auto const& progress) {
                             if (progress.attempts.size() == 2) throw StopAfterControl{};
                           }),
                 StopAfterControl);
    EXPECT_EQ(store.load(0)->attempts.size(), 2u);
  }
  {
    S store(scratch.path, b);
    int controls = 0;
    auto result = store.run(0, [&](auto const&) { ++controls; });
    ASSERT_TRUE(result.accepted()) << result.reason;
    EXPECT_EQ(controls, 2);
    ASSERT_TRUE(result.attempts[0].checked);
    EXPECT_EQ(result.attempts[0].checked->estimate.value.precision(), uni20::Precision::bits(64));
    EXPECT_EQ(result.estimate->value.precision(), uni20::Precision::bits(256));
  }
  S restored(scratch.path, b);
  auto final = restored.run(0);
  S::Driver direct(b.policy);
  EXPECT_EQ(final.estimate->value, direct(b.points[0].beta).estimate->value);
  for (int bits : {80, 128, 160, 256})
    scalars(R(0, uni20::Precision::bits(bits)));
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
#endif
