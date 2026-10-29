// SPDX-License-Identifier: GPL-3.0-or-later
#include "test_support.hpp"
#include <array>
#include <bethe/central_spin.hpp>
#include <bethe/lieb_liniger_open.hpp>
#include <bethe/observable_grid.hpp>
#include <bethe/richardson.hpp>
#include <tuple>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
namespace ll = bethe::lieb_liniger;
namespace pairing = bethe::richardson;
namespace central = bethe::central_spin;
using namespace uni20::async;

template <class Options> Options controls(uni20::Precision p)
{
  if constexpr (std::is_constructible_v<Options, uni20::Precision>)
    return Options(p);
  else
    return Options{};
}
template <typename Real> void finite(Real const& value, uni20::Precision p)
{
  EXPECT_TRUE(uni20::isfinite(value));
  if constexpr (uni20::has_runtime_precision_v<Real>)
  {
    EXPECT_FALSE(value.is_exact());
    EXPECT_EQ(value.precision(), p);
  }
}
template <typename State> void same_ll_state(State const& a, State const& b, uni20::Precision p)
{
  EXPECT_EQ(std::tie(a.particles, a.length, a.interaction, a.quantum_numbers, a.momenta, a.energy, a.residual_norm,
                     a.iterations, a.converged, a.status),
            std::tie(b.particles, b.length, b.interaction, b.quantum_numbers, b.momenta, b.energy, b.residual_norm,
                     b.iterations, b.converged, b.status));
  for (auto const& v : {a.length, a.interaction, a.energy, a.residual_norm})
    finite(v, p);
  for (auto const& v : a.momenta)
    finite(v, p);
  if constexpr (requires { a.momentum; })
  {
    EXPECT_EQ(a.momentum, b.momentum);
    EXPECT_EQ(a.momentum_index, b.momentum_index);
    finite(a.momentum, p);
  }
}
template <typename Scan> void same_scan(Scan const& a, Scan const& b, uni20::Precision p)
{
  same_ll_state(a.ground_state, b.ground_state, p);
  EXPECT_EQ(a.candidate_count, b.candidate_count);
  EXPECT_EQ(a.converged_count, b.converged_count);
  EXPECT_EQ(a.converged(), b.converged());
  ASSERT_EQ(a.levels.size(), b.levels.size());
  for (std::size_t i = 0; i < a.levels.size(); ++i)
  {
    same_ll_state(a.levels[i].state, b.levels[i].state, p);
    EXPECT_EQ(a.levels[i].gap, b.levels[i].gap);
    if (a.levels[i].gap) finite(*a.levels[i].gap, p);
  }
  ASSERT_EQ(bool(a.first_unconverged), bool(b.first_unconverged));
  if (a.first_unconverged) same_ll_state(*a.first_unconverged, *b.first_unconverged, p);
  if constexpr (requires { a.padding; })
  {
    EXPECT_EQ(std::tie(a.particles, a.padding, a.length, a.interaction),
              std::tie(b.particles, b.padding, b.length, b.interaction));
    finite(a.length, p);
    finite(a.interaction, p);
  }
}

template <typename Real> using Scans = std::pair<ll::ExcitationScan<Real>, ll::open::ExcitationScan<Real>>;
template <typename Real>
AsyncTask join_scans(ReadBuffer<ll::ExcitationScan<Real>> periodic, ReadBuffer<ll::open::ExcitationScan<Real>> walls,
                     WriteBuffer<Scans<Real>> output)
{
  auto destination = co_await output;
  auto const& a = co_await periodic;
  auto const& b = co_await walls;
  destination = Scans<Real>{a, b};
}
template <typename Real> auto submit_scans(IAsyncScheduler& scheduler, uni20::Precision p, bool fail)
{
  // All caller-owned model parameters/options and the two intermediate handles
  // disappear before the final consumer is observed (before execution on DebugScheduler).
  auto length = uni20::scalar_at_precision<Real>(3, p), c = uni20::scalar_at_precision<Real>(1, p);
  auto options = controls<ll::SolverOptions<Real>>(p);
  if (fail) options.max_iterations = 0;
  bethe::RealExcitationOptions scan{20, 100, {&scheduler, 128, 2}};
  auto periodic = ll::real_excitations_async(3, length, c, 1, scan, options);
  auto walls = ll::open::real_excitations_async(3, length, c, 1, scan, options);
  Async<Scans<Real>> output;
  scheduler.schedule(join_scans<Real>(periodic.read(), walls.read(), output.write()));
  length = uni20::scalar_at_precision<Real>(0, p);
  c = length;
  options.max_iterations = fail ? 10000 : 0;
  scan.count = 1;
  return output;
}
template <typename Real> auto reference_scans(uni20::Precision p, bool fail)
{
  auto options = controls<ll::SolverOptions<Real>>(p);
  if (fail) options.max_iterations = 0;
  auto length = uni20::scalar_at_precision<Real>(3, p), c = uni20::scalar_at_precision<Real>(1, p);
  Scans<Real> result{ll::real_excitations(3, length, c, 1, {20, 100}, options),
                     ll::open::real_excitations(3, length, c, 1, {20, 100}, options)};
  EXPECT_EQ(result.first.converged(), !fail);
  EXPECT_EQ(result.second.converged(), !fail);
  return result;
}
template <typename Real> void check_scans(Scans<Real> const& a, Scans<Real> const& b, uni20::Precision p)
{
  same_scan(a.first, b.first, p);
  same_scan(a.second, b.second, p);
}

template <typename Real> class SolverComposition : public ::testing::Test {};
TYPED_TEST_SUITE(SolverComposition, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(SolverComposition, FiniteLiebLinigerOwnsInputsAndComposesThroughFinalConsumer)
{
  using Real = TypeParam;
  auto const p = uni20::native_precision<Real>();
  for (bool fail : {false, true})
  {
    auto reference = reference_scans<Real>(p, fail);
    for (auto order : {DebugSchedulerOrder::fifo, DebugSchedulerOrder::reverse, DebugSchedulerOrder::random})
    {
      DebugScheduler scheduler({.order = order, .random_seed = 17});
      auto pending = submit_scans<Real>(scheduler, p, fail);
      check_scans(pending.get_wait(scheduler), reference, p);
    }
    for (int workers : {1, 4})
    {
      TbbScheduler scheduler(workers);
      auto pending = submit_scans<Real>(scheduler, p, fail);
      check_scans(pending.get_wait(scheduler), reference, p);
    }
  }
}

TEST(SolverComposition, FiniteLiebLinigerPreflightAndCancellationAreIsolated)
{
  DebugScheduler scheduler;
  bethe::RealExcitationOptions options{20, 100, {&scheduler, 128, 2}};
  EXPECT_THROW((void)ll::real_excitations_async(3, 0.0, 1.0, 1, options), std::invalid_argument);
  EXPECT_THROW((void)ll::open::real_excitations_async(3, 3.0, -1.0, 1, options), std::invalid_argument);
  options.max_candidates = 1;
  EXPECT_THROW((void)ll::real_excitations_async(3, 3.0, 1.0, 1, options), std::length_error);
  EXPECT_THROW((void)ll::open::real_excitations_async(3, 3.0, 1.0, 1, options), std::length_error);
  options.max_candidates = 100;
  options.execution.max_in_flight = 0;
  EXPECT_THROW((void)ll::real_excitations_async(3, 3.0, 1.0, 1, options), std::invalid_argument);
  EXPECT_THROW((void)ll::open::real_excitations_async(3, 3.0, 1.0, 1, options), std::invalid_argument);
  options.execution.max_in_flight = 2;
  std::stop_source stop;
  options.execution.stop = stop.get_token();
  auto a = ll::real_excitations_async(3, 3.0, 1.0, 1, options);
  auto b = ll::open::real_excitations_async(3, 3.0, 1.0, 1, options);
  options.execution.stop = {};
  auto sibling = ll::real_excitations_async(3, 3.0, 1.0, 1, options);
  stop.request_stop(); // after admission but before either scan runs
  EXPECT_THROW((void)a.get_wait(scheduler), bethe::ScanCancelled);
  EXPECT_THROW((void)b.get_wait(scheduler), bethe::ScanCancelled);
  EXPECT_TRUE(sibling.get_wait(scheduler).converged());
}

template <typename Real>
void same_gaudin(pairing::State<Real> const& a, pairing::State<Real> const& b, uni20::Precision p)
{
  EXPECT_EQ(std::tie(a.levels, a.blocked, a.active, a.eigenvalue_variables, a.pairs, a.iterations, a.stages,
                     a.rejected_stages, a.coupling, a.reached_coupling, a.energy, a.residual_norm,
                     a.target_residual_norm, a.particle_number_error, a.converged, a.status),
            std::tie(b.levels, b.blocked, b.active, b.eigenvalue_variables, b.pairs, b.iterations, b.stages,
                     b.rejected_stages, b.coupling, b.reached_coupling, b.energy, b.residual_norm,
                     b.target_residual_norm, b.particle_number_error, b.converged, b.status));
  for (auto const& v : a.levels)
    finite(v, p);
  for (auto const& v : a.eigenvalue_variables)
    finite(v, p);
  for (auto const& v :
       {a.coupling, a.reached_coupling, a.energy, a.residual_norm, a.target_residual_norm, a.particle_number_error})
    finite(v, p);
}
template <typename Real>
void same_gaudin(central::State<Real> const& a, central::State<Real> const& b, uni20::Precision p)
{
  EXPECT_EQ(std::tie(a.couplings, a.eigenvalue_variables, a.sz, a.up_spins, a.iterations, a.stages, a.rejected_stages,
                     a.field, a.continuation_parameter, a.residual_norm, a.number_error, a.reached_field, a.energy,
                     a.spin_reversed, a.converged, a.status),
            std::tie(b.couplings, b.eigenvalue_variables, b.sz, b.up_spins, b.iterations, b.stages, b.rejected_stages,
                     b.field, b.continuation_parameter, b.residual_norm, b.number_error, b.reached_field, b.energy,
                     b.spin_reversed, b.converged, b.status));
  for (auto const& v : a.couplings)
    finite(v, p);
  for (auto const& v : a.eigenvalue_variables)
    finite(v, p);
  for (auto const& v : {a.field, a.continuation_parameter, a.residual_norm, a.number_error})
    finite(v, p);
  if (a.energy) finite(*a.energy, p);
  if (a.reached_field) finite(*a.reached_field, p);
}

template <typename Real> auto gaudin_evaluator(std::vector<uni20::Precision> const& precisions)
{
  struct Request
  {
      std::vector<Real> levels, couplings;
      Real g, field;
      pairing::SolverOptions<Real> pairing_options;
      central::SolverOptions<Real> central_options;
  };
  std::vector<Request> requests;
  for (std::size_t i = 0; i < precisions.size(); ++i)
  {
    auto const p = precisions[i];
    auto v = [p](int n) { return uni20::scalar_at_precision<Real>(n, p); };
    auto po = controls<pairing::SolverOptions<Real>>(p);
    auto co = controls<central::SolverOptions<Real>>(p);
    if (i % 3 == 1) po.max_stages = co.max_stages = 1;
    if (i % 3 == 2) po.max_iterations = co.max_iterations = 0;
    requests.push_back(
        {{v(0), v(1), v(2), v(3), v(4), v(5)}, {v(3) / v(10), v(1), -v(7) / v(10)}, v(3) / v(7), -v(2) / v(3), po, co});
  }
  // Spans are created only inside each call, from immutable owned vectors.
  // Include blocked levels, unequal signed bath couplings and partial results.
  return [requests = std::move(requests)](std::size_t i) {
    auto const& r = requests.at(i);
    std::array<std::size_t, 2> const blocked{1, 4};
    return std::pair{pairing::ground_state<Real>(r.levels, 2, r.g, blocked, r.pairing_options),
                     central::sector_ground_state<Real>(r.couplings, r.field, uni20::half_int{0}, r.central_options)};
  };
}
template <typename Real> void gaudin_sweep(IAsyncScheduler& scheduler, std::vector<uni20::Precision> const& precisions)
{
  auto serial = gaudin_evaluator<Real>(precisions);
  using Value = decltype(serial(0));
  std::vector<Value> expected, actual;
  for (std::size_t i = 0; i < precisions.size(); ++i)
    expected.push_back(serial(i));
  auto complete = [&] {
    auto evaluate = gaudin_evaluator<Real>(precisions);
    return bethe::evaluate_observable_grid_async(precisions.size(), std::move(evaluate),
                                                 [&actual](std::size_t i, Value state) {
                                                   EXPECT_EQ(i, actual.size());
                                                   actual.push_back(std::move(state));
                                                 },
                                                 {&scheduler, 128, 2});
  }(); // all evaluator-local configuration is gone before the deferred work
  EXPECT_EQ(complete.get_wait(scheduler), precisions.size());
  ASSERT_EQ(actual.size(), expected.size());
  for (std::size_t i = 0; i < actual.size(); ++i)
  {
    same_gaudin(actual[i].first, expected[i].first, precisions[i]);
    same_gaudin(actual[i].second, expected[i].second, precisions[i]);
    EXPECT_EQ(actual[i].first.converged, i % 3 == 0);
    EXPECT_EQ(actual[i].second.converged, i % 3 == 0);
  }
}
TYPED_TEST(SolverComposition, GaudinSweepsPreserveCompleteResultsAndNumericalFailures)
{
  using Real = TypeParam;
  auto p = uni20::native_precision<Real>();
  DebugScheduler deferred({.order = DebugSchedulerOrder::reverse});
  gaudin_sweep<Real>(deferred, {p, p, p});
  for (int workers : {1, 4})
  {
    TbbScheduler scheduler(workers);
    gaudin_sweep<Real>(scheduler, {p, p, p});
  }
}

#if UNI20_HAS_MPREAL
TEST(SolverCompositionMpfr, GaudinSweepsOwnIndependentPrecisionsWithoutChangingDefaults)
{
  auto const ambient = mpfr_get_default_prec();
  std::vector<uni20::Precision> p;
  for (int bits : {128, 256, 333, 333, 128, 256})
    p.push_back(uni20::Precision::bits(bits));
  DebugScheduler deferred({.order = DebugSchedulerOrder::random, .random_seed = 31});
  gaudin_sweep<uni20::mpreal>(deferred, p);
  TbbScheduler scheduler(4);
  gaudin_sweep<uni20::mpreal>(scheduler, p);
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
#endif
#if UNI20_HAS_MPREAL
TEST(SolverCompositionMpfr, FiniteLiebLinigerConcurrentPrecisionsAndFailureDiagnostics)
{
  using Real = uni20::mpreal;
  auto const ambient = mpfr_get_default_prec();
  for (int workers : {1, 4})
    for (bool fail : {false, true})
    {
      TbbScheduler scheduler(workers);
      auto a = uni20::Precision::bits(128), b = uni20::Precision::bits(333);
      auto reference_a = reference_scans<Real>(a, false), reference_b = reference_scans<Real>(b, fail);
      auto pending_a = submit_scans<Real>(scheduler, a, false);
      auto pending_b = submit_scans<Real>(scheduler, b, fail);
      check_scans(pending_b.get_wait(scheduler), reference_b, b);
      check_scans(pending_a.get_wait(scheduler), reference_a, a);
    }
  DebugScheduler deferred({.order = DebugSchedulerOrder::reverse});
  auto p = uni20::Precision::bits(256);
  auto reference = reference_scans<Real>(p, false);
  auto pending = submit_scans<Real>(deferred, p, false);
  check_scans(pending.get_wait(deferred), reference, p);
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
#endif
} // namespace
