// SPDX-License-Identifier: GPL-3.0-or-later
#include "test_support.hpp"
#include <bethe/bose_fermi.hpp>
#include <bethe/lieb_liniger_open.hpp>
#include <bethe/observable_grid.hpp>
#include <tuple>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
namespace model = bethe::bose_fermi;
using R = uni20::mpreal;
using C = uni20::mpcomplex;
class BoseFermiMpfr : public ::testing::TestWithParam<int> {
  protected:
    uni20::Precision p() const { return uni20::Precision::bits(GetParam()); }
    R v(int n) const { return R(n, p()); }
};
INSTANTIATE_TEST_SUITE_P(WorkingPrecisions, BoseFermiMpfr, ::testing::Values(128, 256, 333));

void finite(R const& x, uni20::Precision p)
{
  ASSERT_TRUE(x.initialized());
  EXPECT_FALSE(x.is_exact());
  EXPECT_EQ(x.precision(), p);
  EXPECT_TRUE(uni20::isfinite(x));
}
void finite_state(model::State<R> const& s, uni20::Precision p)
{
  for (auto const& x : {s.length, s.interaction, s.residual_norm, s.momentum})
    finite(x, p);
  if (s.energy) finite(*s.energy, p);
  for (auto const& x : s.momenta)
    finite(x, p);
  for (auto const& x : s.auxiliary)
    finite(x, p);
}
void same_state(model::State<R> const& a, model::State<R> const& b)
{
  EXPECT_EQ(std::tie(a.bosons, a.fermions, a.iterations, a.length, a.interaction, a.residual_norm, a.momentum,
                     a.momentum_index, a.momenta, a.auxiliary, a.charge_numbers, a.auxiliary_numbers,
                     a.free_fermion_modes, a.energy, a.converged, a.status),
            std::tie(b.bosons, b.fermions, b.iterations, b.length, b.interaction, b.residual_norm, b.momentum,
                     b.momentum_index, b.momenta, b.auxiliary, b.charge_numbers, b.auxiliary_numbers,
                     b.free_fermion_modes, b.energy, b.converged, b.status));
  finite_state(a, b.length.precision());
}

#if UNI20_HAS_MPLAPACK_MPREAL
// Independent, unreduced rational equations in physical units. No reflected
// coordinates, analytic rank subtraction or model Jacobian is reused here.
void rational_equations(model::State<R> const& s)
{
  auto p = s.length.precision();
  R const zero(0, p), one(1, p), w = s.interaction / 2;
  R const tolerance = R(8192, p) * R(s.bosons + s.fermions) * uni20::numeric_limits<R>::epsilon(p);
  auto ratio = [w](R d) { return C(d, w) / C(d, -w); };
  for (R k : s.momenta)
  {
    C product(one, zero);
    for (R lambda : s.auxiliary)
      product *= ratio(k - lambda);
    EXPECT_REAL_NEAR(product.real(), uni20::math::cos(k * s.length), tolerance);
    EXPECT_REAL_NEAR(product.imag(), uni20::math::sin(k * s.length), tolerance);
  }
  for (R lambda : s.auxiliary)
  {
    C product(one, zero);
    for (R k : s.momenta)
      product *= ratio(k - lambda);
    EXPECT_REAL_NEAR(product.real(), one, tolerance);
    EXPECT_REAL_NEAR(product.imag(), zero, tolerance);
  }
}

#endif

TEST_P(BoseFermiMpfr, FreeStatesAndPureSingletonNeedOnlyScalars)
{
  model::SolverOptions<R> options(p());
  options.max_iterations = 0;
  R const pi = uni20::pi<R>.at(p()), eps = uni20::numeric_limits<R>::epsilon(p());
  for (std::size_t b = 0; b < 4; ++b)
    for (std::size_t f = 0; f < 6; ++f)
    {
      auto s = model::ground_state(b, f, v(2), v(0), options);
      ASSERT_TRUE(s.converged && s.energy.has_value());
      finite_state(s, p());
      R expected = pi * pi * R(f) * (R(f) * R(f) + (f % 2 ? -v(1) : v(2))) / 12;
      EXPECT_REAL_NEAR(*s.energy, expected, v(64) * eps * (v(1) + expected));
      EXPECT_EQ(s.momentum_index, f % 2 ? 0 : f / 2);
      for (std::size_t j = 0; j < b; ++j)
        EXPECT_EQ(s.momenta[j], v(0));
    }
  auto const pure = model::ground_state(0, 4, v(2), v(8), options);
  auto const free = model::ground_state(0, 4, v(2), v(0), options);
  EXPECT_EQ(pure.energy, free.energy);
  auto const singleton = model::ground_state(1, 0, v(1), v(1), options);
  ASSERT_TRUE(singleton.converged);
  finite_state(singleton, p());
  EXPECT_EQ(*singleton.energy, v(0));
  // A finite energy well beyond native exponent ranges must not be rejected
  // through a type-wide min()/epsilon() bound or a hidden conversion.
  auto tiny = model::ground_state(0, 3, R("1e7000", p()), v(0), options);
  ASSERT_TRUE(tiny.converged && tiny.energy.has_value());
  finite_state(tiny, p());
  EXPECT_GT(*tiny.energy, v(0));
  EXPECT_REAL_NEAR(*tiny.energy * R("1e14000", p()), 8 * pi * pi, v(512) * eps);
}

TEST_P(BoseFermiMpfr, ZeroBudgetAndInputValidation)
{
  model::SolverOptions<R> options(p());
  options.max_iterations = 0;
  for (auto [b, f] : {std::pair{2u, 3u}, std::pair{3u, 0u}})
  {
    auto s = model::ground_state(b, f, v(5), v(1), options);
    EXPECT_FALSE(s.converged);
    EXPECT_FALSE(s.energy);
    EXPECT_EQ(s.iterations, 0u);
    EXPECT_EQ(s.status, model::Status::iteration_limit);
    finite_state(s, p());
    EXPECT_EQ(s.momentum, v(0));
  }
  // The shared driver's zero-update path also remains available to open LL.
  auto open = bethe::lieb_liniger::open::ground_state(2, v(1), v(1), options);
  EXPECT_FALSE(open.converged);
  EXPECT_EQ(open.iterations, 0u);
  finite(open.residual_norm, p());
  for (R bad : {R(), R(1), R(1, uni20::Precision::bits(GetParam() + 1))})
  {
    EXPECT_THROW((void)model::ground_state(0, 0, bad, v(0), options), std::invalid_argument);
    EXPECT_THROW((void)model::ground_state(0, 0, v(1), bad, options), std::invalid_argument);
  }
  EXPECT_THROW((void)model::ground_state(2, 2, v(1), v(1), options), std::invalid_argument);
  EXPECT_THROW((void)model::ground_state(1, 1, v(0), v(1), options), std::invalid_argument);
  EXPECT_THROW((void)model::ground_state(1, 1, v(1), -v(1), options), std::invalid_argument);
  options.residual_tolerance = R(1);
  EXPECT_THROW((void)model::ground_state(0, 0, v(1), v(0), options), std::invalid_argument);
}

#if !UNI20_HAS_MPLAPACK_MPREAL
TEST_P(BoseFermiMpfr, MissingProviderOnlyRejectsActualNewtonUpdates)
{
  for (auto [b, f] : {std::pair{2u, 3u}, std::pair{3u, 0u}})
    try
    {
      (void)model::ground_state(b, f, v(1), v(1), model::SolverOptions<R>(p()));
      FAIL() << "MPFR LU should not be configured";
    }
    catch (std::runtime_error const& error)
    {
      EXPECT_NE(std::string(error.what()).find("require MPLAPACK MPFR"), std::string::npos);
    }
}
#endif

#if UNI20_HAS_MPLAPACK_MPREAL
TEST_P(BoseFermiMpfr, MixedEquationsScalingAndLiebLinigerReductions)
{
  model::SolverOptions<R> options(p());
  R const eps = uni20::numeric_limits<R>::epsilon(p());
  for (auto [b, f] : {std::pair{1u, 1u}, std::pair{2u, 3u}, std::pair{3u, 3u}, std::pair{4u, 5u}})
    for (R c : {v(1) / 100, v(1), v(100)})
    {
      auto s = model::ground_state(b, f, v(1), c, options);
      ASSERT_TRUE(s.converged) << b << ',' << f << ',' << int(s.status);
      finite_state(s, p());
      rational_equations(s);
      auto scaled = model::ground_state(b, f, v(2), c / 2, options);
      ASSERT_TRUE(scaled.converged);
      EXPECT_REAL_NEAR(*scaled.energy * 4, *s.energy, v(4096) * eps * *s.energy);
    }
  for (std::size_t b : {1, 2, 3, 5})
  {
    auto mixed = model::ground_state(b, 1, v(1), v(2), options);
    auto pure = model::ground_state(b + 1, 0, v(1), v(2), options);
    auto ll = bethe::lieb_liniger::ground_state(b + 1, v(1), v(2), options);
    ASSERT_TRUE(mixed.converged && pure.converged && ll.converged);
    EXPECT_EQ(pure.energy, ll.energy);
    EXPECT_REAL_NEAR(*mixed.energy, ll.energy, v(4096) * eps * ll.energy);
    finite_state(mixed, p());
    finite_state(pure, p());
  }
}

TEST(BoseFermiMpfrShared, IndependentTwoBodyAccuracyAndWeakCouplingCorrection)
{
  auto p = uni20::Precision::bits(512);
  R lo(0, p), hi = uni20::pi<R>.at(p), half = R(1, p) / 2;
  // Independent contact jump condition at L=c=1: q*tan(q/2)=1/2.
  for (int i = 0; i < 530; ++i)
  {
    R mid = (lo + hi) / 2;
    if (mid == lo || mid == hi) break;
    if (mid * uni20::math::tan(mid / 2) > half)
      hi = mid;
    else
      lo = mid;
  }
  R const reference = 2 * ((lo + hi) / 2) * ((lo + hi) / 2);
  R previous(1, p);
  for (int bits : {128, 256, 333})
  {
    auto q = uni20::Precision::bits(bits);
    auto s = model::ground_state(1, 1, R(1, q), R(1, q), model::SolverOptions<R>(q));
    ASSERT_TRUE(s.converged && s.energy.has_value());
    R error = uni20::math::abs(s.energy->at(p) - reference);
    EXPECT_LT(error, previous / R("1e15", p));
    EXPECT_LT(error, (R(4096, q) * uni20::numeric_limits<R>::epsilon(q)).at(p));
    previous = error;
  }
  auto q = uni20::Precision::bits(333);
  R const c("1e-60", q);
  auto weak = model::ground_state(1, 1, R(1, q), c, model::SolverOptions<R>(q));
  ASSERT_TRUE(weak.converged && weak.energy.has_value());
  EXPECT_REAL_NEAR((*weak.energy / (2 * c) - 1) / c, -R(1, q) / 12, R("1e-35", q));
  // First-order perturbation theory: c*N_b*(N_b-1+2*N_f)/L.
  // Unlike the two-body zero-energy limit, this resolves an O(c) shift on
  // top of a nonzero Fermi sea, with nontrivial auxiliary roots.
  auto free = model::ground_state(2, 3, R(1, q), R(0, q), model::SolverOptions<R>(q));
  auto mixture = model::ground_state(2, 3, R(1, q), c, model::SolverOptions<R>(q));
  ASSERT_TRUE(free.converged && mixture.converged && free.energy && mixture.energy);
  EXPECT_REAL_NEAR((*mixture.energy - *free.energy) / c, R(14, q), R("1e-35", q));
}
#endif

TEST(BoseFermiMpfrShared, OwningMixedPrecisionGridAndPartialStates)
{
  auto const ambient = mpfr_get_default_prec();
  struct Request
  {
      std::size_t b, f;
      R length, c;
      model::SolverOptions<R> options;
  };
  auto submit = [](uni20::async::IAsyncScheduler& scheduler, auto& collected) {
    std::vector<Request> requests;
    for (int bits : {128, 256, 333})
    {
      auto p = uni20::Precision::bits(bits);
      model::SolverOptions<R> options(p);
      requests.push_back({2, 2, R(1, p), R(0, p), options});
#if UNI20_HAS_MPLAPACK_MPREAL
      requests.push_back({2, 3, R(5, p), R(1, p), options});
#endif
      options.max_iterations = 0;
      requests.push_back({3, 3, R(6, p), R(1, p), options});
    }
    return bethe::evaluate_observable_grid_async(
        requests.size(),
        [requests](std::size_t i) {
          auto const& r = requests[i];
          return model::ground_state(r.b, r.f, r.length, r.c, r.options);
        },
        [&collected](std::size_t i, auto state) {
          EXPECT_EQ(i, collected.size());
          collected.push_back(std::move(state));
        },
        {&scheduler, 128, 2});
  };
  uni20::async::DebugScheduler deferred({.order = uni20::async::DebugSchedulerOrder::reverse});
  std::vector<model::State<R>> reference;
  auto serial = submit(deferred, reference);
  auto const count = serial.get_wait(deferred);
  EXPECT_EQ(count, UNI20_HAS_MPLAPACK_MPREAL ? 9u : 6u);
  EXPECT_TRUE(reference.front().converged);
  EXPECT_FALSE(reference.back().converged);
  for (int workers : {1, 4})
  {
    uni20::async::TbbScheduler scheduler(workers);
    std::vector<model::State<R>> actual;
    auto complete = submit(scheduler, actual);
    EXPECT_EQ(complete.get_wait(scheduler), count);
    ASSERT_EQ(actual.size(), reference.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
      same_state(actual[i], reference[i]);
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
} // namespace
