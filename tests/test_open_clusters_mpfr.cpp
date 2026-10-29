// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/biquadratic_ferromagnetic.hpp>
#include <bethe/observable_grid.hpp>
#include <tuple>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using R = uni20::mpreal;
using C = uni20::complex<R>;
using uni20::Precision;
namespace math = uni20::math;
namespace qg = bethe::xxz::quantum_group;
namespace ferro = bethe::biquadratic::ferromagnetic;

// Independent open TL link-pattern Hamiltonian sum_i e_i at loop weight 3.
// Joining neighbouring paired endpoints contributes 3; otherwise reconnect
// their partners with coefficient 1. The five-dimensional N=6 characteristic
// polynomial is (g-6)(g-7)(g^3-17g^2+80g-106). At N=8 the dimension is 14:
// (g^4-30g^3+329g^2-1560g+2693) times the degree-ten factor below.
// These reference brackets select simple eigenvalues independently of Bethe roots.
R exact_gap(bool six_sites, char const* left, char const* right)
{
  auto p = Precision::bits(768);
  auto polynomial = [&](R const& x) {
    if (six_sites) return ((x - 17) * x + 80) * x - 106;
    R value(1, p);
    for (int coefficient : {-75, 2476, -47340, 579856, -4747746, 26269515, -96768627, 226453872, -302801604, 174717400})
      value = value * x + coefficient;
    return value;
  };
  R a(left, p), b(right, p), fa = polynomial(a);
  EXPECT_LT(fa * polynomial(b), R(0, p));
  for (int i = 0; i < 780; ++i)
  {
    R const m = (a + b) / 2, fm = polynomial(m);
    if (fm == R(0, p)) return m;
    if (fa * fm < 0)
      b = m;
    else
    {
      a = m;
      fa = fm;
    }
  }
  return (a + b) / 2;
}

template <typename S> void check_precision(S const& s, Precision p)
{
  for (auto const* v : {&s.delta, &s.energy, &s.energy_shift, &s.residual_norm, &s.phase_residual, &s.modulus_residual})
    EXPECT_EQ(uni20::precision_of(*v), p);
  if constexpr (requires { s.center; }) EXPECT_EQ(s.center.precision(), p);
  if constexpr (requires { s.log_deviation; })
  {
    EXPECT_EQ(s.log_deviation.precision(), p);
    EXPECT_EQ(s.deviation_phase.precision(), p);
  }
  if constexpr (requires { s.inner_log_deviation; })
    for (auto const* v : {&s.inner_log_deviation, &s.outer_log_deviation, &s.outer_deviation_phase})
      EXPECT_EQ(v->precision(), p);
  if constexpr (requires { s.rapidity; }) EXPECT_EQ(s.rapidity.precision(), p);
  if constexpr (requires { s.centers; })
    for (std::size_t j = 0; j < 2; ++j)
    {
      EXPECT_EQ(s.centers[j].precision(), p);
      EXPECT_EQ(s.log_deviations[j].precision(), p);
    }
}

TEST(OpenClustersMpfr, IndependentSmallModuleEnergiesGainAccuracy)
{
  std::array<R, 4> const exact{exact_gap(true, "2.2", "2.3"), exact_gap(false, "2.2", "2.3"),
                               exact_gap(false, "4.9", "5.0"), exact_gap(false, "3.8", "3.9")};
  R previous(1, Precision::bits(768));
  auto const ambient = mpfr_get_default_prec();
  for (int bits : {80, 160, 256})
  {
    SCOPED_TRACE(bits);
    auto p = Precision::bits(bits);
    bethe::SolverOptions<R> options(p);
    R error(0, Precision::bits(768));
    auto check = [&](auto const& s, std::size_t i) {
      EXPECT_TRUE(s.reference.converged) << i << ' ' << int(s.reference.status);
      check_precision(s.reference, p);
      EXPECT_EQ(s.tl_energy.precision(), p);
      EXPECT_EQ(s.energy.precision(), p);
      EXPECT_EQ(s.tl_energy, -2 * s.reference.energy_shift);
      error = std::max(error, math::abs(s.tl_energy.at(Precision::bits(768)) - exact[i]));
    };
    check(ferro::bound_triple<R>(6, 1, options), 0);
    check(ferro::bound_quartet<R>(8, 1, options), 1);
    check(ferro::two_bound_pairs<R>(8, {1, 2}, options), 2);
    check(ferro::triple_defect<R>(8, 5, 1, options), 3);
    EXPECT_LT(error, 65536 * uni20::numeric_limits<R>::epsilon(p).at(Precision::bits(768)));
    EXPECT_LT(error, previous * R("1e-15", Precision::bits(768)));
    if (bits == 256) EXPECT_LT(error, R("1e-70", Precision::bits(768)));
    previous = error;
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

template <typename System> void check_jacobian(System const& system, std::vector<R> x)
{
  auto p = x.front().precision();
  R const one(1, p), h = math::ldexp(one, -50);
  for (bool ideal : {false, true})
  {
    std::vector<R> jac;
    auto const f = system.evaluate(x, &jac, ideal);
    for (auto const* v : {&f.norm, &f.phase_norm, &f.modulus_norm})
      EXPECT_EQ(v->precision(), p);
    for (auto const& v : jac)
      EXPECT_EQ(v.precision(), p);
    for (std::size_t j = 0; j < x.size(); ++j)
    {
      auto a = x, b = x;
      a[j] += h;
      b[j] -= h;
      auto const fa = system.evaluate(a, nullptr, ideal), fb = system.evaluate(b, nullptr, ideal);
      for (std::size_t i = 0; i < x.size(); ++i)
        EXPECT_LT(math::abs((fa.residual[i] - fb.residual[i]) / (2 * h) - jac[i * x.size() + j]),
                  1000000 * h * h * std::max(one, math::abs(jac[i * x.size() + j])));
    }
  }
}

TEST(OpenClustersMpfr, ComplexAndCoupledJacobians)
{
  auto p = Precision::bits(160);
  for (auto const& d : {R(5, p) / 4, R(3, p) / 2, R(3, p)})
  {
    check_jacobian(qg::three_string::detail::System<R>(16, d, 11), {R(2, p), R(4, p), R(1, p)});
    check_jacobian(qg::four_string::detail::System<R>(16, d, 9), {R(2, p), R(4, p), R(5, p), R(1, p)});
    check_jacobian(qg::two_pairs::detail::System<R>(16, d, {9, 10}), {R(1, p), R(3, p), R(2, p), R(4, p)});
    check_jacobian(qg::triple_defect::detail::System<R>(16, d, 7, 8), {R(2, p), R(4, p), R(1, p), R(1, p)});
  }
}

// Separate real-component implementation checks the shared MPC scattering
// amplitudes and their derivatives without calling complex sinh/cosh.
TEST(OpenClustersMpfr, ScatteringMatchesRealComponentReference)
{
  auto oracle = Precision::bits(768);
  auto sinh = [](C const& z) -> C {
    return {math::sinh(z.real()) * math::cos(z.imag()), math::cosh(z.real()) * math::sin(z.imag())};
  };
  auto cosh = [](C const& z) -> C {
    return {math::cosh(z.real()) * math::cos(z.imag()), math::sinh(z.real()) * math::sin(z.imag())};
  };
  for (int bits : {80, 160, 256})
  {
    auto p = Precision::bits(bits);
    C const u{R(5, p) / 4, R(7, p) / 8}, v{R(3, p) / 8, -R(1, p) / 4};
    R const eta = R(3, p) / 4;
    auto const got = qg::detail::open_scattering(u, v, eta);
    C const a = u.at(oracle), b = v.at(oracle);
    R const e = eta.at(oracle);
    C const ap = a - b + e, am = a - b - e, bp = a + b + e, bm = a + b - e;
    C const value = math::log(sinh(ap)) - math::log(sinh(am)) + math::log(sinh(bp)) - math::log(sinh(bm));
    C const first = cosh(ap) / sinh(ap) - cosh(am) / sinh(am) + cosh(bp) / sinh(bp) - cosh(bm) / sinh(bm);
    C const second = -cosh(ap) / sinh(ap) + cosh(am) / sinh(am) + cosh(bp) / sinh(bp) - cosh(bm) / sinh(bm);
    for (auto const& [actual, expected] : {std::pair{got.value, value}, {got.first, first}, {got.second, second}})
    {
      EXPECT_EQ(uni20::precision_of(actual), p);
      EXPECT_LT(math::abs(actual.at(oracle) - expected), 4096 * uni20::numeric_limits<R>::epsilon(p).at(oracle));
    }
  }
}

TEST(OpenClustersMpfr, LongDeviationsAndFailureStates)
{
  for (int bits : {80, 160, 256})
  {
    auto p = Precision::bits(bits);
    bethe::SolverOptions<R> options(p);
    auto triple = ferro::bound_triple<R>(256, 1, options).reference;
    auto quartet = ferro::bound_quartet<R>(256, 1, options).reference;
    auto pairs = ferro::two_bound_pairs<R>(256, {249, 250}, options).reference;
    auto defect = ferro::triple_defect<R>(32, 29, 25, options).reference;
    EXPECT_TRUE(triple.converged);
    EXPECT_TRUE(quartet.converged);
    EXPECT_TRUE(pairs.converged);
    EXPECT_TRUE(defect.converged);
    check_precision(triple, p);
    check_precision(quartet, p);
    check_precision(pairs, p);
    check_precision(defect, p);
    EXPECT_GT(triple.log_deviation, R(100, p));
    EXPECT_GT(quartet.outer_log_deviation, R(100, p));
    for (auto const& L : pairs.log_deviations)
      EXPECT_GT(L, R(100, p));
    for (std::size_t budget : {0, 1})
    {
      options.max_iterations = budget;
      R const d = R(3, p) / 2;
      auto check = [&](auto const& s, auto const& system, std::vector<R> x) {
        EXPECT_FALSE(s.converged);
        EXPECT_EQ(s.status, qg::SolveStatus::iteration_limit);
        EXPECT_EQ(s.iterations, budget);
        check_precision(s, p);
        EXPECT_EQ(s.residual_norm, system.evaluate(x).norm);
      };
      auto const t = qg::three_string::bound_triple<R>(16, d, 1, options);
      check(t, qg::three_string::detail::System<R>(16, d, t.string_label),
            {t.center, t.log_deviation, t.deviation_phase});
      auto const q = qg::four_string::bound_quartet<R>(16, d, 1, options);
      check(q, qg::four_string::detail::System<R>(16, d, q.string_label),
            {q.center, q.inner_log_deviation, q.outer_log_deviation, q.outer_deviation_phase});
      auto const b = qg::two_pairs::solve<R>(16, d, {9, 10}, options);
      check(b, qg::two_pairs::detail::System<R>(16, d, b.string_labels),
            {b.centers[0], b.log_deviations[0], b.centers[1], b.log_deviations[1]});
      auto const f = qg::triple_defect::solve<R>(16, d, 13, 9, options);
      check(f, qg::triple_defect::detail::System<R>(16, d, f.real_label, f.string_label),
            {f.center, f.log_deviation, f.deviation_phase, f.rapidity});
    }
  }
}

TEST(OpenClustersMpfr, RejectsInvalidPrecisionAndControls)
{
  auto p = Precision::bits(160);
  bethe::SolverOptions<R> options(p);
  auto reject = [&](R const& d) {
    EXPECT_THROW(qg::three_string::bound_triple<R>(8, d, 1, options), std::invalid_argument);
    EXPECT_THROW(qg::four_string::bound_quartet<R>(8, d, 1, options), std::invalid_argument);
    EXPECT_THROW(qg::two_pairs::solve<R>(8, d, {1, 2}, options), std::invalid_argument);
    EXPECT_THROW(qg::triple_defect::solve<R>(8, d, 1, 1, options), std::invalid_argument);
  };
  for (auto const& bad : {R(uni20::uninitialized), R{}, R(2), R(2, Precision::bits(80)), R("nan", p), R(1, p)})
    reject(bad);
  for (auto const& bad : {R(uni20::uninitialized), R{}, R(1), R(1, Precision::bits(80)), R("nan", p), R(0, p)})
  {
    options.residual_tolerance = bad;
    reject(R(3, p) / 2);
  }
}

TEST(OpenClustersMpfr, OwningMixedPrecisionQueriesPreserveAllFields)
{
  struct Request
  {
      bethe::SolverOptions<R> options;
  };
  std::vector<Request> requests;
  for (int bits : {80, 160})
    for (bool fail : {false, true})
    {
      bethe::SolverOptions<R> options(Precision::bits(bits));
      if (fail) options.max_iterations = 0;
      requests.push_back({options});
    }
  auto solve = [](Request const& r) {
    return std::tuple{ferro::bound_triple<R>(8, 1, r.options), ferro::bound_quartet<R>(8, 1, r.options),
                      ferro::two_bound_pairs<R>(8, {1, 2}, r.options), ferro::triple_defect<R>(8, 5, 1, r.options)};
  };
  auto fields = [](auto const& s) {
    auto const& r = s.reference;
    auto common = std::tuple{s.sites,
                             s.through_lines,
                             s.mode,
                             s.energy,
                             s.tl_energy,
                             s.multiplicity,
                             r.sites,
                             r.delta,
                             r.energy,
                             r.energy_shift,
                             r.residual_norm,
                             r.phase_residual,
                             r.modulus_residual,
                             r.iterations,
                             r.converged,
                             r.status};
    if constexpr (requires { r.centers; })
      return std::tuple_cat(common, std::tuple{r.string_labels, r.centers, r.log_deviations, r.deviation_signs});
    else if constexpr (requires { r.inner_log_deviation; })
      return std::tuple_cat(common, std::tuple{r.string_label, r.center, r.inner_log_deviation, r.outer_log_deviation,
                                               r.outer_deviation_phase, r.inner_deviation_sign});
    else if constexpr (requires { r.rapidity; })
      return std::tuple_cat(
          common, std::tuple{r.string_label, r.center, r.log_deviation, r.deviation_phase, r.real_label, r.rapidity});
    else
      return std::tuple_cat(common, std::tuple{r.string_label, r.center, r.log_deviation, r.deviation_phase});
  };
  using Result = decltype(solve(requests.front()));
  std::vector<Result> expected;
  for (auto const& r : requests)
    expected.push_back(solve(r));
  for (int workers : {1, 2})
  {
    uni20::async::TbbScheduler scheduler(workers);
    std::vector<Result> actual;
    auto done = [&] {
      auto owned = requests;
      return bethe::evaluate_observable_grid_async(
          owned.size(), [owned, solve](std::size_t i) { return solve(owned[i]); },
          [&actual](std::size_t, Result value) { actual.push_back(std::move(value)); },
          {.scheduler = &scheduler, .max_in_flight = 2});
    }();
    (void)done.get_wait(scheduler);
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
      [&]<std::size_t... J>(std::index_sequence<J...>) {
        auto equal = [&](auto const& a, auto const& b) { EXPECT_EQ(fields(a), fields(b)); };
        (equal(std::get<J>(actual[i]), std::get<J>(expected[i])), ...);
      }(std::make_index_sequence<4>{});
  }
}
} // namespace
