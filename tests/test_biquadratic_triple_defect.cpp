// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include "tl_ed.hpp"
#include <bethe/biquadratic_ferromagnetic.hpp>
#include <bethe/xxz_open_triple_defect.hpp>

namespace
{
namespace qg = bethe::xxz::quantum_group;
TEST(TripleDefectED, SmallModules)
{
  for (unsigned n = 8; n <= 10; ++n)
  {
    auto ed = bethe::test::quantum_group_module_ed(n, n - 8, 1.5);
    auto remove = [&](double energy) {
      auto nearest = std::min_element(ed.begin(), ed.end(),
                                      [&](auto a, auto b) { return std::abs(a - energy) < std::abs(b - energy); });
      ASSERT_NE(nearest, ed.end());
      ASSERT_NEAR(*nearest, energy, 2e-11);
      ed.erase(nearest);
    };
    namespace ferro = bethe::biquadratic::ferromagnetic;
    auto const real = ferro::real_excitations<double>(n, n - 8, {.count = 1000});
    ASSERT_TRUE(real.converged());
    for (auto const& s : real.levels)
      ASSERT_NO_FATAL_FAILURE(remove(s.state.reference.energy));
    for (std::size_t i = 1; i < n - 4; ++i)
      for (std::size_t k = i + 1; k <= n - 4; ++k)
        for (std::size_t j = 1; j <= n - 7; ++j)
        {
          auto const s = ferro::pair_with_real_roots<double>(n, std::array<std::size_t, 2>{i, k}, j);
          ASSERT_TRUE(s.reference.converged);
          ASSERT_NO_FATAL_FAILURE(remove(s.reference.energy));
        }
    for (std::size_t i = 1; i < n - 6; ++i)
      for (std::size_t j = i + 1; j <= n - 6; ++j)
      {
        auto const s = ferro::two_bound_pairs<double>(n, {i, j});
        ASSERT_TRUE(s.reference.converged);
        ASSERT_NO_FATAL_FAILURE(remove(s.reference.energy));
      }
    for (std::size_t i = 1; i <= n - 3; ++i)
      for (std::size_t j = 1; j <= n - 7; ++j)
      {
        qg::triple_defect::detail::System<double> system(n, 1.5, i, j);
        auto const s = qg::detail::solve_log_string<double>(system, {});
        EXPECT_TRUE(s.converged) << n << " " << i << " " << j << " status " << int(s.status) << " x=" << s.x[0] << ","
                                 << s.x[1] << "," << s.x[2] << "," << s.x[3] << " norm=" << system.evaluate(s.x).norm;
        if (!s.converged) continue;
        double const energy = (n - 1) * 1.5 / 4 + system.energy_shift(s.x);
        ASSERT_NO_FATAL_FAILURE(remove(energy));
      }
    EXPECT_EQ(ed.size(), n - 7);
    for (std::size_t mode = 1; mode <= n - 7; ++mode)
    {
      auto const s = ferro::bound_quartet<double>(n, mode);
      ASSERT_TRUE(s.reference.converged);
      ASSERT_NO_FATAL_FAILURE(remove(s.reference.energy));
    }
    EXPECT_TRUE(ed.empty()); // All five four-defect string topologies, as a multiset.
  }
}
template <typename R> class TripleDefect : public ::testing::Test {};
TYPED_TEST_SUITE(TripleDefect, test_support::RealTypes, test_support::PrecisionNames);
TEST(TripleDefectED, OtherAnisotropies)
{
  for (double delta : {1.25, 2.0, 3.0})
    for (unsigned n = 8; n <= 10; ++n)
    {
      auto ed = bethe::test::quantum_group_module_ed(n, n - 8, delta);
      for (std::size_t i = 1; i <= n - 3; ++i)
        for (std::size_t j = 1; j <= n - 7; ++j)
        {
          auto const s = qg::triple_defect::solve<double>(n, delta, i, j);
          EXPECT_TRUE(s.converged) << n << " " << delta << " " << i << " " << j << " x=" << s.center << ","
                                   << s.log_deviation << "," << s.deviation_phase << "," << s.rapidity
                                   << " residual=" << s.residual_norm << " status=" << int(s.status);
          if (!s.converged) continue;
          auto nearest = std::min_element(
              ed.begin(), ed.end(), [&](auto a, auto b) { return std::abs(a - s.energy) < std::abs(b - s.energy); });
          ASSERT_NE(nearest, ed.end());
          ASSERT_NEAR(*nearest, s.energy, 2e-11);
          ed.erase(nearest);
        }
    }
}
TYPED_TEST(TripleDefect, CoupledJacobian)
{
  using R = TypeParam;
  R const h = std::cbrt(uni20::numeric_limits<R>::epsilon());
  qg::triple_defect::detail::System<R> system(16, R{3} / R{2}, 4, 3);
  for (bool ideal : {false, true})
    for (R alpha : {R{1} / R{2}, R{1}, R{2}})
    {
      if (!ideal && alpha == R{1}) continue; // Coincident roots are not physical.
      std::vector<R> x{R{1}, R{4}, R{1} / R{3}, alpha}, jac;
      if (alpha == R{1}) EXPECT_FALSE(system.physical(x));
      system.evaluate(x, &jac, ideal);
      for (std::size_t j = 0; j < 4; ++j)
      {
        auto a = x, b = x;
        a[j] += h;
        b[j] -= h;
        auto const fa = system.evaluate(a, nullptr, ideal), fb = system.evaluate(b, nullptr, ideal);
        for (std::size_t i = 0; i < 4; ++i)
          EXPECT_REAL_NEAR((fa.residual[i] - fb.residual[i]) / (R{2} * h), jac[4 * i + j], R{4096} * h * h);
      }
    }
}

TYPED_TEST(TripleDefect, OriginalEquationsAndLongChains)
{
  using R = TypeParam;
  using C = std::complex<R>;
  R const eps = uni20::numeric_limits<R>::epsilon();
  for (std::size_t n : {8, 9, 10, 16, 65, 128, 1024, 100000})
    for (auto labels : {std::pair<std::size_t, std::size_t>{1, 1},
                        {n - 3, n - 7},
                        {1, n - 7},
                        {n - 3, 1},
                        {n - 3, n > 8 ? n - 8 : 1}})
    {
      auto const state = bethe::biquadratic::ferromagnetic::triple_defect<R>(n, labels.first, labels.second);
      auto const& s = state.reference;
      ASSERT_TRUE(s.converged) << n << " " << labels.first << " " << labels.second << " status " << int(s.status);
      EXPECT_EQ(state.through_lines, n - 8);
      EXPECT_EQ(state.tl_energy, -R{2} * s.energy_shift);
      R const eta = std::acosh(s.delta);
      C const z = std::exp(-s.log_deviation) * C{std::cos(s.deviation_phase), std::sin(s.deviation_phase)};
      C const u = C{eta, s.center / R{2}} + z, v{0, s.center / R{2}}, w{0, s.rapidity / R{2}};
      auto scattering = [&](C a, C b) {
        return std::sinh(a - b + eta) / std::sinh(a - b - eta) * std::sinh(a + b + eta) / std::sinh(a + b - eta);
      };
      C const drive = std::sinh(u + eta / R{2}) / std::sinh(u - eta / R{2});
      C const regular = std::sinh(R{2} * eta + z) * std::sinh(u + v + eta) / std::sinh(u + v - eta) *
                        scattering(u, std::conj(u)) * scattering(u, w);
      C const corr = std::abs(z) < std::sqrt(eps) ? z * z / R{6} : std::log(std::sinh(z) / z);
      R const phase = R{2} * R(n) * std::arg(drive) - std::arg(regular) + s.deviation_phase + corr.imag();
      EXPECT_REAL_NEAR(std::sin(phase / R{2}), R{0}, R{1024} * R(n) * eps);
      EXPECT_REAL_NEAR(std::log(std::abs(drive)) -
                           (std::log(std::abs(regular)) + s.log_deviation - corr.real()) / (R{2} * R(n)),
                       R{0}, R{512} * eps);
      C const central_drive = std::sinh(v + eta / R{2}) / std::sinh(v - eta / R{2});
      C const b = std::sinh(u + v + eta) / std::sinh(u + v - eta);
      R const central_phase =
          R{2} * R(n) * std::arg(central_drive) -
          R{2} * (s.deviation_phase + corr.imag() - std::arg(std::sinh(R{2} * eta + z)) + std::arg(b)) -
          std::arg(scattering(v, w));
      EXPECT_REAL_NEAR(std::sin(central_phase / R{2}), R{0}, R{2048} * R(n) * eps);
      C const real_drive = std::sinh(w + eta / R{2}) / std::sinh(w - eta / R{2});
      C const real_rhs = scattering(w, v) * scattering(w, u) * scattering(w, std::conj(u));
      EXPECT_REAL_NEAR(std::sin((R{2} * R(n) * std::arg(real_drive) - std::arg(real_rhs)) / R{2}), R{0},
                       R{1024} * R(n) * eps);
      EXPECT_REAL_NEAR(std::abs(real_rhs), R{1}, R{512} * eps);
      if (labels.first == n - 3 && labels.second == n - 7)
      {
        EXPECT_GT(state.tl_energy, R{3});
        if (n >= 16) EXPECT_LT(state.tl_energy - R{3}, R{50} / (R(n) * R(n)));
        if (n == 100000) EXPECT_EQ(std::exp(-s.log_deviation), R{0});
      }
    }
}

TYPED_TEST(TripleDefect, BudgetsAndValidation)
{
  using R = TypeParam;
  for (std::size_t budget : {0, 1})
  {
    auto const s = qg::triple_defect::solve<R>(16, R{3} / R{2}, 13, 9, {.max_iterations = budget});
    EXPECT_FALSE(s.converged);
    EXPECT_EQ(s.iterations, budget);
    EXPECT_EQ(s.status, qg::SolveStatus::iteration_limit);
    qg::triple_defect::detail::System<R> system(16, s.delta, 13, 9);
    std::vector<R> x{s.center, s.log_deviation, s.deviation_phase, s.rapidity};
    EXPECT_EQ(s.residual_norm, system.evaluate(x).norm);
    EXPECT_EQ(s.energy_shift, system.energy_shift(x));
  }
  for (std::size_t n : {0, 6, 7})
    EXPECT_THROW((void)qg::triple_defect::solve<R>(n, R{3} / R{2}, 1, 1), std::invalid_argument);
  for (auto labels : {std::pair<std::size_t, std::size_t>{0, 1}, {1, 0}, {6, 1}, {1, 2}})
    EXPECT_THROW((void)qg::triple_defect::solve<R>(8, R{3} / R{2}, labels.first, labels.second), std::invalid_argument);
  for (R delta : {R{1}, R{0}, uni20::numeric_limits<R>::infinity(), uni20::numeric_limits<R>::quiet_NaN()})
    EXPECT_THROW((void)qg::triple_defect::solve<R>(8, delta, 1, 1), std::invalid_argument);
  EXPECT_THROW((void)qg::triple_defect::solve<R>(8, R{3} / R{2}, 1, 1, {.residual_tolerance = R{0}}),
               std::invalid_argument);
  EXPECT_THROW((void)qg::triple_defect::solve<R>(std::numeric_limits<std::size_t>::max(), R{3} / R{2}, 1, 1),
               std::invalid_argument);
}

TYPED_TEST(TripleDefect, RoundedAngleUnderflowPolish)
{
  using R = TypeParam;
  for (std::size_t i : {99996, 99997})
    for (std::size_t j : {99992, 99993})
    {
      auto const s = qg::triple_defect::solve<R>(100000, R{3} / R{2}, i, j);
      ASSERT_TRUE(s.converged) << i << " " << j;
      EXPECT_LE(s.residual_norm, R{32} * uni20::numeric_limits<R>::epsilon());
      EXPECT_EQ(std::exp(-s.log_deviation), R{0});
      qg::triple_defect::detail::System<R> system(100000, s.delta, i, j);
      std::vector<R> x{s.center, s.log_deviation + R{1}, s.deviation_phase + R{1}, s.rapidity};
      auto ideal = x;
      system.normalize(ideal, true);
      EXPECT_EQ(ideal[1], x[1]); // Never apply the finite-row correction in initialization.
      system.normalize(x, false);
      EXPECT_EQ(x[0], s.center);
      EXPECT_EQ(x[3], s.rapidity);
      auto const f = system.evaluate(x);
      EXPECT_LE(f.modulus_norm, R{32} * uni20::numeric_limits<R>::epsilon());
      EXPECT_LE(f.phase_norm, R{32} * uni20::numeric_limits<R>::epsilon());
    }
}
} // namespace
