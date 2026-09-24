// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include "tl_ed.hpp"
#include <bethe/biquadratic_ferromagnetic.hpp>
#include <bethe/xxz_open_four_string.hpp>

namespace
{
namespace qg = bethe::xxz::quantum_group;
TEST(FourStringED, DropletBranches)
{
  for (double delta : {1.25, 1.5, 2.0, 3.0})
    for (unsigned n = 8; n <= 10; ++n)
    {
      auto ed = bethe::test::quantum_group_module_ed(n, n - 8, delta);
      for (std::size_t mode = 1; mode <= n - 7; ++mode)
      {
        qg::four_string::detail::System<double> system(n, delta, n - 6 - mode);
        auto const s = qg::detail::solve_log_string<double>(system, {});
        ASSERT_TRUE(s.converged) << n << " " << delta << " " << mode << " status " << int(s.status);
        double const energy = (n - 1) * delta / 4 + system.energy_shift(s.x);
        auto nearest = std::min_element(ed.begin(), ed.end(),
                                        [&](auto a, auto b) { return std::abs(a - energy) < std::abs(b - energy); });
        ASSERT_NE(nearest, ed.end());
        EXPECT_NEAR(*nearest, energy, 2e-11) << n << " " << delta << " " << mode;
        if (mode == 1) EXPECT_NEAR(ed.back(), energy, 2e-11); // Largest reference energy = smallest ferro gap.
        ed.erase(nearest);
      }
    }
}
template <typename R> class FourString : public ::testing::Test {};
TYPED_TEST_SUITE(FourString, test_support::RealTypes, test_support::PrecisionNames);
TYPED_TEST(FourString, Jacobian)
{
  using R = TypeParam;
  R const h = std::cbrt(uni20::numeric_limits<R>::epsilon());
  for (std::size_t label : {9, 10})
    for (bool ideal : {false, true})
    {
      qg::four_string::detail::System<R> system(16, R{3} / R{2}, label);
      std::vector<R> x{R{2}, R{4}, R{5}, R{1}}, jac;
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

TYPED_TEST(FourString, OriginalEquationsAndLongChains)
{
  using R = TypeParam;
  using C = std::complex<R>;
  R const eps = uni20::numeric_limits<R>::epsilon(), pi = R{4} * std::atan(R{1});
  for (std::size_t n : {8, 9, 10, 16, 65, 128, 1024, 100000})
    for (std::size_t mode : {std::size_t{1}, std::min(std::size_t{2}, n - 7), n - 7})
    {
      auto const state = bethe::biquadratic::ferromagnetic::bound_quartet<R>(n, mode);
      auto const& s = state.reference;
      ASSERT_TRUE(s.converged) << n << " " << mode << " status=" << int(s.status);
      EXPECT_EQ(s.inner_deviation_sign, mode % 2 ? 1 : -1);
      EXPECT_EQ(state.through_lines, n - 8);
      EXPECT_EQ(state.mode, mode);
      EXPECT_EQ(state.tl_energy, -R{2} * s.energy_shift);
      R const eta = std::acosh(s.delta), d = R(s.inner_deviation_sign) * std::exp(-s.inner_log_deviation);
      C const z =
          std::exp(-s.outer_log_deviation) * C{std::cos(s.outer_deviation_phase), std::sin(s.outer_deviation_phase)};
      C const u{(eta + d) / R{2}, s.center / R{2}}, w = u + eta + z;
      auto scattering = [&](C a, C b) {
        return std::sinh(a - b + eta) / std::sinh(a - b - eta) * std::sinh(a + b + eta) / std::sinh(a + b - eta);
      };
      auto correction = [&](C v) { return std::abs(v) < std::sqrt(eps) ? v * v / R{6} : std::log(std::sinh(v) / v); };
      C const cz = correction(z), cd = correction(C{d});
      C const sum_ratio = std::sinh(w + u + eta) / std::sinh(w + u - eta);
      C const outer_regular =
          std::sinh(R{2} * eta + z) * sum_ratio * scattering(w, std::conj(u)) * scattering(w, std::conj(w));
      C const outer_drive = std::sinh(w + eta / R{2}) / std::sinh(w - eta / R{2});
      EXPECT_REAL_NEAR(std::log(std::abs(outer_drive)) -
                           (std::log(std::abs(outer_regular)) + s.outer_log_deviation - cz.real()) / (R{2} * R(n)),
                       R{0}, R{512} * eps);
      R const op = R{2} * R(n) * std::arg(outer_drive) - std::arg(outer_regular) + s.outer_deviation_phase + cz.imag();
      EXPECT_REAL_NEAR(std::sin(op / R{2}), R{0}, R{1024} * R(n) * eps);
      C const inner_regular = std::sinh(u - std::conj(u) + eta) / std::sinh(u - std::conj(u) - eta) *
                              std::sinh(R{2} * eta + d) / std::sinh(R{2} * eta + z) * sum_ratio *
                              scattering(u, std::conj(w));
      C const inner_drive = std::sinh(u + eta / R{2}) / std::sinh(u - eta / R{2});
      EXPECT_REAL_NEAR(std::log(std::abs(inner_drive)) - (std::log(std::abs(inner_regular)) - s.outer_log_deviation +
                                                          s.inner_log_deviation + cz.real() - cd.real()) /
                                                             (R{2} * R(n)),
                       R{0}, R{512} * eps);
      R const ip = R{2} * R(n) * std::arg(inner_drive) - std::arg(inner_regular) - s.outer_deviation_phase - cz.imag() +
                   (s.inner_deviation_sign < 0 ? pi : R{0});
      EXPECT_REAL_NEAR(std::sin(ip / R{2}), R{0}, R{2048} * R(n) * eps);
      if (mode == 1)
      {
        EXPECT_GT(state.tl_energy, R{15} / R{7});
        if (n >= 16) EXPECT_LT(state.tl_energy - R{15} / R{7}, R{2} / (R(n) * R(n)));
      }
      if (n == 100000)
      {
        EXPECT_EQ(std::exp(-s.inner_log_deviation), R{0});
        EXPECT_EQ(std::exp(-s.outer_log_deviation), R{0});
      }
    }
}

TYPED_TEST(FourString, BudgetsAndValidation)
{
  using R = TypeParam;
  for (std::size_t budget : {0, 1})
  {
    auto const s = qg::four_string::bound_quartet<R>(16, R{3} / R{2}, 1, {.max_iterations = budget});
    EXPECT_FALSE(s.converged);
    EXPECT_EQ(s.iterations, budget);
    EXPECT_EQ(s.status, qg::SolveStatus::iteration_limit);
    qg::four_string::detail::System<R> system(16, s.delta, s.string_label);
    std::vector<R> x{s.center, s.inner_log_deviation, s.outer_log_deviation, s.outer_deviation_phase};
    EXPECT_EQ(s.residual_norm, system.evaluate(x).norm);
    EXPECT_EQ(s.energy_shift, system.energy_shift(x));
  }
  for (std::size_t n : {0, 6, 7})
    EXPECT_THROW((void)qg::four_string::bound_quartet<R>(n, R{3} / R{2}), std::invalid_argument);
  for (std::size_t mode : {0, 2})
    EXPECT_THROW((void)qg::four_string::bound_quartet<R>(8, R{3} / R{2}, mode), std::invalid_argument);
  for (R delta : {R{0}, R{1}, uni20::numeric_limits<R>::infinity(), uni20::numeric_limits<R>::quiet_NaN()})
    EXPECT_THROW((void)qg::four_string::bound_quartet<R>(8, delta), std::invalid_argument);
  EXPECT_THROW((void)qg::four_string::bound_quartet<R>(8, R{3} / R{2}, 1, {.residual_tolerance = R{0}}),
               std::invalid_argument);
  EXPECT_THROW((void)qg::four_string::bound_quartet<R>(std::numeric_limits<std::size_t>::max(), R{3} / R{2}),
               std::invalid_argument);
}
} // namespace
