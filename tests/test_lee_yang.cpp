// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/lee_yang.hpp>
#include <bethe/lee_yang_excited.hpp>

namespace
{
namespace model = bethe::lee_yang;
template <typename Real> class LeeYang : public ::testing::Test {};
TYPED_TEST_SUITE(LeeYang, test_support::RealTypes, test_support::PrecisionNames);

TEST(LeeYangScan, RegularOneParticleDomain)
{
  double previous_gap = 2, previous_displacement = 1;
  for (int i = 0; i <= 50; ++i)
  {
    double const r = 5 + i / 2.0;
    auto const excited = model::one_particle(1.0, r);
    auto const vacuum = model::ground_state(1.0, r);
    ASSERT_TRUE(excited.converged) << r << " status=" << static_cast<int>(excited.status);
    ASSERT_TRUE(vacuum.converged);
    double const gap = *excited.casimir_energy - *vacuum.casimir_energy;
    EXPECT_GT(gap, 1);
    EXPECT_LT(gap, previous_gap);
    EXPECT_LT(*excited.pole_displacement, previous_displacement);
    previous_gap = gap;
    previous_displacement = *excited.pole_displacement;
  }
}

TYPED_TEST(LeeYang, OneParticleSourceFunctions)
{
  using R = TypeParam;
  R const pi = R{4} * std::atan(R{1}), beta = pi / R{5}, a = std::sqrt(R{3}) / R{2};
  R const eps = uni20::numeric_limits<R>::epsilon();
  EXPECT_REAL_NEAR(model::detail::particle_source(R{0}, beta),
                   -R{2} * std::log((a + std::sin(beta)) / (a - std::sin(beta))), R{128} * eps);
  EXPECT_REAL_NEAR(model::detail::continued_kernel(R{0}, beta),
                   a / pi * std::cos(beta) / (a * a - std::sin(beta) * std::sin(beta)), R{32} * eps);
  for (R x : {R{0}, R{1}, R{5}, R{30}})
  {
    EXPECT_REAL_NEAR(model::detail::continued_kernel(x, R{0}), model::kernel(x), R{32} * eps);
    EXPECT_EQ(model::detail::particle_source(x, R{0}), R{0});
    EXPECT_LE(model::detail::particle_source(x, beta), R{0});
  }
  EXPECT_EQ(model::detail::particle_source(R{1000000}, beta), R{0});
  EXPECT_EQ(model::detail::continued_kernel(R{1000000}, beta), R{0});
}

TYPED_TEST(LeeYang, OneParticleOracleAndNativePrecision)
{
  using R = TypeParam;
  model::OneParticleOptions<R> options;
  options.tolerance = R{1e-10};
  for (auto const& ref : {std::pair{5, "5.146781165267865"},
                          {8, "8.019437400300701"},
                          {10, "10.004545725250384"},
                          {20, "20.00000175700884"},
                          {30, "30.00000000046612"}})
  {
    auto const s = model::one_particle(R{1}, R(ref.first), options);
    ASSERT_TRUE(s.converged) << ref.first << " status=" << static_cast<int>(s.status)
                             << " products=" << s.kernel_products;
    EXPECT_REAL_NEAR(*s.scaling_function, uni20::parse_real<R>(ref.second), R{2e-12});
    EXPECT_LE(s.quantization_residual, options.tolerance);
    EXPECT_LE(s.source_error, options.tolerance / R{24});
    EXPECT_LE(s.kernel_products, options.max_kernel_products);
  }
  options = {};
  auto const a = model::one_particle(R{1}, R{5}, options);
  ASSERT_TRUE(a.converged) << static_cast<int>(a.status) << " work=" << a.kernel_products;
  options.initial_intervals = 48;
  options.initial_cutoff = a.cutoff + R{1} / R{2};
  auto const b = model::one_particle(R{2}, R{5} / R{2}, options);
  ASSERT_TRUE(b.converged) << static_cast<int>(b.status) << " work=" << b.kernel_products << " n=" << b.intervals
                           << " roots=" << b.root_iterations << " source_error=" << uni20::format_real(b.source_error)
                           << " quantization=" << uni20::format_real(b.quantization_residual)
                           << " nonlinear=" << uni20::format_real(b.nonlinear_error);
  EXPECT_REAL_NEAR(*a.scaling_function, *b.scaling_function, options.tolerance);
  EXPECT_REAL_NEAR(*a.casimir_energy * R{2}, *b.casimir_energy, options.tolerance);
  EXPECT_REAL_NEAR(*a.pole_displacement, *b.pole_displacement, options.tolerance);
}

TYPED_TEST(LeeYang, OneParticleFailureContracts)
{
  using R = TypeParam;
  auto check = [](auto const& s, model::Status status) {
    EXPECT_FALSE(s.converged);
    EXPECT_EQ(s.status, status);
    EXPECT_FALSE(s.scaling_function);
    EXPECT_FALSE(s.casimir_energy);
    EXPECT_FALSE(s.beta);
    EXPECT_FALSE(s.pole_displacement);
  };
  model::OneParticleOptions<R> o;
  o.max_root_iterations = 0;
  check(model::one_particle(R{1}, R{5}, o), model::Status::iteration_limit);
  o = {};
  o.max_iterations = 0;
  check(model::one_particle(R{1}, R{5}, o), model::Status::iteration_limit);
  o = {};
  o.max_kernel_products = 1100;
  auto const budget = model::one_particle(R{1}, R{5}, o);
  check(budget, model::Status::work_limit);
  EXPECT_LE(budget.kernel_products, 1100u);
  o = {};
  o.max_intervals = 32;
  check(model::one_particle(R{1}, R{5}, o), model::Status::mesh_limit);
  o = {};
  o.max_cutoffs = 0;
  check(model::one_particle(R{1}, R{5}, o), model::Status::cutoff_limit);
  o = {};
  o.max_cutoffs = 1;
  o.tolerance = R{1e-9};
  check(model::one_particle(R{1}, R{5}, o), model::Status::cutoff_limit);
  o.max_cutoffs = 2;
  o.initial_cutoff = R{1} / R{4};
  o.tolerance = R{1e-4};
  check(model::one_particle(R{1}, R{5}, o), model::Status::cutoff_limit);
  o = {};
  o.initial_cutoff = uni20::numeric_limits<R>::max();
  check(model::one_particle(R{1}, R{5}, o), model::Status::precision_limit);
  o = {};
  o.tolerance = uni20::numeric_limits<R>::epsilon();
  check(model::one_particle(R{1}, R{5}, o), model::Status::precision_limit);
  EXPECT_THROW(model::one_particle(R{1}, R{4}), std::invalid_argument);
  EXPECT_THROW(model::one_particle(R{1}, R{31}), std::invalid_argument);
}

TYPED_TEST(LeeYang, NativeKernel)
{
  using R = TypeParam;
  R const pi = R{4} * std::atan(R{1}), eps = uni20::numeric_limits<R>::epsilon();
  for (R x : {R{0}, R{1} / R{10}, R{1}, R{5}, R{30}})
  {
    R const sinh = std::sinh(x);
    R const direct = std::sqrt(R{3}) / (R{2} * pi) * std::cosh(x) / (sinh * sinh + R{3} / R{4});
    EXPECT_REAL_NEAR(model::kernel(x), direct, R{32} * eps * direct);
    EXPECT_EQ(model::kernel(x), model::kernel(-x));
    EXPECT_GT(model::kernel(x), R{0});
  }
  EXPECT_EQ(model::kernel(R{1000000}), R{0});
  EXPECT_THROW(model::kernel(uni20::numeric_limits<R>::infinity()), std::invalid_argument);
}

TYPED_TEST(LeeYang, IndependentGaussOracle)
{
  using R = TypeParam;
  struct Ref
  {
      int numerator, denominator;
      char const* y;
  };
  for (auto ref :
       {Ref{1, 1000, "-0.20943937151285794"}, Ref{1, 100, "-0.20942648592833846"}, Ref{1, 10, "-0.20835015785667577"},
        Ref{1, 2, "-0.19017376405211875"}, Ref{1, 1, "-0.15320688011013006"}, Ref{2, 1, "-0.0812547352032044"},
        Ref{5, 1, "-0.006408441082412908"}, Ref{10, 1, "-0.000059359269102363974"}})
  {
    auto const s = model::ground_state(R{1}, R(ref.numerator) / R(ref.denominator), {.tolerance = R{1e-11}});
    ASSERT_TRUE(s.converged) << ref.y << " status=" << static_cast<int>(s.status);
    EXPECT_REAL_NEAR(*s.scaling_function, uni20::parse_real<R>(ref.y), R{1e-11});
    EXPECT_GE(s.cutoffs, 2u);
    EXPECT_LE(s.mesh_error, R{1e-11} / R{8});
    EXPECT_LE(s.cutoff_error, R{1e-11} / R{2});
    EXPECT_LE(s.kernel_products, 200000000u);
    EXPECT_LT(*s.scaling_function, R{0});
    EXPECT_GT(*s.effective_central_charge, R{0});
    EXPECT_LT(*s.effective_central_charge, R{2} / R{5});
  }
}

TYPED_TEST(LeeYang, NativeDefaultAndScaling)
{
  using R = TypeParam;
  auto const a = model::ground_state(R{1}, R{1});
  ASSERT_TRUE(a.converged) << static_cast<int>(a.status) << " products=" << a.kernel_products;
  auto const b = model::ground_state(R{2}, R{1} / R{2});
  ASSERT_TRUE(b.converged);
  EXPECT_EQ(a.scaling_function, b.scaling_function);
  EXPECT_EQ(*b.casimir_energy, R{2} * *a.casimir_energy);
  EXPECT_EQ(a.effective_central_charge, b.effective_central_charge);
  test_support::expect_exact(uni20::parse_real<R>(uni20::format_real(*a.scaling_function)), *a.scaling_function,
                             "Lee-Yang Y");
  // A second calculation with a different mesh/cutoff must agree at the
  // selected precision; a double-only solver cannot satisfy fp128 checks.
  model::Options<R> options;
  options.initial_intervals = 48;
  options.initial_cutoff = a.cutoff + R{1} / R{2};
  auto const c = model::ground_state(R{1}, R{1}, options);
  ASSERT_TRUE(c.converged) << static_cast<int>(c.status);
  EXPECT_REAL_NEAR(*a.scaling_function, *c.scaling_function, R{2} * options.tolerance);
}

TYPED_TEST(LeeYang, LimitsAndFailureContracts)
{
  using R = TypeParam;
  auto const uv = model::ground_state(R{1}, R{1} / R{100000}, {.tolerance = R{1e-10}});
  ASSERT_TRUE(uv.converged);
  EXPECT_REAL_NEAR(*uv.effective_central_charge, R{2} / R{5}, R{1e-8});
  auto const ir = model::ground_state(R{1}, R{20}, {.tolerance = R{1e-15}});
  ASSERT_TRUE(ir.converged);
  R const asymptotic = uni20::parse_real<R>("-3.745271025404687e-9");
  EXPECT_REAL_NEAR(*ir.scaling_function / asymptotic, R{1}, R{1e-7});
  auto check = [](auto const& s, model::Status status) {
    EXPECT_FALSE(s.converged);
    EXPECT_EQ(s.status, status);
    EXPECT_FALSE(s.scaling_function);
    EXPECT_FALSE(s.casimir_energy);
    EXPECT_FALSE(s.effective_central_charge);
  };
  check(model::ground_state(R{1}, R{1}, {.max_iterations = 0}), model::Status::iteration_limit);
  check(model::ground_state(R{1}, R{1}, {.max_intervals = 32}), model::Status::mesh_limit);
  check(model::ground_state(R{1}, R{1}, {.max_cutoffs = 0}), model::Status::cutoff_limit);
  check(model::ground_state(R{1}, R{1}, {.tolerance = R{1e-9}, .max_cutoffs = 1}), model::Status::cutoff_limit);
  check(model::ground_state(R{1}, R{1}, {.max_kernel_products = 0}), model::Status::work_limit);
  auto const work = model::ground_state(R{1}, R{1}, {.max_kernel_products = 1100});
  check(work, model::Status::work_limit);
  EXPECT_EQ(work.kernel_products, 1089u);
  check(model::ground_state(R{1}, R{1}, {.tolerance = R{1e-4}, .max_cutoffs = 2, .initial_cutoff = R{1} / R{4}}),
        model::Status::cutoff_limit);
  check(model::ground_state(R{1}, R{1}, {.initial_cutoff = uni20::numeric_limits<R>::max()}),
        model::Status::precision_limit);
  check(model::ground_state(R{1}, R{1}, {.tolerance = uni20::numeric_limits<R>::epsilon()}),
        model::Status::precision_limit);
  EXPECT_THROW(model::ground_state(R{0}, R{1}), std::invalid_argument);
  EXPECT_THROW(model::ground_state(R{1}, R{-1}), std::invalid_argument);
  EXPECT_THROW(model::ground_state(R{1}, uni20::numeric_limits<R>::infinity()), std::invalid_argument);
  EXPECT_THROW(model::ground_state(uni20::numeric_limits<R>::max(), R{2}), std::invalid_argument);
  EXPECT_THROW(model::ground_state(uni20::numeric_limits<R>::min(), uni20::numeric_limits<R>::min()),
               std::invalid_argument);
  EXPECT_THROW(model::ground_state(R{1}, R{1}, {.tolerance = R{0}}), std::invalid_argument);
  EXPECT_THROW(model::ground_state(R{1}, R{1}, {.initial_intervals = 0}), std::invalid_argument);
  EXPECT_THROW(model::ground_state(R{1}, R{1}, {.max_intervals = 16384}), std::invalid_argument);
  EXPECT_THROW(model::ground_state(R{1}, R{1}, {.initial_cutoff = R{0}}), std::invalid_argument);
}
} // namespace
