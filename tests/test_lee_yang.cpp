// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/lee_yang.hpp>

namespace
{
namespace model = bethe::lee_yang;
template <typename Real> class LeeYang : public ::testing::Test {};
TYPED_TEST_SUITE(LeeYang, test_support::RealTypes, test_support::PrecisionNames);

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
