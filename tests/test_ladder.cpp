// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/heisenberg.hpp>
#include <bethe/ladder.hpp>
#include <bethe/su3.hpp>
#include <map>
#include <uni20/linalg/ops/self_adjoint_eigh.hpp>

namespace
{
namespace model = bethe::ladder;
template <typename Real> class Ladder : public ::testing::Test {};
TYPED_TEST_SUITE(Ladder, test_support::RealTypes, test_support::PrecisionNames);

struct Exact
{
    double energy, translation, gap;
};
Exact permutation_ground(unsigned l, model::detail::Shape const& counts)
{
  // Independent color-word matrix, including both bonds at L=2.
  static std::map<model::detail::Shape, Exact> cache;
  if (auto it = cache.find(counts); it != cache.end()) return it->second;
  std::vector<unsigned> basis;
  std::vector<std::size_t> index(1u << (2 * l));
  for (unsigned word = 0; word < index.size(); ++word)
  {
    model::detail::Shape found{};
    for (unsigned j = 0; j < l; ++j)
      ++found[(word >> (2 * j)) & 3u];
    if (found == counts)
    {
      index[word] = basis.size();
      basis.push_back(word);
    }
  }
  uni20::DenseMatrix<double> h(basis.size(), basis.size());
  for (std::size_t i = 0; i < basis.size(); ++i)
    for (std::size_t j = 0; j < basis.size(); ++j)
      h[i, j] = 0;
  for (std::size_t col = 0; col < basis.size(); ++col)
    for (unsigned j = 0; j < l; ++j)
    {
      unsigned const k = (j + 1) % l, diff = ((basis[col] >> (2 * j)) ^ (basis[col] >> (2 * k))) & 3u;
      h[index[basis[col] ^ (diff << (2 * j)) ^ (diff << (2 * k))], col] += 1;
    }
  auto const eig = uni20::linalg::eigh(std::move(h));
  double translation = 0;
  for (std::size_t j = 0; j < basis.size(); ++j)
  {
    unsigned const word = ((basis[j] << 2) & (index.size() - 1)) | (basis[j] >> (2 * (l - 1)));
    translation += eig.eigenvectors[index[word], 0] * eig.eigenvectors[j, 0];
  }
  return cache[counts] = {eig.eigenvalues[0], translation,
                          basis.size() > 1 ? eig.eigenvalues[1] - eig.eigenvalues[0] : 1};
}

double spin_ground(unsigned l, double rung)
{
  // Literal spin-1/2 dot products, not the permutation identity used by BA.
  unsigned const size = 1u << (2 * l);
  uni20::DenseMatrix<double> h(size, size);
  for (unsigned row = 0; row < size; ++row)
    for (unsigned col = 0; col < size; ++col)
      h[row, col] = 0;
  auto dot = [](unsigned word, unsigned a, unsigned b) {
    bool const unlike = ((word >> a) ^ (word >> b)) & 1u;
    return std::array<std::pair<unsigned, double>, 2>{
        {{word, unlike ? -.25 : .25}, {word ^ (1u << a) ^ (1u << b), unlike ? .5 : 0}}};
  };
  for (unsigned col = 0; col < size; ++col)
    for (unsigned j = 0; j < l; ++j)
    {
      unsigned const k = (j + 1) % l;
      for (auto [row, v] : dot(col, 2 * j, 2 * k))
        h[row, col] += v;
      for (auto [row, v] : dot(col, 2 * j + 1, 2 * k + 1))
        h[row, col] += v;
      for (auto [row, v] : dot(col, 2 * j, 2 * j + 1))
        h[row, col] += rung * v;
      for (auto [mid, v] : dot(col, 2 * j, 2 * k))
        for (auto [row, w] : dot(mid, 2 * j + 1, 2 * k + 1))
          h[row, col] += 4 * v * w;
    }
  return uni20::linalg::eigh(std::move(h)).eigenvalues[0];
}

template <uni20::Real Real> void check_branch(model::State<Real> const& state)
{
  auto const& b = state.highest_weight;
  ASSERT_TRUE(b.converged);
  EXPECT_TRUE(model::detail::dominates(b.shape, state.populations));
  using Complex = uni20::complex<Real>;
  auto e = [](Real x, Real w) { return Complex{x, w / Real{2}} / Complex{x, -w / Real{2}}; };
  Real const pi = Real{4} * std::atan(Real{1}),
             tol = Real{1024} * Real(state.rungs) * uni20::numeric_limits<Real>::epsilon();
  Real energy = Real(state.rungs), phase{};
  std::size_t count = state.rungs;
  for (std::size_t a = 0; a < b.rapidities.size(); ++a)
  {
    count -= b.shape[a];
    EXPECT_EQ(b.rapidities[a].size(), count);
    ASSERT_EQ(b.rapidities[a].size(), b.labels[a].size());
    for (std::size_t j = 0; j < count; ++j)
    {
      Real const root = b.rapidities[a][j];
      if (j) EXPECT_GT(root, b.rapidities[a][j - 1]);
      Complex product{Real{-1}, Real{0}};
      if (a == 0)
      {
        for (std::size_t site = 0; site < state.rungs; ++site)
          product *= e(root, Real{-1});
        energy -= Real{1} / (root * root + Real{.25});
        phase += pi - Real{2} * std::atan(Real{2} * root);
      }
      for (Real other : b.rapidities[a])
        product *= e(root - other, Real{2});
      if (a)
        for (Real other : b.rapidities[a - 1])
          product *= e(root - other, Real{-1});
      if (a + 1 < b.rapidities.size())
        for (Real other : b.rapidities[a + 1])
          product *= e(root - other, Real{-1});
      EXPECT_REAL_NEAR(product.real(), Real{1}, tol);
      EXPECT_REAL_NEAR(product.imag(), Real{0}, tol);
    }
  }
  EXPECT_REAL_NEAR(energy, b.energy, tol);
  Real const p = Real{2} * pi * Real(b.momentum_index) / Real(state.rungs);
  EXPECT_REAL_NEAR(std::cos(phase), std::cos(p), tol);
  EXPECT_REAL_NEAR(std::sin(phase), std::sin(p), tol);
  ASSERT_TRUE(state.energy);
  EXPECT_REAL_NEAR(*state.energy,
                   energy - Real(state.rungs) / Real{4} +
                       state.rung_coupling * (Real(state.rungs) / Real{4} - Real(state.singlets)),
                   tol);
}

TYPED_TEST(Ladder, AllSingletSectorsAgainstIndependentED)
{
  using Real = TypeParam;
  for (unsigned l = 2; l <= 7; ++l)
  {
    auto const scan = model::sector_ground_states<Real>(l, Real{.75});
    ASSERT_TRUE(scan.complete) << l << " " << int(scan.status);
    ASSERT_EQ(scan.sectors.size(), l + 1);
    for (auto const& s : scan.sectors)
    {
      SCOPED_TRACE(::testing::Message() << l << "," << s.singlets);
      ASSERT_TRUE(s.converged);
      auto const exact = permutation_ground(l, s.populations);
      EXPECT_REAL_NEAR(s.highest_weight.energy, Real(exact.energy), Real{1} / Real{10000000000});
      if (exact.gap > 1e-7)
        EXPECT_NEAR(std::cos(2 * std::acos(-1.) * double(s.highest_weight.momentum_index) / l), exact.translation,
                    1e-10);
      ASSERT_NO_FATAL_FAILURE(check_branch(s));
    }
  }
}

TEST(LadderOracles, LiteralSpinHamiltonian)
{
  for (unsigned l : {2, 3, 4})
    for (double rung : {-8., -1., 0., 1., 3., 4., 5.})
    {
      SCOPED_TRACE(::testing::Message() << l << "," << rung);
      auto const s = model::ground_state(l, rung);
      ASSERT_TRUE(s.converged);
      EXPECT_NEAR(*s.energy, spin_ground(l, rung), 1e-11);
    }
}

TYPED_TEST(Ladder, DescendantsBeatOwnWeightSeas)
{
  using Real = TypeParam;
  Real const tol = Real{1024} * uni20::numeric_limits<Real>::epsilon();
  for (std::size_t ns : {2, 4})
  {
    auto const s = model::sector_ground_state<Real>(6, ns, Real{0});
    ASSERT_TRUE(s.converged);
    ASSERT_TRUE(s.descendant);
    model::detail::Shape const expected = ns == 4 ? model::detail::Shape{4, 2, 0, 0} : model::detail::Shape{2, 2, 2, 0};
    EXPECT_EQ(s.highest_weight.shape, expected);
    Real const exact = ns == 4 ? Real{1} - std::sqrt(Real{5}) : -Real{1} - std::sqrt(Real{13});
    EXPECT_REAL_NEAR(s.highest_weight.energy, exact, tol);
    if constexpr (uni20::numeric_limits<Real>::digits > 53)
      EXPECT_GT(std::abs(Real(double(exact)) - exact), Real{32} * uni20::numeric_limits<Real>::epsilon());
    auto const own = model::detail::populations(6, ns);
    Real own_best = Real{6};
    for (auto shift : model::detail::shifts(own))
    {
      auto const b = model::detail::solve<Real>(own, shift, Real{32} * uni20::numeric_limits<Real>::epsilon(), 1000);
      ASSERT_TRUE(b.converged);
      own_best = std::min(own_best, b.energy);
    }
    EXPECT_GT(own_best - s.highest_weight.energy, Real{1} / Real{100});
    EXPECT_REAL_NEAR(own_best, ns == 4 ? Real{-1} : -(Real{5} + std::sqrt(Real{17})) / Real{2}, tol);
  }
}

TYPED_TEST(Ladder, ExactLimitsAndLargerScans)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), pi = Real{4} * std::atan(Real{1});
  for (std::size_t l : {2, 3, 4, 7, 12, 24})
  {
    auto const one = model::sector_ground_state<Real>(l, l - 1, Real{2});
    ASSERT_TRUE(one.converged);
    Real const sine = std::sin(pi * Real(l / 2) / Real(l));
    EXPECT_REAL_NEAR(*one.energy, -Real{3} * Real(l) / Real{4} + Real{2} - Real{4} * sine * sine,
                     Real{1024} * Real(l) * eps);
    model::SolverOptions<Real> zero;
    zero.max_iterations = zero.max_branches = 0;
    for (Real rung : {Real{4}, Real{7}})
    {
      auto const product = model::ground_state<Real>(l, rung, zero);
      ASSERT_TRUE(product.converged && product.analytic);
      EXPECT_EQ(product.singlets, l);
      EXPECT_EQ(product.iterations, 0u);
      EXPECT_EQ(product.branches, 0u);
      EXPECT_REAL_NEAR(*product.energy, Real{3} * Real(l) * (Real{1} - rung) / Real{4}, Real{32} * Real(l) * eps);
    }
  }
  for (std::size_t l : {3, 6, 12})
  {
    auto const s = model::sector_ground_state<Real>(l, 0, Real{-2});
    auto const su3 = bethe::su3::ground_state<Real>(l);
    ASSERT_TRUE(s.converged && su3.converged);
    EXPECT_REAL_NEAR(s.highest_weight.energy, su3.energy, Real{1024} * Real(l) * eps);
  }
  for (std::size_t l : {4, 12, 16})
  {
    auto const scan = model::sector_ground_states<Real>(l, Real{0});
    ASSERT_TRUE(scan.complete) << l << " " << int(scan.status);
    for (auto const& s : scan.sectors)
      ASSERT_NO_FATAL_FAILURE(check_branch(s));
    if (l == 4) EXPECT_REAL_NEAR(*scan.sectors[1].energy, Real{-5}, Real{1024} * eps);
  }
}

TYPED_TEST(Ladder, SU2ReductionAndAnalyticJacobian)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (std::size_t l : {4, 6, 12, 24})
  {
    auto const b = model::detail::solve<Real>({l / 2, l / 2, 0, 0}, {}, Real{32} * eps, 1000);
    auto const xxx = bethe::heisenberg::ground_state<Real>(l);
    ASSERT_TRUE(b.converged && xxx.converged);
    EXPECT_REAL_NEAR(b.energy, Real{2} * xxx.energy + Real(l) / Real{2}, Real{1024} * Real(l) * eps);
  }
  model::detail::Shape const shape{4, 3, 2, 1};
  for (auto shift : model::detail::shifts(shape))
  {
    model::detail::System<Real> const sys(shape, shift);
    auto x = sys.seed();
    std::vector<Real> jac;
    (void)sys.evaluate(x, &jac);
    Real const h = std::cbrt(eps), tol = Real{64} * h * h;
    for (std::size_t col = 0; col < x.size(); ++col)
    {
      auto up = x, down = x;
      up[col] += h;
      down[col] -= h;
      auto const p = sys.evaluate(up), m = sys.evaluate(down);
      for (std::size_t row = 0; row < x.size(); ++row)
        EXPECT_REAL_NEAR((p.residual[row] - m.residual[row]) / (Real{2} * h), jac[row * x.size() + col], tol);
    }
  }
}

TYPED_TEST(Ladder, DisplacedSeasAndReflection)
{
  using Real = TypeParam;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  model::detail::Shape const shape{4, 3, 2, 1};
  auto const shifts = model::detail::shifts(shape);
  ASSERT_EQ(shifts.size(), 4u);
  Real lowest = Real{10}, highest = Real{-10};
  for (auto shift : shifts)
  {
    EXPECT_EQ(shift[0], -1); // only the simultaneous-reflection redundancy is removed
    auto const b = model::detail::solve<Real>(shape, shift, Real{32} * eps, 1000);
    ASSERT_TRUE(b.converged);
    for (auto& value : shift)
      value = -value;
    auto const reflected = model::detail::solve<Real>(shape, shift, Real{32} * eps, 1000);
    ASSERT_TRUE(reflected.converged);
    EXPECT_REAL_NEAR(b.energy, reflected.energy, Real{4096} * eps);
    EXPECT_EQ((10 - b.momentum_index) % 10, reflected.momentum_index);
    for (std::size_t a = 0; a < b.rapidities.size(); ++a)
    {
      std::size_t const prev = a == 0 ? 10 : b.rapidities[a - 1].size();
      std::size_t const next = a + 1 == b.rapidities.size() ? 0 : b.rapidities[a + 1].size();
      for (std::size_t j = 0; j < b.rapidities[a].size(); ++j)
      {
        // 2I has parity M_prev+M_next-M_current+1 in the original rational equations.
        auto const parity = std::int64_t(prev + next + 1) - std::int64_t(b.rapidities[a].size());
        EXPECT_EQ((b.labels[a][j].twice() - parity) % 2, 0);
        EXPECT_REAL_NEAR(b.rapidities[a][j], -reflected.rapidities[a][b.rapidities[a].size() - 1 - j],
                         Real{4096} * eps);
      }
    }
    model::State<Real> s;
    s.rungs = 10;
    s.singlets = 4;
    s.populations = shape;
    model::detail::set_energy(s, b);
    ASSERT_NO_FATAL_FAILURE(check_branch(s));
    lowest = std::min(lowest, b.energy);
    highest = std::max(highest, b.energy);
  }
  EXPECT_GT(highest - lowest, Real{1} / Real{100}); // independent displacements are not all equivalent
}

TYPED_TEST(Ladder, BudgetsAndValidation)
{
  using Real = TypeParam;
  model::SolverOptions<Real> options;
  options.max_branches = 0;
  auto none = model::ground_state<Real>(6, Real{0}, options);
  EXPECT_FALSE(none.converged);
  EXPECT_FALSE(none.energy);
  EXPECT_EQ(none.status, model::SolveStatus::branch_limit);
  options.max_branches = 10000;
  options.max_iterations = 0;
  auto candidate = model::ground_state<Real>(6, Real{0}, options);
  EXPECT_FALSE(candidate.converged);
  ASSERT_TRUE(candidate.energy);
  EXPECT_TRUE(candidate.highest_weight.converged); // exact zero-rapidity multiplet, not an unsolved seed
  EXPECT_EQ(candidate.status, model::SolveStatus::iteration_limit);
  EXPECT_EQ(candidate.iterations, 0u);
  options.max_iterations = 1;
  auto one = model::ground_state<Real>(6, Real{0}, options);
  EXPECT_FALSE(one.converged);
  EXPECT_EQ(one.iterations, 1u);
  options.max_branches = 1;
  EXPECT_EQ(model::ground_state<Real>(6, Real{0}, options).status, model::SolveStatus::branch_limit);
  EXPECT_TRUE(model::sector_ground_state<Real>(6, 6, Real{-2}, options).converged);
  EXPECT_THROW((void)model::ground_state<Real>(1, Real{0}), std::invalid_argument);
  EXPECT_THROW((void)model::sector_ground_state<Real>(6, 7, Real{0}), std::invalid_argument);
  EXPECT_THROW((void)model::ground_state<Real>(6, uni20::numeric_limits<Real>::infinity()), std::invalid_argument);
  EXPECT_THROW((void)model::ground_state<Real>(6, uni20::numeric_limits<Real>::quiet_NaN()), std::invalid_argument);
  options.residual_tolerance = Real{0};
  EXPECT_THROW((void)model::ground_state<Real>(6, Real{4}, options), std::invalid_argument);
}
} // namespace
