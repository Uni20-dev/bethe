// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/haldane_shastry.hpp>
#include <bit>
#include <uni20/linalg/ops/self_adjoint_eigh.hpp>

namespace
{
namespace model = bethe::haldane_shastry;
template <typename Real> class HaldaneShastry : public ::testing::Test {};
TYPED_TEST_SUITE(HaldaneShastry, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(HaldaneShastry, NativeGroundFormulasAndChiralPartners)
{
  using Real = TypeParam;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  for (std::size_t n = 2; n <= 51; ++n)
  {
    auto const levels = model::ground_levels<Real>(n);
    ASSERT_EQ(levels.size(), n % 2 ? 2 : 1);
    Real const expected = -pi * pi / Real{24} * (Real(n) + Real(n % 2 ? -1 : 5) / Real(n));
    for (auto const& s : levels)
    {
      EXPECT_REAL_NEAR(s.energy, expected, Real{64} * eps * Real(n));
      EXPECT_EQ(s.gap, Real{0});
      EXPECT_EQ(s.spinons, n % 2);
      EXPECT_EQ(s.maximum_spin.twice(), n % 2);
      ASSERT_TRUE(s.degeneracy);
      EXPECT_EQ(*s.degeneracy, n % 2 ? 2 : 1);
    }
    if (n % 2) EXPECT_EQ((levels[0].momentum_index + levels[1].momentum_index) % n, 0);
  }
  auto const ferro = model::evaluate<Real>(4, {});
  EXPECT_EQ(ferro.energy_numerator, 60);
  EXPECT_EQ(ferro.degeneracy, 5);
  auto const reducible = model::evaluate<Real>(4, {2});
  EXPECT_EQ(reducible.degeneracy, 4); // A singlet plus a triplet, not spin 3/2.
  EXPECT_EQ(reducible.maximum_spin, uni20::half_int{1});
}

TYPED_TEST(HaldaneShastry, MotifCompletenessAndExactOrdering)
{
  using Real = TypeParam;
  for (std::size_t n = 2; n <= 16; ++n)
  {
    auto const all = model::spectrum<Real>(n);
    ASSERT_TRUE(all.complete);
    ASSERT_EQ(all.levels.size(), all.total_motifs);
    std::uint64_t states = 0;
    for (std::size_t i = 0; i < all.levels.size(); ++i)
    {
      auto const& s = all.levels[i];
      ASSERT_TRUE(s.degeneracy);
      states += *s.degeneracy;
      EXPECT_GE(s.gap, Real{0});
      if (i) EXPECT_LE(all.levels[i - 1].energy_numerator, s.energy_numerator);
    }
    EXPECT_EQ(states, std::uint64_t{1} << n);
    auto const few = model::spectrum<Real>(n, 2);
    ASSERT_TRUE(few.complete);
    ASSERT_EQ(few.levels.size(), 2);
    EXPECT_EQ(few.levels[1].motif, all.levels[1].motif);
  }
}

TYPED_TEST(HaldaneShastry, BudgetsInputsAndCountOverflow)
{
  using Real = TypeParam;
  for (auto n : {std::size_t{0}, std::size_t{1}, model::max_sites + 1})
    EXPECT_THROW(model::ground_levels<Real>(n), std::invalid_argument);
  for (auto motif : {std::vector<std::size_t>{0}, {6}, {3, 2}, {1, 2}, {2, 2}})
    EXPECT_THROW(model::evaluate<Real>(6, motif), std::invalid_argument);
  EXPECT_THROW(model::sector_ground_levels<Real>(6, uni20::from_twice(1)), std::invalid_argument);
  EXPECT_THROW(model::sector_ground_levels<Real>(6, uni20::half_int{4}), std::invalid_argument);
  EXPECT_THROW(model::spectrum<Real>(6, 0), std::invalid_argument);
  auto const no = model::spectrum<Real>(6, {}, 12);
  EXPECT_FALSE(no.complete);
  EXPECT_TRUE(no.levels.empty());
  EXPECT_TRUE(model::spectrum<Real>(6, {}, 13).complete);
  EXPECT_FALSE(model::spectrum<Real>(model::max_sites, 1).complete);
  std::vector<std::size_t> motif;
  for (std::size_t m = 2; m < 200; m += 3)
    motif.push_back(m);
  auto const large = model::evaluate<Real>(200, motif);
  EXPECT_FALSE(large.degeneracy);
  EXPECT_TRUE(uni20::isfinite(large.energy));
  EXPECT_GT(large.gap, Real{0});
  auto const top = model::evaluate<Real>(model::max_sites, {});
  EXPECT_EQ(top.energy_numerator, 999999999999000000LL);
}

std::vector<double> exact_energies(unsigned n, std::optional<unsigned> down = {}, double translation_weight = 0)
{
  unsigned const size = 1u << n;
  std::vector<unsigned> basis, index(size);
  for (unsigned word = 0; word < size; ++word)
    if (!down || std::popcount(word) == int(*down))
    {
      index[word] = basis.size();
      basis.push_back(word);
    }
  uni20::DenseMatrix<double> h(basis.size(), basis.size());
  for (std::size_t i = 0; i < basis.size(); ++i)
    for (std::size_t j = 0; j < basis.size(); ++j)
      h[i, j] = 0;
  double const pi = 4 * std::atan(1.);
  for (std::size_t col = 0; col < basis.size(); ++col)
  {
    auto const word = basis[col];
    for (unsigned i = 0; i < n; ++i)
      for (unsigned j = i + 1; j < n; ++j)
      {
        double const chord = std::sin(pi * double(j - i) / n);
        double const coupling = (pi / n) * (pi / n) / (chord * chord);
        bool const unlike = ((word >> i) ^ (word >> j)) & 1u;
        h[col, col] += coupling * (unlike ? -.25 : .25);
        if (unlike) h[index[word ^ (1u << i) ^ (1u << j)], col] += coupling * .5;
      }
    auto const translated = ((word << 1) & (size - 1)) | (word >> (n - 1));
    h[index[translated], col] += translation_weight / 2;
    h[col, index[translated]] += translation_weight / 2;
  }
  auto const eig = uni20::linalg::eigh(std::move(h));
  std::vector<double> result(basis.size());
  for (std::size_t j = 0; j < basis.size(); ++j)
    result[j] = eig.eigenvalues[j];
  return result;
}

TEST(HaldaneShastryOracles, FullSpinHamiltonianAndTranslation)
{
  for (unsigned n = 2; n <= 8; ++n)
    for (double shift : {0., .137})
    {
      SCOPED_TRACE(::testing::Message() << n << ',' << shift);
      auto const spectrum = model::spectrum<double>(n);
      std::vector<double> expected;
      for (auto const& s : spectrum.levels)
        for (std::uint64_t j = 0; j < *s.degeneracy; ++j)
          expected.push_back(s.energy + shift * std::cos(s.momentum));
      std::sort(expected.begin(), expected.end());
      auto const exact = exact_energies(n, {}, shift);
      ASSERT_EQ(expected.size(), exact.size());
      for (std::size_t i = 0; i < exact.size(); ++i)
        EXPECT_NEAR(expected[i], exact[i], 2.e-11);
    }
}

TEST(HaldaneShastryOracles, EveryMagnetizationSector)
{
  for (unsigned n = 2; n <= 9; ++n)
    for (unsigned down = 0; down <= n; ++down)
    {
      auto const sz = uni20::from_twice(int(n) - 2 * int(down));
      auto const levels = model::sector_ground_levels<double>(n, sz);
      auto const exact = exact_energies(n, down);
      for (auto const& s : levels)
        EXPECT_NEAR(s.energy, exact.front(), 2.e-11);
    }
}

TEST(SpinMultiplets, ClebschGordanAndBudgets)
{
  namespace spin = bethe::spin;
  auto const empty = spin::tensor_product({});
  ASSERT_TRUE(empty.complete);
  ASSERT_EQ(empty.multiplets.size(), 1);
  EXPECT_EQ(empty.multiplets[0].spin, uni20::half_int{0});
  EXPECT_EQ(empty.dimension, 1);
  auto const half = uni20::from_twice(1);
  auto const triple = spin::tensor_product({half, half, half});
  ASSERT_TRUE(triple.complete);
  ASSERT_EQ(triple.multiplets.size(), 2);
  EXPECT_EQ(triple.multiplets[0].spin, half);
  EXPECT_EQ(triple.multiplets[0].multiplicity, 2);
  EXPECT_EQ(triple.multiplets[1].spin, uni20::from_twice(3));
  EXPECT_EQ(triple.multiplets[1].multiplicity, 1);
  EXPECT_EQ(triple.dimension, 8);
  auto const limited = spin::tensor_product({half, half, half}, {.max_updates = triple.updates - 1});
  EXPECT_FALSE(limited.complete);
  EXPECT_EQ(limited.status, spin::DecompositionStatus::work_limit);
  EXPECT_EQ(limited.updates, triple.updates - 1);
  EXPECT_TRUE(limited.multiplets.empty());
  EXPECT_FALSE(limited.dimension);
  EXPECT_TRUE(spin::tensor_product({half, half, half}, {.max_updates = triple.updates}).complete);
  EXPECT_THROW(spin::tensor_product({uni20::half_int{-1}}), std::invalid_argument);
  auto const huge = uni20::from_twice(std::numeric_limits<std::int64_t>::max());
  EXPECT_THROW(spin::tensor_product({huge, half}), std::invalid_argument);
  auto const singleton = spin::tensor_product({huge});
  ASSERT_TRUE(singleton.complete);
  EXPECT_EQ(singleton.multiplets[0].spin, huge);
  EXPECT_EQ(singleton.updates, 1);
  // Total dimension can overflow before the individual irrep counts do.
  auto const many = spin::tensor_product(std::vector<uni20::half_int>(64, half));
  EXPECT_TRUE(many.complete);
  EXPECT_FALSE(many.dimension);
  auto const overflow = spin::tensor_product(std::vector<uni20::half_int>(100, half));
  EXPECT_FALSE(overflow.complete);
  EXPECT_EQ(overflow.status, spin::DecompositionStatus::count_overflow);
  EXPECT_TRUE(overflow.multiplets.empty());
}

TEST(HaldaneShastrySpinContent, MotifsAndHilbertSpace)
{
  auto const reducible = model::spin_decomposition(4, {2});
  ASSERT_TRUE(reducible.complete);
  ASSERT_EQ(reducible.multiplets.size(), 2);
  EXPECT_EQ(reducible.multiplets[0].spin, uni20::half_int{0});
  EXPECT_EQ(reducible.multiplets[1].spin, uni20::half_int{1});
  EXPECT_EQ(reducible.multiplets[0].multiplicity, 1);
  EXPECT_EQ(reducible.multiplets[1].multiplicity, 1);
  EXPECT_EQ(reducible.dimension, 4);
  for (unsigned n = 2; n <= 16; ++n)
  {
    std::vector<std::uint64_t> counts(n + 1);
    for (auto const& level : model::spectrum<double>(n).levels)
    {
      auto const content = model::spin_decomposition(n, level.motif);
      ASSERT_TRUE(content.complete);
      EXPECT_EQ(content.dimension, level.degeneracy);
      ASSERT_FALSE(content.multiplets.empty());
      EXPECT_EQ(content.multiplets.back().spin, level.maximum_spin);
      for (auto const& term : content.multiplets)
        counts[term.spin.twice()] += term.multiplicity;
    }
    // Independently count spin words with fixed Sz, then subtract adjacent weights.
    std::vector<std::uint64_t> weights(n + 3);
    for (unsigned word = 0; word < (1u << n); ++word)
    {
      int const twice = int(n) - 2 * std::popcount(word);
      if (twice >= 0) ++weights[twice];
    }
    for (unsigned twice = n % 2; twice <= n; twice += 2)
      EXPECT_EQ(counts[twice], weights[twice] - weights[twice + 2]);
  }
  for (auto motif : {std::vector<std::size_t>{0}, {6}, {1, 2}, {3, 2}})
    EXPECT_THROW(model::spin_decomposition(6, motif), std::invalid_argument);
  auto const refused = model::spin_decomposition(4, {2}, {.max_updates = 0});
  EXPECT_FALSE(refused.complete);
  EXPECT_TRUE(refused.multiplets.empty());
  auto const large = model::spin_decomposition(model::max_sites, {});
  ASSERT_TRUE(large.complete);
  EXPECT_EQ(large.dimension, model::max_sites + 1);
  EXPECT_EQ(large.updates, 1);
}

TEST(HaldaneShastryOracles, FullMagnetizationResolvedSpectrumAndTranslation)
{
  for (unsigned n = 2; n <= 8; ++n)
    for (unsigned down = 0; down <= n / 2; ++down)
      for (double shift : {0., .137})
      {
        SCOPED_TRACE(::testing::Message() << n << ',' << down << ',' << shift);
        std::vector<double> expected;
        for (auto const& level : model::spectrum<double>(n).levels)
        {
          auto const content = model::spin_decomposition(n, level.motif);
          ASSERT_TRUE(content.complete);
          for (auto const& term : content.multiplets)
            if (term.spin.twice() >= n - 2 * down)
              for (std::uint64_t i = 0; i < term.multiplicity; ++i)
                expected.push_back(level.energy + shift * std::cos(level.momentum));
        }
        std::sort(expected.begin(), expected.end());
        auto const exact = exact_energies(n, down, shift);
        ASSERT_EQ(expected.size(), exact.size());
        for (std::size_t i = 0; i < exact.size(); ++i)
          EXPECT_NEAR(expected[i], exact[i], 2.e-11);
      }
}
} // namespace
