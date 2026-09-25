// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/xxz_quantum_group_free.hpp>
#include <bit>
#include <complex>

namespace
{
namespace model = bethe::xxz::quantum_group::free;
template <typename Real> class FreeQuantumGroup : public ::testing::Test {};
TYPED_TEST_SUITE(FreeQuantumGroup, test_support::RealTypes, test_support::PrecisionNames);

TYPED_TEST(FreeQuantumGroup, NativeEnergiesAndBlocks)
{
  using R = TypeParam;
  R const eps = uni20::numeric_limits<R>::epsilon();
  auto const four = model::block<R>(4, std::vector<std::size_t>{3}, 1);
  EXPECT_REAL_NEAR(four.energy, -std::sqrt(R{2}) / R{2}, R{16} * eps);
  EXPECT_EQ(four.block_size, 2u);
  EXPECT_EQ(four.down, 2u);
  EXPECT_EQ(model::block<R>(4, std::vector<std::size_t>{1, 3}, 0).energy, R{0});
  for (std::size_t n = 2; n <= 16; ++n)
    for (std::size_t down = 0; down <= n; ++down)
    {
      auto const rows = model::sector<R>(n, down);
      std::size_t dimension = 0, pairs = 0;
      for (auto const& row : rows)
      {
        EXPECT_EQ(row.down, down);
        EXPECT_EQ(row.modes.size() + row.zero_occupation, down);
        dimension += row.block_size;
        pairs += row.block_size == 2;
        test_support::expect_exact(uni20::parse_real<R>(uni20::format_real(row.energy)), row.energy, "free QG I/O");
      }
      EXPECT_EQ(dimension, bethe::detail::bounded_binomial(n, down, 1000000));
      auto const expected_pairs = n % 2 || down == 0 ? 0 : bethe::detail::bounded_binomial(n - 2, down - 1, 1000000);
      EXPECT_EQ(pairs, expected_pairs);
    }
}

TYPED_TEST(FreeQuantumGroup, SectorGroundFormulas)
{
  using R = TypeParam;
  R const pi = R{4} * std::atan(R{1}), eps = uni20::numeric_limits<R>::epsilon();
  for (std::size_t n = 2; n <= 16; ++n)
  {
    auto const rows = model::sector<R>(n, n / 2);
    auto const ground =
        std::min_element(rows.begin(), rows.end(), [](auto const& a, auto const& b) { return a.energy < b.energy; });
    R const angle = pi / (R{2} * R(n));
    R const exact = (R{1} - (n % 2 ? R{1} / std::sin(angle) : R{1} / std::tan(angle))) / R{2};
    EXPECT_REAL_NEAR(ground->energy, exact, R{64} * R(n) * eps);
    EXPECT_EQ(ground->block_size, n % 2 ? 1u : 2u);
  }
}

TEST(FreeQuantumGroup, ValidationAndPreflightBudgets)
{
  EXPECT_THROW(model::sector<>(1, 0), std::invalid_argument);
  EXPECT_THROW(model::sector<>(4, 5), std::invalid_argument);
  EXPECT_THROW(model::block<>(4, std::vector<std::size_t>{2}, 0), std::invalid_argument);
  EXPECT_THROW(model::block<>(4, std::vector<std::size_t>{1, 1}, 0), std::invalid_argument);
  EXPECT_THROW(model::block<>(4, std::vector<std::size_t>{3, 1}, 0), std::invalid_argument);
  EXPECT_THROW(model::block<>(4, std::vector<std::size_t>{0}, 0), std::invalid_argument);
  EXPECT_THROW(model::block<>(4, std::vector<std::size_t>{4}, 0), std::invalid_argument);
  EXPECT_THROW(model::block<>(3, {}, 2), std::invalid_argument);
  EXPECT_THROW(model::block<>(4, {}, 3), std::invalid_argument);
  EXPECT_THROW(model::sector<>(4, 2, {.max_blocks = 3}), std::length_error);
  EXPECT_EQ(model::sector<>(4, 2, {.max_blocks = 4, .max_mode_entries = 4}).size(), 4u);
  EXPECT_THROW(model::sector<>(4, 2, {.max_mode_entries = 3}), std::length_error);
  EXPECT_THROW(model::sector<>(2, 0, {.max_blocks = 0}), std::length_error);
  EXPECT_EQ(model::sector<>(2, 0, {.max_mode_entries = 0}).size(), 1u);
  EXPECT_THROW(model::sector<>(1000000000, 500000000), std::length_error);
  EXPECT_THROW(model::sector<>(1000000000, 1000000000), std::length_error);
}

// Independent spin-basis rank oracle: no free-mode diagonalization and no
// Hermitian eigensolver. At each predicted energy, nullity(H-E) counts blocks,
// nullity((H-E)^2) counts algebraic states; the third power must add none.
using C = std::complex<double>;
using Matrix = std::vector<std::vector<C>>;
Matrix product(Matrix const& a, Matrix const& b)
{
  auto const n = a.size();
  Matrix out(n, std::vector<C>(n));
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t k = 0; k < n; ++k)
      for (std::size_t j = 0; j < n; ++j)
        out[i][j] += a[i][k] * b[k][j];
  return out;
}
std::size_t nullity(Matrix a)
{
  auto const n = a.size();
  for (std::size_t r = 0; r < n; ++r)
  {
    std::size_t pi = r, pj = r;
    for (std::size_t i = r; i < n; ++i)
      for (std::size_t j = r; j < n; ++j)
        if (std::abs(a[i][j]) > std::abs(a[pi][pj]))
        {
          pi = i;
          pj = j;
        }
    if (std::abs(a[pi][pj]) < 1e-9) return n - r;
    std::swap(a[r], a[pi]);
    for (auto& row : a)
      std::swap(row[r], row[pj]);
    for (std::size_t i = r + 1; i < n; ++i)
    {
      auto const factor = a[i][r] / a[r][r];
      for (std::size_t j = r + 1; j < n; ++j)
        a[i][j] -= factor * a[r][j];
      a[i][r] = 0;
    }
  }
  return 0;
}
TEST(FreeQuantumGroup, IndependentSpinBasisJordanNullities)
{
  for (unsigned n = 2; n <= 8; ++n)
    for (unsigned down = 0; down <= n; ++down)
    {
      std::vector<unsigned> basis;
      for (unsigned word = 0; word < (1u << n); ++word)
        if (std::popcount(word) == static_cast<int>(down)) basis.push_back(word);
      Matrix h(basis.size(), std::vector<C>(basis.size()));
      for (std::size_t col = 0; col < basis.size(); ++col)
      {
        auto const word = basis[col];
        h[col][col] = C(0, .5 * (static_cast<int>((word >> (n - 1)) & 1) - static_cast<int>(word & 1)));
        for (unsigned j = 0; j + 1 < n; ++j)
          if (((word >> j) ^ (word >> (j + 1))) & 1)
          {
            auto const other = word ^ (3u << j);
            auto const row = std::lower_bound(basis.begin(), basis.end(), other) - basis.begin();
            h[row][col] += .5;
          }
      }
      auto rows = model::sector<>(n, down);
      std::sort(rows.begin(), rows.end(), [](auto const& a, auto const& b) { return a.energy < b.energy; });
      for (std::size_t begin = 0; begin < rows.size();)
      {
        auto end = begin + 1;
        auto algebraic = rows[begin].block_size;
        while (end < rows.size() && std::abs(rows[end].energy - rows[begin].energy) < 1e-10)
          algebraic += rows[end++].block_size;
        SCOPED_TRACE(::testing::Message() << "N=" << n << " down=" << down << " E=" << rows[begin].energy);
        auto shifted = h;
        for (std::size_t i = 0; i < h.size(); ++i)
          shifted[i][i] -= rows[begin].energy;
        auto const square = product(shifted, shifted);
        EXPECT_EQ(nullity(shifted), end - begin);
        EXPECT_EQ(nullity(square), algebraic);
        EXPECT_EQ(nullity(product(square, shifted)), algebraic);
        begin = end;
      }
    }
}
} // namespace
