// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include <bethe/xxz_phantom_wave.hpp>

#include "test_support.hpp"
#include <array>
#include <bethe/polynomial_roots.hpp>
#include <bethe/xxz_coordinate_wave.hpp>
#include <bethe/xxz_odd_continuation.hpp>
#include <bethe/xxz_spin_helix.hpp>
#include <bit>

namespace
{
namespace engine = bethe::xxz::detail;
template <typename Real> class XXZPhantomWave : public ::testing::Test {};
TYPED_TEST_SUITE(XXZPhantomWave, test_support::RealTypes, test_support::PrecisionNames);

std::vector<std::size_t> occupied_sites(unsigned n, unsigned bits)
{
  std::vector<std::size_t> out;
  for (unsigned j = 0; j < n; ++j)
    if (bits & (1U << j)) out.push_back(j);
  return out;
}

unsigned bits_of(std::span<std::size_t const> occupied)
{
  unsigned out = 0;
  for (auto j : occupied)
    out |= 1U << j;
  return out;
}

// Direct spin Hamiltonian, not a Bethe formula. On the boundary bond,
// row 0 <- column N-1 carries exp(-i*phi); exp(i*N*k)=exp(i*phi).
template <typename Real>
std::vector<std::complex<Real>> hamiltonian(unsigned n, unsigned m, Real delta, std::complex<Real> twist,
                                            std::vector<std::complex<Real>> const& wave)
{
  using C = std::complex<Real>;
  std::vector<C> result(1U << n);
  for (unsigned bits = 0; bits < (1U << n); ++bits)
    if (std::popcount(bits) == int(m))
      for (unsigned j = 0; j < n; ++j)
      {
        auto const k = (j + 1) % n;
        bool const anti = ((bits >> j) & 1U) != ((bits >> k) & 1U);
        result[bits] += delta * (anti ? -Real{1} : Real{1}) * wave[bits] / Real{4};
        if (anti)
        {
          C const factor = j + 1 == n ? ((bits & 1U) ? std::conj(twist) : twist) : C{1};
          result[bits] += factor * wave[bits ^ (1U << j) ^ (1U << k)] / Real{2};
        }
      }
  return result;
}

template <typename Real> Real max_magnitude(std::vector<std::complex<Real>> const& v)
{
  Real result{0};
  for (auto x : v)
    result = std::max(result, std::abs(x));
  return result;
}

TYPED_TEST(XXZPhantomWave, VacuumAndIdentityLimits)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  for (unsigned n : {3, 5, 7})
    for (unsigned m = 0; m <= n; ++m)
      for (unsigned winding : {n / 2, n / 2 + 1})
      {
        auto const helix =
            bethe::xxz::periodic_spin_helix<Real>(n, winding, uni20::from_twice(std::int64_t(n) - 2 * std::int64_t(m)));
        auto const phase = engine::helix_phase<Real>(n, winding);
        engine::PhantomDressing<Real> const dressing(n, 0, m, phase), identity(n, m, 0, phase);
        EXPECT_EQ(dressing.terms_per_amplitude(), 1U);
        EXPECT_EQ(identity.terms_per_amplitude(), 1U);
        for (unsigned bits = 0; bits < (1U << n); ++bits)
          if (std::popcount(bits) == int(m))
          {
            auto const occupied = occupied_sites(n, bits);
            auto const result = dressing.amplitude(occupied, [](auto selected) {
              EXPECT_TRUE(selected.empty());
              return C{1};
            });
            EXPECT_LT(std::abs(result - helix.unnormalized_amplitude(occupied)), Real{256} * Real(n) * eps);
            C const arbitrary{Real(bits) / Real{7}, -Real(bits + 1) / Real{11}};
            auto const unchanged = identity.amplitude(occupied, [&](auto selected) {
              EXPECT_EQ(bits_of(selected), bits);
              return arbitrary;
            });
            EXPECT_EQ(unchanged.real(), arbitrary.real());
            EXPECT_EQ(unchanged.imag(), arbitrary.imag());
          }
      }
}

// Check every column of H_(r+p)*D - D*H_r^twist, not selected eigenvectors.
template <typename Real>
Real intertwining_error(unsigned n, unsigned r, unsigned p, std::complex<Real> phase, std::complex<Real> twist)
{
  using C = std::complex<Real>;
  engine::PhantomDressing<Real> const dressing(n, r, p, phase);
  std::vector<unsigned> basis;
  std::vector<std::size_t> index(1U << n);
  for (unsigned bits = 0; bits < (1U << n); ++bits)
    if (std::popcount(bits) == int(r))
    {
      index[bits] = basis.size();
      basis.push_back(bits);
    }
  auto const dim = basis.size();
  std::vector<C> matrix((1U << n) * dim);
  for (unsigned bits = 0; bits < (1U << n); ++bits)
    if (std::popcount(bits) == int(r + p))
      dressing.for_each_term(occupied_sites(n, bits), [&](auto selected, C coefficient) {
        matrix[bits * dim + index[bits_of(selected)]] = coefficient;
      });
  Real error{0};
  for (std::size_t col = 0; col < dim; ++col)
  {
    std::vector<C> source(1U << n), image(1U << n);
    source[basis[col]] = C{1};
    for (unsigned bits = 0; bits < (1U << n); ++bits)
      image[bits] = matrix[bits * dim + col];
    auto const lhs = hamiltonian(n, r + p, phase.real(), C{1}, image);
    auto const acted = hamiltonian(n, r, phase.real(), twist, source);
    for (unsigned bits = 0; bits < (1U << n); ++bits)
      if (std::popcount(bits) == int(r + p))
      {
        C rhs{};
        for (std::size_t j = 0; j < dim; ++j)
          if (acted[basis[j]] != C{}) rhs += matrix[bits * dim + j] * acted[basis[j]];
        error = std::max(error, std::abs(lhs[bits] - rhs));
      }
  }
  return error;
}

TYPED_TEST(XXZPhantomWave, FullHamiltonianIntertwiningIdentity)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  for (auto const sizes : {std::array<unsigned, 3>{5, 1, 1}, {7, 2, 1}, {9, 3, 1}, {9, 2, 2}, {8, 1, 2}, {9, 6, 1}})
    for (int chirality : {-1, 1})
    {
      auto const [n, r, p] = sizes;
      auto const period = n > 2 * r ? n - 2 * r : 2 * r - n;
      Real const angle = Real(chirality) * Real{2} * pi * (Real((period - 1) / 2) / Real(period));
      C const phase{std::cos(angle), std::sin(angle)};
      C const twist{std::cos(-Real{2} * Real(p) * angle), std::sin(-Real{2} * Real(p) * angle)};
      EXPECT_LT(intertwining_error(n, r, p, phase, twist), Real{1024} * Real(n) * eps);
    }
  Real const angle = Real{2} * pi / Real{3};
  C const phase{std::cos(angle), std::sin(angle)};
  C const wrong_twist{std::cos(Real{2} * angle), std::sin(Real{2} * angle)};
  EXPECT_GT(intertwining_error(7, 2, 1, phase, wrong_twist), Real{1} / Real{10});
  C const noncommensurate{std::cos(Real{1}), std::sin(Real{1})};
  C const twist{std::cos(Real{2}), -std::sin(Real{2})};
  engine::PhantomDressing<Real> const invalid_identity(7, 2, 1, noncommensurate);
  EXPECT_GT(invalid_identity.commensurability_error(), Real{1});
  EXPECT_GT(intertwining_error(7, 2, 1, noncommensurate, twist), Real{1} / Real{10});
}

TYPED_TEST(XXZPhantomWave, NonzeroReducedEigenstateCanBeAnnihilated)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon(), d = -Real{1} / Real{2};
  C const q{d, std::sqrt(Real{3}) / Real{2}}; // q^3=1
  engine::PhantomDressing<Real> const dressing(7, 2, 1, q);
  auto finite = [&](std::span<std::size_t const> selected) {
    return engine::phantom_phase_power(q, selected[0] + selected[1]);
  };
  std::vector<C> input(128), output(128);
  for (unsigned bits = 0; bits < 128; ++bits)
  {
    auto const occupied = occupied_sites(7, bits);
    if (occupied.size() == 2) input[bits] = finite(occupied);
    if (occupied.size() == 3) output[bits] = dressing.amplitude(occupied, finite);
  }
  auto action = hamiltonian(7, 2, d, std::conj(q * q), input);
  for (unsigned bits = 0; bits < 128; ++bits)
    action[bits] -= Real{7} * d * input[bits] / Real{4};
  EXPECT_LT(max_magnitude(action), Real{256} * eps);
  EXPECT_GT(max_magnitude(input), Real{9} / Real{10});
  // D maps this state to q^(sum J)*(1+q^2+q^4)=0.
  EXPECT_LT(max_magnitude(output), Real{128} * eps);
}

TYPED_TEST(XXZPhantomWave, ContinuedTwoFiniteRootStatesLiftNontrivially)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  for (unsigned n : {7, 9, 11, 13})
  {
    SCOPED_TRACE(::testing::Message() << "N=" << n);
    unsigned const m = n / 2, p = m - 2;
    Real const d = n == 7 ? -Real{1} / Real{2} : -std::cos(pi / Real(n - 4));
    auto const branch = engine::continue_odd_polynomial(n, d, uni20::from_twice(std::int64_t{1}));
    ASSERT_TRUE(branch.equations_converged);
    engine::PolynomialBetheSystem<Real> const system(n, m, branch.center, branch.coordinate_scale);
    auto const reduced = engine::reduce_phantom_polynomial<Real>(system, branch.coefficients, d, p, -1);
    ASSERT_EQ(reduced.status, engine::PhantomReductionStatus::regular_reduced_equations);
    auto const& c = reduced.finite_coefficients;
    ASSERT_EQ(c.size(), 2U);
    Real const discriminant = c[1] * c[1] - Real{4} * c[0];
    C const root = discriminant >= Real{0} ? C{std::sqrt(discriminant), 0} : C{0, std::sqrt(-discriminant)};
    C const z1 = branch.center + branch.coordinate_scale * (-c[1] + root) / Real{2};
    C const z2 = branch.center + branch.coordinate_scale * (-c[1] - root) / Real{2}, imaginary{0, 1};
    C const v1 = -(Real{1} - imaginary * z1) / (Real{1} + imaginary * z1);
    C const v2 = -(Real{1} - imaginary * z2) / (Real{1} + imaginary * z2);
    C const a = Real{1} + v1 * v2 - Real{2} * d * v1;
    C const b = -(Real{1} + v1 * v2 - Real{2} * d * v2);
    Real const normalization = std::max({Real{1}, std::abs(a), std::abs(b)});
    ASSERT_GT(std::max(std::abs(a), std::abs(b)), Real{64} * eps);
    engine::CoordinateBetheWave<Real> const finite_wave(n, d, std::array<C, 2>{v1, v2});
    auto finite = [&](std::span<std::size_t const> selected) {
      auto const result = finite_wave.evaluate(selected);
      C const explicit_wave =
          (a * engine::phantom_phase_power(v1, selected[0]) * engine::phantom_phase_power(v2, selected[1]) +
           b * engine::phantom_phase_power(v2, selected[0]) * engine::phantom_phase_power(v1, selected[1])) /
          normalization;
      EXPECT_LT(std::abs(result.value - explicit_wave), Real{32} * eps * std::max(Real{1}, result.absolute_term_sum));
      return result.value;
    };
    C const q{d, -std::sqrt(Real{1} + d) * std::sqrt(Real{1} - d)};
    engine::PhantomDressing<Real> const dressing(n, 2, p, q);
    std::vector<C> input(1U << n), output(1U << n);
    for (unsigned bits = 0; bits < (1U << n); ++bits)
    {
      auto const occupied = occupied_sites(n, bits);
      if (occupied.size() == 2) input[bits] = finite(occupied);
      if (occupied.size() == m) output[bits] = dressing.amplitude(occupied, finite);
    }
    C const twist = reduced.equation_rotation * reduced.equation_rotation;
    auto source_action = hamiltonian(n, 2, d, twist, input);
    auto target_action = hamiltonian(n, m, d, C{1}, output);
    for (unsigned bits = 0; bits < (1U << n); ++bits)
    {
      source_action[bits] -= branch.energy * input[bits];
      target_action[bits] -= branch.energy * output[bits];
    }
    ASSERT_GT(max_magnitude(input), Real{1} / Real{1000});
    ASSERT_GT(max_magnitude(output), Real{1} / Real{1000});
    EXPECT_LT(max_magnitude(source_action) / max_magnitude(input), Real{8192} * Real(n) * eps);
    EXPECT_LT(max_magnitude(target_action) / max_magnitude(output), Real{8192} * Real(n) * eps);
  }
}

TYPED_TEST(XXZPhantomWave, GeneralRecoveredFiniteRootsLiftNontrivially)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  Real const pi = Real{4} * std::atan(Real{1}), eps = uni20::numeric_limits<Real>::epsilon();
  for (auto const sizes : {std::array<unsigned, 3>{9, 3, 1}, {11, 4, 1}, {11, 3, 2}, {13, 4, 2}, {13, 3, 3}})
  {
    auto const [n, r, p] = sizes;
    auto const m = r + p;
    SCOPED_TRACE(::testing::Message() << "N=" << n << " r=" << r << " p=" << p);
    Real const d = p == 1 ? -Real{1} / Real{2} : -std::cos(pi / Real(n - 2 * r));
    auto const branch = engine::continue_odd_polynomial(n, d, uni20::from_twice(std::int64_t(n - 2 * m)));
    ASSERT_TRUE(branch.equations_converged);
    engine::PolynomialBetheSystem<Real> const system(n, m, branch.center, branch.coordinate_scale);
    auto const reduced = engine::reduce_phantom_polynomial<Real>(system, branch.coefficients, d, p, -1);
    ASSERT_EQ(reduced.status, engine::PhantomReductionStatus::regular_reduced_equations);
    auto const recovered = bethe::detail::recover_polynomial_roots<Real>(reduced.finite_coefficients);
    ASSERT_EQ(recovered.status, bethe::detail::PolynomialRootStatus::resolved);
    std::vector<C> v;
    for (C x : recovered.roots)
    {
      C const z = branch.center + branch.coordinate_scale * x, imaginary{0, 1};
      v.push_back(-(Real{1} - imaginary * z) / (Real{1} + imaginary * z));
    }
    engine::CoordinateBetheWave<Real> const wave(n, d, v);
    C const q{d, -std::sqrt(Real{1} + d) * std::sqrt(Real{1} - d)};
    engine::PhantomDressing<Real> const dressing(n, r, p, q);
    std::vector<C> input(1U << n), output(1U << n);
    Real absolute_terms{0};
    for (unsigned bits = 0; bits < (1U << n); ++bits)
      if (std::popcount(bits) == int(r))
      {
        auto const amplitude = wave.evaluate(occupied_sites(n, bits));
        input[bits] = amplitude.value;
        absolute_terms = std::max(absolute_terms, amplitude.absolute_term_sum);
      }
    for (unsigned bits = 0; bits < (1U << n); ++bits)
      if (std::popcount(bits) == int(m))
        output[bits] =
            dressing.amplitude(occupied_sites(n, bits), [&](auto selected) { return input[bits_of(selected)]; });
    ASSERT_GT(max_magnitude(input), Real{64} * eps * absolute_terms);
    ASSERT_GT(max_magnitude(output), Real{64} * eps * Real(dressing.terms_per_amplitude()) * max_magnitude(input));
    auto source_action = hamiltonian(n, r, d, reduced.equation_rotation * reduced.equation_rotation, input);
    auto target_action = hamiltonian(n, m, d, C{1}, output);
    for (unsigned bits = 0; bits < (1U << n); ++bits)
    {
      source_action[bits] -= branch.energy * input[bits];
      target_action[bits] -= branch.energy * output[bits];
    }
    EXPECT_LT(max_magnitude(source_action) / max_magnitude(input), Real{32768} * Real(n) * eps);
    EXPECT_LT(max_magnitude(target_action) / max_magnitude(output), Real{32768} * Real(n) * eps);
  }
}

TYPED_TEST(XXZPhantomWave, BudgetAndValidationBeforeCallbacks)
{
  using Real = TypeParam;
  using C = std::complex<Real>;
  EXPECT_THROW((engine::PhantomDressing<Real>(10, 2, 3, C{1}, 9)), std::length_error);
  EXPECT_THROW((engine::PhantomDressing<Real>(100, 34, 34, C{1}, std::numeric_limits<std::size_t>::max())),
               std::length_error);
  EXPECT_THROW((engine::PhantomDressing<Real>(5, 0, 0, C{1}, 0)), std::length_error);
  EXPECT_THROW((engine::PhantomDressing<Real>(1, 0, 0, C{1})), std::invalid_argument);
  EXPECT_THROW((engine::PhantomDressing<Real>(5, 3, 3, C{1})), std::invalid_argument);
  EXPECT_THROW((engine::PhantomDressing<Real>(5, 1, 1, C{0})), std::invalid_argument);
  EXPECT_THROW((engine::PhantomDressing<Real>(5, 1, 1, C{uni20::numeric_limits<Real>::quiet_NaN(), 0})),
               std::invalid_argument);
  engine::PhantomDressing<Real> const dressing(10, 2, 3, C{1}, 10);
  EXPECT_EQ(dressing.terms_per_amplitude(), 10U);
  std::size_t calls = 0;
  auto finite = [&](auto) {
    ++calls;
    return C{1};
  };
  for (auto const& occupied : {std::vector<std::size_t>{0, 1}, {0, 1, 2, 2, 4}, {0, 1, 2, 3, 10}, {0, 2, 1, 3, 4}})
    EXPECT_THROW(dressing.amplitude(occupied, finite), std::invalid_argument);
  EXPECT_EQ(calls, 0U);
  std::array<std::size_t, 5> const occupied{0, 1, 2, 3, 4};
  auto const sum = dressing.amplitude(occupied, finite);
  EXPECT_EQ(sum.real(), Real{10});
  EXPECT_EQ(sum.imag(), Real{0});
  EXPECT_EQ(calls, 10U);
  // Integer phase exponents remain exact even near the supported N limit.
  auto const n = std::size_t(std::numeric_limits<std::int64_t>::max() / 4);
  engine::PhantomDressing<Real> const large(n, 1, 2, C{0, 1});
  std::array<std::size_t, 3> const far{n - 3, n - 2, n - 1};
  auto const exact = large.amplitude(far, [](auto) { return C{1}; });
  EXPECT_EQ(exact.real(), Real{1});
  EXPECT_EQ(exact.imag(), Real{0});
  EXPECT_THROW(dressing.amplitude(occupied, [](auto) { return C{uni20::numeric_limits<Real>::infinity(), 0}; }),
               std::overflow_error);
}
} // namespace
