// SPDX-License-Identifier: GPL-3.0-or-later
#include "test_support.hpp"
#include <bethe/observable_grid.hpp>
#include <bethe/sine_gordon_particles.hpp>
#include <uni20/async/tbb_scheduler.hpp>

namespace
{
using R = uni20::mpreal;
namespace sg = bethe::sine_gordon;
using Spectrum = sg::ParticleSpectrum<R>;
constexpr sg::Particle soliton{sg::ParticleKind::soliton}, anti{sg::ParticleKind::antisoliton};
constexpr sg::Particle b1{sg::ParticleKind::breather, 1}, b2{sg::ParticleKind::breather, 2};
class SineGordonParticlesMpfr : public ::testing::TestWithParam<int> {
  protected:
    uni20::Precision p() const { return uni20::Precision::bits(GetParam()); }
    R v(int n) const { return R(n, p()); }
    R eps() const { return uni20::numeric_limits<R>::epsilon(p()); }
};
INSTANTIATE_TEST_SUITE_P(WorkingPrecisions, SineGordonParticlesMpfr, ::testing::Values(128, 256, 333));

void finite(R const& x, uni20::Precision p)
{
  EXPECT_FALSE(x.is_exact());
  EXPECT_EQ(x.precision(), p);
  EXPECT_TRUE(uni20::isfinite(x));
}

// At p=2/5 the two squared breather masses are (5 +/- sqrt(5))/2.
// This algebraic reference contains no trigonometric evaluation.
R reference(unsigned n, R momentum)
{
  auto p = momentum.precision();
  R root = uni20::math::sqrt(R(5, p));
  R square = n == 1 ? (5 - root) / 2 : (5 + root) / 2;
  return uni20::math::sqrt(square + momentum * momentum);
}

TEST_P(SineGordonParticlesMpfr, AlgebraicMassesAndRelativisticThresholds)
{
  Spectrum spectrum(v(1), v(2) / 5);
  EXPECT_EQ(spectrum.particles().size(), 4u);
  for (auto particle : {soliton, anti, b1, b2})
  {
    R mass = spectrum.mass(particle);
    finite(mass, p());
    EXPECT_EQ(spectrum.energy(particle, v(0)), mass);
    for (R k : {v(0), v(3) / 7, -v(3), v(12)})
    {
      R energy = spectrum.energy(particle, k);
      finite(energy, p());
      EXPECT_EQ(energy, spectrum.energy(particle, -k));
      EXPECT_REAL_NEAR(energy * energy - k * k, mass * mass, 1024 * eps());
      if (particle.kind == sg::ParticleKind::breather)
        EXPECT_REAL_NEAR(energy, reference(particle.index, k), 32 * eps());
    }
  }
  // Stationarity at a common velocity for unequal masses, checked against
  // direct energy partitions on both sides of the minimizing momenta.
  std::array content{soliton, b1, b2};
  R rest = spectrum.mass(soliton) + spectrum.mass(b1) + spectrum.mass(b2);
  R total = v(7) / 3;
  R first = total * spectrum.mass(soliton) / rest;
  R second = total * spectrum.mass(b1) / rest;
  R third = total - first - second;
  R edge = spectrum.threshold(content, total);
  finite(edge, p());
  EXPECT_REAL_NEAR(edge, spectrum.energy(soliton, first) + spectrum.energy(b1, second) + spectrum.energy(b2, third),
                   64 * eps());
  for (R shift : {-v(1) / 10, v(1) / 10})
    EXPECT_GT(spectrum.energy(soliton, first + shift) + spectrum.energy(b1, second - shift) +
                  spectrum.energy(b2, third),
              edge);
  EXPECT_EQ(spectrum.threshold(content, -total), edge);
}

TEST_P(SineGordonParticlesMpfr, ExponentRangeAndStableSpecies)
{
  Spectrum base(v(1), v(2) / 5);
  for (R scale : {R("1e-7000", p()), R("1e7000", p())})
  {
    Spectrum scaled(scale, v(2) / 5);
    R e = scaled.energy(b1, scale * 3 / 7);
    finite(e, p());
    EXPECT_GT(e, 0);
    EXPECT_REAL_NEAR(e / scale, base.energy(b1, v(3) / 7), 64 * eps());
    std::array content{b1, b2};
    EXPECT_REAL_NEAR(scaled.threshold(content, -scale) / scale, base.threshold(content, -v(1)), 64 * eps());
  }
  R tiny("1e-7000", p());
  Spectrum weak(v(1), tiny);
  R light = weak.mass(b1);
  finite(light, p());
  EXPECT_GT(light, 0);
  EXPECT_REAL_NEAR(light / tiny, uni20::pi<R>.at(p()), 32 * eps());
  EXPECT_THROW(weak.particles(256), std::length_error);
  R half = v(1) / 2, step = 8 * eps();
  EXPECT_TRUE(Spectrum(v(1), half - step).exists(b2));
  EXPECT_FALSE(Spectrum(v(1), half).exists(b2));
  EXPECT_FALSE(Spectrum(v(1), half + step).exists(b2));
  EXPECT_EQ(Spectrum(v(1), v(1)).particles().size(), 2u);
  EXPECT_EQ(Spectrum(v(1), v(2)).particles().size(), 2u);
}

TEST_P(SineGordonParticlesMpfr, ExplicitPrecisionAndDomainChecks)
{
  Spectrum spectrum(v(1), v(2) / 5);
  std::array content{soliton, anti};
  for (R bad : {R(uni20::uninitialized), R(1), R(1, uni20::Precision::bits(GetParam() + 1))})
  {
    EXPECT_THROW((Spectrum{v(1), bad}), std::invalid_argument);
    if (!bad.initialized() || bad.is_exact()) EXPECT_THROW((Spectrum{bad, v(1)}), std::invalid_argument);
    EXPECT_THROW(spectrum.energy(b1, bad), std::invalid_argument);
    EXPECT_THROW(spectrum.threshold(content, bad), std::invalid_argument);
  }
  for (R bad : {v(0), -v(1), R("nan", p()), R("inf", p())})
  {
    EXPECT_THROW((Spectrum{bad, v(1)}), std::invalid_argument);
    EXPECT_THROW((Spectrum{v(1), bad}), std::invalid_argument);
  }
  EXPECT_THROW(spectrum.energy(b1, R("nan", p())), std::invalid_argument);
  EXPECT_THROW(spectrum.energy(b1, R("inf", p())), std::invalid_argument);
  EXPECT_THROW(spectrum.mass({sg::ParticleKind::breather, 3}), std::invalid_argument);
  EXPECT_THROW(spectrum.mass({sg::ParticleKind::soliton, 1}), std::invalid_argument);
  EXPECT_THROW(spectrum.threshold({}, v(0)), std::invalid_argument);
  EXPECT_THROW(spectrum.particles(1), std::length_error);
  mpfr_t raw;
  mpfr_init2(raw, p().bit_count());
  mpfr_set_ui_2exp(raw, 1, mpfr_get_emin() - 1, MPFR_RNDN);
  R smallest(raw, p());
  mpfr_clear(raw);
  EXPECT_THROW(Spectrum(smallest, smallest).mass(b1), std::underflow_error);
}

TEST(SineGordonParticlesMpfrShared, AccuracyGrowsBeyondBinary128)
{
  auto p = uni20::Precision::bits(512);
  for (auto particle : {b1, b2})
  {
    R expected = reference(particle.index, R(7, p) / 10), previous(1, p);
    for (int bits : {128, 256, 333})
    {
      auto q = uni20::Precision::bits(bits);
      Spectrum spectrum(R(1, q), R(2, q) / 5);
      R value = spectrum.energy(particle, R(7, q) / 10);
      R error = uni20::math::abs(value.at(p) - expected);
      EXPECT_LT(error, previous / R("1e15", p));
      EXPECT_LT(error, (32 * uni20::numeric_limits<R>::epsilon(q)).at(p));
      previous = error;
    }
  }
  // The stable species changes within 1e-60 of the pole: binary128 cannot
  // distinguish this from p=1/2, but explicit 333-bit input can.
  auto q = uni20::Precision::bits(333);
  EXPECT_EQ(Spectrum(R(1, q), R(1, q) / 2 - R("1e-60", q)).particles().size(), 4u);
  EXPECT_EQ(Spectrum(R(1, q), R(1, q) / 2 + R("1e-60", q)).particles().size(), 3u);
}

TEST(SineGordonParticlesMpfrShared, OwningMixedPrecisionQueriesAndFailureIsolation)
{
  auto const ambient = mpfr_get_default_prec();
  using Row = std::vector<R>;
  auto evaluate = [](Spectrum const& spectrum, R const& k) {
    Row row;
    auto particles = spectrum.particles();
    for (std::size_t i = 0; i < particles.size(); ++i)
    {
      row.push_back(spectrum.mass(particles[i]));
      row.push_back(spectrum.energy(particles[i], k));
      for (std::size_t j = i; j < particles.size(); ++j)
      {
        std::array content{particles[i], particles[j]};
        row.push_back(spectrum.threshold(content, k));
      }
    }
    return row;
  };
  std::vector<Row> reference;
  for (int bits : {128, 256, 333})
  {
    auto p = uni20::Precision::bits(bits);
    reference.push_back(evaluate(Spectrum(R(7, p), R(2, p) / 5), R(3, p) / 7));
  }
  auto submit = [&](uni20::async::IAsyncScheduler& scheduler, std::vector<Row>& rows, bool fail) {
    struct Request
    {
        Spectrum spectrum;
        R momentum;
    };
    std::vector<Request> requests;
    for (int bits : {128, 256, 333})
    {
      auto p = uni20::Precision::bits(bits);
      requests.push_back({Spectrum(R(7, p), R(2, p) / 5), R(3, p) / 7});
    }
    if (fail) requests.back().momentum = R(1);
    return bethe::evaluate_observable_grid_async(
        requests.size(),
        [requests, evaluate](std::size_t i) { return evaluate(requests[i].spectrum, requests[i].momentum); },
        [&rows](std::size_t i, Row row) {
          EXPECT_EQ(i, rows.size());
          rows.push_back(std::move(row));
        },
        {&scheduler, 128, 2});
  };
  auto check = [&](auto& scheduler) {
    std::vector<Row> rows, failed;
    auto invalid = submit(scheduler, failed, true);
    auto operation = submit(scheduler, rows, false);
    EXPECT_THROW((void)invalid.get_wait(scheduler), std::invalid_argument);
    EXPECT_EQ(operation.get_wait(scheduler), reference.size());
    ASSERT_EQ(rows.size(), reference.size());
    for (std::size_t i = 0; i < rows.size(); ++i)
    {
      ASSERT_EQ(rows[i].size(), reference[i].size());
      for (std::size_t j = 0; j < rows[i].size(); ++j)
      {
        finite(rows[i][j], reference[i][j].precision());
        EXPECT_EQ(rows[i][j], reference[i][j]);
      }
    }
  };
  uni20::async::DebugScheduler deferred({.order = uni20::async::DebugSchedulerOrder::reverse});
  check(deferred);
  for (int workers : {1, 4})
  {
    uni20::async::TbbScheduler scheduler(workers);
    check(scheduler);
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}
} // namespace
