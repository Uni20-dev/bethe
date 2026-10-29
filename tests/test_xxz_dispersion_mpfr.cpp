// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "test_support.hpp"
#include <bethe/xxz_bulk.hpp>
#include <bethe/xxz_dispersion.hpp>

namespace
{
using R = uni20::mpreal;
using uni20::Precision;
namespace model = bethe::xxz;
namespace math = uni20::math;
// Independent mpmath theta constants and original Fourier bulk expressions.
// scripts/reference_xxz_dispersion_mpfr.py: 128/160 digits agree beyond 120.
struct Fixture
{
    char const *delta, *maximum, *gap, *energy, *bulk;
};
Fixture const fixtures[] = {
    {"-0.75",
     "0."
     "4295349028190754789138207216423490913358250178105059686416772841148708848067108071599124530185009435273427907458"
     "3",
     "0.0",
     "0."
     "1269362432493515963814335712872541155041697835628607506088190673314028191907263547288432527341742153720945667380"
     "6",
     "-0."
     "258815227610764176177266609981691875265553453477375827131353522390754832608311998727852409836455285408274498551"
     "2"},
    {"0.5",
     "1."
     "2990381056766579701455847561294042752071039403577854710418552345889497626816000278108596400679364317567196060615"
     "5",
     "0.0",
     "0."
     "3838920094505210415827867884813960331766108015760679786250323112604334647544936358215022668974797000431130389532"
     "2",
     "-0.375"},
    {"1.0078125",
     "1."
     "5748848156539544539484581936286597756506865484420657905392643118229948007458053373257927561116718054541490623043"
     "4",
     "0."
     "00000000000000004394403365495107781556670501508773572180778665838746595180932219401351039734898965153746943785873"
     "8983806493674803337536189273086",
     "0."
     "4654102861898622997698735075598630859007071798687388431349776651290327295096695489007728531730930926707656711683"
     "0",
     "-0."
     "444302644254333433552569624296847819038663213571795387869335226226313347991977210634085688662684897485650098362"
     "2"},
    {"1.125",
     "1."
     "6357164087606463741937888990516701607318323509017464057657935709603139572136287040024358610418714493838630549534"
     "5",
     "0."
     "00030589359195605845759489886791366745112076920288044351145150645919111380091018259969091071938439622954459430474"
     "713649656457060577",
     "0."
     "4833873394903648515587980651567851902084938065816586415499720319125805190113075028062297766269059144833559226303"
     "9",
     "-0."
     "461984157634790441030839697958010732438226735259084766039838672090467237459725107627935584418297456751833792652"
     "2"},
    {"2",
     "2."
     "0704961548071344346591600215573356585753030750010460849619498191053763046322402556131075592201835688657715376215"
     "6",
     "0."
     "19490113571001470727701163515916987988002921827177622392600310333178972920817909479715176261830939674621882141798"
     "702897810621407",
     "0."
     "6395765264796260783381100949669387120888956782479856587113999331502892409226089964965500376219563880234738733584"
     "8",
     "-0."
     "617222045975865643852780656734321406704332851676111761438426802683835507328062013265006748112753675187489636782"
     "2"},
    {"100",
     "51."
     "0024750625000155089746771518275259570445846487037805310423228804276435850618383970435298257822059362473737323586",
     "49."
     "00252506249998425663079042624131989869578411961931394527826925875047927326962648118594903118949048878318490787609"
     "0144139674069",
     "49."
     "1804263673154190878316270259864683636413280312440806444403880489982701918972385145487799346591955588253156820576",
     "-25."
     "002499937500000039064453173826904037454222354833125090736768933342202050475572909212878544502931401909449451456"},
};

TEST(XXZDispersionMpfr, IndependentBandAndBulkAccuracyBeyondBinary128)
{
  auto const oracle = Precision::bits(768);
  R previous(1, oracle);
  auto const ambient = mpfr_get_default_prec();
  for (int bits : {80, 160, 256})
  {
    SCOPED_TRACE(bits);
    auto const p = Precision::bits(bits);
    model::BulkOptions<R> controls(p);
    R error(0, oracle);
    for (auto const& f : fixtures)
    {
      SCOPED_TRACE(f.delta);
      R const delta(f.delta, p);
      model::SpinonDispersion<R> const band(delta);
      auto const bulk = model::bulk_energy_density(delta, R(1, p), controls);
      ASSERT_TRUE(bulk.converged) << int(bulk.status) << ", " << bulk.evaluations;
      std::array<R, 4> const actual{band.maximum_energy(), band.gap(), band.energy(R("0.3", p)), *bulk.energy};
      std::array<char const*, 4> const expected{f.maximum, f.gap, f.energy, f.bulk};
      for (std::size_t i = 0; i < actual.size(); ++i)
      {
        EXPECT_EQ(uni20::precision_of(actual[i]), p);
        R const reference(expected[i], oracle);
        if (reference == 0)
          EXPECT_EQ(actual[i], 0);
        else
          error = std::max(error, math::abs(actual[i].at(oracle) / reference - 1));
      }
      EXPECT_EQ(uni20::precision_of(*bulk.error), p);
      EXPECT_LT(math::abs(bulk.energy->at(oracle) - R(f.bulk, oracle)), bulk.error->at(oracle));
    }
    EXPECT_LT(error, 4096 * uni20::numeric_limits<R>::epsilon(p).at(oracle));
    EXPECT_LT(error, previous * R("1e-15", oracle));
    if (bits == 256) EXPECT_LT(error, R("1e-70", oracle));
    previous = error;
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

TEST(XXZDispersionMpfr, EndpointsFoldingAndPrecisionOwnership)
{
  auto const ambient = mpfr_get_default_prec();
  for (int bits : {80, 160, 256})
  {
    auto const p = Precision::bits(bits);
    R const zero(0, p), one(1, p), pi = bethe::detail::spinon_pi(zero);
    model::SpinonDispersion<R> const band(R(2, p));
    EXPECT_EQ(band.energy(zero), band.gap());
    EXPECT_EQ(band.energy(pi), band.gap());
    EXPECT_EQ(band.energy(pi / 2), band.maximum_energy());
    auto const unfolded = band.continuum(zero);
    auto const folded = band.continuum(zero, bethe::SpinonMomentum::folded);
    EXPECT_EQ(unfolded.lower, 2 * band.gap());
    EXPECT_EQ(unfolded.upper, 2 * band.gap());
    EXPECT_EQ(folded.lower, 2 * band.gap());
    EXPECT_EQ(folded.upper, 2 * band.maximum_energy());
    EXPECT_EQ(uni20::precision_of(folded.lower), p);
    model::BulkOptions<R> controls(p);
    auto const xxx = model::bulk_energy_density(one, one, controls);
    auto const xx = model::bulk_energy_density(zero, one, controls);
    ASSERT_TRUE(xxx.converged);
    ASSERT_TRUE(xx.converged);
    EXPECT_EQ(*xxx.energy, bethe::heisenberg::bulk_energy_density(one));
    EXPECT_EQ(*xx.energy, -one / pi);
    model::SpinonDispersion<R> const gapless(zero);
    EXPECT_EQ(gapless.gap(), zero);
    EXPECT_EQ(uni20::precision_of(gapless.gap()), p);
    model::SpinonDispersion<R> const near(R("1.0001", p));
    EXPECT_GT(near.gap(), zero);
    EXPECT_LT(near.gap(), R("1e-100", p));
    EXPECT_EQ(near.energy(pi), near.gap());
    model::SpinonDispersion<R> const other(R(2, Precision::bits(bits + 32)));
    EXPECT_EQ(uni20::precision_of(band.gap()), p);
    EXPECT_THROW(band.energy(R(1, Precision::bits(bits + 32))), std::invalid_argument);
    EXPECT_GT(other.gap(), 0);
  }
  EXPECT_EQ(mpfr_get_default_prec(), ambient);
}

TEST(XXZDispersionMpfr, BudgetsAndInvalidPrecisionRemainExplicit)
{
  auto const p = Precision::bits(160);
  R const one(1, p);
  model::BulkOptions<R> controls(p);
  controls.max_evaluations = 0;
  for (char const* text : {"-0.75", "1.125", "2"})
  {
    auto const s = model::bulk_energy_density(R(text, p), one, controls);
    EXPECT_FALSE(s.converged);
    EXPECT_FALSE(s.energy);
    EXPECT_FALSE(s.error);
    EXPECT_EQ(s.evaluations, 0u);
    EXPECT_EQ(s.status, model::BulkStatus::evaluation_limit);
  }
  controls.max_evaluations = 200000;
  controls.max_levels = 0;
  EXPECT_EQ(model::bulk_energy_density(R("1.125", p), one, controls).status, model::BulkStatus::mesh_limit);
  for (R const& invalid : {R{}, R(1), R(1, Precision::bits(80)), R("nan", p)})
    EXPECT_THROW(model::bulk_energy_density(invalid, one, controls), std::invalid_argument);
  EXPECT_THROW((void)model::BulkOptions<R>(Precision::exact()), std::invalid_argument);
  EXPECT_THROW((void)model::SpinonDispersion<R>(R(2)), std::invalid_argument);
  EXPECT_THROW((void)model::SpinonDispersion<R>(one, R(1, Precision::bits(80))), std::invalid_argument);
  controls.tolerance = R(1, Precision::bits(80));
  EXPECT_THROW(model::bulk_energy_density(one, one, controls), std::invalid_argument);
}
} // namespace
