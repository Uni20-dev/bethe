// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
// Opt-in arithmetic probe. NDJSON outputs and per-point inputs belong in store.
#include <bethe/detail/xxx_spinon_proposal.hpp>
#include <bethe/detail/xxx_spinon_surface.hpp>
#include <bethe/sampling.hpp>
#include <bethe/xxx_six_spinon_structure_factor.hpp>
#include <bethe/xxx_spinon_density.hpp>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace math = uni20::math;
namespace model = bethe::heisenberg;
using Clock = std::chrono::steady_clock;

template <class R> R parse(std::string const& text, R const& one)
{
  if constexpr (uni20::has_runtime_precision_v<R>)
    return R(text, one.precision());
  else
    return uni20::parse_real<R>(text);
}
template <class R> void field(std::ostream& out, char const* name, R const& value)
{
  out << ',' << std::quoted(name) << ':' << std::quoted(uni20::format_scalar(value));
}

template <class R> void geometry(R const& one)
{
  using Proposal = model::detail::SpinonJointProposal<R>;
  using Kind = model::detail::SpinonProposalKind;
  R const pi = 4 * math::atan(one);
  Proposal proposal(1040 * one, one / 2,
                    {{Kind::background, one / 2, one, one},
                     {Kind::one_pair, one / 4, 4 * one, one / 2},
                     {Kind::two_pairs, one / 4, 8 * one, one}});
  R const q = parse("1.5707963267948966192313216916397514420985846996875529104874722961539082031431045", one);
  struct Query
  {
      std::string name;
      R q, omega;
  };
  std::vector<Query> queries{{"pi2_w2", q, 2 * one},    {"pi2_w1p6", q, parse("1.6", one)},
                             {"pi2_w3", q, 3 * one},    {"pi8_w0p7", pi / 8, parse("0.7", one)},
                             {"pi_w0p25", pi, one / 4}, {"pi_w2", pi, 2 * one},
                             {"pi_w4", pi, 4 * one}};
  auto emit = [&](std::string const& id, Query const& query, std::array<R, 4> const& free) {
    for (int section = 0; section < 3; ++section)
    {
      std::ostringstream out;
      out << "{\"id\":" << std::quoted(id + "_" + query.name + "_s" + std::to_string(section));
      field(out, "q", query.q);
      field(out, "omega", query.omega);
      try
      {
        auto const point =
            model::detail::SpinonEnergySurface<R>(6, 1040 * one, query.q, section, query.omega).solve(free);
        out << ",\"status\":" << std::quoted(point.inside ? "inside" : "outside")
            << ",\"working_bits\":" << point.working_bits << ",\"refinements\":" << point.refinements;
        if (point.inside)
        {
          out << ",\"branches\":[";
          for (int branch = 0; branch < 2; ++branch)
          {
            if (branch) out << ',';
            std::array<R, 6> beta;
            std::copy(point.beta[branch].begin(), point.beta[branch].end(), beta.begin());
            auto const weight = proposal.balanced_log_jacobian(beta);
            out << "{\"beta\":[";
            for (int j = 0; j < 6; ++j)
            {
              if (j) out << ',';
              out << std::quoted(uni20::format_scalar(beta[j]));
            }
            out << ']';
            field(out, "log_weight", weight);
            out << '}';
          }
          out << ']';
        }
      }
      catch (std::exception const& error)
      {
        // Do not retain a partial object if weighting failed after a solve.
        out.str("");
        out.clear();
        out << "{\"id\":" << std::quoted(id + "_" + query.name + "_s" + std::to_string(section))
            << ",\"status\":\"error\",\"reason\":" << std::quoted(error.what());
      }
      std::cout << out.str() << "}\n";
    }
  };
  for (unsigned replicate = 0; replicate < 4; ++replicate)
  {
    bethe::MonteCarloStream<5> stream(7303, replicate, 52);
    for (unsigned sample = 0; sample < 16; ++sample)
    {
      auto const free = proposal.draw(stream.next().at_precision(one)).beta;
      auto const id = "draw_" + std::to_string(replicate) + "_" + std::to_string(sample);
      for (auto const& query : queries)
        emit(id, query, free);
    }
  }
  std::array<std::array<int, 4>, 4> const stress{{{-1039, -1037, -1033, -1029},
                                                  {1029, 1033, 1037, 1039},
                                                  {-1039, -1037, 1033, 1039},
                                                  {-1039, -1037, -1033, 1039}}};
  for (std::size_t i = 0; i < stress.size(); ++i)
  {
    std::array<R, 4> free;
    for (int j = 0; j < 4; ++j)
      free[j] = stress[i][j] * one;
    emit("stress_" + std::to_string(i), queries[0], free);
  }
  std::array<R, 4> const free{-1039 * one, -1037 * one, -1033 * one, -1029 * one};
  emit("stress_edge", {"pi2_edge", q, pi * math::sin(q / 2) - math::ldexp(one, -40)}, free);
}

template <class R> void density(R const& one, int terms, int order)
{
  model::XXXSpinonDensity<R> context(uni20::precision_of(one));
  std::string id, log_weight;
  while (std::cin >> id >> log_weight)
  {
    std::array<std::string, 6> texts;
    for (auto& text : texts)
      if (!(std::cin >> text)) throw std::invalid_argument("incomplete density input");
    for (int control = 0; control < 2; ++control)
    {
      auto const start = Clock::now();
      std::ostringstream out;
      out << "{\"id\":" << std::quoted(id) << ",\"control\":" << control;
      try
      {
        std::array<R, 6> beta;
        for (int j = 0; j < 6; ++j)
          beta[j] = context.parse(texts[j]);
        auto value = context.evaluate(beta, {terms + 16 * control, order + 4 * control});
        out << ",\"status\":\"value\"";
        field(out, "value", value.value);
        field(out, "contour_tail", value.contour_tail_error);
        field(out, "arithmetic_indicator", value.arithmetic_indicator);
        field(out, "log_scale", value.log_scale);
      }
      catch (model::SpinonDensityError const& error)
      {
        out << ",\"status\":\"error\",\"failure\":" << int(error.failure())
            << ",\"reason\":" << std::quoted(error.what());
      }
      out << ",\"seconds\":" << std::chrono::duration<double>(Clock::now() - start).count() << '}';
      std::cout << out.str() << std::endl;
    }
  }
}

template <class R> void checked(R const& one, int low_bits)
{
  using Kind = model::detail::SpinonProposalKind;
  model::detail::SpinonJointProposal<R> proposal(1040 * one, one / 2,
                                                 {{Kind::background, one / 2, one, one},
                                                  {Kind::one_pair, one / 4, 4 * one, one / 2},
                                                  {Kind::two_pairs, one / 4, 8 * one, one}});
  model::SpinonRecoveryPolicy policy{
      "1e-6", "0", {{low_bits, uni20::numeric_limits<R>::digits(one), {{{1040, 32}, {1056, 36}}}}}, 60};
  model::ThermodynamicSixSpinonStructureFactor<R, R> spectrum(std::move(proposal), policy, one);
  R const q = parse("1.5707963267948966192313216916397514420985846996875529104874722961539082031431045", one);
  auto const point = spectrum.at(q, 2 * one);
  for (unsigned replicate = 0; replicate < 4; ++replicate)
  {
    bethe::MonteCarloStream<5> stream(7303, replicate, 52);
    for (unsigned sample = 0; sample < 16; ++sample)
    {
      auto const value = point.observation(stream.next());
      std::ostringstream out;
      out << "{\"id\":" << std::quoted("draw_" + std::to_string(replicate) + "_" + std::to_string(sample))
          << ",\"supported_roots\":" << value.supported_roots << ",\"density_attempts\":" << value.density_attempts;
      if (value.values)
      {
        out << ",\"status\":\"value\"";
        field(out, "value", (*value.values)[0]);
        field(out, "contour_diagnostic", (*value.values)[1]);
        field(out, "arithmetic_diagnostic", (*value.values)[2]);
      }
      else
        out << ",\"status\":\"error\",\"reason\":" << std::quoted(value.reason);
      std::cout << out.str() << '}' << std::endl;
    }
  }
}

int main(int argc, char** argv)
try
{
  if (argc < 3)
    throw std::invalid_argument(
        "usage: probe geometry|density|checked bits-N|fp64|extended80|fp128 [terms order | low_bits]");
  std::string const mode = argv[1], kind = argv[2];
  auto run = [&](auto const& one) {
    if (mode == "geometry")
      geometry(one);
    else if (mode == "density")
      density(one, argc > 3 ? std::stoi(argv[3]) : 1040, argc > 4 ? std::stoi(argv[4]) : 32);
    else if (mode == "checked")
    {
      if constexpr (uni20::has_runtime_precision_v<std::remove_cvref_t<decltype(one)>>)
        checked(one, argc > 3 ? std::stoi(argv[3]) : 80);
      else
        throw std::invalid_argument("checked mode requires MPFR");
    }
    else
      throw std::invalid_argument("unknown mode");
  };
  if (kind == "fp64")
    run(1.);
  else if (kind == "extended80")
    run(1.L);
#if UNI20_HAS_FLOAT128
  else if (kind == "fp128")
    run(uni20::float128{1});
#endif
#if UNI20_HAS_MPREAL
  else if (kind.starts_with("bits-"))
    run(uni20::mpreal(1, uni20::Precision::bits(std::stoi(kind.substr(5)))));
#endif
  else
    throw std::invalid_argument("unsupported precision");
}
catch (std::exception const& error)
{
  std::cerr << error.what() << '\n';
  return 1;
}
