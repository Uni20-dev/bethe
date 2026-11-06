// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include "program-options.hpp"
#include <array>
#include <bethe/momentum_input.hpp>
#include <iostream>

namespace bethe::cli
{
struct MomentumArguments
{
    std::optional<std::string> radians, from_pi, beta;
    std::string beta_tolerance = "1e-12";
    explicit operator bool() const { return radians || from_pi || beta; }
    MomentumInput input() const
    {
      if (int(bool(radians)) + int(bool(from_pi)) + int(bool(beta)) != 1)
        throw std::invalid_argument("choose exactly one of --momentum, --momentum-from-pi and --momentum-beta");
      if (radians) return {MomentumInputKind::radians, *radians};
      if (from_pi) return {MomentumInputKind::distance_from_pi, *from_pi};
      return {MomentumInputKind::beta, *beta};
    }
    template <uni20::Real R> ResolvedMomentum<R> resolve(R const& exemplar) const
    {
      return input().resolve(exemplar, beta_tolerance);
    }
};

inline std::array<CLI::Option*, 3> add_momentum_options(CLI::App& app, MomentumArguments& a)
{
  std::array options{
      text_option(app, "--momentum", a.radians, "Momentum in [0,2*pi], radians/site; also pi or 2*pi"),
      text_option(app, "--momentum-from-pi", a.from_pi, "Nonnegative distance d in q=pi-d, radians/site"),
      text_option(app, "--momentum-beta", a.beta, "Logarithmic external momentum b=log(tan(q/2)), 0<q<pi")};
  for (auto* option : options)
    for (auto* other : options)
      if (option != other) option->excludes(other);
  option(app, "--momentum-beta-tolerance", a.beta_tolerance,
         "Warning threshold for estimated absolute beta conversion error; not a kernel tolerance")
      ->capture_default_str();
  return options;
}

template <uni20::Real R> R momentum_radians(ResolvedMomentum<R> const& value)
{
  R const result = value.require_radians();
  if (value.beta_accuracy_warning)
    std::cerr << "warning: momentum input to working beta loses approximately " << value.beta_absolute_error
              << " in beta; increase working precision\n";
  if (value.scalar_accuracy_warning)
    std::cerr << "warning: scalar momentum geometry loses approximately " << value.scalar_beta_absolute_error
              << " in beta; increase geometry precision\n";
  return result;
}

template <class Report, uni20::Real R>
void report_momentum_input(Report& report, MomentumArguments const& args, ResolvedMomentum<R> const& value)
{
  auto const input = args.input();
  auto const kind = input.kind == MomentumInputKind::radians            ? "radians"
                    : input.kind == MomentumInputKind::distance_from_pi ? "distance_from_pi"
                                                                        : "beta";
  report.field("momentum_input_kind", "Momentum input kind", std::string(kind))
      .field("momentum_input_text", "Momentum input", input.text)
      .field("momentum_beta", "Folded momentum beta", value.coordinate.beta())
      .field("momentum_reflected", "Reflected from [0,pi]", value.coordinate.reflected())
      .field("momentum_input_bits", "Momentum input precision (bits)", value.input_bits)
      .field("momentum_input_mpfr", "MPFR input conversion", value.mpfr_input)
      .field("momentum_beta_error", "Estimated beta conversion error", value.beta_absolute_error)
      .field("momentum_scalar_beta_error", "Estimated scalar-q beta error", value.scalar_beta_absolute_error)
      .field("momentum_input_warning", "Momentum conversion warning",
             value.beta_accuracy_warning || value.scalar_accuracy_warning);
}
} // namespace bethe::cli
