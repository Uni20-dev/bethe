// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include "program-options.hpp"
#include "result-output.hpp"
#if UNI20_HAS_MPREAL
#include <mpc.h>
#endif

namespace
{
namespace cli = bethe::cli;
auto program_info()
{
  auto info = cli::program_info("bethe-info", "Inspect this Bethe build and its available arithmetic.");
  info.examples = {{"bethe-info", "Show arithmetic formats and configured providers"},
                   {"bethe-info --format json", "Machine-readable build and precision information"}};
  info.notes = {
      "No arguments prints information; no solver or benchmark is run.",
      "fp80 requires native extended80: 64 significand bits and the extended exponent range.",
      "long-double always means this build's native type; storage bytes may include padding.",
      "bits-N selects N significand bits, not a total storage width: bits-128 differs from fp128.",
      "digits-N requests at least N decimal digits of working precision, not guaranteed result accuracy.",
      "Exponent columns follow numeric_limits: min normal=2^(min_exponent-1), max finite<2^max_exponent.",
      "Configured providers are capabilities; individual operations use Uni20 dispatch. --build-info gives details."};
  return info;
}
int run(cli::DataOutputOptions const& options, int argc, char** argv)
{
  uni20::run_context context(program_info(), {.invocation = std::vector<std::string>(argv, argv + argc)});
  cli::RunReport report(context, "Bethe build and arithmetic");
  report.field("bethe_version", "Version", cli::build_info::version)
      .field("bethe_source", "Bethe source", cli::build_info::revision)
      .field("uni20_source", "Uni20 source", cli::build_info::uni20_revision)
      .field("blas_enabled", "BLAS configured", bool(UNI20_BACKEND_BLAS))
      .field("float80_enabled", "Uni20 float80 configured", bool(UNI20_HAS_FLOAT80))
      .field("float128_enabled", "Uni20 float128 configured", bool(UNI20_HAS_FLOAT128))
      .field("mpreal_enabled", "Uni20 arbitrary-precision real/complex configured", bool(UNI20_HAS_MPREAL))
      .field("mplapack_float128_enabled", "MPLAPACK binary128 configured", bool(UNI20_HAS_MPLAPACK_FLOAT128))
      .field("mplapack_float80_enabled", "MPLAPACK binary80 configured", bool(UNI20_HAS_MPLAPACK_FLOAT80))
      .field("mplapack_mpfr_enabled", "MPLAPACK MPFR configured", bool(UNI20_HAS_MPLAPACK_MPREAL))
      .field("arbitrary_frontends", "Arbitrary-precision frontends",
             cli::has_arbitrary_structure_factor
                 ? "bethe-xxx-pbc, bethe-xxx-obc, bethe-lieb-liniger-pbc, bethe-lieb-liniger-obc, "
                   "bethe-lieb-liniger-dispersion, bethe-richardson, bethe-central-spin, "
                   "bethe-xxx-structure-factor, bethe-su3-pbc, bethe-xxz-pbc, bethe-xxz-obc, "
                   "bethe-hubbard-pbc, bethe-hubbard-obc, bethe-gaudin-yang-pbc, bethe-bose-fermi-pbc, bethe-q-boson-pbc, bethe-su3-dispersion, "
                   "bethe-haldane-shastry-pbc, bethe-sutherland-pbc, bethe-tb-dispersion, bethe-sine-gordon-dispersion, bethe-xxz-qg-obc, bethe-tb-pbc"
             : UNI20_HAS_MPREAL ? "bethe-xxx-pbc, bethe-xxx-obc, bethe-richardson, bethe-central-spin, bethe-su3-dispersion, "
                   "bethe-haldane-shastry-pbc, bethe-sutherland-pbc, bethe-tb-dispersion, bethe-sine-gordon-dispersion, bethe-xxz-qg-obc"
                                : "none enabled in this build")
      .field("precision_scope", "Native precision scope", "all numerical frontends, subject to availability")
      .field("precision_note", "Precision units", "significand bits; fp64/fp80/fp128 name complete formats")
      .field("storage_note", "Storage", "native object bytes, including any padding")
      .result(true, "available_build_information");
#if UNI20_HAS_MPREAL
  report.field("mpfr_version", "MPFR version", mpfr_get_version());
  report.field("mpc_version", "MPC version", mpc_get_version());
#endif
  cli::ResultOutput output(report, options, {"precisions"});
  using cli::column;
  using Optional = std::optional<std::int64_t>;
  output.table(
      "precisions", "Arithmetic choices",
      [&](auto& table) {
        auto native = [&]<typename Real>(std::string name, bool available = true) {
          auto const info = cli::native_arithmetic<Real>();
          table.append(name, available, available ? info.format : std::string("unavailable"),
                       available ? Optional(info.bits) : Optional{},
                       available ? Optional(info.decimal_digits) : Optional{},
                       available ? Optional(info.min_exponent) : Optional{},
                       available ? Optional(info.max_exponent) : Optional{},
                       available ? Optional(info.storage_bytes) : Optional{});
        };
        native.template operator()<double>("fp64");
        native.template operator()<long double>("long-double");
        native.template operator()<long double>("fp80", cli::has_fp80);
#if UNI20_HAS_FLOAT128
        native.template operator()<uni20::float128>("fp128");
#else
        table.append("fp128", false, "unavailable", Optional{}, Optional{}, Optional{}, Optional{}, Optional{});
#endif
        for (std::string const choice : {"bits-N", "digits-N"})
          table.append(choice, bool(UNI20_HAS_MPREAL), "arbitrary binary", Optional{}, Optional{}, Optional{},
                       Optional{}, Optional{});
      },
      column<std::string>("precision"), column<bool>("available"), column<std::string>("format"),
      column<Optional>("significand_bits"), column<Optional>("decimal_digits"), column<Optional>("min_exponent"),
      column<Optional>("max_exponent"), column<Optional>("storage_bytes"));
  output.finish();
  return 0;
}
} // namespace

int main(int argc, char** argv)
{
  cli::DataOutputOptions options;
  return cli::program_main(
      argc, argv, program_info(),
      [&](auto& app) {
        cli::add_data_output_options(
            app, options,
            {{.name = "precisions", .description = "Arithmetic formats and availability", .primary = true}});
      },
      [&](auto&) { return run(options, argc, argv); }, {.empty = uni20::cli::no_arguments::run});
}
