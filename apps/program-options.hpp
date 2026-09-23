// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include "bethe-build-info.hpp"
#include "data-output.hpp"
#include <bethe/citations.hpp>
#include <limits>
#include <uni20/cli/cli.hpp>

namespace bethe::cli
{
// A separate informational request: let it escape Uni20's ordinary parse-error
// handling, just as help/version stop parsing before model validation and I/O.
struct ReferencesRequested
{};

inline void add_references_option(CLI::App& app)
{
  app.get_help_ptr()->description("Show command-line help");
  app.add_flag_callback(
         "--references", [] { throw ReferencesRequested{}; }, "Show literature references and their applicability")
      ->group("Program")
      ->configurable(false)
      ->callback_priority(CLI::CallbackPriority::First);
}

inline uni20::presentation::report_builder references_report(uni20::presentation::program_info const& program)
{
  // Reuse Uni20's bibliography rendering without the option listing, examples,
  // or model notes. The normal help document has no reference provider.
  auto references = program;
  references.examples.clear();
  references.notes.clear();
  return uni20::presentation::help_report(references, program.name + " --references", {});
}

// Identity and literature stay application-owned; Uni20 owns parsing and rendering.
inline uni20::presentation::program_info program_info(std::string name, std::string description, citations::Tool tool)
{
  return {.name = std::move(name),
          .description = std::move(description),
          .version = build_info::version,
          .revision = build_info::revision,
          .copyright = "Copyright (C) 2026 Ian McCulloch",
          .project_url = "https://github.com/Uni20-dev/bethe",
          .license = "GPL-3.0-or-later; see COPYING",
          .references = [tool] {
            std::vector<uni20::presentation::program_reference> references;
            for (auto const& use : citations::for_tool(tool))
            {
              auto const& ref = *use.reference;
              std::string citation =
                  std::string(ref.authors) + '\n' + std::string(ref.title) + '\n' + std::string(ref.publication);
              if (ref.year) citation += " (" + std::to_string(ref.year) + ')';
              std::string links;
              for (auto const& link : ref.links)
              {
                if (!links.empty()) links += '\n';
                links += std::string(link.label) + ": " + std::string(link.url);
              }
              references.push_back({.key = "[" + std::string(ref.id) + "]",
                                    .citation = std::move(citation),
                                    .link = std::move(links),
                                    .applicability = "Used for: " + std::string(use.context)});
            }
            return references;
          }};
}

// One informational/error lifecycle for every executable. Registration callbacks
// only bind arguments; validation, solvers and output files belong in run().
template <typename Register, typename Run>
int program_main(int argc, char** argv, uni20::presentation::program_info const& program, Register&& register_options,
                 Run&& run)
{
  namespace options = uni20::cli;
  auto help_program = program;
  help_program.references = {};
  CLI::App app;
  try
  {
    options::configure(app, help_program);
    add_references_option(app);
    register_options(app);
    options::parse_result result;
    try
    {
      result = options::parse(app, argc, argv);
    }
    catch (ReferencesRequested const&)
    {
      uni20::display::emit(references_report(program), uni20::display::stream::out);
      return 0;
    }
    if (result.requested != options::action::run)
    {
      uni20::display::emit(options::result_report(app, help_program, result), result.destination);
      return result.exit_code;
    }
    return run(app);
  }
  catch (data::data_delivery_error const& error)
  {
    print_output_error(std::cerr, error);
  }
  catch (std::exception const& error)
  {
    uni20::display::emit(
        options::result_report(app, help_program, {.requested = options::action::error, .message = error.what()}),
        uni20::display::stream::err);
  }
  return 1;
}

inline CLI::Option* count_option(CLI::App& app, std::string names, std::optional<std::size_t>& value,
                                 std::string description)
{
  return app
      .add_option_function<std::string>(
          names, [&value](std::string const& token) { value = parse_size(token); }, std::move(description))
      ->type_name("COUNT");
}
// CLI11's generic optional binding treats an empty token as disengaged. Here
// presence is meaningful (an empty motif/quantum-number list is a valid state).
inline CLI::Option* text_option(CLI::App& app, std::string names, std::optional<std::string>& value,
                                std::string description)
{
  return app.add_option_function<std::string>(
      names, [&value](std::string const& token) { value = token; }, std::move(description));
}
inline CLI::Option* count_option(CLI::App& app, std::string names, std::size_t& value, std::string description)
{
  return uni20::cli::add_count_option(app, std::move(names), value, std::move(description));
}
inline void precision_option(CLI::App& app, std::string& precision)
{
  app.add_option("--precision", precision, "Real scalar type; fp128 requires MPLAPACK")
      ->check(CLI::IsMember({"fp64", "long-double", "fp128"}))
      ->capture_default_str();
}

inline CLI::Option* option(CLI::App& app, std::string names, bool& value, std::string description)
{
  return app.add_flag(std::move(names), value, std::move(description));
}
inline CLI::Option* option(CLI::App& app, std::string names, std::string& value, std::string description)
{
  return app.add_option(std::move(names), value, std::move(description));
}
inline CLI::Option* option(CLI::App& app, std::string names, std::optional<std::string>& value, std::string description)
{
  return text_option(app, std::move(names), value, std::move(description));
}
inline CLI::Option* option(CLI::App& app, std::string names, std::size_t& value, std::string description)
{
  return count_option(app, std::move(names), value, std::move(description));
}
inline CLI::Option* option(CLI::App& app, std::string names, std::optional<std::size_t>& value, std::string description)
{
  return count_option(app, std::move(names), value, std::move(description));
}
inline CLI::Option* option(CLI::App& app, std::string names, uni20::half_int& value, std::string description)
{
  return uni20::cli::add_half_int_option(app, std::move(names), value, std::move(description));
}
inline CLI::Option* option(CLI::App& app, std::string names, std::optional<uni20::half_int>& value,
                           std::string description)
{
  return app
      .add_option_function<std::string>(
          names, [&value](std::string const& token) { value = uni20::half_int::parse(token); }, std::move(description))
      ->type_name("HALF_INT");
}
inline CLI::Option* option(CLI::App& app, std::string names, std::optional<std::vector<uni20::half_int>>& value,
                           std::string description)
{
  return app
      .add_option_function<std::string>(
          names,
          [&value](std::string const& token) {
            value = token == "none" ? std::vector<uni20::half_int>{} : parse_quantum_numbers(token);
          },
          std::move(description))
      ->type_name("LIST");
}
inline CLI::Option* all_count_option(CLI::App& app, std::string names, std::optional<std::size_t>& value,
                                     std::string description)
{
  return app
      .add_option_function<std::string>(
          names,
          [&value](std::string const& token) {
            // The numerical scan checks its candidate budget before clamping this sentinel.
            value = token == "all" ? std::numeric_limits<std::size_t>::max() : parse_size(token);
          },
          std::move(description))
      ->type_name("COUNT|all");
}
} // namespace bethe::cli
