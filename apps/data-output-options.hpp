// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include "data-output.hpp"
#include <uni20/cli/cli.hpp>

namespace bethe::cli
{
// Register configuration only. Files are opened by DataOutput after validation.
inline void add_data_output_options(CLI::App& app, DataOutputOptions& options)
{
  auto* output = app.add_option_group("Output");
  output->add_option("--format", options.format, "Stdout format; file exports are independent")
      ->check(CLI::IsMember({"auto", "pretty", "plain", "csv", "tsv", "json"}))
      ->capture_default_str();
  for (std::string const format : {"csv", "tsv", "json"})
    output
        ->add_option_function<std::vector<std::string>>(
            "--" + format,
            [&options, format](std::vector<std::string> const& paths) {
              for (auto const& path : paths)
                options.files.push_back({format, path});
            },
            "Additional " + format + " export; repeat for multiple files")
        ->type_name("FILE")
        ->expected(1)
        ->take_all();
  output->add_flag("--quiet", options.quiet, "Suppress stdout, not files or warnings");
  output->add_flag("!--no-preamble", options.preamble, "Omit CSV/TSV metadata and summary comments");
  output->add_flag("--force", options.force, "Permit replacing existing regular files");
  output->add_flag("--stream", options.stream, "Display human-readable rows as computed");
  output->add_flag("!--no-retain", options.retain, "Discard rows after delivery; needs live output or --quiet");
}
} // namespace bethe::cli
