// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <fmt/format.h>
#include <uni20/common/run_context.hpp>

namespace bethe::cli
{
// Temporary consumer projection for the existing table metadata format. Values
// stay native in the run document; labels are not used to discover export keys.
inline uni20::presentation::table_metadata run_metadata(uni20::metadata_document const& document,
                                                        std::map<std::string, std::string> const& keys = {})
{
  auto result = document.strings(keys);
  // Existing Bethe metadata uses this spelling rather than an empty string.
  // The source document still distinguishes a missing value from empty text.
  for (auto const& group : document.groups())
    for (auto const& field : group.fields)
      if (field.value.missing()) result.at(keys.contains(field.id) ? keys.at(field.id) : field.id) = "unavailable";
  return result;
}

inline uni20::presentation::table_metadata run_summary(uni20::metadata_document const& summary)
{
  auto result = run_metadata(summary, {{"status", "Status"},
                                       {"rows", "Rows"},
                                       {"outcome", "Outcome"},
                                       {"compute_cpu_seconds", "CPU time"},
                                       {"run_cpu_seconds", "Run CPU seconds"},
                                       {"elapsed_seconds", "Elapsed seconds"}});
  // Preserve the established compute-only CPU field, including its units and
  // six fractional digits. New timing fields export native seconds directly.
  auto const& seconds = summary.find("compute_cpu_seconds")->value.get<std::optional<long double>>();
  result.at("CPU time") = seconds ? fmt::format("{:.6f} s", *seconds) : "unavailable";
  return result;
}
} // namespace bethe::cli
