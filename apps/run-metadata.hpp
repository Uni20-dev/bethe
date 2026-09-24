// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <fmt/format.h>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <uni20/common/run_context.hpp>
#include <utility>
#include <vector>

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
  std::map<std::string, std::string> keys{{"status", "Status"},
                                          {"outcome", "Outcome"},
                                          {"compute_cpu_seconds", "CPU time"},
                                          {"run_cpu_seconds", "Run CPU seconds"},
                                          {"elapsed_seconds", "Elapsed seconds"}};
  if (summary.find("rows")) keys.emplace("rows", "Rows");
  auto result = run_metadata(summary, keys);
  // Preserve the established compute-only CPU field, including its units and
  // six fractional digits. New timing fields export native seconds directly.
  auto const& seconds = summary.find("compute_cpu_seconds")->value.get<std::optional<long double>>();
  result.at("CPU time") = seconds ? fmt::format("{:.6f} s", *seconds) : "unavailable";
  return result;
}

inline std::map<std::string, std::string> provenance_keys(uni20::run_context const& context)
{
  std::map<std::string, std::string> keys{{"program", "Program"},         {"version", "Bethe version"},
                                          {"revision", "Bethe revision"}, {"uni20_revision", "Uni20 revision"},
                                          {"started_utc", "Date"},        {"compiler", "Compiler"},
                                          {"build_type", "Build type"},   {"platform", "Platform"}};
  if (!context.invocation().empty()) keys.emplace("invocation", "Command");
  return keys;
}

// Bethe's presentation policy for a completed batch calculation. Uni20 owns the
// native values and clocks; explicit stable IDs and legacy export keys live at
// the model call sites. No formatted report is ever parsed back into metadata.
class RunReport {
  public:
    RunReport(uni20::run_context& context, std::string title)
        : context_(context), title_(std::move(title)), keys_(provenance_keys(context))
    {
      context_.metadata().group("calculation", "Calculation");
    }
    RunReport(RunReport const&) = delete;
    RunReport(RunReport&&) = default;
    RunReport& field(std::string id, std::string key, uni20::metadata_value value,
                     uni20::presentation::data_column_display display = {.missing = "unavailable"})
    {
      context_.metadata().add("calculation", id, std::move(value), {.label = key, .display = std::move(display)});
      keys_.emplace(std::move(id), std::move(key));
      return *this;
    }
    RunReport& status(uni20::presentation::semantic_glyph glyph, std::string label)
    {
      statuses_.emplace_back(glyph, std::move(label));
      return *this;
    }
    RunReport& result(bool complete, std::string status)
    {
      if (outcome_) throw std::logic_error("scientific outcome already specified");
      outcome_ = complete ? uni20::run_outcome::success : uni20::run_outcome::partial;
      summary_.group("results", "Results");
      summary_.add("results", "status", std::move(status), {.label = "Status"});
      return *this;
    }
    uni20::run_context& context() const { return context_; }
    uni20::presentation::table_metadata metadata() const
    {
      auto snapshot = context_.snapshot();
      auto values = run_metadata(snapshot, keys_);
      // Preserve model-specific missing-value explanations (overflow, seeds, etc.).
      for (auto const& group : snapshot.groups())
        if (group.id == "calculation")
          for (auto const& field : group.fields)
            if (field.value.missing()) values.at(keys_.at(field.id)) = field.options.display.missing;
      return values;
    }
    uni20::presentation::table_metadata finish()
    {
      if (!outcome_) throw std::logic_error("model must specify its scientific outcome");
      return run_summary(context_.finish(*outcome_, summary_));
    }
    uni20::presentation::report_builder overview(uni20::presentation::table_metadata const& summary) const
    {
      uni20::presentation::report_builder report(title_);
      for (auto const& [glyph, label] : statuses_)
        report.status(glyph, label);
      auto snapshot = context_.snapshot();
      for (auto const& group : snapshot.groups())
        if (group.id == "calculation")
          for (auto const& field : group.fields)
            report.field(field.options.label, field.value.text(field.options.display));
      for (auto const& [key, value] : summary)
        report.field(key, value);
      return report;
    }

  private:
    uni20::run_context& context_;
    std::string title_;
    std::map<std::string, std::string> keys_;
    std::vector<std::pair<uni20::presentation::semantic_glyph, std::string>> statuses_;
    uni20::metadata_document summary_;
    std::optional<uni20::run_outcome> outcome_;
};
} // namespace bethe::cli
