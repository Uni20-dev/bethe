// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#pragma once

#include "cli-common.hpp"
#include <uni20/common/presentation.hpp>

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

// Presentation belongs to the front end, not the numerical library. Typed data
// tables supply their native-precision cell text; scalar overviews are projected
// from Uni20's typed run metadata without narrowing through double.
namespace bethe::cli
{
namespace presentation = uni20::presentation;
using presentation::report_builder;
using presentation::semantic_glyph;
using presentation::table_alignment;

// Our tables have one labeled column per cell, with no spans or separator rows.
// Prefer vertical records when a table would exceed the terminal width. Never
// ask the table renderer to hard-wrap a high-precision numeric token.
inline void print_report(report_builder const& report)
{
  auto policy = presentation::terminal_policy(stdout);
  policy.wrap_width = std::nullopt;
  auto const columns = static_cast<std::size_t>(terminal::columns(stdout));
  report_builder header(report.title());
  for (auto const& [glyph, label] : report.statuses())
    header.status(glyph, label);
  for (auto const& [key, value] : report.fields())
    header.field(key, value);
  std::cout << presentation::render_terminal(header, policy);

  for (auto const& table : report.tables())
  {
    report_builder block;
    block.table("") = table;
    auto const plain = presentation::render_plain(block, policy);
    bool fits = true;
    for (std::string_view remaining = plain; !remaining.empty();)
    {
      auto const newline = remaining.find('\n');
      fits = fits && presentation::display_width(remaining.substr(0, newline), policy) <= columns;
      if (newline == std::string_view::npos) break;
      remaining.remove_prefix(newline + 1);
    }
    std::cout << '\n';
    if (fits)
    {
      std::cout << presentation::render_terminal(block, policy);
      continue;
    }

    presentation::styled_text records;
    records.append(presentation::style("Cyan")(table.title())).append("\n");
    std::size_t label_width = 0;
    for (auto const& column : table.columns())
      label_width = std::max(label_width, presentation::display_width(column.heading, policy));
    bool first = true;
    for (auto const& entry : table.entries())
    {
      auto const& row = std::get<std::vector<presentation::table_cell>>(entry);
      if (!first) records.append("\n");
      first = false;
      for (std::size_t i = 0; i < row.size(); ++i)
        records.append("  ")
            .append(presentation::style("LightGray")(
                presentation::pad_right(table.columns()[i].heading, label_width + 2, policy)))
            .append(row[i].content)
            .append("\n");
    }
    std::cout << presentation::render_terminal(records, policy);
  }
}

// Explicit-format front ends share stable key/value metadata and unwrapped
// numeric tables. Plain output is independent of terminal width and color.
// All real-valued cells must already have been formatted in their own type.
inline void print_report(report_builder const& report, std::string_view format)
{
  if (format == "pretty" || (format == "auto" && terminal::is_a_terminal(stdout)))
    print_report(report);
  else
  {
    std::cout << "# " << report.title() << '\n';
    for (auto const& [key, value] : report.fields())
      std::cout << key << ": " << value << '\n';
    auto policy = presentation::plain_policy();
    policy.wrap_width = std::nullopt;
    for (auto const& table : report.tables())
    {
      report_builder block;
      block.table("") = table;
      std::cout << '\n' << presentation::render_plain(block, policy);
    }
  }
}

template <typename Report> void add_scan_status(Report& report, std::size_t converged, std::size_t total)
{
  bool const all = converged == total;
  report.status(all ? semantic_glyph::success : semantic_glyph::warning,
                std::to_string(converged) + "/" + std::to_string(total) + " states converged" +
                    (all ? "" : "; remaining energies are unconverged estimates"));
}

} // namespace bethe::cli
