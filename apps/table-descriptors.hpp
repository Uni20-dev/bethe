// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <algorithm>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace bethe::cli
{
// Frontend-owned contract, independent of CLI11 and scalar/column types.
// A screen flag requests presentation, never changes a table's scientific domain.
struct TableDescriptor
{
    std::string name, description, availability = "always", screen_option = {};
    bool primary = false;
    std::string required_option = {};
};
using TableCatalogue = std::vector<TableDescriptor>;

inline void validate_table_catalogue(TableCatalogue const& catalogue)
{
  bool primary = false;
  std::vector<std::string> names;
  for (auto const& table : catalogue)
  {
    if (table.name.empty() ||
        table.name.find_first_not_of("abcdefghijklmnopqrstuvwxyz0123456789_") != std::string::npos)
      throw std::logic_error("invalid output table identifier: " + table.name);
    if (std::ranges::find(names, table.name) != names.end())
      throw std::logic_error("duplicate output table descriptor: " + table.name);
    names.push_back(table.name);
    primary |= table.primary;
  }
  if (!catalogue.empty() && !primary) throw std::logic_error("table catalogue has no declared primary table");
}

inline TableDescriptor const& table_descriptor(TableCatalogue const& catalogue, std::string_view name)
{
  auto found = std::ranges::find(catalogue, name, &TableDescriptor::name);
  if (found != catalogue.end()) return *found;
  std::string valid;
  for (auto const& table : catalogue)
    valid += (valid.empty() ? "" : ", ") + table.name;
  throw std::invalid_argument("unknown output table '" + std::string(name) + "'; valid tables: " + valid);
}

} // namespace bethe::cli
