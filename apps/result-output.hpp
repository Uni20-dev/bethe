// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include "data-output-options.hpp"

namespace bethe::cli
{
template <data::DataTableValue T> auto column(std::string id, std::string label = {})
{
  data::data_column<T> result(std::move(id));
  if (!label.empty()) result.label(std::move(label));
  if constexpr (requires { result.round_trip(); }) result.round_trip();
  return result;
}

// Application metadata stays with the model. This adapter only adds provenance
// and gives a collection of differently typed tables one output lifetime.
class ResultOutput {
  public:
    ResultOutput(report_builder const& report, DataOutputOptions const& options, std::string program, int argc,
                 char** argv, std::vector<std::string> names)
        : options_{.retain = options.retain ? data::retention::all : data::retention::none,
                   .metadata = provenance(std::move(program), argc, argv)},
          output_(options, std::move(names))
    {
      output_.overview(report);
      for (auto const& [key, value] : report.fields())
        if (key == "CPU time" || key == "Status")
          summary_[key] = value;
        else
          options_.metadata[key] = value;
    }
    template <typename Fill, data::DataTableValue... Ts>
    void table(std::string name, std::string title, Fill&& fill, data::data_column<Ts>... columns)
    {
      auto table = data::make_data_table(std::move(title), options_, std::move(columns)...);
      output_.write_table(std::move(name), table, std::forward<Fill>(fill), summary_);
    }
    void finish() { output_.finish_document(); }

  private:
    data::data_table_options options_;
    data::table_metadata summary_;
    DataOutput output_;
};
} // namespace bethe::cli
