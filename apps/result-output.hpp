// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include "data-output-options.hpp"
#include "run-metadata.hpp"

namespace bethe::cli
{
template <data::DataTableValue T> auto column(std::string id, std::string label = {})
{
  data::data_column<T> result(std::move(id));
  if (!label.empty()) result.label(std::move(label));
  if constexpr (requires { result.round_trip(); }) result.round_trip();
  return result;
}

// Models supply typed metadata and an explicit scientific outcome. Freeze one
// numerical summary before emitting the batch, and share it across all tables.
// Output completion remains an independent document-level property.
class ResultOutput {
  public:
    ResultOutput(RunReport& report, DataOutputOptions const& options, std::vector<std::string> names,
                 bool overview = true)
        : options_{.retain = options.retain ? data::retention::all : data::retention::none,
                   .metadata = report.metadata()},
          summary_(report.finish()), output_(options, std::move(names))
    {
      if (overview) output_.overview(report.overview(summary_));
    }
    template <typename Fill, data::DataTableValue... Ts>
    void table(std::string name, std::string title, Fill&& fill, data::data_column<Ts>... columns)
    {
      if (!output_.accepts(name)) return;
      auto table = data::make_data_table(std::move(title), options_, std::move(columns)...);
      output_.write_table(std::move(name), table, std::forward<Fill>(fill), summary_);
    }
    void finish() { output_.finish_document(); }

  private:
    data::data_table_options options_;
    data::table_metadata summary_;
    DataOutput output_;
};

// Streaming calculations finish their numerical report only after the last
// row. Keep table/document completion and exception cleanup identical across
// frontends. The fill callback returns the model's scientific success flag.
template <typename Table, typename Fill>
bool stream_result_table(RunReport& report, uni20::run_context& context, DataOutputOptions const& options,
                         std::string name, Table& table, Fill&& fill,
                         std::string failure_status = "precision_limit; unavailable quantities omitted")
{
  DataOutput output(options, {name});
  try
  {
    output.attach(table, name);
    bool const complete = std::forward<Fill>(fill)(table);
    report.result(complete, complete ? "converged" : std::move(failure_status));
    auto const summary = report.finish();
    output.overview(report.overview(summary));
    output.finish(table, summary);
    output.finish_document();
    return complete;
  }
  catch (...)
  {
    data::table_metadata aborted{{"Status", "aborted"}};
    if (!context.finished()) try
      {
        aborted = run_summary(context.finish(uni20::run_outcome::failed));
      }
      catch (...)
      {}
    output.abort(table, std::move(aborted));
    throw;
  }
}
} // namespace bethe::cli
