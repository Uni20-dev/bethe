// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include "bethe-build-info.hpp"
#include "report-common.hpp"
#include <ctime>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sys/stat.h>
#include <uni20/common/data_table_sinks.hpp>

namespace bethe::cli
{
namespace data = uni20::presentation;
struct DataFile
{
    std::string format, path;
    std::string table = {};
};
struct DataOutputOptions
{
    std::string format = "auto";
    std::vector<DataFile> files;
    std::string table;
    bool quiet = false, preamble = true, force = false, stream = false, retain = true;

    bool human() const { return format == "auto" || format == "pretty" || format == "plain"; }
    void validate() const
    {
      if (!human() && format != "csv" && format != "tsv" && format != "json")
        throw std::invalid_argument("unknown output format: " + std::string(format));
      if (!retain && !quiet && human() && !stream)
        throw std::invalid_argument("--no-retain needs --stream, a machine --format, or --quiet");
      for (auto const& file : files)
        if (file.format != "csv" && file.format != "tsv" && file.format != "json")
          throw std::invalid_argument("unknown export format: " + file.format);
    }
};
// Quote argv for POSIX shells, including empty arguments and literal apostrophes.
inline std::string quote_argument(std::string_view value)
{
  std::string result = "'";
  for (char c : value)
    result += c == '\'' ? "'\\''" : std::string(1, c);
  return result + "'";
}
inline std::string command_line(int argc, char** argv)
{
  std::string result;
  for (int i = 0; i < argc; ++i)
  {
    if (i) result += ' ';
    result += quote_argument(argv[i]);
  }
  return result;
}
inline data::table_metadata provenance(std::string program, int argc, char** argv)
{
  auto const now = std::time(nullptr);
  std::tm utc{};
  char date[32]{};
  if (now == std::time_t{-1} || !gmtime_r(&now, &utc) || !std::strftime(date, sizeof(date), "%Y-%m-%dT%H:%M:%SZ", &utc))
    throw std::runtime_error("could not determine output timestamp");
  return {{"Program", std::move(program)},
          {"Bethe version", build_info::version},
          {"Bethe revision", build_info::revision},
          {"Uni20 revision", build_info::uni20_revision},
          {"Date", date},
          {"Command", command_line(argc, argv)}};
}
// Keep every metadata field on one comment line; argv may contain real newlines.
inline std::string comment_text(std::string_view value)
{
  std::string result;
  for (unsigned char c : value)
  {
    if (c == '\\')
      result += "\\\\";
    else if (c == '\n')
      result += "\\n";
    else if (c == '\r')
      result += "\\r";
    else if (c == '\t')
      result += "\\t";
    else if (c < 32 || c == 127)
      result += fmt::format("\\x{:02x}", c);
    else
      result += static_cast<char>(c);
  }
  return result;
}
inline void write_comments(std::ostream& out, data::table_metadata const& values)
{
  for (auto const& [key, value] : values)
    out << "# " << comment_text(key) << ": " << comment_text(value) << '\n';
  if (!out) throw std::ios_base::failure("metadata output failed");
}

// Add Bethe's comment dialect without changing Uni20's rectangular CSV/TSV writers.
template <typename Sink> class CommentedSink {
  public:
    CommentedSink(std::ostream& out, Sink sink, bool preamble) : out_(out), sink_(std::move(sink)), preamble_(preamble)
    {}
    template <typename Schema>
    void begin(std::string const& title, Schema const& schema, data::table_metadata const& metadata,
               data::data_sink_start start)
    {
      if (preamble_)
      {
        write_comments(out_, metadata);
        write_comments(out_, {{"First row", std::to_string(start.first_row)}});
      }
      sink_.begin(title, schema, metadata, start);
    }
    template <typename Schema, typename Row> void row(Schema const& schema, Row const& row) { sink_.row(schema, row); }
    void finish(data::table_metadata const& summary)
    {
      if (preamble_) write_comments(out_, summary);
      sink_.finish(summary);
    }

  private:
    std::ostream& out_;
    Sink sink_;
    bool preamble_;
};

template <typename Sink> class NamedSink {
  public:
    NamedSink(std::string name, Sink sink) : name_(std::move(name)), sink_(std::move(sink)) {}
    template <typename... Args> void begin(Args const&... args)
    {
      this->call([&] { sink_.begin(args...); });
    }
    template <typename... Args> void row(Args const&... args)
    {
      this->call([&] { sink_.row(args...); });
    }
    void finish(data::table_metadata const& summary)
    {
      this->call([&] { sink_.finish(summary); });
    }

  private:
    template <typename Function> void call(Function&& function)
    {
      try
      {
        function();
      }
      catch (std::exception const& e)
      {
        throw std::runtime_error(name_ + ": " + e.what());
      }
    }
    std::string name_;
    Sink sink_;
};
inline void print_output_error(std::ostream& out, std::exception const& error)
{
  out << "Error: " << error.what() << '\n';
  if (auto delivery = dynamic_cast<data::data_delivery_error const*>(&error))
    for (auto const& failure : delivery->report().failures)
      try
      {
        std::rethrow_exception(failure.exception);
      }
      catch (std::exception const& cause)
      {
        out << "  " << cause.what() << '\n';
      }
      catch (...)
      {
        out << "  non-standard output exception\n";
      }
}

// Own files before attaching their borrowed streams. Call finish or abort while
// the table is still alive. No global output state is changed except the scoped
// plain-screen router used for explicitly requested human streaming.
class DataOutput {
  public:
    explicit DataOutput(DataOutputOptions options, std::vector<std::string> tables = {})
        : options_(std::move(options)), tables_(std::move(tables)), json_files_(options_.files.size())
    {
      options_.validate();
      if (!tables_.empty())
      {
        for (auto const& name : tables_)
        {
          if (name.empty() || name.find_first_not_of("abcdefghijklmnopqrstuvwxyz0123456789_") != std::string::npos)
            throw std::logic_error("invalid output table identifier: " + name);
          if (std::count(tables_.begin(), tables_.end(), name) != 1)
            throw std::logic_error("duplicate output table identifier: " + name);
        }
        if (options_.table.empty()) options_.table = tables_.front();
        auto check = [&](std::string const& name) {
          if (std::find(tables_.begin(), tables_.end(), name) == tables_.end())
            throw std::invalid_argument("output table is unavailable in this calculation: " + name);
        };
        check(options_.table);
        for (auto const& file : options_.files)
          if (!file.table.empty()) check(file.table);
      }
      this->check_paths(); // Reject collisions/existing files before opening any target.
      for (auto const& file : options_.files)
      {
        auto stream = std::make_unique<std::ofstream>();
        auto mode = std::ios::out | (options_.force ? std::ios::trunc : std::ios::noreplace);
        stream->open(file.path, mode);
        if (!*stream) throw std::runtime_error("cannot open output file: " + file.path);
        files_.push_back(std::move(stream));
      }
    }
    ~DataOutput()
    {
      // Also cover failures between tables (e.g. allocating an auxiliary schema).
      if (!tables_.empty() && !document_finished_) try
        {
          finish_document(true);
        }
        catch (std::exception const& error)
        {
          print_output_error(std::cerr, error);
        }
    }
    // Keep a compact, model-authored human overview. Export metadata still
    // carries the complete provenance on every independently readable table.
    void overview(report_builder report) { overview_ = std::move(report); }
    template <typename Table> void attach(Table& table, std::string name = {})
    {
      if (!tables_.empty())
      {
        if (document_finished_) throw std::logic_error("output document is already finished");
        if (!active_table_.empty()) throw std::logic_error("finish the current output table before attaching another");
        if (std::find(tables_.begin(), tables_.end(), name) == tables_.end() ||
            std::find(delivered_.begin(), delivered_.end(), name) != delivered_.end())
          throw std::logic_error("unknown or duplicate output table: " + name);
        active_table_ = name;
        delivered_.push_back(name);
      }
      // Files first, so they receive the accepted row even if the screen fails.
      for (std::size_t i = 0; i < files_.size(); ++i)
      {
        auto const& file = options_.files[i];
        if (tables_.empty() || file.format == "json" || name == (file.table.empty() ? options_.table : file.table))
          this->attach_stream(table, *files_[i], file.format, file.path, name, json_files_[i]);
      }
      if (options_.quiet) return;
      if (!options_.human())
      {
        if (tables_.empty() || options_.format == "json" || name == options_.table)
          this->attach_stream(table, std::cout, options_.format, "stdout", name, json_stdout_);
      }
      else if (options_.stream)
      {
        if (options_.format == "plain" && !plain_router_)
          plain_router_.emplace([](uni20::display::event const& event) {
            auto policy = data::plain_policy();
            policy.wrap_width = std::nullopt;
            std::cout << std::visit([&](auto const& content) { return data::render_plain(content, policy); },
                                    event.content);
            if (event.newline) std::cout << '\n';
            std::cout.flush();
            if (!std::cout) throw std::ios_base::failure("stdout write/flush failed");
          });
        table.attach(NamedSink("stdout", data::terminal_sink()));
      }
    }
    template <typename Table> void finish(Table& table, data::table_metadata summary)
    {
      std::exception_ptr failure;
      try
      {
        table.finish(std::move(summary));
      }
      catch (...)
      {
        failure = std::current_exception();
      }
      active_table_.clear();
      if (tables_.empty()) this->close_files(failure);
      if (!failure && !options_.quiet && options_.human() && !options_.stream)
      {
        report_builder report(table.title());
        if (overview_)
        {
          if (!overview_printed_) report = *overview_;
          overview_printed_ = true;
        }
        else
        {
          for (auto const& [key, value] : table.metadata())
            report.field(key, value);
          for (auto const& [key, value] : *table.summary())
            report.field(key, value);
        }
        report.table("") = data::to_report_table(table);
        print_report(report, options_.format);
        std::cout.flush();
        if (!std::cout) throw std::ios_base::failure("stdout write/flush failed");
      }
      if (failure) std::rethrow_exception(failure);
    }
    void finish_document(bool aborted = false)
    {
      if (tables_.empty()) throw std::logic_error("not a multi-table output document");
      if (document_finished_) return;
      if (!aborted && (!active_table_.empty() || delivered_.size() != tables_.size()))
        throw std::logic_error("not all output tables were finished");
      document_finished_ = true;
      std::exception_ptr failure;
      auto close_json = [&](std::ostream& stream, JsonDocument& state, std::string const& destination) {
        if (!state.started || state.closed || state.failed) return;
        state.closed = true;
        stream << "},\"status\":\"" << (aborted ? "aborted" : "complete") << "\"}\n";
        stream.flush();
        if (!stream && !failure)
          failure = std::make_exception_ptr(std::runtime_error(destination + ": JSON document finalization failed"));
      };
      for (std::size_t i = 0; i < files_.size(); ++i)
        close_json(*files_[i], json_files_[i], options_.files[i].path);
      close_json(std::cout, json_stdout_, "stdout");
      this->close_files(failure);
      if (failure) std::rethrow_exception(failure);
    }
    template <typename Table, typename Fill>
    void write_table(std::string name, Table& table, Fill&& fill, data::table_metadata summary = {})
    {
      try
      {
        attach(table, std::move(name));
        fill(table);
        finish(table, std::move(summary));
      }
      catch (...)
      {
        abort(table, {{"Status", "aborted"}});
        try
        {
          finish_document(true);
        }
        catch (std::exception const& error)
        {
          print_output_error(std::cerr, error);
        }
        throw;
      }
    }
    template <typename Table> void abort(Table& table, data::table_metadata summary) noexcept
    {
      // Best effort: close healthy JSON documents with an aborted status, never
      // retry the accepted row or replace a summary already finalized.
      try
      {
        if (!table.finished()) table.finish(std::move(summary));
      }
      catch (std::exception const& e)
      {
        print_output_error(std::cerr, e);
      }
      catch (...)
      {}
      std::exception_ptr failure;
      active_table_.clear();
      if (tables_.empty()) this->close_files(failure);
      if (failure)
      {
        try
        {
          std::rethrow_exception(failure);
        }
        catch (std::exception const& e)
        {
          print_output_error(std::cerr, e);
        }
        catch (...)
        {}
      }
    }

  private:
    struct JsonDocument
    {
        bool started = false, failed = false, closed = false;
    };
    // Each named table is an ordinary Uni20 JSON table, enclosed in a single
    // document. A failed destination is quarantined, not retried or repaired.
    class DocumentSink {
      public:
        DocumentSink(std::ostream& out, std::string name, JsonDocument& document)
            : out_(out), name_(std::move(name)), document_(document), sink_(data::json_sink(out))
        {}
        template <typename Schema>
        void begin(std::string const& title, Schema const& schema, data::table_metadata const& metadata,
                   data::data_sink_start start)
        {
          guarded([&] {
            out_ << (document_.started ? "," : "{\"tables\":{") << '"' << name_ << "\":";
            document_.started = true;
            sink_.begin(title, schema, metadata, start);
          });
        }
        template <typename Schema, typename Row> void row(Schema const& schema, Row const& row)
        {
          guarded([&] { sink_.row(schema, row); });
        }
        void finish(data::table_metadata const& summary)
        {
          guarded([&] { sink_.finish(summary); });
        }

      private:
        template <typename Function> void guarded(Function&& function)
        {
          if (document_.failed) throw std::runtime_error("JSON destination previously failed");
          try
          {
            function();
          }
          catch (...)
          {
            document_.failed = true;
            throw;
          }
        }
        std::ostream& out_;
        std::string name_;
        JsonDocument& document_;
        decltype(data::json_sink(std::declval<std::ostream&>())) sink_;
    };
    template <typename Table>
    void attach_stream(Table& table, std::ostream& out, std::string_view format, std::string const& name,
                       std::string const& table_name, JsonDocument& document)
    {
      if (format == "json" && !tables_.empty())
        table.attach(NamedSink(name, DocumentSink(out, table_name, document)));
      else if (format == "json")
        table.attach(NamedSink(name, data::json_sink(out)));
      else if (format == "csv")
        table.attach(NamedSink(name, CommentedSink(out, data::csv_sink(out), options_.preamble)));
      else
        table.attach(NamedSink(name, CommentedSink(out, data::tsv_sink(out), options_.preamble)));
    }
    void check_paths() const
    {
      namespace fs = std::filesystem;
      std::vector<fs::path> paths;
      struct stat stdout_stat
      {};
      bool const regular_stdout =
          !options_.quiet && fstat(fileno(stdout), &stdout_stat) == 0 && S_ISREG(stdout_stat.st_mode);
      for (auto const& file : options_.files)
      {
        if (file.path.empty() || file.path == "-")
          throw std::invalid_argument("export needs a file path, not '-' or an empty name");
        if (fs::is_symlink(file.path) && !fs::exists(file.path))
          throw std::invalid_argument("output destination is a dangling symlink: " + file.path);
        auto path = fs::weakly_canonical(file.path);
        for (auto const& previous : paths)
          if (path == previous || (fs::exists(path) && fs::exists(previous) && fs::equivalent(path, previous)))
            throw std::invalid_argument("output destinations refer to the same file: " + file.path);
        if (fs::exists(path))
        {
          if (!fs::is_regular_file(path))
            throw std::invalid_argument("output destination is not a regular file: " + file.path);
          if (!options_.force) throw std::invalid_argument("output file exists (use --force to replace): " + file.path);
          struct stat target
          {};
          if (regular_stdout && stat(file.path.c_str(), &target) == 0 && target.st_dev == stdout_stat.st_dev &&
              target.st_ino == stdout_stat.st_ino)
            throw std::invalid_argument("output file aliases stdout: " + file.path);
        }
        paths.push_back(std::move(path));
      }
    }
    void close_files(std::exception_ptr& failure) noexcept
    {
      for (std::size_t i = 0; i < files_.size(); ++i)
      {
        if (!files_[i]->is_open()) continue;
        files_[i]->close();
        if (!*files_[i] && !failure)
        {
          try
          {
            throw std::runtime_error("output close failed: " + options_.files[i].path);
          }
          catch (...)
          {
            failure = std::current_exception();
          }
        }
      }
    }
    DataOutputOptions options_;
    std::vector<std::string> tables_, delivered_;
    std::string active_table_;
    std::vector<JsonDocument> json_files_;
    JsonDocument json_stdout_;
    bool document_finished_ = false;
    std::optional<report_builder> overview_;
    bool overview_printed_ = false;
    std::vector<std::unique_ptr<std::ofstream>> files_;
    std::optional<uni20::display::scoped_sink> plain_router_;
};
} // namespace bethe::cli
