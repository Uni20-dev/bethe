// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include "result-output.hpp"
#include <tuple>

namespace bethe::cli
{
template <typename State> std::string spin_status(State const& state)
{
  if (state.converged) return "converged";
  if constexpr (requires { state.status; })
  {
    using Status = decltype(state.status);
    if constexpr (requires { Status::stalled; })
      if (state.status == Status::stalled) return "line search stalled; unconverged estimate";
  }
  return "iteration limit reached; unconverged estimate";
}

template <uni20::Real Real, typename State>
void spin_rows(ResultOutput& output, std::string name, std::string title, std::vector<State const*> const& states,
               std::vector<std::optional<Real>> const& gaps, std::size_t offset = 0)
{
  auto const extra_columns = [] {
    auto momentum = [] {
      if constexpr (requires(State s) { s.momentum; })
        return std::tuple{column<std::size_t>("momentum_index", "Momentum index"), column<Real>("p", "P")};
      else
        return std::tuple{};
    }();
    auto root_delta = [] {
      if constexpr (requires(State s) { s.root_delta; })
        return std::tuple{column<Real>("root_delta", "Root Delta")};
      else
        return std::tuple{};
    }();
    return std::tuple_cat(momentum, root_delta);
  }();
  std::apply(
      [&](auto... extra) {
        output.table(
            std::move(name), std::move(title),
            [&](auto& table) {
              for (std::size_t i = 0; i < states.size(); ++i)
              {
                auto const& s = *states[i];
                auto momentum = [&] {
                  if constexpr (requires { s.momentum; })
                    return std::tuple{s.momentum_index, s.momentum};
                  else
                    return std::tuple{};
                }();
                auto root_delta = [&] {
                  if constexpr (requires { s.root_delta; })
                    return std::tuple{s.root_delta};
                  else
                    return std::tuple{};
                }();
                std::apply(
                    [&](auto... values) {
                      table.append(offset + i, s.sz, s.spin_reversed, s.energy, gaps.empty() ? std::nullopt : gaps[i],
                                   s.residual_norm, s.iterations, s.converged, spin_status(s), values...);
                    },
                    std::tuple_cat(momentum, root_delta));
              }
            },
            column<std::size_t>("state_id"), column<uni20::half_int>("sz", "Sz"), column<bool>("spin_reversed"),
            column<Real>("energy", "Energy"), column<std::optional<Real>>("gap", "E-E0"),
            column<Real>("residual", "Residual"), column<std::size_t>("iterations", "Iterations"),
            column<bool>("converged"), column<std::string>("status", "Status"), extra...);
      },
      extra_columns);
}

template <uni20::Real Real, typename State>
void spin_roots(ResultOutput& output, std::vector<State const*> const& states, bool roots, bool labels)
{
  if (labels)
    output.table(
        "quantum_numbers", "Quantum numbers",
        [&](auto& table) {
          for (std::size_t i = 0; i < states.size(); ++i)
            for (std::size_t j = 0; j < states[i]->quantum_numbers.size(); ++j)
              table.append(i, j, states[i]->quantum_numbers[j]);
        },
        column<std::size_t>("state_id"), column<std::size_t>("index", "Index"),
        column<uni20::half_int>("quantum_number", "I"));
  if (!roots) return;
  std::string title = "Rapidities";
  if (states.size() == 1 && states.front()->rapidities.empty())
  {
    bool boundary = false;
    if constexpr (requires(State s) { s.boundary_root; }) boundary = states.front()->boundary_root.has_value();
    title += boundary ? " (no bulk roots)" : " (none; fully polarized)";
  }
  auto const extra_columns = [] {
    if constexpr (requires(State s) { s.log_rapidities; })
      return std::tuple{column<std::optional<Real>>("lambda", "Lambda")};
    else
      return std::tuple{};
  }();
  std::apply(
      [&](auto... extra) {
        output.table(
            "roots", title,
            [&](auto& table) {
              for (std::size_t i = 0; i < states.size(); ++i)
              {
                auto const& s = *states[i];
                for (std::size_t j = 0; j < s.rapidities.size(); ++j)
                {
                  auto values = [&] {
                    if constexpr (requires { s.log_rapidities; })
                      return std::tuple{s.log_rapidities.empty() ? std::nullopt
                                                                 : std::optional<Real>{s.log_rapidities[j]}};
                    else
                      return std::tuple{};
                  }();
                  std::apply([&](auto... tail) { table.append(i, j, s.quantum_numbers[j], s.rapidities[j], tail...); },
                             values);
                }
              }
            },
            column<std::size_t>("state_id"), column<std::size_t>("index", "Index"),
            column<uni20::half_int>("quantum_number", "I"), column<Real>("rapidity", "Rapidity"), extra...);
      },
      extra_columns);
  if constexpr (requires(State s) { s.boundary_root; })
    output.table(
        "boundary_roots", "Boundary roots (not bulk rapidities)",
        [&](auto& table) {
          for (std::size_t i = 0; i < states.size(); ++i)
            if (states[i]->boundary_root)
            {
              auto const& root = *states[i]->boundary_root;
              std::string kind = root.inverse_square > Real{0}   ? "real"
                                 : root.inverse_square < Real{0} ? "imaginary"
                                                                 : "infinity";
              table.append(i, root.quantum_number, root.inverse_square, root.log_distance, kind);
            }
        },
        column<std::size_t>("state_id"), column<uni20::half_int>("quantum_number", "I"),
        column<Real>("inverse_square", "Inverse square y"), column<Real>("log_distance", "Log distance w"),
        column<std::string>("kind", "Kind of z"));
}

template <typename State>
std::vector<std::string> spin_tables(bool roots, bool labels, bool reference = false, bool failed = false)
{
  std::vector<std::string> names{"states"};
  if (reference) names.push_back("reference");
  if (failed) names.push_back("failed");
  if (labels) names.push_back("quantum_numbers");
  if (roots)
  {
    names.push_back("roots");
    if constexpr (requires(State s) { s.boundary_root; }) names.push_back("boundary_roots");
  }
  return names;
}

template <uni20::Real Real, typename State>
bool spin_state(report_builder report, std::size_t sites, State const& state, bool roots,
                DataOutputOptions const& options, std::string program, int argc, char** argv)
{
  report.status(state.converged ? semantic_glyph::success : semantic_glyph::warning, spin_status(state))
      .field("Sz", uni20::to_string(state.sz))
      .field("Reference vacuum", state.spin_reversed ? "all down (spin reversed)" : "all up")
      .field("Spin-reversed reference", state.spin_reversed ? 1 : 0)
      .field("Status", spin_status(state))
      .field("Iterations", state.iterations)
      .field("Residual norm", uni20::format_real(state.residual_norm))
      .field("Total energy", uni20::format_real(state.energy))
      .field("Energy per site", uni20::format_real(state.energy / Real(sites)));
  if constexpr (requires { state.momentum; })
    report.field("Momentum index", state.momentum_index).field("Momentum P", uni20::format_real(state.momentum));
  if constexpr (requires { state.boundary_root; })
  {
    auto const kind = !state.boundary_root                            ? "none"
                      : state.boundary_root->inverse_square > Real{0} ? "real"
                      : state.boundary_root->inverse_square < Real{0} ? "imaginary"
                                                                      : "infinity";
    report.field("Root Delta", uni20::format_real(state.root_delta)).field("Boundary root", kind);
  }
  ResultOutput output(report, options, std::move(program), argc, argv, spin_tables<State>(roots, false));
  spin_rows<Real>(output, "states", "State", std::vector{&state}, {});
  spin_roots<Real>(output, std::vector{&state}, roots, false);
  output.finish();
  return state.converged;
}

template <uni20::Real Real, typename State>
bool spin_sectors(report_builder report, std::vector<State> const& states, bool roots, DataOutputOptions const& options,
                  std::string program, int argc, char** argv)
{
  std::vector<State const*> rows;
  std::size_t converged = 0;
  for (auto const& s : states)
  {
    rows.push_back(&s);
    converged += s.converged;
  }
  add_scan_status(report, converged, states.size());
  report.field("Status", converged == states.size() ? "converged" : "incomplete; unconverged estimates");
  if constexpr (requires(State s) { s.momentum; })
    report.field("Momentum convention", "P = 2*pi*momentum_index/N (mod 2*pi)");
  ResultOutput output(report, options, std::move(program), argc, argv, spin_tables<State>(roots, false));
  spin_rows<Real>(output, "states", "Sector energies and convergence", rows, {});
  spin_roots<Real>(output, rows, roots, false);
  output.finish();
  return converged == states.size();
}
} // namespace bethe::cli
