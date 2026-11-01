// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include "spinon-precision-stages.hpp"
#include <bethe/xxx_six_spinon_precision_recovery.hpp>

namespace bethe::cli::spinon
{
inline char const* precision_stage_name(PrecisionStage stage)
{
  switch (stage)
  {
    case PrecisionStage::fp80:
      return "fp80";
    case PrecisionStage::fp128:
      return "fp128";
    case PrecisionStage::reference:
      return "reference";
  }
  throw std::logic_error("unknown precision stage");
}
inline char const* observation_failure_name(heisenberg::SpinonObservationFailure failure)
{
  using F = heisenberg::SpinonObservationFailure;
  switch (failure)
  {
    case F::none:
      return "none";
    case F::geometry:
      return "geometry";
    case F::density:
      return "density";
    case F::exponent_range:
      return "exponent_range";
    case F::cancelled:
      return "cancelled";
  }
  throw std::logic_error("unknown observation failure");
}
inline char const* recovery_status_name(heisenberg::SpinonRecoveryStatus status)
{
  using S = heisenberg::SpinonRecoveryStatus;
  switch (status)
  {
    case S::accepted:
      return "accepted";
    case S::exhausted:
      return "exhausted";
    case S::cancelled:
      return "cancelled";
    case S::time_limit:
      return "time_limit";
    case S::invalid_input:
      return "invalid_input";
    case S::unsupported_limit:
      return "unsupported_limit";
    case S::inconsistent_controls:
      return "inconsistent_controls";
  }
  throw std::logic_error("unknown recovery status");
}
inline char const* density_failure_name(heisenberg::SpinonDensityFailure failure)
{
  using F = heisenberg::SpinonDensityFailure;
  switch (failure)
  {
    case F::none:
      return "none";
    case F::contour_tail:
      return "contour_tail";
    case F::arithmetic_precision:
      return "arithmetic_precision";
    case F::tail_domain:
      return "tail_domain";
    case F::exponent_range:
      return "exponent_range";
    case F::special_function_budget:
      return "special_function_budget";
    case F::rounded_collision:
      return "rounded_collision";
  }
  throw std::logic_error("unknown density failure");
}

// Worst case: each attempted observation has its own failure reason. Check the
// bound before opening any destination, rather than losing rows after sampling.
inline std::size_t recovery_detail_bound(std::size_t replicates, std::size_t samples, std::size_t stages,
                                         std::size_t remaining)
{
  if (!replicates || !stages) return 0;
  auto const per_stage = remaining / replicates / stages;
  if (samples >= per_stage)
    throw std::invalid_argument(
        "requested recovery table can exceed --detail-row-limit; split the grid or raise the limit");
  return replicates * stages * (samples + 1);
}

template <class Table, class R>
void append_recovery_records(Table& table, R const& q, R const& omega, std::uint64_t replicate,
                             std::vector<PrecisionStage> const& policy,
                             std::vector<heisenberg::SixSpinonPrecisionStageWork> const& work)
{
  if (policy.size() != work.size()) throw std::logic_error("recovery stage accounting does not match its policy");
  using Count = std::optional<std::size_t>;
  // Empty enum tags mean absent; valid enum spellings are always nonempty.
  // Uni20 deliberately excludes optional strings until CSV null encoding is specified.
  using Text = std::string;
  for (std::size_t i = 0; i < policy.size(); ++i)
  {
    auto const& stage = work[i];
    std::string const name = precision_stage_name(policy[i]);
    table.append(q, omega, replicate, name, std::string("stage"), stage.attempts, Count(stage.accepted),
                 Count(stage.retried), stage.evaluated_roots, stage.density_attempts, Text{}, Text{}, Text{},
                 std::string{});
    for (auto const& [key, counts] : stage.failures)
    {
      auto const& [failure, status, density_failure, reason] = key;
      table.append(q, omega, replicate, name, std::string("failure"), counts.attempts, Count{}, Count{},
                   counts.evaluated_roots, counts.density_attempts, Text(observation_failure_name(failure)),
                   status ? Text(recovery_status_name(*status)) : Text{},
                   density_failure ? Text(density_failure_name(*density_failure)) : Text{}, reason);
    }
  }
}
} // namespace bethe::cli::spinon
