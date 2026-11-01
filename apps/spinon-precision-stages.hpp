// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include "precision.hpp"
#include <optional>
#include <vector>

namespace bethe::cli::spinon
{
enum class PrecisionStage
{
  fp80,
  fp128,
  reference
};

inline std::vector<PrecisionStage> precision_stages(std::optional<std::string> const& text)
{
  if (!text) return {};
  std::vector<PrecisionStage> result;
  std::string_view remaining = *text;
  while (true)
  {
    auto const comma = remaining.find(',');
    auto const token = remaining.substr(0, comma);
    PrecisionStage stage;
    if (token == "fp80")
      stage = PrecisionStage::fp80;
    else if (token == "fp128")
      stage = PrecisionStage::fp128;
    else if (token == "reference")
      stage = PrecisionStage::reference;
    else
      throw std::invalid_argument("--six-precision-stages accepts fp80,fp128,reference in that order");
    if (!result.empty() && stage <= result.back())
      throw std::invalid_argument("--six-precision-stages requires distinct, increasing stages");
    result.push_back(stage);
    if (comma == std::string_view::npos) break;
    remaining.remove_prefix(comma + 1);
  }
  if (result.back() != PrecisionStage::reference)
    throw std::invalid_argument("--six-precision-stages must end with reference");
  return result;
}

inline void validate_precision_stages(std::vector<PrecisionStage> const& stages, PrecisionRequest reference,
                                      bool has_six, std::optional<std::size_t> guard_bits)
{
  if (stages.empty()) return;
  if (!has_six) throw std::invalid_argument("--six-precision-stages requires sector 6");
  if (!reference.arbitrary())
    throw std::invalid_argument("--six-precision-stages requires an MPFR reference: --precision bits-N or digits-N");
  auto const output = arithmetic_info(reference);
  auto const guard = guard_bits.value_or(std::size_t(output.bits));
  if (guard <= 24 || guard > std::size_t(output.bits))
    throw std::invalid_argument("--density-guard-bits must be in 25..reference geometry bits");
  for (auto stage : stages)
  {
    if (stage == PrecisionStage::reference) continue;
    PrecisionRequest const request{stage == PrecisionStage::fp80 ? PrecisionKind::fp80 : PrecisionKind::fp128};
    // Never silently omit a requested native format on a different build/host.
    auto const native = arithmetic_info(request);
    if (guard < std::size_t(native.bits))
      throw std::invalid_argument("reference density guard must retain every preceding native precision stage");
  }
}
} // namespace bethe::cli::spinon
