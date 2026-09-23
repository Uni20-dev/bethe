// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include "bethe-build-info.hpp"
#include <bethe/citations.hpp>
#include <uni20/cli/cli.hpp>

namespace bethe::cli
{
// Identity and literature stay application-owned; Uni20 owns parsing and rendering.
inline uni20::presentation::program_info program_info(std::string name, std::string description, citations::Tool tool)
{
  return {.name = std::move(name),
          .description = std::move(description),
          .version = build_info::version,
          .revision = build_info::revision,
          .copyright = "Copyright (C) 2026 Ian McCulloch",
          .project_url = "https://github.com/Uni20-dev/bethe",
          .license = "GPL-3.0-or-later; see COPYING",
          .references = [tool] {
            std::vector<uni20::presentation::program_reference> references;
            for (auto const& use : citations::for_tool(tool))
            {
              auto const& ref = *use.reference;
              std::string citation =
                  std::string(ref.authors) + '\n' + std::string(ref.title) + '\n' + std::string(ref.publication);
              if (ref.year) citation += " (" + std::to_string(ref.year) + ')';
              std::string links;
              for (auto const& link : ref.links)
              {
                if (!links.empty()) links += '\n';
                links += std::string(link.label) + ": " + std::string(link.url);
              }
              references.push_back({.key = "[" + std::string(ref.id) + "]",
                                    .citation = std::move(citation),
                                    .link = std::move(links),
                                    .applicability = "Used for: " + std::string(use.context)});
            }
            return references;
          }};
}
} // namespace bethe::cli
