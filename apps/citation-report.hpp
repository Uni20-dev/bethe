// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/citations.hpp>
#include <ostream>

namespace bethe::cli
{
inline void print_citations(std::ostream& out, citations::Tool tool)
{
  out << "\nReferences:\n"
      << "Cite the references relevant to the modes used; see CITATIONS.md for conventions and provenance.\n";
  for (auto const& use : citations::for_tool(tool))
  {
    auto const& ref = *use.reference;
    out << "\n  [" << ref.id << "] " << ref.authors << '\n' << "    " << ref.title << '\n' << "    " << ref.publication;
    if (ref.year) out << " (" << ref.year << ')';
    out << '\n';
    for (auto const& link : ref.links)
      out << "    " << link.label << ": " << link.url << '\n';
    out << "    Used for: " << use.context << '\n';
  }
}
} // namespace bethe::cli
