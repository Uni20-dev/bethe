// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include <bethe/citations.hpp>

#include <iostream>
#include <set>
#include <stdexcept>
#include <string_view>

namespace refs = bethe::citations;
static_assert(refs::find("lieb-wu-2003")->year == 2003);
static_assert(refs::find("missing") == nullptr);
static_assert(refs::for_tool(refs::Tool::hubbard_pbc).front().reference == refs::find("lieb-wu-2003"));

int main()
{
  try
  {
    std::set<std::string_view> ids;
    for (auto const& ref : refs::references)
    {
      if (!ids.insert(ref.id).second || refs::find(ref.id) != &ref || ref.links.empty())
        throw std::runtime_error("invalid citation registry");
    }
    for (auto tool : {refs::Tool::xxx_pbc, refs::Tool::xxx_obc, refs::Tool::xxz_pbc, refs::Tool::xxz_obc,
                      refs::Tool::hubbard_pbc, refs::Tool::hubbard_obc})
    {
      ids.clear();
      auto const uses = refs::for_tool(tool);
      if (uses.empty()) throw std::runtime_error("empty tool references");
      for (auto const& use : uses)
        if (!use.reference || !ids.insert(use.reference->id).second || refs::find(use.reference->id) != use.reference ||
            use.context.empty())
          throw std::runtime_error("invalid tool reference selection");
    }
    try
    {
      (void)refs::for_tool(static_cast<refs::Tool>(-1));
      throw std::runtime_error("invalid tool silently accepted");
    }
    catch (std::invalid_argument const&)
    {}
    return 0;
  }
  catch (std::exception const& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
