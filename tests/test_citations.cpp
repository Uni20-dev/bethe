// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include <bethe/citations.hpp>
#include <gtest/gtest.h>

#include <set>
#include <stdexcept>
#include <string>
#include <string_view>

namespace refs = bethe::citations;
static_assert(refs::find("lieb-wu-2003")->year == 2003);
static_assert(refs::find("missing") == nullptr);
static_assert(refs::for_tool(refs::Tool::hubbard_pbc).front().reference == refs::find("lieb-wu-2003"));

TEST(Citations, RegistryUniqueIdsAndLookup)
{
  std::set<std::string_view> ids;
  for (auto const& ref : refs::references)
  {
    SCOPED_TRACE(std::string(ref.id));
    EXPECT_TRUE(ids.insert(ref.id).second);
    EXPECT_EQ(refs::find(ref.id), &ref);
    EXPECT_FALSE(ref.links.empty());
  }
}

TEST(Citations, ToolReferenceSelections)
{
  for (auto tool : {refs::Tool::xxx_pbc,
                    refs::Tool::xxx_obc,
                    refs::Tool::xxz_pbc,
                    refs::Tool::xxz_obc,
                    refs::Tool::hubbard_pbc,
                    refs::Tool::hubbard_obc,
                    refs::Tool::lieb_liniger_pbc,
                    refs::Tool::su3_pbc,
                    refs::Tool::gaudin_yang_pbc,
                    refs::Tool::tj_pbc,
                    refs::Tool::tb_pbc,
                    refs::Tool::richardson,
                    refs::Tool::central_spin,
                    refs::Tool::sun_fermions_pbc,
                    refs::Tool::ladder_pbc,
                    refs::Tool::hubbard_dispersion,
                    refs::Tool::haldane_shastry_pbc,
                    refs::Tool::sutherland_pbc,
                    refs::Tool::biquadratic_obc,
                    refs::Tool::lieb_liniger_obc,
                    refs::Tool::lieb_liniger_dispersion,
                    refs::Tool::lieb_liniger_thermal})
  {
    SCOPED_TRACE(static_cast<int>(tool));
    std::set<std::string_view> ids;
    auto const uses = refs::for_tool(tool);
    ASSERT_FALSE(uses.empty());
    for (auto const& use : uses)
    {
      ASSERT_NE(use.reference, nullptr);
      SCOPED_TRACE(std::string(use.reference->id));
      EXPECT_TRUE(ids.insert(use.reference->id).second);
      EXPECT_EQ(refs::find(use.reference->id), use.reference);
      EXPECT_FALSE(use.context.empty());
    }
  }
}

TEST(Citations, InvalidTool) { EXPECT_THROW((void)refs::for_tool(static_cast<refs::Tool>(-1)), std::invalid_argument); }
