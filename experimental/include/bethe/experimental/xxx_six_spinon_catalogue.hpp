// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/experimental/xxx_six_spinon.hpp>
#include <cstdint>
#include <limits>
#include <tuple>

namespace bethe::experimental::xxx::six_spinon
{
namespace detail
{
inline int modulo(std::int64_t value, int sites)
{
  value %= sites;
  return int(value < 0 ? value + sites : value);
}

inline Configuration canonical_request(Configuration c)
{
  (void)labels(c);
  std::sort(c.holes.begin(), c.holes.end());
  std::sort(c.branches.begin(), c.branches.end());
  return c;
}

inline bool request_less(Configuration const& a, Configuration const& b)
{
  return std::tie(a.sites, a.family, a.branches, a.holes) < std::tie(b.sites, b.family, b.branches, b.holes);
}

inline std::size_t checked_count_sum(std::size_t a, std::size_t b)
{
  if (b > std::numeric_limits<std::size_t>::max() - a) throw std::overflow_error("six-spinon catalogue count overflow");
  return a + b;
}
} // namespace detail

/// Physical label identity: all hole/branch numbers after h/b are TWICE the quantum numbers.
inline std::string candidate_id(Configuration c)
{
  c = detail::canonical_request(std::move(c));
  std::string id = "N" + std::to_string(c.sites) + (c.family == Family::two_pairs ? "-22-h" : "-3-h");
  for (std::size_t i = 0; i < c.holes.size(); ++i)
    id += (i ? "," : "") + std::to_string(c.holes[i].twice());
  id += "-b";
  for (std::size_t i = 0; i < c.branches.size(); ++i)
    id += (i ? "," : "") + std::to_string(c.branches[i].twice());
  return id;
}

/// q=P0-Pn, in [0,N); the supported N divisible by four has P0=0.
inline int transfer_index(Configuration const& c)
{
  return detail::modulo(-std::int64_t(labels(c).momentum_index), c.sites);
}

/// Count and unrank nominal Takahashi candidates at one momentum transfer.
/// This enumerates labels, not a guarantee of regular/resolved Bethe states.
/// Suffix counts avoid scanning other momenta. No reflection/parity folding:
/// each physical label tuple has multiplicity one.
class MomentumCatalogue {
  public:
    /// The budget counts stored size_t cells, not bytes. It is checked before allocation.
    explicit MomentumCatalogue(int sites, int transfer, std::size_t max_table_entries = 4'000'000)
        : sites_(sites), transfer_(transfer)
    {
      if (sites < 12 || sites % 4 || transfer < 0 || transfer >= sites)
        throw std::invalid_argument("require N>=12 divisible by four and transfer index in [0,N)");
      auto const n = static_cast<std::size_t>(sites);
      // Two tables: 7*N*(N/2+2 + N/2+3) cells. Divide before multiplying.
      if (n > max_table_entries / 7 / (n + 5))
        throw std::length_error("six-spinon catalogue exceeds max_table_entries");
      tables_[0] = make_table(Family::two_pairs);
      tables_[1] = make_table(Family::triple);
      for (int a : {-3, -1, 1, 3})
        for (int b : {-3, -1, 1, 3})
          if (a < b) add_group(0, {uni20::from_twice(a), uni20::from_twice(b)});
      for (int b : {-2, 0, 2})
        add_group(1, {uni20::from_twice(b)});
    }

    int sites() const { return sites_; }
    int momentum_transfer_index() const { return transfer_; }
    std::size_t size() const { return size_; }
    std::size_t family_size(Family family) const
    {
      if (family != Family::two_pairs && family != Family::triple)
        throw std::invalid_argument("invalid six-spinon family");
      return family_counts_[family == Family::two_pairs ? 0 : 1];
    }

    /// Stable order: family 22, then 3; increasing branches, then lexicographic holes.
    Configuration candidate(std::size_t ordinal) const
    {
      if (ordinal >= size_) throw std::out_of_range("six-spinon candidate ordinal outside catalogue");
      for (auto const& group : groups_)
      {
        if (ordinal >= group.count)
        {
          ordinal -= group.count;
          continue;
        }
        auto const& table = tables_[group.table];
        Configuration result{sites_, group.table == 0 ? Family::two_pairs : Family::triple, {}, group.branches};
        int start = 0, residue = group.residue;
        for (int remaining = 6; remaining > 0; --remaining)
        {
          bool found = false;
          for (int i = start; i <= table.length - remaining; ++i)
          {
            int const next = detail::modulo(residue - i, sites_);
            auto const count = table.at(i + 1, remaining - 1, next, sites_);
            if (ordinal < count)
            {
              result.holes.push_back(uni20::from_twice(table.base + 2 * i));
              start = i + 1;
              residue = next;
              found = true;
              break;
            }
            ordinal -= count;
          }
          if (!found) throw std::logic_error("six-spinon suffix unranking exhausted");
        }
        if (residue || ordinal) throw std::logic_error("six-spinon suffix unranking did not terminate");
        return result;
      }
      throw std::logic_error("six-spinon group unranking exhausted");
    }

    /// Materialization budget is independent of the suffix-table storage budget.
    std::vector<Configuration> requests(std::size_t max_candidates = 10000) const
    {
      if (!max_candidates) throw std::invalid_argument("max_candidates must be positive");
      if (size_ > max_candidates) throw std::length_error("six-spinon catalogue exceeds max_candidates");
      std::vector<Configuration> result;
      result.reserve(size_);
      for (std::size_t i = 0; i < size_; ++i)
        result.push_back(candidate(i));
      return result;
    }

  private:
    struct Table
    {
        int base{}, length{};
        std::vector<std::size_t> values;
        std::size_t& at(int start, int count, int residue, int n)
        {
          return values[(static_cast<std::size_t>(start) * 7 + count) * n + residue];
        }
        std::size_t at(int start, int count, int residue, int n) const
        {
          return values[(static_cast<std::size_t>(start) * 7 + count) * n + residue];
        }
    };
    struct Group
    {
        std::size_t table;
        std::vector<uni20::half_int> branches;
        int residue;
        std::size_t count;
    };
    Table make_table(Family family) const
    {
      int const edge = sites_ / 2 + (family == Family::triple ? 1 : 0);
      Table table{-edge, edge + 1, {}};
      table.values.resize(static_cast<std::size_t>(table.length + 1) * 7 * sites_);
      table.at(table.length, 0, 0, sites_) = 1;
      for (int start = table.length - 1; start >= 0; --start)
      {
        table.at(start, 0, 0, sites_) = 1;
        for (int k = 1; k <= std::min(6, table.length - start); ++k)
          for (int r = 0; r < sites_; ++r)
            table.at(start, k, r, sites_) =
                detail::checked_count_sum(table.at(start + 1, k, r, sites_),
                                          table.at(start + 1, k - 1, detail::modulo(r - start, sites_), sites_));
      }
      return table;
    }
    void add_group(std::size_t index, std::vector<uni20::half_int> branches)
    {
      auto const& table = tables_[index];
      int branch_sum = 0;
      for (auto b : branches)
        branch_sum += b.twice();
      // Six holes sum to 3*base + sum(indices). Pn = pi + 2*pi*(holes-branches)/N.
      int const residue =
          detail::modulo(-std::int64_t(transfer_) - sites_ / 2 - 3LL * table.base + branch_sum / 2, sites_);
      auto const count = table.at(0, 6, residue, sites_);
      groups_.push_back({index, std::move(branches), residue, count});
      size_ = detail::checked_count_sum(size_, count);
      family_counts_[index] = detail::checked_count_sum(family_counts_[index], count);
    }
    int sites_, transfer_;
    std::array<Table, 2> tables_;
    std::vector<Group> groups_;
    std::size_t size_{};
    std::array<std::size_t, 2> family_counts_{};
};
} // namespace bethe::experimental::xxx::six_spinon
