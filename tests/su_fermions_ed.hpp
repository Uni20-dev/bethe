// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <bit>
#include <cstdint>
#include <functional>
#include <map>
#include <numbers>
#include <span>
#include <stdexcept>
#include <uni20/linalg/ops/self_adjoint_eigh.hpp>
#include <vector>

namespace bethe::test
{
// Independent variational continuum Hamiltonian in a truncated plane-wave
// Fock basis, with total momentum zero. No Bethe equations or BA roots used.
inline double su_fermions_cutoff_ground(std::span<unsigned const> populations, int cutoff, double length, double c)
{
  if (cutoff < 1 || cutoff > 4 || populations.empty() || populations.size() > 3)
    throw std::invalid_argument("small SU fermion ED oracle domain");
  unsigned const modes = 2 * cutoff + 1, mask = (1U << modes) - 1;
  std::vector<std::vector<unsigned>> choices(populations.size());
  for (std::size_t a = 0; a < populations.size(); ++a)
    for (unsigned b = 0; b <= mask; ++b)
      if (std::popcount(b) == int(populations[a])) choices[a].push_back(b);
  std::vector<std::uint64_t> basis;
  std::map<std::uint64_t, std::size_t> index;
  std::function<void(std::size_t, std::uint64_t, int)> enumerate = [&](std::size_t a, std::uint64_t key, int momentum) {
    if (a == populations.size())
    {
      if (!momentum)
      {
        index[key] = basis.size();
        basis.push_back(key);
      }
      return;
    }
    for (unsigned b : choices[a])
    {
      int p = momentum;
      for (unsigned j = 0; j < modes; ++j)
        if (b & (1U << j)) p += int(j) - cutoff;
      enumerate(a + 1, key | (std::uint64_t(b) << (a * modes)), p);
    }
  };
  enumerate(0, 0, 0);
  if (basis.empty()) throw std::invalid_argument("empty cutoff ED sector");
  auto const dim = basis.size();
  uni20::DenseMatrix<double> h(dim, dim);
  for (std::size_t i = 0; i < dim; ++i)
    for (std::size_t j = 0; j < dim; ++j)
      h[i, j] = 0;
  auto move = [](unsigned bits, unsigned from, unsigned to) {
    unsigned const lo = std::min(from, to), hi = std::max(from, to);
    unsigned const between = hi == lo ? 0 : ((1U << hi) - 1) ^ ((1U << (lo + 1)) - 1);
    return std::popcount(bits & between) % 2 ? -1. : 1.;
  };
  double const unit = 2 * std::numbers::pi / length;
  for (std::size_t col = 0; col < dim; ++col)
  {
    auto const key = basis[col];
    for (std::size_t a = 0; a < populations.size(); ++a)
    {
      unsigned const first = (key >> (a * modes)) & mask;
      for (unsigned j = 0; j < modes; ++j)
        if (first & (1U << j)) h[col, col] += unit * unit * (int(j) - cutoff) * (int(j) - cutoff);
      for (std::size_t b = a + 1; b < populations.size(); ++b)
      {
        unsigned const second = (key >> (b * modes)) & mask;
        for (unsigned p = 0; p < modes; ++p)
          if (first & (1U << p))
            for (unsigned q = 0; q < modes; ++q)
              if (second & (1U << q))
                for (unsigned target = 0; target < modes; ++target)
                {
                  int const other = int(p) + int(q) - int(target);
                  if (other < 0 || other >= int(modes)) continue;
                  if (((first ^ (1U << p)) & (1U << target)) || ((second ^ (1U << q)) & (1U << other))) continue;
                  auto const newkey = key ^ (std::uint64_t((1U << p) ^ (1U << target)) << (a * modes)) ^
                                      (std::uint64_t((1U << q) ^ (1U << other)) << (b * modes));
                  h[index.at(newkey), col] +=
                      2 * c / length * move(first, p, target) * move(second, q, unsigned(other));
                }
      }
    }
  }
  return uni20::linalg::eigh(std::move(h)).eigenvalues[0];
}
} // namespace bethe::test
