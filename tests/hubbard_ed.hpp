// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bit>
#include <cmath>
#include <uni20/linalg/ops/self_adjoint_eigh.hpp>
#include <vector>

namespace bethe::test
{
// Independent Fock-space Hamiltonian, ordered up orbitals then down orbitals.
// Both hopping and translation include the fermionic reordering signs.
struct ExactGround
{
    double energy, translation;
    std::vector<double> momentum_weights;
};
inline ExactGround exact_ground(unsigned n, double u, unsigned num_up, unsigned num_down, bool periodic = true)
{
  unsigned const mask = (1U << n) - 1;
  std::vector<unsigned> basis;
  std::vector<std::size_t> index(1U << (2 * n));
  for (unsigned up = 0; up <= mask; ++up)
    for (unsigned down = 0; down <= mask; ++down)
      if (std::popcount(up) == int(num_up) && std::popcount(down) == int(num_down))
      {
        unsigned const state = up | (down << n);
        index[state] = basis.size();
        basis.push_back(state);
      }
  uni20::DenseMatrix<double> h(basis.size(), basis.size());
  for (std::size_t row = 0; row < basis.size(); ++row)
    for (std::size_t col = 0; col < basis.size(); ++col)
      h[row, col] = 0;
  for (std::size_t col = 0; col < basis.size(); ++col)
  {
    unsigned const state = basis[col];
    h[col, col] = u * std::popcount((state & mask) & (state >> n));
    for (unsigned spin = 0; spin < 2; ++spin)
      for (unsigned j = 0; j < (periodic ? n : n - 1); ++j)
        for (unsigned direction = 0; direction < 2; ++direction)
        {
          unsigned const from = spin * n + (direction ? j : (j + 1) % n);
          unsigned const to = spin * n + (direction ? (j + 1) % n : j);
          if (!(state & (1U << from)) || (state & (1U << to))) continue;
          unsigned const removed = state ^ (1U << from), added = removed | (1U << to);
          int const parity = std::popcount(state & ((1U << from) - 1)) + std::popcount(removed & ((1U << to) - 1));
          h[index[added], col] -= parity % 2 ? -1.0 : 1.0;
        }
  }
  auto const eig = uni20::linalg::eigh(std::move(h));
  if (!periodic) return {eig.eigenvalues[0], 0, {}};
  double translation = 0;
  std::vector<std::size_t> translated_index(basis.size());
  std::vector<double> translation_sign(basis.size());
  for (std::size_t col = 0; col < basis.size(); ++col)
  {
    unsigned translated = 0;
    int parity = 0;
    for (unsigned spin = 0; spin < 2; ++spin)
    {
      unsigned const bits = (basis[col] >> (spin * n)) & mask;
      translated |= (((bits << 1) & mask) | (bits >> (n - 1))) << (spin * n);
      if (bits & (1U << (n - 1))) parity += std::popcount(bits) - 1;
    }
    translation += (parity % 2 ? -1.0 : 1.0) * eig.eigenvectors[index[translated], 0] * eig.eigenvectors[col, 0];
    translated_index[col] = index[translated];
    translation_sign[col] = parity % 2 ? -1.0 : 1.0;
  }
  // Fourier-project the entire degenerate ground eigenspace, not a single
  // arbitrary real eigenvector (which need not have definite momentum).
  std::vector<double> weights(n, 0);
  double const pi = 4 * std::atan(1.0);
  for (std::size_t v = 0; v < basis.size() && eig.eigenvalues[v] - eig.eigenvalues[0] < 1e-9; ++v)
    for (std::size_t col = 0; col < basis.size(); ++col)
    {
      std::size_t row = col;
      double sign = 1;
      for (unsigned power = 0; power < n; ++power)
      {
        double const contribution = sign * eig.eigenvectors[row, v] * eig.eigenvectors[col, v] / n;
        for (unsigned q = 0; q < n; ++q)
          weights[q] += contribution * std::cos(2 * pi * q * power / n);
        sign *= translation_sign[row];
        row = translated_index[row];
      }
    }
  return {eig.eigenvalues[0], translation, weights};
}

inline ExactGround exact_ground(unsigned n, double u) { return exact_ground(n, u, n / 2, n / 2); }

} // namespace bethe::test
