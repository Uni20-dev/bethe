// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once
#include <array>
#include <stdexcept>
#include <uni20/linalg/ops/self_adjoint_eigh.hpp>
#include <vector>

namespace bethe::test
{
struct TBExactGround
{
    double energy, translation;
};
// Independent spin-basis H=sum(X-X^2), X=S.S. Work in Sz=0, which contains
// every integer-spin SU(2) multiplet, so its minimum is the absolute minimum.
inline TBExactGround tb_exact_ground(unsigned n)
{
  if (n < 3 || n > 8) throw std::invalid_argument("small TB ED oracle domain");
  std::vector<unsigned> powers(n + 1, 1);
  for (unsigned j = 1; j <= n; ++j)
    powers[j] = 3 * powers[j - 1];
  std::vector<unsigned> basis, index(powers[n]);
  for (unsigned s = 0; s < powers[n]; ++s)
  {
    unsigned sum = 0;
    for (unsigned j = 0; j < n; ++j)
      sum += (s / powers[j]) % 3;
    if (sum == n)
    {
      index[s] = basis.size();
      basis.push_back(s);
    }
  }
  std::array<std::array<double, 9>, 9> dot{}, local{};
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b)
    {
      dot[3 * a + b][3 * a + b] = (a - 1) * (b - 1);
      if (a < 2 && b > 0) dot[3 * (a + 1) + b - 1][3 * a + b] = 1;
      if (a > 0 && b < 2) dot[3 * (a - 1) + b + 1][3 * a + b] = 1;
    }
  for (unsigned a = 0; a < 9; ++a)
    for (unsigned b = 0; b < 9; ++b)
    {
      local[a][b] = dot[a][b];
      for (unsigned k = 0; k < 9; ++k)
        local[a][b] -= dot[a][k] * dot[k][b];
    }
  auto const dim = basis.size();
  uni20::DenseMatrix<double> h(dim, dim);
  for (std::size_t i = 0; i < dim; ++i)
    for (std::size_t j = 0; j < dim; ++j)
      h[i, j] = 0;
  for (std::size_t c = 0; c < dim; ++c)
    for (unsigned j = 0; j < n; ++j)
    {
      unsigned const k = (j + 1) % n, s = basis[c], a = (s / powers[j]) % 3, b = (s / powers[k]) % 3;
      for (unsigned aa = 0; aa < 3; ++aa)
        for (unsigned bb = 0; bb < 3; ++bb)
        {
          double const v = local[3 * aa + bb][3 * a + b];
          if (v != 0) h[index[s - a * powers[j] - b * powers[k] + aa * powers[j] + bb * powers[k]], c] += v;
        }
    }
  auto const eig = uni20::linalg::eigh(std::move(h));
  double translation = 0;
  for (std::size_t c = 0; c < dim; ++c)
  {
    unsigned const s = basis[c], translated = (s * 3) % powers[n] + s / powers[n - 1];
    translation += eig.eigenvectors[index[translated], 0] * eig.eigenvectors[c, 0];
  }
  return {eig.eigenvalues[0], translation};
}
} // namespace bethe::test
