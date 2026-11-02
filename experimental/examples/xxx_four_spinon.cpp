// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include <bethe/experimental/xxx_four_spinon.hpp>
#include <iostream>

int main()
try
{
  namespace xxx = bethe::experimental::xxx;
#if UNI20_HAS_MPREAL
  using R = uni20::mpreal;
  R const zero(0, uni20::Precision::bits(128));
  bethe::SolverOptions<R> options(zero.precision());
#else
  using R = double;
  R const zero = 0;
  bethe::SolverOptions<R> options;
#endif
  xxx::FourSpinonRequest const request{12, {-3, -1, 0, 2}, 1};
  auto const gs = bethe::heisenberg::sector_ground_state<R>(request.sites, uni20::half_int{0}, options);
  xxx::PreparedSinglet<R> const ground(request.sites, gs);
  if (!ground.ready())
  {
    std::cerr << "ground state unresolved\n";
    return 1;
  }
  auto const selected = xxx::solve_four_spinon<R>(request, zero);
  if (!selected.converged())
  {
    std::cerr << selected.diagnostic << '\n';
    return 1;
  }
  auto const weight = xxx::four_spinon_raising_weight(ground, *selected.state);
  if (!weight.converged())
  {
    std::cerr << "form factor unresolved\n";
    return 1;
  }
  std::cout << "# Selected finite N=12 four-spinon triplet, holes -3,-1,0,2, branch 1\n"
            << "# J=1; Fourier raising weight; no 2*pi factor; one state, not a sector sum\n"
            << "momentum_index,gap,raising_weight,deviation,rational_residual\n"
            << selected.state->momentum_index() << ','
            << uni20::format_scalar(selected.state->energy() - *ground.energy()) << ','
            << uni20::format_scalar(*weight.weight) << ',' << uni20::format_scalar(selected.state->deviation()) << ','
            << uni20::format_scalar(selected.state->rational_residual()) << '\n';
}
catch (std::exception const& error)
{
  std::cerr << error.what() << '\n';
  return 1;
}
