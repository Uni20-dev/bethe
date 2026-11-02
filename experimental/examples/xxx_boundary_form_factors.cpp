// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include <bethe/experimental/xxx_boundary_form_factors.hpp>
#include <iomanip>
#include <iostream>

int main()
{
  namespace model = bethe::heisenberg::open;
  namespace boundary = bethe::experimental::xxx::open;
  auto ground = model::ground_state<double>(8);
  auto excited = model::sector_ground_state<double>(8, uni20::half_int{1});
  boundary::BoundaryFormFactors<double> prepared(8, ground);
  auto factor = prepared(excited);
  if (!factor.converged())
  {
    std::cerr << "boundary form factor rejected at working precision\n";
    return 2;
  }
  std::cout << std::setprecision(17) << "sites=8\nground_energy=" << *prepared.ground_energy()
            << "\ngap=" << excited.energy - *prepared.ground_energy() << "\nraising_weight=" << *factor.weight
            << "\nzz_weight=" << *factor.weight / 2 << '\n';
}
