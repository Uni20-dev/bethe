// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include <bethe/xxx_four_spinon_thermodynamics.hpp>
#include <iostream>
#include <uni20/async/tbb_scheduler.hpp>

int main()
try
{
  using namespace bethe::heisenberg;
#if UNI20_HAS_MPREAL
  using R = uni20::mpreal;
  using Low = R;
  R const one(1, uni20::Precision::bits(192));
  SpinonRecoveryPolicy recovery{"1e-6", "1e-30", {{128, 192, {{{64, 16}, {80, 20}}}}}, 60};
#else
  using R = long double;
  using Low = double;
  R const one = 1;
  SpinonRecoveryPolicy recovery{
      "1e-4", "1e-24", {{53, uni20::numeric_limits<R>::digits(one), {{{64, 16}, {80, 20}}}}}, 60};
#endif
  ThermodynamicFourSpinon<R, Low, R> model(recovery, one);
  uni20::async::TbbScheduler scheduler(2);
  auto pending = model.refined_global_moments(2, 3, {one, one, one}, {2 * one, 2 * one, 2 * one})
                     .evaluate_async({&scheduler, 8, 2});
  auto result = pending.get_wait(scheduler);
  if (!result.estimates)
  {
    std::cerr << result.reason << '\n';
    return 1;
  }
  std::cout << "Finite gap box (0,2)^3; unrestricted common shift. This is not a converged thermodynamic moment.\n";
  std::cout << "moment,value,fraction,order_difference,cutoff_change\n";
  for (int moment = 0; moment < 2; ++moment)
    std::cout << moment << ',' << uni20::format_scalar((*result.estimates)[moment].value) << ','
              << uni20::format_scalar(*result.sum_rule_fractions[moment]) << ','
              << uni20::format_scalar((*result.quadrature_difference)[moment]) << ','
              << uni20::format_scalar((*result.cutoff_change)[moment]) << '\n';
}
catch (std::exception const& error)
{
  std::cerr << error.what() << '\n';
  return 1;
}
