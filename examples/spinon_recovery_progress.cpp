// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include <bethe/xxx_spinon_density_recovery.hpp>
#include <iostream>

int main()
{
  using namespace bethe::heisenberg;
  using Driver = XXXSpinonDensityRecovery<double, long double>;
  int const low = uni20::numeric_limits<double>::digits(1.);
  int const guard = uni20::numeric_limits<long double>::digits(1.L);
  if (guard <= low)
  {
    std::cerr << "This example requires a distinct native long-double guard.\n";
    return 77;
  }
  Driver driver({"1e-6", "0", {{low, guard, {{{32, 12}, {64, 20}}}}}, 60});
  std::vector<std::string> const beta{"-2", "1"};
  std::stop_source stop;
  // Stop between complete controls; keep the full returned history in memory.
  auto partial = driver.resume(
      beta, {},
      [&](auto const& progress) {
        std::cout << "completed controls: " << progress.attempts.size() << '\n';
        stop.request_stop();
      },
      stop.get_token());
  if (partial.status != SpinonRecoveryStatus::cancelled || partial.attempts.size() != 1) return 2;
  std::size_t new_controls = 0;
  auto result = driver.resume(beta, partial, [&](auto const&) { ++new_controls; });
  if (!result.accepted() || new_controls != 1) return 2;
  std::cout << "new controls on resume: " << new_controls << '\n'
            << "longitudinal W2: " << uni20::format_real(result.estimate->value) << '\n';
  auto restored = driver.restore(result);
  if (!restored.accepted() || restored.estimate->value != result.estimate->value) return 2;
}
