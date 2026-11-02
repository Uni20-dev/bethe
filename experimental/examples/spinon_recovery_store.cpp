// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#include <bethe/experimental/spinon_recovery_store.hpp>
#include <iostream>

int main(int argc, char** argv)
try
{
  namespace x = bethe::experimental::xxx;
  if (argc < 2 || argc > 3 || (argc == 3 && std::string_view(argv[2]) != "--stop-after-control"))
  {
    std::cerr << "Usage: bethe-experimental-recovery-store-example DIRECTORY [--stop-after-control]\n";
    return 1;
  }
  int const bits = uni20::numeric_limits<long double>::digits();
  if (bits <= 53)
  {
    std::cerr << "This example requires long double wider than binary64.\n";
    return 77;
  }
  x::RecoveryChunkBinding binding{
      "bethe-recovery-store-example-v1",
      "fresh W2/W4/W6 at fixed decimal rapidities; longitudinal labelled dp",
      {"1e-6", "0", {{53, bits, {{{32, 12}, {64, 20}}}}}, 60},
      {{"two", {"-2", "1"}}, {"four", {"-2.1", "-.4", ".3", "2"}}, {"six", {"-2.1", "-1.3", "-.4", ".3", "1.2", "2"}}}};
  x::SpinonRecoveryStore<double, long double> store(argv[1], std::move(binding));
  std::stop_source stop;
  for (std::size_t i = 0; i < store.size(); ++i)
  {
    std::size_t controls = 0;
    auto result = store.run(
        i,
        [&](auto const&) {
          ++controls;
          if (argc == 3) stop.request_stop();
        },
        stop.get_token());
    std::cout << store.point(i).id << ": " << controls << " new controls; ";
    if (!result.accepted())
    {
      std::cout << "no accepted value; " << result.reason << '\n';
      return 2;
    }
    std::cout << "W=" << uni20::format_real(result.estimate->value) << '\n';
  }
}
catch (std::exception const& error)
{
  std::cerr << error.what() << '\n';
  return 1;
}
