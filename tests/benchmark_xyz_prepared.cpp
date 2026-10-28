// SPDX-License-Identifier: GPL-3.0-or-later
// Library parameter-sweep benchmark; complete results checked before reporting.
#include <bethe/observable_grid.hpp>
#include <bethe/xyz_dispersion.hpp>
#include <chrono>
#include <ctime>
#include <iostream>
#include <string>
#include <uni20/async/tbb_scheduler.hpp>
#include <vector>

template <uni20::Real R> int run(std::size_t spectra, std::size_t points)
{
  auto evaluate = [=](std::size_t i) {
    bethe::xyz::SpinonDispersion<R> band(R{0.75} + R(i % 31) / R{256}, R{0.5} + R(i % 7) / R{4});
    std::vector<R> result;
    result.reserve(2 * points);
    R const end = 2 * bethe::detail::pi<R>();
    for (unsigned s : {1u, 2u})
    {
      auto prepared = band.prepare_bound(s);
      for (std::size_t j = 0; j < points; ++j)
        result.push_back(prepared.energy(j == points - 1 ? end : end * (R(j) / R(points - 1))));
    }
    return result;
  };
  std::vector<std::vector<R>> reference(spectra);
  auto wall = std::chrono::steady_clock::now();
  auto cpu = std::clock();
  for (std::size_t i = 0; i < spectra; ++i)
    reference[i] = evaluate(i);
  auto elapsed = [&] { return std::chrono::duration<double>(std::chrono::steady_clock::now() - wall).count(); };
  std::cout << "workers,wall_seconds,cpu_seconds,values\n0," << elapsed() << ','
            << double(std::clock() - cpu) / CLOCKS_PER_SEC << ',' << 2 * spectra * points << '\n';
  for (int workers : {1, 4, 8})
  {
    uni20::async::TbbScheduler scheduler(workers);
    std::vector<std::vector<R>> actual(spectra);
    wall = std::chrono::steady_clock::now();
    cpu = std::clock();
    auto operation = bethe::evaluate_observable_grid_async(
        spectra, evaluate, [&](std::size_t i, std::vector<R> result) { actual[i] = std::move(result); },
        {&scheduler, 128, 16});
    auto const completed = operation.get_wait(scheduler);
    double const seconds = elapsed(), cpu_seconds = double(std::clock() - cpu) / CLOCKS_PER_SEC;
    if (completed != spectra || actual != reference)
      throw std::runtime_error("prepared XYZ sweep differs from complete serial results");
    std::cout << workers << ',' << seconds << ',' << cpu_seconds << ',' << 2 * spectra * points << '\n';
  }
  return 0;
}

int main(int argc, char** argv)
{
  try
  {
    if (argc != 4) throw std::invalid_argument("usage: benchmark fp64|long-double|fp128 SPECTRA POINTS");
    std::string precision = argv[1];
    auto spectra = std::stoull(argv[2]), points = std::stoull(argv[3]);
    if (!spectra || spectra > 1024 || points < 2 || points > 32769 || spectra * points > 5000000)
      throw std::invalid_argument("benchmark size exceeds bounds");
    if (precision == "fp64") return run<double>(spectra, points);
    if (precision == "long-double") return run<long double>(spectra, points);
#if UNI20_HAS_FLOAT128
    if (precision == "fp128") return run<uni20::float128>(spectra, points);
#endif
    throw std::invalid_argument("unavailable benchmark precision");
  }
  catch (std::exception const& e)
  {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
