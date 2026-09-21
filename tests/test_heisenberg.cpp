// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch

#include <bethe/heisenberg.hpp>
#include <uni20/core/scalar_io.hpp>

#include <iostream>
#include <string>
#include <string_view>

namespace
{
// These checks remain active in Release builds (unlike assert).
void require(bool condition, std::string_view message)
{
  if (!condition)
    throw std::runtime_error(std::string(message));
}

template <typename Function> void invalid_argument(Function&& function)
{
  try
  {
    function();
  }
  catch (std::invalid_argument const&)
  {
    return;
  }
  throw std::runtime_error("expected invalid_argument");
}

template <uni20::Real Real> void check_diagnostics(std::size_t sites, bethe::heisenberg::GroundState<Real> const& result)
{
  using std::abs;
  using std::atan;
  Real const pi = Real{4} * atan(Real{1});
  Real residual = Real{0};
  Real energy = static_cast<Real>(sites) / Real{4};
  for (std::size_t i = 0; i < result.rapidities.size(); ++i)
  {
    Real const z = result.rapidities[i];
    Real const quantum_number = static_cast<Real>(i) - static_cast<Real>(sites) / Real{4} + Real{1} / Real{2};
    // Direct evaluation of N*phi - 2*pi*I - sum(phi), independently of the
    // solver's fixed-point angle calculation and compensated accumulation.
    Real equation = static_cast<Real>(sites) * Real{2} * atan(z) - Real{2} * pi * quantum_number;
    for (std::size_t j = 0; j < result.rapidities.size(); ++j)
      if (i != j)
        equation -= Real{2} * atan((z - result.rapidities[j]) / Real{2});
    residual = std::max(residual, abs(equation) / static_cast<Real>(sites));
    energy -= Real{2} / (Real{1} + z * z);
  }
  Real const rounding = Real{4} * static_cast<Real>(sites) * uni20::numeric_limits<Real>::epsilon();
  require(abs(residual - result.residual_norm) <= rounding, "residual does not describe the returned roots");
  require(abs(energy - result.energy) <= rounding, "energy does not describe the returned roots");
}

template <uni20::Real Real> void tests()
{
  using std::abs;
  using std::sqrt;
  using bethe::heisenberg::ground_state;
  using Options = bethe::heisenberg::SolverOptions<Real>;
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  Real const tolerance = Real{2048} * eps;

  auto const two = ground_state<Real>(2);
  require(two.converged && two.iterations == 0, "N=2 initial roots should already be exact");
  require(two.rapidities.size() == 1 && two.rapidities[0] == Real{0}, "N=2 rapidity");
  require(two.energy == -Real{3} / Real{2} && two.residual_norm == Real{0}, "N=2 double-bond normalization");

  auto const four = ground_state<Real>(4);
  require(four.converged, "N=4 convergence");
  require(abs(four.energy + Real{2}) < tolerance, "N=4 exact energy");
  require(abs(four.rapidities[0] + Real{1} / sqrt(Real{3})) < tolerance, "N=4 exact negative root");
  require(abs(four.rapidities[1] - Real{1} / sqrt(Real{3})) < tolerance, "N=4 exact positive root");

  // At N=6 the roots are {-a,0,a}; a^2 solves 3*a^4+10*a^2-9=0.
  // This gives E=-(2+sqrt(13))/2, an irrational precision-sensitive oracle.
  auto const six = ground_state<Real>(6);
  Real const exact_six = -(Real{2} + sqrt(Real{13})) / Real{2};
  Real const a = sqrt((Real{2} * sqrt(Real{13}) - Real{5}) / Real{3});
  Real const precision_tolerance = Real{128} * eps;
  require(six.converged && abs(six.energy - exact_six) < precision_tolerance, "N=6 exact energy at selected precision");
  require(abs(six.rapidities[0] + a) < tolerance && abs(six.rapidities[2] - a) < tolerance,
          "N=6 exact outer roots");
  require(abs(six.rapidities[1]) < tolerance, "N=6 central root");
  if constexpr (uni20::numeric_limits<Real>::digits > uni20::numeric_limits<double>::digits)
    require(abs(static_cast<Real>(static_cast<double>(exact_six)) - exact_six) > precision_tolerance,
            "high-precision oracle must distinguish a double-only implementation");
  require(uni20::parse_real<Real>(uni20::format_real(six.energy)) == six.energy, "precision-preserving text round trip");

  for (std::size_t const sites : {4, 6, 16, 64})
  {
    auto const result = ground_state<Real>(sites);
    require(result.converged && result.residual_norm <= Options{}.residual_tolerance, "equation convergence");
    require(result.iterations <= Options{}.max_iterations, "iteration budget");
    require(result.rapidities.size() == sites / 2, "ground-state root count");
    check_diagnostics(sites, result);
    for (std::size_t i = 0; i < result.rapidities.size(); ++i)
    {
      if (i != 0)
        require(result.rapidities[i - 1] < result.rapidities[i], "roots must be ordered");
      require(abs(result.rapidities[i] + result.rapidities[result.rapidities.size() - 1 - i]) < tolerance,
              "ground-state reflection symmetry");
    }
  }

  // Independent published finite-size value: Table I of cond-mat/9809163,
  // relative to the ferromagnetic energy. Its decimal precision limits the test.
  auto const sixteen = ground_state<Real>(16);
  Real const reference = uni20::parse_real<Real>("-0.696393522538549");
  require(abs(sixteen.energy / Real{16} - Real{1} / Real{4} - reference) < uni20::parse_real<Real>("2e-14"),
          "N=16 published energy convention");

  auto const baseline = ground_state<double>(32);
  auto const selected = ground_state<Real>(32);
  require(baseline.converged && selected.converged, "cross-precision convergence");
  require(abs(selected.energy - static_cast<Real>(baseline.energy)) <
              Real{2048} * static_cast<Real>(uni20::numeric_limits<double>::epsilon()),
          "cross-precision energy agreement");

  auto const no_updates = ground_state<Real>(4, Options{.max_iterations = 0});
  require(!no_updates.converged && no_updates.iterations == 0, "zero update budget");
  require(no_updates.energy == -Real{3}, "initial-guess energy");
  check_diagnostics(4, no_updates);
  auto const one_update = ground_state<Real>(4, Options{.max_iterations = 1});
  require(!one_update.converged && one_update.iterations == 1, "iteration exhaustion must not look converged");
  require(abs(one_update.energy + Real{1} + sqrt(Real{2})) < tolerance, "one simultaneous update, not two");
  check_diagnostics(4, one_update);

  for (std::size_t const sites : {0, 1, 3, 7})
    invalid_argument([&] { (void)ground_state<Real>(sites); });
  for (Real const bad_tolerance : {Real{0}, -Real{1}, uni20::numeric_limits<Real>::infinity(),
                                  uni20::numeric_limits<Real>::quiet_NaN()})
    invalid_argument([&] { (void)ground_state<Real>(4, Options{.residual_tolerance = bad_tolerance}); });
}
} // namespace

int main(int argc, char** argv)
{
  try
  {
    require(argc == 2, "pass a precision");
    std::string_view const precision = argv[1];
    if (precision == "fp64")
      tests<double>();
    else if (precision == "long-double")
      tests<long double>();
#if UNI20_HAS_FLOAT128
    else if (precision == "fp128")
      tests<uni20::float128>();
#endif
    else
      throw std::invalid_argument("unavailable test precision");
    std::cout << "Heisenberg checks passed for " << precision << '\n';
    return 0;
  }
  catch (std::exception const& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
