// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/polynomial_roots.hpp>
#include <bethe/real_excitations.hpp>
#include <bethe/xxz_quantum_group.hpp>

namespace bethe::xxz::quantum_group::qsystem
{
/// Open quantum-group XXZ, NOT the zero-field open chain. Bajnok et al.
/// arXiv:1910.07805, Sec. 5, especially the Wronskian (5.17).
/// Q(x)=x^M+sum c[k]x^k, x=cosh(2u)=cos(alpha). Real coefficients can
/// represent conjugate complex roots; no ideal-string approximation is used.
template <uni20::Real Real> class System {
  public:
    System(std::size_t sites, std::size_t roots, Real delta) : sites(sites), order(roots), delta(delta)
    {
      (void)quantum_group::detail::sector_roots(sites, 0);
      if (roots > sites / 2) throw std::invalid_argument("open Q-system requires M <= N/2");
      if (!uni20::isfinite(delta) || delta <= Real{1})
        throw std::invalid_argument("open Q-system requires finite Delta > 1");
      // This coefficient formulation is a small-system/selected-state engine.
      // Bound its O(N^3) cached polynomial basis before allocating anything.
      if (sites > 32) throw std::length_error("open Q-system currently supports N <= 32");
      degree = sites - order + 1;
      Real const s2 = (delta - Real{1}) * (delta + Real{1});
      // B_ij=[x_+^i x_-^j-x_-^i x_+^j]/(2 sinh(eta) sqrt(x^2-1)),
      // x_+ + x_-=2 Delta x, x_+ x_-=x^2+sinh(eta)^2.
      // D_k=(x_+^k-x_-^k)/(x_+-x_-) has a division-free recurrence.
      std::vector<std::vector<Real>> d(degree + 1, std::vector<Real>(sites + 1));
      d[1][0] = Real{1};
      for (std::size_t k = 1; k < degree; ++k)
        for (std::size_t j = 0; j <= sites; ++j)
          d[k + 1][j] =
              (j ? Real{2} * delta * d[k][j - 1] : Real{0}) - s2 * d[k - 1][j] - (j >= 2 ? d[k - 1][j - 2] : Real{0});
      basis.resize((degree + 1) * (order + 1), std::vector<Real>(sites + 1));
      for (std::size_t i = 0; i <= degree; ++i)
        for (std::size_t j = 0; j <= order; ++j)
        {
          if (i == j) continue;
          auto b = d[i > j ? i - j : j - i];
          for (std::size_t k = 0; k < std::min(i, j); ++k)
            for (std::size_t a = sites + 1; a-- > 0;)
              b[a] = s2 * b[a] + (a >= 2 ? b[a - 2] : Real{0});
          if (i < j)
            for (auto& value : b)
              value = -value;
          for (Real value : b)
            if (!uni20::isfinite(value)) throw std::overflow_error("open Q-system basis overflow");
          basis[i * (order + 1) + j] = std::move(b);
        }
      rhs.resize(sites + 1);
      rhs[0] = Real{1}; // N is even: (x-1)^N.
      for (std::size_t k = 1; k <= sites; ++k)
        rhs[k] = -rhs[k - 1] * Real(sites - k + 1) / Real(k);
    }

    struct Evaluation
    {
        std::vector<Real> residual, partner;
        Real norm{}; // Full Wronskian coefficient backward residual.
    };

    /// Eliminate P linearly, fixing its gauge P_M=0. There are M remaining
    /// equations (coefficients 0,...,M-2 and 2M-1). Only constant nonzero
    /// q-number pivots are divided by; there is no division by unknown Q.
    /// Jacobian differentiates the raw residual, including the elimination.
    Evaluation evaluate(std::span<Real const> c, std::vector<Real>* jacobian = nullptr) const
    {
      validate(c);
      if (jacobian) jacobian->assign(order * order, Real{0});
      Evaluation result;
      for (std::size_t direction = 0; direction < (jacobian && order ? order : 1); ++direction)
      {
        std::vector<Real> rem = rhs, tangent(sites + 1), p(degree + 1);
        for (std::size_t i = degree + 1; i-- > 0;)
        {
          if (i == order) continue; // P -> P+aQ gauge.
          auto const k = i + order - 1;
          auto const& top = basis[i * (order + 1) + order];
          Real const pivot = top[k];
          p[i] = rem[k] / pivot;
          Real const dp = tangent[k] / pivot;
          for (std::size_t a = 0; a <= sites; ++a)
          {
            bethe::detail::CompensatedSum<Real> sum;
            sum.add(top[a]);
            for (std::size_t j = 0; j < order; ++j)
              sum.add(c[j] * basis[i * (order + 1) + j][a]);
            Real const value = sum.value();
            rem[a] -= p[i] * value;
            if (jacobian && order) tangent[a] -= dp * value + p[i] * basis[i * (order + 1) + direction][a];
          }
        }
        for (Real value : rem)
          if (!uni20::isfinite(value)) throw std::overflow_error("nonfinite open Q-system residual");
        if (direction == 0)
        {
          result.partner = p;
          result.residual.resize(order);
          for (std::size_t a = 0; a < order; ++a)
            result.residual[a] = rem[a + 1 == order ? 2 * order - 1 : a];
          // Independent coefficient-wise reconstruction, with a componentwise
          // absolute-product scale (not a rootwise or energy error bound).
          for (std::size_t a = 0; a <= sites; ++a)
          {
            bethe::detail::CompensatedSum<Real> sum;
            Real scale = std::abs(rhs[a]);
            sum.add(-rhs[a]);
            for (std::size_t i = 0; i <= degree; ++i)
              for (std::size_t j = 0; j <= order; ++j)
              {
                Real const term = p[i] * (j == order ? Real{1} : c[j]) * basis[i * (order + 1) + j][a];
                sum.add(term);
                scale += std::abs(term);
              }
            if (!uni20::isfinite(scale) || !uni20::isfinite(sum.value()))
              throw std::overflow_error("nonfinite open Q-system reconstruction");
            result.norm = std::max(result.norm, std::abs(sum.value()) / std::max(Real{1}, scale));
          }
        }
        if (jacobian)
          for (std::size_t a = 0; a < order; ++a)
          {
            auto const value = tangent[a + 1 == order ? 2 * order - 1 : a];
            if (!uni20::isfinite(value)) throw std::overflow_error("nonfinite open Q-system Jacobian");
            (*jacobian)[a * order + direction] = value;
          }
      }
      return result;
    }

    /// E=(N-1)Delta/4-(Delta^2-1) Q'(Delta)/Q(Delta).
    std::optional<Real> energy(std::span<Real const> c) const
    {
      validate(c);
      Real q{1}, derivative{};
      for (std::size_t k = order; k-- > 0;)
      {
        derivative = derivative * delta + q;
        q = q * delta + c[k];
      }
      if (q == Real{0}) return std::nullopt;
      Real const energy = Real(sites - 1) * delta / Real{4} - (delta - Real{1}) * (delta + Real{1}) * (derivative / q);
      return uni20::isfinite(energy) ? std::optional<Real>{energy} : std::nullopt;
    }

    std::size_t sites, order, degree;
    Real delta;

  private:
    void validate(std::span<Real const> c) const
    {
      if (c.size() != order) throw std::invalid_argument("open Q-system seed size differs from M");
      for (Real value : c)
        if (!uni20::isfinite(value)) throw std::invalid_argument("nonfinite open Q-system coefficient");
    }
    std::vector<std::vector<Real>> basis;
    std::vector<Real> rhs;
};

enum class Status
{
  converged,
  iteration_limit,
  singular_jacobian,
  stalled,
  unresolved_roots,
  inadmissible
};

template <uni20::Real Real> struct State
{
    std::size_t sites = 0, through_lines = 0, iterations = 0;
    Real delta{}, residual_norm{}, bethe_residual = uni20::numeric_limits<Real>::infinity();
    Real bethe_residual_bound = uni20::numeric_limits<Real>::infinity();
    // Root-recovery uncertainty propagated through sinh factors; not an interval enclosure.
    std::vector<Real> coefficients, partner;
    /// x=cosh(2u), NOT the real solver's alpha coordinate.
    bethe::detail::PolynomialRoots<Real> roots;
    std::vector<std::complex<Real>> rapidities; // Bajnok u, alpha=-2iu.
    std::optional<Real> energy;
    Status status = Status::iteration_limit;
    bool converged = false;
};

namespace detail
{
template <uni20::Real Real> void check_roots(State<Real>& state)
{
  using C = std::complex<Real>;
  auto const m = state.coefficients.size();
  Real const eps = uni20::numeric_limits<Real>::epsilon();
  state.roots = bethe::detail::recover_polynomial_roots<Real>(state.coefficients);
  if (state.roots.status != bethe::detail::PolynomialRootStatus::resolved)
  {
    state.status = Status::unresolved_roots;
    return;
  }
  // Numerical admissibility, not an interval certificate. Do not accept roots
  // indistinguishable from 0,i*pi/2,eta/2 in the u coordinate.
  Real const margin = Real{16} * state.roots.max_root_uncertainty + Real{128} * eps;
  for (C x : state.roots.roots)
  {
    if (std::abs(x - C{1}) <= margin || std::abs(x + C{1}) <= margin || std::abs(x - C{state.delta}) <= margin)
    {
      state.status = Status::inadmissible;
      return;
    }
    C u = std::acosh(x) / Real{2};
    if (u.real() == Real{0} && u.imag() < Real{0}) u = -u;
    state.rapidities.push_back(u);
  }
  Real const eta = std::acosh(state.delta);
  state.bethe_residual = Real{0};
  state.bethe_residual_bound = Real{0};
  std::vector<Real> radii;
  for (C x : state.roots.roots)
    radii.push_back(state.roots.max_root_uncertainty / (Real{2} * std::abs(std::sqrt(x * x - C{1}))));
  bool verified = true;
  for (std::size_t i = 0; i < m; ++i)
  {
    C a{1}, b{1};
    bool regular = true;
    Real error{};
    auto multiply = [&](C z, C w, Real radius) {
      C const f = std::sinh(z), g = std::sinh(w);
      // First-order relative product uncertainty. Require it small enough
      // that a linear propagation is meaningful; unresolved strings FAIL.
      auto relative_error = [&](C argument, C factor) {
        Real const magnitude = std::abs(factor);
        if (magnitude == Real{0}) return uni20::numeric_limits<Real>::infinity();
        return (std::abs(std::cosh(argument)) * radius + Real{32} * eps * std::max(Real{1}, magnitude)) / magnitude;
      };
      error += relative_error(z, f) + relative_error(w, g);
      a *= f;
      b *= g;
      Real const scale = std::max(std::abs(a), std::abs(b));
      if (!uni20::isfinite(scale) || scale == Real{0})
      {
        regular = false;
        return;
      }
      a /= scale;
      b /= scale;
    };
    C const u = state.rapidities[i];
    for (std::size_t k = 0; k < 2 * state.sites; ++k)
      multiply(u + eta / Real{2}, u - eta / Real{2}, radii[i]);
    for (std::size_t j = 0; j < m; ++j)
      if (i != j)
      {
        C const v = state.rapidities[j];
        multiply(u - v - eta, u - v + eta, radii[i] + radii[j]);
        multiply(u + v - eta, u + v + eta, radii[i] + radii[j]);
      }
    Real const residual = regular ? std::abs(a - b) : uni20::numeric_limits<Real>::infinity();
    Real const bound = std::max(Real{8192} * Real(state.sites) * eps, Real{8} * error);
    state.bethe_residual = std::max(state.bethe_residual, residual);
    state.bethe_residual_bound = std::max(state.bethe_residual_bound, bound);
    verified = verified && uni20::isfinite(error) && error < Real{1} / Real{10000} && residual <= bound;
  }
  state.converged = state.energy.has_value() && verified;
  state.status = state.converged ? Status::converged : Status::inadmissible;
}

template <uni20::Real Real>
State<Real> solve(System<Real> const& system, std::span<Real const> seed, SolverOptions<Real> const& options)
{
  if (!uni20::isfinite(options.residual_tolerance) || options.residual_tolerance <= Real{0})
    throw std::invalid_argument("open Q-system tolerance must be finite and positive");
  State<Real> state;
  state.sites = system.sites;
  state.through_lines = system.sites - 2 * system.order;
  state.delta = system.delta;
  state.coefficients.assign(seed.begin(), seed.end());
  std::vector<Real> jacobian;
  for (;;)
  {
    auto const f = system.evaluate(state.coefficients, &jacobian);
    state.residual_norm = f.norm;
    state.partner = f.partner;
    state.energy = system.energy(state.coefficients);
    if (f.norm <= options.residual_tolerance)
    {
      check_roots(state);
      return state;
    }
    if (state.iterations == options.max_iterations) return state;
    auto step = f.residual;
    for (auto& value : step)
      value = -value;
    if (!bethe::detail::newton_step(jacobian, step))
    {
      state.status = Status::singular_jacobian;
      return state;
    }
    bool accepted = false;
    Real damping{1};
    auto trial = state.coefficients;
    for (int backtrack = 0; backtrack <= uni20::numeric_limits<Real>::digits; ++backtrack)
    {
      for (std::size_t i = 0; i < trial.size(); ++i)
        trial[i] = state.coefficients[i] + damping * step[i];
      if (!std::all_of(trial.begin(), trial.end(), [](Real x) { return uni20::isfinite(x); }))
      {
        damping /= Real{2};
        continue;
      }
      try
      {
        auto const value = system.evaluate(trial);
        // Raw equations determine the descent; the backward-error normalization
        // alone can decrease spuriously as coefficients grow without bound.
        Real before{}, after{};
        for (std::size_t i = 0; i < trial.size(); ++i)
        {
          before = std::max(before, std::abs(f.residual[i]));
          after = std::max(after, std::abs(value.residual[i]));
        }
        if (after < (Real{1} - damping / Real{10000}) * before)
        {
          accepted = true;
          break;
        }
      }
      catch (std::overflow_error const&)
      {} // Backtrack a nonrepresentable trial only.
      damping /= Real{2};
    }
    if (!accepted)
    {
      state.status = Status::stalled;
      return state;
    }
    state.coefficients = std::move(trial);
    ++state.iterations;
  }
}
} // namespace detail

/// Solve one selected polynomial seed, not necessarily a lowest-energy state.
/// The same coefficients can seed continuation in Delta. No precision fallback.
template <uni20::Real Real = double>
[[nodiscard]] State<Real> solve(std::size_t sites, Real delta, std::span<Real const> seed,
                                SolverOptions<Real> const& options = {})
{
  return detail::solve(System<Real>(sites, seed.size(), delta), seed, options);
}

struct SearchOptions
{
    std::size_t max_attempts = 4000;
};

template <uni20::Real Real> struct Spectrum
{
    std::vector<State<Real>> states;
    std::size_t expected_count = 0, attempts = 0, failed_attempts = 0;
    /// Count-matched NUMERICAL completeness, not a rigorous completeness proof.
    bool complete() const { return states.size() == expected_count; }
};

/// Bounded deterministic multistart search of one small TL module, N<=8.
/// Failed attempts and duplicates consume the budget. Returned levels are
/// sorted discoveries, NOT guaranteed lowest levels unless complete().
/// No ED seeds, random-device state, hidden precision changes or energy merging.
template <uni20::Real Real = double>
[[nodiscard]] Spectrum<Real> spectrum(std::size_t sites, Real delta, std::size_t through_lines,
                                      SearchOptions const& search = {}, SolverOptions<Real> const& options = {})
{
  auto const m = quantum_group::detail::sector_roots(sites, through_lines);
  if (sites > 8) throw std::length_error("open Q-system spectrum search currently requires N <= 8");
  if (!uni20::isfinite(options.residual_tolerance) || options.residual_tolerance <= Real{0})
    throw std::invalid_argument("open Q-system tolerance must be finite and positive");
  System<Real> const system(sites, m, delta);
  Spectrum<Real> out;
  out.expected_count =
      bethe::detail::bounded_binomial(sites, m, 256) - (m ? bethe::detail::bounded_binomial(sites, m - 1, 256) : 0);
  // Deliberately conservative polynomial deduplication, independent of energy
  // ties. Close/ill-conditioned roots may leave the search incomplete.
  Real const merge = Real{32} * std::sqrt(uni20::numeric_limits<Real>::epsilon());
  std::uint64_t random = 0x243f6a8885a308d3ULL;
  auto uniform = [&] {
    random ^= random << 13;
    random ^= random >> 7;
    random ^= random << 17;
    return Real(random & 0xffffffULL) / Real{8388608} - Real{1};
  };
  std::vector<Real> seed(m);
  while (out.attempts < search.max_attempts && !out.complete())
  {
    auto const attempt = out.attempts++;
    Real const scale = std::ldexp(Real{1}, int((attempt / 32) % 7) - 1);
    for (auto& value : seed)
    {
      // Exact 24-bit rational seed generation in every supported precision.
      value = scale * uniform();
    }
    if (attempt % 2)
    {
      // Also sample polynomials with explicit real roots / conjugate pairs.
      // Coefficient boxes alone poorly cover multi-string basins.
      auto const pairs = (attempt / 2) % (m / 2 + 1);
      std::vector<Real> polynomial{Real{1}};
      auto multiply = [&](std::vector<Real> const& factor) {
        std::vector<Real> next(polynomial.size() + factor.size() - 1);
        for (std::size_t i = 0; i < polynomial.size(); ++i)
          for (std::size_t j = 0; j < factor.size(); ++j)
            next[i + j] += polynomial[i] * factor[j];
        polynomial = std::move(next);
      };
      for (std::size_t j = 0; j < pairs; ++j)
      {
        Real const center = Real{2} * uniform();
        Real const radius = std::ldexp(Real{1} + uniform() / Real{2}, int((attempt / 16 + j) % 4) - 1);
        multiply({center * center + radius * radius, -Real{2} * center, Real{1}});
      }
      for (std::size_t j = 2 * pairs; j < m; ++j)
        multiply({-uniform(), Real{1}});
      polynomial.pop_back();
      seed = std::move(polynomial);
    }
    State<Real> state;
    try
    {
      state = detail::solve<Real>(system, seed, options);
    }
    catch (std::overflow_error const&)
    {
      ++out.failed_attempts;
      continue;
    }
    if (!state.converged)
    {
      ++out.failed_attempts;
      continue;
    }
    bool duplicate = false;
    for (auto const& previous : out.states)
    {
      Real distance{};
      for (std::size_t j = 0; j < m; ++j)
        distance = std::max(
            distance, std::abs(state.coefficients[j] - previous.coefficients[j]) /
                          std::max({Real{1}, std::abs(state.coefficients[j]), std::abs(previous.coefficients[j])}));
      if (distance <= merge)
      {
        duplicate = true;
        break;
      }
    }
    if (!duplicate) out.states.push_back(std::move(state));
  }
  std::sort(out.states.begin(), out.states.end(), [](auto const& a, auto const& b) {
    if (*a.energy != *b.energy) return *a.energy < *b.energy;
    return a.coefficients < b.coefficients;
  });
  return out;
}
} // namespace bethe::xxz::quantum_group::qsystem
