// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/detail/xxx_spinon_proposal.hpp>
#include <bethe/detail/xxx_spinon_surface.hpp>
#include <bethe/xxx_spinon_observation.hpp>
#include <memory>
#include <set>

namespace bethe::experimental::xxx
{
template <uni20::Real R> struct SpinonSectionRequest
{
    std::string id;
    R momentum, omega, inner_cutoff, outer_cutoff;
    int section = 0;
};

enum class SpinonSectionSupport
{
  selected,
  support_zero,
  masked_zero,
  geometry_failure
};

/// A bounded set of fixed-momentum/frequency/section observables sharing one
/// proposal draw. Masks do not modify the complete proposal or its denominator.
/// This numerical adapter does not itself persist observations or admit tasks.
template <uni20::Real G> class SixSpinonObservationGroup {
  public:
    using Proposal = heisenberg::detail::SpinonJointProposal<G>;
    using Surface = heisenberg::detail::SpinonEnergySurface<G>;
    using Request = SpinonSectionRequest<G>;
    using Observation = bethe::MonteCarloObservation<G, heisenberg::SpinonSampleMetadata>;
    struct Section
    {
        std::optional<typename Surface::Point> geometry;
        std::optional<std::array<G, 2>> log_weights;
        std::string reason;
    };
    struct Row
    {
        std::size_t section;
        SpinonSectionSupport support;
        std::string reason;
    };

  private:
    struct Kernel
    {
        Proposal proposal;
        G exchange, zero, one, pi, log_exchange;
        std::vector<Request> requests;
        std::vector<Surface> surfaces;
        std::vector<bool> boundary_zero;
        std::vector<std::size_t> sections;

        Kernel(Proposal p, G j, std::vector<Request> r)
            : proposal(std::move(p)), exchange(std::move(j)), zero(uni20::scalar_like(proposal.cutoff(), 0)),
              one(uni20::scalar_like(zero, 1)), pi(4 * uni20::math::atan(one)), requests(std::move(r))
        {
          bethe::detail::require_same_precision(exchange, one);
          if (!uni20::isfinite(exchange) || exchange <= zero)
            throw std::invalid_argument("observation group requires finite J>0");
          log_exchange = uni20::math::log(exchange);
          if (requests.empty() || requests.size() > 32)
            throw std::invalid_argument("observation group requires 1 to 32 requests");
          std::set<std::string> ids;
          std::vector<std::size_t> first;
          for (std::size_t i = 0; i < requests.size(); ++i)
          {
            auto const& request = requests[i];
            for (auto const* x : {&request.momentum, &request.omega, &request.inner_cutoff, &request.outer_cutoff})
            {
              bethe::detail::require_same_precision(*x, one);
              if (!uni20::isfinite(*x)) throw std::invalid_argument("nonfinite observation request");
            }
            if (request.id.empty() || !ids.insert(request.id).second || request.momentum <= zero ||
                request.momentum >= 2 * pi || request.omega <= zero || request.inner_cutoff < zero ||
                request.inner_cutoff >= request.outer_cutoff || request.outer_cutoff > proposal.cutoff() ||
                request.section < 0 || request.section >= 3)
              throw std::invalid_argument(
                  "invalid observation identity, interior momentum, frequency, mask or section");
            G const folded = std::min(request.momentum, G(2 * pi - request.momentum));
            G const lower = folded == pi ? zero : G(exchange * pi / 2 * uni20::math::sin(folded));
            // Concavity fixes each section's maximum at six equal momenta.
            // Another section may extend higher, but cannot give this request
            // support above its own edge or remove its vanishing edge value.
            G const upper = exchange * 3 * pi * uni20::math::sin((request.momentum + 2 * request.section * pi) / 6);
            if (!uni20::isfinite(lower) || !uni20::isfinite(upper))
              throw std::overflow_error("observation support exceeds exponent range");
            if (folded != pi && lower <= zero) throw std::underflow_error("observation lower boundary lost scale");
            if (request.omega == lower)
              throw std::invalid_argument("lower spectral threshold requires a limiting calculation");
            // Equal geometric requests share their solved roots and density,
            // even when their masks and output identities differ.
            auto found = std::find_if(first.begin(), first.end(), [&](std::size_t k) {
              auto const& prior = requests[k];
              return request.momentum == prior.momentum && request.omega == prior.omega &&
                     request.section == prior.section;
            });
            if (found != first.end())
            {
              sections.push_back(std::size_t(found - first.begin()));
              continue;
            }
            sections.push_back(first.size());
            first.push_back(i);
            G const normalized = request.omega / exchange;
            if (!uni20::isfinite(normalized) || normalized <= zero)
              throw std::underflow_error("observation frequency lost exchange scale");
            surfaces.emplace_back(6, proposal.cutoff(), request.momentum, request.section, normalized);
            boundary_zero.push_back(request.omega < lower || request.omega >= upper);
          }
        }
    };
    std::shared_ptr<Kernel const> kernel_;

  public:
    class Prepared {
        friend class SixSpinonObservationGroup;
        struct Data
        {
            std::shared_ptr<Kernel const> kernel;
            std::string id;
            bethe::UniformSample<5> identity;
            std::optional<typename Proposal::Sample> draw;
            std::vector<Section> sections;
            std::vector<Row> rows;
            Data(std::shared_ptr<Kernel const> k, std::string name, bethe::UniformSample<5> sample)
                : kernel(std::move(k)), id(std::move(name)), identity(std::move(sample))
            {}
        };
        std::shared_ptr<Data const> data_;
        explicit Prepared(std::shared_ptr<Data const> data) : data_(std::move(data)) {}
        Data const& data() const
        {
          if (!data_) throw std::logic_error("moved observation preparation");
          return *data_;
        }

      public:
        std::string const& id() const { return data().id; }
        bethe::UniformSample<5> const& identity() const { return data().identity; }
        std::optional<typename Proposal::Sample> const& draw() const { return data().draw; }
        std::span<Section const> sections() const { return data().sections; }
        std::span<Row const> rows() const { return data().rows; }
    };

    SixSpinonObservationGroup(Proposal proposal, G exchange, std::vector<Request> requests)
        : kernel_(std::make_shared<Kernel>(std::move(proposal), std::move(exchange), std::move(requests)))
    {}
    std::span<Request const> requests() const { return kernel().requests; }
    Proposal const& proposal() const { return kernel().proposal; }
    G const& exchange() const { return kernel().exchange; }

    Prepared prepare(std::string id, bethe::UniformSample<5> identity) const
    {
      auto const& k = kernel();
      if (id.empty()) throw std::invalid_argument("observation identity is empty");
      auto data = std::make_shared<typename Prepared::Data>(kernel_, std::move(id), std::move(identity));
      // Invalid bit counts/precision are request errors, not random zero samples.
      auto const uniforms = data->identity.at_precision(k.one);
      data->sections.resize(k.surfaces.size());
      try
      {
        data->draw = k.proposal.draw(uniforms);
      }
      catch (std::runtime_error const& error)
      {
        for (auto& section : data->sections)
          section.reason = error.what();
      }
      for (std::size_t s = 0; s < k.surfaces.size(); ++s)
      {
        auto& section = data->sections[s];
        if (k.boundary_zero[s])
        {
          section.geometry = typename Surface::Point{false, {}, k.zero, uni20::numeric_limits<G>::digits(k.one)};
          section.reason.clear();
          continue;
        }
        if (!data->draw) continue;
        try
        {
          section.geometry = k.surfaces[s].solve(data->draw->beta);
          if (!section.geometry->inside) continue;
          std::array<G, 2> weights;
          for (std::size_t root = 0; root < 2; ++root)
          {
            std::array<G, 6> beta;
            std::copy(section.geometry->beta[root].begin(), section.geometry->beta[root].end(), beta.begin());
            weights[root] = k.proposal.balanced_log_jacobian(beta) - k.log_exchange;
            if (!uni20::isfinite(weights[root])) throw std::runtime_error("nonfinite balanced observation weight");
          }
          section.log_weights = std::move(weights);
        }
        catch (std::runtime_error const& error)
        {
          section.reason = error.what();
        }
      }
      for (std::size_t i = 0; i < k.requests.size(); ++i)
      {
        auto const s = k.sections[i];
        auto const& section = data->sections[s];
        Row row{s, SpinonSectionSupport::geometry_failure, section.reason};
        if (section.geometry && section.reason.empty())
        {
          if (!section.geometry->inside)
            row.support = SpinonSectionSupport::support_zero;
          else
          {
            G depth = k.zero;
            for (auto const& beta : section.geometry->beta[0])
              depth = std::max(depth, uni20::math::abs(beta));
            auto const& request = k.requests[i];
            G const tolerance = uni20::math::sqrt(uni20::numeric_limits<G>::epsilon(k.one)) * std::max(k.one, depth);
            bool const near_inner =
                request.inner_cutoff > k.zero && uni20::math::abs(depth - request.inner_cutoff) <= tolerance;
            bool const near_outer = request.outer_cutoff < k.proposal.cutoff() &&
                                    uni20::math::abs(depth - request.outer_cutoff) <= tolerance;
            if (near_inner || near_outer)
              row.reason = "annular mask boundary unresolved at geometry precision";
            else
              row.support = request.inner_cutoff <= depth && depth < request.outer_cutoff
                                ? SpinonSectionSupport::selected
                                : SpinonSectionSupport::masked_zero;
          }
        }
        data->rows.push_back(std::move(row));
      }
      return Prepared(std::move(data));
    }

    /// Root callback seam for bound persistence. Called once per selected
    /// geometric root, irrespective of how many masks include that root.
    /// Returns one physical value and its two diagnostics, never a survivor sum.
    template <class Evaluate>
    Observation observe(Prepared const& prepared, Evaluate&& evaluate, std::stop_token stop = {}) const
    {
      auto const& k = kernel();
      if (prepared.data().kernel != kernel_) throw std::invalid_argument("observation belongs to another group");
      Observation result;
      std::vector<std::optional<std::array<G, 3>>> values(prepared.sections().size());
      for (auto const& row : prepared.rows())
      {
        if (stop.stop_requested())
        {
          result.cancelled = true;
          result.metadata.failure = heisenberg::SpinonObservationFailure::cancelled;
          result.reason = "grouped observation cancelled";
          return result;
        }
        if (row.support == SpinonSectionSupport::geometry_failure)
        {
          result.metadata.failure = heisenberg::SpinonObservationFailure::geometry;
          result.reason = row.reason;
          return result;
        }
        if (row.support != SpinonSectionSupport::selected || values[row.section]) continue;
        auto const& section = prepared.sections()[row.section];
        std::array<bethe::detail::SpinonSum<G>, 3> sum{bethe::detail::SpinonSum<G>(k.zero),
                                                       bethe::detail::SpinonSum<G>(k.zero),
                                                       bethe::detail::SpinonSum<G>(k.zero)};
        for (std::size_t root = 0; root < 2; ++root)
        {
          auto value = evaluate(row.section, root, std::span<G const>(section.geometry->beta[root]),
                                (*section.log_weights)[root], stop);
          result.metadata += value.metadata;
          if (stop.stop_requested() || value.cancelled)
          {
            result.cancelled = true;
            result.metadata.failure = heisenberg::SpinonObservationFailure::cancelled;
            result.reason = "grouped observation cancelled";
            return result;
          }
          if (!value.values || value.cancelled)
          {
            result.reason = value.reason;
            result.cancelled = value.cancelled;
            return result;
          }
          if (value.values->size() != 3) throw std::invalid_argument("root observation must have three fields");
          for (std::size_t field = 0; field < 3; ++field)
          {
            auto const& v = (*value.values)[field];
            bethe::detail::require_same_precision(v, k.one);
            if (!uni20::isfinite(v) || v < k.zero) throw std::invalid_argument("invalid nonnegative root observation");
            sum[field].add(v);
          }
        }
        values[row.section] = std::array{sum[0].value(), sum[1].value(), sum[2].value()};
      }
      std::vector<G> output(3 * prepared.rows().size(), k.zero);
      for (std::size_t i = 0; i < prepared.rows().size(); ++i)
        if (auto const& row = prepared.rows()[i]; row.support == SpinonSectionSupport::selected)
          for (std::size_t field = 0; field < 3; ++field)
          {
            auto const& v = (*values[row.section])[field];
            if (!uni20::isfinite(v))
            {
              result.metadata.failure = heisenberg::SpinonObservationFailure::exponent_range;
              result.reason = "grouped observation sum exceeded exponent range";
              return result;
            }
            output[field * prepared.rows().size() + i] = v;
          }
      result.values = std::move(output);
      return result;
    }

    template <uni20::Real Low, uni20::Real High>
    Observation observation(Prepared const& prepared, heisenberg::XXXSpinonDensityRecovery<Low, High> const& recovery,
                            std::stop_token stop = {}) const
    {
      if (uni20::numeric_limits<G>::digits(kernel().one) < recovery.policy().stages.back().guard_bits)
        throw std::invalid_argument("observation geometry must retain the largest density guard precision");
      return observe(
          prepared,
          [&](auto, auto, std::span<G const> beta, G const& weight, std::stop_token token) {
            return heisenberg::weighted_spinon_observation<G>(beta, std::span<G const>(&weight, 1), recovery, token);
          },
          stop);
    }

  private:
    Kernel const& kernel() const
    {
      if (!kernel_) throw std::logic_error("moved observation group");
      return *kernel_;
    }
};
} // namespace bethe::experimental::xxx
