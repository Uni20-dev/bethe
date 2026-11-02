// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/experimental/detail/spinon_exact_decimal.hpp>
#include <bethe/experimental/spinon_observation_group.hpp>
#include <bethe/experimental/spinon_recovery_store.hpp>
#include <bethe/observable_grid.hpp>

namespace bethe::experimental::xxx
{
struct SpinonObservationIdentity
{
    std::string id;
    bethe::UniformSample<5> uniforms;
};

/// Bounded experimental restart unit. Preparation is reconstructed and compared
/// byte-for-byte with the frozen geometry binding. Density controls are restored
/// by the existing driver; saved weighted observations are rederived, not trusted.
template <uni20::Real G, uni20::Real Low, uni20::Real High = Low> class SpinonObservationStore {
  public:
    using Group = SixSpinonObservationGroup<G>;
    using Prepared = typename Group::Prepared;
    using Observation = typename Group::Observation;
    using Controls = SpinonRecoveryStore<Low, High>;
    using Result = typename Controls::Result;

    SpinonObservationStore(std::filesystem::path directory, Group group, std::vector<SpinonObservationIdentity> samples,
                           heisenberg::SpinonRecoveryPolicy policy, std::string producer, std::string provenance,
                           checkpoint_detail::PublicationObserver publication = {})
        : directory_(std::move(directory)), group_(std::move(group)), publication_(std::move(publication))
    {
      // Preflight and prepare before creating or locking any filesystem object.
      heisenberg::XXXSpinonDensityRecovery<Low, High> validate(policy);
      checkpoint_detail::require(!producer.empty() && !provenance.empty(), "observation provenance is required");
      checkpoint_detail::require(!samples.empty() && samples.size() <= checkpoint_detail::max_points,
                                 "observation chunk requires 1 to 64 samples");
      checkpoint_detail::require(uni20::precision_of(group_.exchange()).bit_count() >= policy.stages.back().guard_bits,
                                 "observation geometry must retain the largest density guard precision");
      RecoveryChunkBinding binding{std::move(producer), std::move(provenance), std::move(policy), {}};
      std::set<std::string> ids;
      checkpoint_detail::Encoder e;
      e.text("bethe-six-spinon-observation-identity-v1");
      e.text(binding.producer);
      e.text(binding.provenance);
      e.text("physical longitudinal W6; labelled roots; fifteen-chart balance; 1/J; no additional factorial");
      encode_group(e);
      e.integer(samples.size());
      for (auto& sample : samples)
      {
        checkpoint_detail::require(!sample.id.empty() && ids.insert(sample.id).second,
                                   "empty or duplicate observation identity");
        prepared_.push_back(group_.prepare(std::move(sample.id), std::move(sample.uniforms)));
        auto const& p = prepared_.back();
        encode_prepared(e, p);
        roots_.emplace_back(p.sections().size());
        for (auto const& row : p.rows())
          if (row.support == SpinonSectionSupport::selected && !roots_.back()[row.section][0])
            for (std::size_t root = 0; root < 2; ++root)
            {
              checkpoint_detail::require(binding.points.size() < checkpoint_detail::max_points,
                                         "observation chunk exceeds 64 selected labelled roots; split the chunk");
              auto const index = binding.points.size();
              roots_.back()[row.section][root] = index;
              RecoveryPoint point{std::to_string(prepared_.size() - 1) + "/" + std::to_string(row.section) + "/" +
                                      std::to_string(root),
                                  {}};
              for (auto const& beta : p.sections()[row.section].geometry->beta[root])
                point.beta.push_back(checkpoint_detail::exact_decimal(beta));
              binding.points.push_back(std::move(point));
            }
      }
      // The W1 codec also binds numerical policy and native precision/range.
      // An all-zero/geometry-failure chunk has no selected root; use a dummy
      // binding only for policy encoding, without creating a density store.
      auto identity_binding = binding;
      if (identity_binding.points.empty()) identity_binding.points.push_back({"policy-only", {"-1", "1"}});
      auto control_identity = checkpoint_detail::encode_binding<Low, High>(identity_binding);
      e.integer(control_identity.size());
      e.bytes(control_identity);
      identity_ = std::move(e.data);
      records_.resize(prepared_.size());
      active_.resize(prepared_.size(), false);
      owner_ = std::make_unique<checkpoint_detail::Ownership>(directory_);
      if (auto bytes = checkpoint_detail::read_snapshot(directory_ / "snapshot.bin"))
        decode(*bytes);
      else
        save_locked();
      if (!binding.points.empty()) controls_ = std::make_unique<Controls>(directory_ / "controls", std::move(binding));
      for (std::size_t i = 0; i < size(); ++i)
        if (records_[i]) (void)replay_locked(i);
    }
    SpinonObservationStore(SpinonObservationStore const&) = delete;
    SpinonObservationStore& operator=(SpinonObservationStore const&) = delete;
    std::size_t size() const { return prepared_.size(); }
    Prepared const& prepared(std::size_t index) const { return prepared_.at(index); }

    /// Absent means this observation has not been atomically committed. A
    /// failed committed observation remains present with absent values.
    std::optional<Observation> load(std::size_t index)
    {
      std::lock_guard lock(mutex_);
      healthy();
      checkpoint_detail::require(!active_.at(index), "observation is active");
      if (!records_.at(index)) return {};
      return replay_locked(index);
    }

    Observation run(std::size_t index, std::stop_token stop = {})
    {
      return run(index, [](std::size_t, std::size_t, Result const&) {}, stop);
    }

    /// The callback runs after each durable density control, outside the store
    /// mutex. It may throw to interrupt work; its exception is never a zero.
    template <class Progress> Observation run(std::size_t index, Progress&& progress, std::stop_token stop = {})
    {
      {
        std::lock_guard lock(mutex_);
        healthy();
        checkpoint_detail::require(!active_.at(index), "observation already active");
        if (stop.stop_requested()) return group_.observe(prepared_.at(index), unreachable(), stop);
        if (records_.at(index)) return replay_locked(index);
        active_[index] = true;
      }
      struct Release
      {
          SpinonObservationStore& store;
          std::size_t index;
          ~Release()
          {
            std::lock_guard lock(store.mutex_);
            store.active_[index] = false;
          }
      } release{*this, index};
      auto result = group_.observe(
          prepared_[index],
          [&](std::size_t section, std::size_t root, auto, G const& weight, std::stop_token token) {
            auto recovered = controls_->run(
                *roots_[index][section][root],
                [&](Result const& value) { std::invoke(progress, section, root, value); }, token);
            return heisenberg::detail::weight_recovered_spinon_observation(std::span<G const>(&weight, 1), recovered);
          },
          stop);
      if (!result.cancelled)
      {
        std::lock_guard lock(mutex_);
        healthy();
        records_[index] = encode_observation(result);
        save_locked();
      }
      return result;
    }

  private:
    using RootMap = std::vector<std::array<std::optional<std::size_t>, 2>>;
    std::filesystem::path directory_;
    Group group_;
    std::vector<Prepared> prepared_;
    std::vector<RootMap> roots_;
    std::string identity_;
    checkpoint_detail::PublicationObserver publication_;
    std::unique_ptr<checkpoint_detail::Ownership> owner_;
    std::unique_ptr<Controls> controls_;
    std::vector<std::optional<std::string>> records_;
    std::vector<bool> active_;
    std::mutex mutex_;
    bool poisoned_ = false;

    static auto unreachable()
    {
      return [](auto, auto, auto, auto, auto) -> Observation {
        throw std::logic_error("cancelled observation entered a density kernel");
      };
    }
    void healthy() const
    {
      if (poisoned_) throw std::runtime_error("observation publication failed; close and reopen the chunk");
    }
    void encode_group(checkpoint_detail::Encoder& e) const
    {
      e.real(group_.exchange());
      e.real(group_.proposal().cutoff());
      e.real(group_.proposal().soft_fraction());
      e.integer(group_.proposal().components().size());
      for (auto const& c : group_.proposal().components())
      {
        e.integer(static_cast<int>(c.kind));
        e.real(c.weight);
        e.real(c.scale);
        e.real(c.gap_scale);
      }
      e.integer(group_.requests().size());
      for (auto const& r : group_.requests())
      {
        e.text(r.id);
        e.real(r.momentum);
        e.real(r.omega);
        e.real(r.inner_cutoff);
        e.real(r.outer_cutoff);
        e.integer(r.section);
      }
    }
    static void encode_prepared(checkpoint_detail::Encoder& e, Prepared const& p)
    {
      e.text(p.id());
      e.integer(p.identity().bits());
      for (auto const& coordinate : p.identity().words())
      {
        e.integer(coordinate.size());
        for (auto word : coordinate)
          e.integer(word);
      }
      e.integer(p.draw().has_value());
      if (p.draw())
      {
        for (auto const& b : p.draw()->beta)
          e.real(b);
        e.real(p.draw()->log_density);
        e.integer(p.draw()->component);
        e.integer(p.draw()->pairing);
      }
      e.integer(p.sections().size());
      for (auto const& section : p.sections())
      {
        e.text(section.reason);
        e.integer(section.geometry.has_value());
        if (section.geometry)
        {
          auto const& g = *section.geometry;
          e.integer(g.inside);
          e.integer(g.working_bits);
          e.integer(g.refinements);
          e.real(g.log_coarea_jacobian);
          for (auto const& root : g.beta)
          {
            e.integer(root.size());
            for (auto const& b : root)
              e.real(b);
          }
        }
        e.integer(section.log_weights.has_value());
        if (section.log_weights)
          for (auto const& w : *section.log_weights)
            e.real(w);
      }
      e.integer(p.rows().size());
      for (auto const& row : p.rows())
      {
        e.integer(row.section);
        e.integer(static_cast<int>(row.support));
        e.text(row.reason);
      }
    }
    static std::string encode_observation(Observation const& result)
    {
      checkpoint_detail::Encoder e;
      e.integer(result.values.has_value());
      if (result.values)
      {
        e.integer(result.values->size());
        for (auto const& value : *result.values)
          e.real(value);
      }
      e.integer(result.cancelled);
      e.text(result.reason);
      e.integer(result.metadata.supported_configurations);
      e.integer(result.metadata.density_attempts);
      e.integer(static_cast<int>(result.metadata.failure));
      e.integer(result.metadata.density_status.has_value());
      if (result.metadata.density_status) e.integer(static_cast<int>(*result.metadata.density_status));
      e.integer(result.metadata.density_failure.has_value());
      if (result.metadata.density_failure) e.integer(static_cast<int>(*result.metadata.density_failure));
      return std::move(e.data);
    }
    Observation replay_locked(std::size_t index)
    {
      auto value =
          group_.observe(prepared_[index], [&](std::size_t section, std::size_t root, auto, G const& weight, auto) {
            checkpoint_detail::require(controls_ && roots_[index][section][root],
                                       "committed observation lacks root identity");
            auto recovered = controls_->load(*roots_[index][section][root]);
            checkpoint_detail::require(recovered.has_value(), "committed observation lacks density history");
            return heisenberg::detail::weight_recovered_spinon_observation(std::span<G const>(&weight, 1), *recovered);
          });
      checkpoint_detail::require(encode_observation(value) == *records_[index],
                                 "committed observation disagrees with its validated density histories");
      return value;
    }
    void save_locked()
    {
      try
      {
        checkpoint_detail::Encoder e;
        e.text("bethe-six-spinon-observation-chunk-v1");
        e.integer(identity_.size());
        e.bytes(identity_);
        e.integer(records_.size());
        for (auto const& record : records_)
        {
          e.integer(record.has_value());
          if (record)
          {
            e.integer(record->size());
            e.bytes(*record);
          }
        }
        checkpoint_detail::publish(directory_, e.finish(), publication_);
      }
      catch (...)
      {
        poisoned_ = true;
        throw;
      }
    }
    void decode(std::string_view data)
    {
      using checkpoint_detail::require;
      require(data.size() >= 8, "truncated observation snapshot");
      checkpoint_detail::Decoder trailer{data.substr(data.size() - 8)};
      data.remove_suffix(8);
      require(checkpoint_detail::checksum(data) == trailer.integer(), "observation checksum mismatch");
      checkpoint_detail::Decoder d{data};
      require(d.text() == "bethe-six-spinon-observation-chunk-v1", "unknown observation snapshot version");
      require(d.bytes(d.count(checkpoint_detail::max_bytes)) == identity_,
              "observation producer, policy, proposal, requests, samples or prepared geometry changed");
      require(d.count(checkpoint_detail::max_points) == size(), "observation count changed");
      for (auto& record : records_)
        if (d.flag()) record = std::string(d.bytes(d.count(checkpoint_detail::max_bytes)));
      require(d.at == d.data.size(), "trailing observation payload");
    }
};

/// Own the chunk through every admitted task, including exception/cancellation
/// drainage. Collection is ordered and failed observations remain explicit.
/// The existing scheduler/window bounds admission; no extra pool or batch join.
template <uni20::Real G, uni20::Real Low, uni20::Real High, class Collect>
[[nodiscard]] uni20::async::Async<std::size_t>
evaluate_spinon_observation_chunk_async(std::shared_ptr<SpinonObservationStore<G, Low, High>> store, Collect collect,
                                        bethe::ExecutionOptions execution = {})
{
  if (!store) throw std::invalid_argument("observation chunk is null");
  auto const count = store->size();
  return bethe::evaluate_observable_grid_async(
      count, [store = std::move(store), stop = execution.stop](std::size_t i) { return store->run(i, stop); },
      std::move(collect), execution);
}
} // namespace bethe::experimental::xxx
