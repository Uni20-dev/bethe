// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Ian McCulloch
#pragma once

#include <bethe/experimental/detail/spinon_checkpoint_codec.hpp>
#include <filesystem>
#include <functional>
#include <mutex>
#include <system_error>

#if !defined(__unix__)
#error "Experimental recovery storage currently requires POSIX file locking and atomic rename."
#endif
#include <fcntl.h>
#include <sys/file.h>
#include <sys/stat.h>
#include <unistd.h>

namespace bethe::experimental::xxx
{
namespace checkpoint_detail
{
inline void system_failure(char const* operation)
{
  throw std::system_error(errno, std::generic_category(), operation);
}

class File {
  public:
    explicit File(int fd = -1) : fd_(fd) {}
    ~File()
    {
      if (fd_ >= 0) ::close(fd_);
    }
    File(File const&) = delete;
    File& operator=(File const&) = delete;
    int get() const { return fd_; }
    void close_checked()
    {
      int fd = std::exchange(fd_, -1);
      if (::close(fd) != 0) system_failure("closing checkpoint");
    }

  private:
    int fd_;
};

class Ownership {
  public:
    explicit Ownership(std::filesystem::path const& directory) : lock_(open_lock(directory))
    {
      if (::flock(lock_.get(), LOCK_EX | LOCK_NB) != 0)
        system_failure("recovery chunk already owned or lock unavailable");
    }

  private:
    static int open_lock(std::filesystem::path const& path)
    {
      std::filesystem::create_directories(path);
      int fd = ::open((path / "owner.lock").c_str(), O_CREAT | O_RDWR | O_CLOEXEC | O_NOFOLLOW, 0600);
      if (fd < 0) system_failure("opening recovery ownership lock");
      return fd;
    }
    File lock_; // Never unlink the lock inode, including after process death.
};

enum class PublicationPhase
{
  before_replace,
  after_replace
};
/// Fault-injection/observation seam; runs under the store mutex, never reenter.
using PublicationObserver = std::function<void(PublicationPhase)>;

inline std::optional<std::string> read_snapshot(std::filesystem::path const& path)
{
  int raw = ::open(path.c_str(), O_RDONLY | O_CLOEXEC | O_NOFOLLOW);
  if (raw < 0)
  {
    if (errno == ENOENT) return {};
    system_failure("opening recovery snapshot");
  }
  File fd(raw);
  struct stat info
  {};
  if (::fstat(raw, &info) != 0) system_failure("reading recovery snapshot size");
  require(S_ISREG(info.st_mode) && info.st_size >= 0 && static_cast<std::uint64_t>(info.st_size) <= max_bytes,
          "invalid recovery snapshot type or size");
  std::string data(static_cast<std::size_t>(info.st_size), '\0');
  std::size_t at = 0;
  while (at < data.size())
  {
    auto n = ::read(raw, data.data() + at, data.size() - at);
    if (n < 0 && errno == EINTR) continue;
    if (n < 0) system_failure("reading recovery snapshot");
    require(n != 0, "truncated recovery snapshot");
    at += static_cast<std::size_t>(n);
  }
  return data;
}

inline void publish(std::filesystem::path const& directory, std::string const& data,
                    PublicationObserver const& observer)
{
  auto const temporary = directory / "snapshot.pending", final = directory / "snapshot.bin";
  File fd(::open(temporary.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC | O_NOFOLLOW, 0600));
  if (fd.get() < 0) system_failure("opening pending recovery snapshot");
  std::size_t at = 0;
  while (at < data.size())
  {
    auto n = ::write(fd.get(), data.data() + at, data.size() - at);
    if (n < 0 && errno == EINTR) continue;
    if (n < 0) system_failure("writing recovery snapshot");
    if (n == 0) throw std::runtime_error("checkpoint write made no progress");
    at += static_cast<std::size_t>(n);
  }
  fd.close_checked();
  if (observer) observer(PublicationPhase::before_replace);
  if (::rename(temporary.c_str(), final.c_str()) != 0) system_failure("replacing recovery snapshot");
  if (observer) observer(PublicationPhase::after_replace);
  // Rename gives process-restart atomicity, not fsync-based power-loss durability.
}

inline double epoch_seconds()
{
  return std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count();
}
} // namespace checkpoint_detail

/// Bounded POSIX single-writer chunk store for the existing density recovery.
/// Concurrent distinct-point calls are allowed; publication is serialized.
/// The store must outlive its callers. It owns no scheduler or cluster policy.
template <uni20::Real Low, uni20::Real High = Low> class SpinonRecoveryStore {
  public:
    using Driver = heisenberg::XXXSpinonDensityRecovery<Low, High>;
    using Result = typename Driver::Result;

    explicit SpinonRecoveryStore(std::filesystem::path directory, RecoveryChunkBinding binding,
                                 checkpoint_detail::PublicationObserver publication = {})
        : directory_(std::move(directory)), binding_(std::move(binding)),
          identity_(checkpoint_detail::encode_binding<Low, High>(binding_)), driver_(binding_.policy),
          owner_(directory_), publication_(std::move(publication)), entries_(binding_.points.size()),
          active_(binding_.points.size(), false)
    {
      if (auto data = checkpoint_detail::read_snapshot(directory_ / "snapshot.bin"))
        decode(*data);
      else
        save_locked();
    }
    SpinonRecoveryStore(SpinonRecoveryStore const&) = delete;
    SpinonRecoveryStore& operator=(SpinonRecoveryStore const&) = delete;
    std::size_t size() const { return binding_.points.size(); }
    RecoveryPoint const& point(std::size_t index) const { return binding_.points.at(index); }

    /// Revalidate a snapshot. A running point's time includes downtime. A
    /// finalized result is reusable without a new numerical deadline.
    std::optional<Result> load(std::size_t index)
    {
      std::lock_guard lock(mutex_);
      healthy();
      auto const& e = entries_.at(index);
      if (!e.started) return {};
      return reassess(e, e.finished ? last_epoch_ : now_locked());
    }

    Result run(std::size_t index, std::stop_token stop = {})
    {
      return run(index, [](Result const&) {}, stop);
    }

    /// Observer runs after each atomically published control, outside the mutex.
    /// Copy any retained view. Its exception propagates without becoming a
    /// numerical failure. The latest published prefix remains restartable.
    template <class Observer> Result run(std::size_t index, Observer&& observer, std::stop_token stop = {})
    {
      Result prior;
      {
        std::lock_guard lock(mutex_);
        healthy();
        auto& e = entries_.at(index);
        checkpoint_detail::require(!active_[index], "recovery point already active in this store");
        if (e.finished) return driver_.restore(e.result);
        double const now = now_locked();
        if (!e.started)
        {
          e.started = now;
          e.deadline = now + binding_.policy.max_wall_seconds;
          checkpoint_detail::require(std::isfinite(e.deadline), "recovery deadline overflow");
          e.result.policy = binding_.policy;
          e.result.decimal_inputs = binding_.points[index].beta;
          save_locked(); // Bind the original deadline before admitting kernels.
        }
        prior = reassess(e, now_locked());
        active_[index] = true;
      }
      struct Release
      {
          SpinonRecoveryStore& store;
          std::size_t index;
          ~Release()
          {
            std::lock_guard lock(store.mutex_);
            store.active_[index] = false;
          }
      } release{*this, index};
      auto result = driver_.resume(
          binding_.points[index].beta, prior,
          [&](Result const& progress) {
            {
              std::lock_guard lock(mutex_);
              healthy();
              auto& e = entries_[index];
              e.result = progress;
              e.result.wall_seconds = std::max(e.result.wall_seconds, now_locked() - e.started);
              save_locked();
            }
            observer(progress);
          },
          stop);
      {
        std::lock_guard lock(mutex_);
        healthy();
        auto& e = entries_[index];
        bool const cancelled = result.status == heisenberg::SpinonRecoveryStatus::cancelled;
        e.result = result;
        e.result.wall_seconds = std::max(e.result.wall_seconds, now_locked() - e.started);
        e.result = driver_.restore(e.result);
        e.finished = terminal(e.result);
        save_locked();
        result = e.result;
        // Preserve the caller's cancellation without trusting/storing a top-level
        // flag. A later resume reconstructs the prefix's numerical status.
        if (cancelled && !e.finished)
        {
          result.status = heisenberg::SpinonRecoveryStatus::cancelled;
          result.reason = "recovery cancelled; completed controls saved";
        }
      }
      return result;
    }

  private:
    struct Entry
    {
        double started = 0, deadline = 0;
        bool finished = false;
        Result result;
    };
    std::filesystem::path directory_;
    RecoveryChunkBinding binding_;
    std::string identity_;
    Driver driver_;
    checkpoint_detail::Ownership owner_;
    checkpoint_detail::PublicationObserver publication_;
    std::vector<Entry> entries_;
    std::vector<bool> active_;
    std::mutex mutex_;
    bool poisoned_ = false;
    double last_epoch_ = 0;

    void healthy() const
    {
      if (poisoned_) throw std::runtime_error("recovery publication failed; close and reopen the chunk");
    }
    double now_locked()
    {
      double const now = checkpoint_detail::epoch_seconds();
      checkpoint_detail::require(std::isfinite(now) && now > 0 && now >= last_epoch_,
                                 "checkpoint clock moved backwards; refusing to reset time");
      last_epoch_ = now;
      return now;
    }
    bool terminal(Result const& r) const
    {
      using S = heisenberg::SpinonRecoveryStatus;
      return r.accepted() || r.status == S::inconsistent_controls || r.status == S::time_limit ||
             r.attempts.size() == 2 * binding_.policy.stages.size();
    }
    Result reassess(Entry const& e, double now) const
    {
      auto saved = e.result;
      if (!e.finished) saved.wall_seconds = std::max(saved.wall_seconds, now - e.started);
      return driver_.restore(saved);
    }
    void save_locked()
    {
      try
      {
        checkpoint_detail::Encoder e;
        e.text("bethe-spinon-recovery-chunk-v1");
        e.integer(identity_.size());
        e.bytes(identity_);
        e.clock(now_locked());
        e.integer(entries_.size());
        for (std::size_t i = 0; i < entries_.size(); ++i)
        {
          auto const& point = entries_[i];
          e.integer(i);
          e.clock(point.started);
          e.clock(point.deadline);
          e.integer(point.finished);
          if (!point.started) continue;
          e.clock(point.result.wall_seconds);
          e.integer(point.result.attempts.size());
          for (auto const& attempt : point.result.attempts)
            checkpoint_detail::encode_attempt(e, attempt);
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
      require(data.size() >= 8, "truncated recovery checksum");
      checkpoint_detail::Decoder trailer{data.substr(data.size() - 8)};
      data.remove_suffix(8);
      require(checkpoint_detail::checksum(data) == trailer.integer(), "recovery checksum mismatch");
      checkpoint_detail::Decoder d{data};
      require(d.text() == "bethe-spinon-recovery-chunk-v1", "unknown recovery snapshot version");
      require(d.bytes(d.count(checkpoint_detail::max_bytes)) == identity_,
              "recovery producer, policy or requests changed");
      last_epoch_ = d.clock();
      require(last_epoch_ > 0, "invalid recovery publication time");
      require(d.count(checkpoint_detail::max_points) == entries_.size(), "recovery point count changed");
      for (std::size_t i = 0; i < entries_.size(); ++i)
      {
        require(d.integer() == i, "recovery point ordering changed");
        auto& e = entries_[i];
        e.started = d.clock();
        e.deadline = d.clock();
        e.finished = d.flag();
        if (!e.started)
        {
          require(e.deadline == 0 && !e.finished, "unstarted recovery point has a deadline/result");
          continue;
        }
        require(e.started <= last_epoch_ && e.deadline == e.started + binding_.policy.max_wall_seconds,
                "original recovery deadline changed");
        e.result.policy = binding_.policy;
        e.result.decimal_inputs = binding_.points[i].beta;
        e.result.wall_seconds = d.clock();
        auto count = d.count(2 * binding_.policy.stages.size());
        e.result.attempts.reserve(count);
        for (std::size_t j = 0; j < count; ++j)
          e.result.attempts.push_back(checkpoint_detail::decode_attempt<High>(d, binding_.policy));
        e.result = driver_.restore(e.result);
        require(!e.finished || terminal(e.result), "finalized recovery point lacks terminal evidence");
      }
      require(d.at == d.data.size(), "trailing recovery payload");
      (void)now_locked();
    }
};
} // namespace bethe::experimental::xxx
