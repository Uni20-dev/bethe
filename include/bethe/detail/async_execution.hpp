// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once

#include <bethe/execution.hpp>
#include <exception>
#include <mutex>
#include <stop_token>
#include <uni20/async/async.hpp>

namespace bethe::detail
{
/// Operation-local cancellation and first execution error. The owner must join
/// every admitted operation before propagating the error or destroying this.
class ExecutionControl {
    struct ForwardStop
    {
        std::stop_source* source;
        void operator()() const noexcept { source->request_stop(); }
    };
    std::stop_source stop_;
    std::stop_callback<ForwardStop> external_;
    mutable std::mutex mutex_;
    std::exception_ptr error_;

  public:
    explicit ExecutionControl(std::stop_token external) : external_(external, ForwardStop{&stop_}) {}
    bool stopped() const noexcept { return stop_.stop_requested(); }
    std::stop_token token() const noexcept { return stop_.get_token(); }
    void fail(std::exception_ptr error = std::current_exception())
    {
      stop_.request_stop();
      std::lock_guard lock(mutex_);
      if (!error_) error_ = error;
    }
    void rethrow_error() const
    {
      std::lock_guard lock(mutex_);
      if (error_) std::rethrow_exception(error_);
    }
};

/// A rejected, still-owned task must be explicitly cancelled before destruction.
/// Callers submitting children additionally retain and join accepted children.
inline void schedule_owned(uni20::async::IAsyncScheduler& scheduler, uni20::async::AsyncTask task)
{
  try
  {
    scheduler.schedule(std::move(task));
  }
  catch (...)
  {
    if (task) task.set_cancel_on_resume();
    throw;
  }
}

template <typename T> T take_result(uni20::async::Async<T>& result, uni20::async::IAsyncScheduler& scheduler)
{
  (void)result.get_wait(scheduler);
  return result.move_from_wait();
}
} // namespace bethe::detail
