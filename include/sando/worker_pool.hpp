/* ----------------------------------------------------------------------------
 * Copyright 2026, Kota Kondo, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Kota Kondo, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#pragma once

#include <condition_variable>
#include <deque>
#include <fstream>
#include <functional>
#include <mutex>
#include <set>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace sando {

/** @brief Number of physical CPU cores (hyperthreads not counted): unique (package, core) pairs in
 *  /proc/cpuinfo, else std::thread::hardware_concurrency().
 */
inline std::size_t physicalCoreCount() {
  std::set<std::pair<int, int>> cores;
  std::ifstream f("/proc/cpuinfo");
  std::string line;
  int package = -1;
  while (std::getline(f, line)) {
    const auto colon = line.find(':');
    if (colon == std::string::npos) continue;
    try {
      if (line.rfind("physical id", 0) == 0) package = std::stoi(line.substr(colon + 1));
      if (line.rfind("core id", 0) == 0) cores.emplace(package, std::stoi(line.substr(colon + 1)));
    } catch (...) {
    }
  }
  if (!cores.empty()) return cores.size();
  const unsigned n = std::thread::hardware_concurrency();
  return n > 0 ? n : 1;
}

/** @brief Fixed set of threads running submitted tasks in FIFO order.
 *
 *  Used for the parallel local-trajectory attempts of every replan, so they reuse the same
 *  threads instead of starting new ones each time. Tasks must not throw.
 */
class WorkerPool {
 public:
  explicit WorkerPool(std::size_t num_threads) {
    if (num_threads == 0) num_threads = 1;
    threads_.reserve(num_threads);
    for (std::size_t i = 0; i < num_threads; ++i) threads_.emplace_back([this] { run(); });
  }

  ~WorkerPool() {
    {
      std::lock_guard<std::mutex> lk(mtx_);
      stop_ = true;
    }
    cv_.notify_all();
    for (auto& t : threads_) t.join();
  }

  WorkerPool(const WorkerPool&) = delete;
  WorkerPool& operator=(const WorkerPool&) = delete;

  void submit(std::function<void()> task) {
    {
      std::lock_guard<std::mutex> lk(mtx_);
      tasks_.push_back(std::move(task));
    }
    cv_.notify_one();
  }

  std::size_t size() const { return threads_.size(); }

 private:
  void run() {
    for (;;) {
      std::function<void()> task;
      {
        std::unique_lock<std::mutex> lk(mtx_);
        cv_.wait(lk, [this] { return stop_ || !tasks_.empty(); });
        if (tasks_.empty()) return;  // stopping and nothing left to run
        task = std::move(tasks_.front());
        tasks_.pop_front();
      }
      task();
    }
  }

  std::vector<std::thread> threads_;
  std::deque<std::function<void()>> tasks_;
  std::mutex mtx_;
  std::condition_variable cv_;
  bool stop_ = false;
};

}  // namespace sando
