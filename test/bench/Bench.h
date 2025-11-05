// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_BENCH_H_
#define SFCGAL_BENCH_H_

#include <iostream>
#include <stack>

#include <chrono>
#include <iomanip>
#include <sstream>

#include <boost/format.hpp>

namespace SFCGAL {

/**
 * @brief Simple timer
 */
class CpuTimer {

public:
  CpuTimer() = default;

  /**
   * @brief Start the timer
   */
  void
  start()
  {
    _start_time = std::chrono::high_resolution_clock::now();
    _is_running = true;
  }

  /**
   * @brief Stop the timer
   */
  void
  stop()
  {
    _stop_time  = std::chrono::high_resolution_clock::now();
    _is_running = false;
  }

  /**
   * @brief Returns elapsed time in nanoseconds
   *
   * @return Elapsed time in nanoseconds since start() was called
   */
  [[nodiscard]] auto
  elapsed() const -> long long
  {
    auto end =
        _is_running ? std::chrono::high_resolution_clock::now() : _stop_time;
    return std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                _start_time)
        .count();
  }

private:
  std::chrono::high_resolution_clock::time_point _start_time;
  std::chrono::high_resolution_clock::time_point _stop_time;
  bool                                           _is_running = false;
};

/**
 * @brief helper class to write formated benchs
 */
class Bench {
public:
  /**
   * @brief Type alias for the timer used in benchmarks
   */
  using timer_t = CpuTimer;

  /**
   * destructor
   */
  ~Bench();

  /**
   * @brief Start a bench
   * @param description Description of the benchmark
   */
  void
  start(const std::string &description);
  /**
   * @brief Start a bench
   * @param description Description of the benchmark as formatted string
   */
  void
  start(const boost::basic_format<char> &description);
  /**
   * stop a bench
   */
  void
  stop();

  /**
   * @brief Get bench instance
   * @return Reference to the singleton bench instance
   */
  static Bench &
  instance();

  /**
   * @brief Get output stream
   * @return Reference to the output stream
   */
  std::ostream &
  s();

private:
  /**
   * default constructor
   */
  Bench();
  /**
   * copy constructor
   */
  Bench(const Bench &bench);

  /**
   * output stream to write bench result (default is std::cout)
   */
  std::ostream *_s;
  /**
   * timer stack with description
   */
  std::stack<std::pair<std::string, timer_t>> _timers;
};

/**
 * @brief Get bench instance
 * @return Reference to the singleton bench instance
 */
inline Bench &
bench()
{
  return Bench::instance();
}

} // namespace SFCGAL

#endif
