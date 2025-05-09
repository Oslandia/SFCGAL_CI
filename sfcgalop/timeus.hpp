// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

/**
 * timeus.hpp - High-resolution timing utilities
 */

#ifndef SFCGALOP_TIMEUS_HPP
#define SFCGALOP_TIMEUS_HPP

#include <chrono>
#include <string>

/**
 * @brief Return the current time as microseconds since the epoch.
 *
 * Retrieves the current time from std::chrono::high_resolution_clock and
 * returns the duration since the clock's epoch expressed in microseconds as a
 * double.
 *
 * @return double Current time in microseconds since the clock epoch.
 */
inline auto
get_time_us() -> double
{
  using namespace std::chrono;
  auto now      = high_resolution_clock::now();
  auto duration = now.time_since_epoch();
  return static_cast<double>(duration_cast<microseconds>(duration).count());
}

/**
 * @brief Compute the elapsed time between two timestamps in microseconds.
 *
 * Returns end - start in microseconds. The result may be negative if `end` is
 * earlier than `start`.
 *
 * @param start Start timestamp in microseconds.
 * @param end End timestamp in microseconds.
 * @return Elapsed time in microseconds (double).
 */
inline auto
elapsed_time_us(double start, double end) -> double
{
  return end - start;
}

/**
 * @brief Convert a duration in microseconds to a human-readable string with an
 * appropriate unit.
 *
 * The function selects the unit based on magnitude:
 * - "<1000" -> microseconds ("μs")
 * - ">=1000 and <1,000,000" -> milliseconds ("ms")
 * - ">=1,000,000" -> seconds ("s")
 *
 * @param microseconds Duration expressed in microseconds.
 * @return std::string Numeric value and unit (e.g., "123.456 ms", "1.234 s",
 * "500.0 μs").
 */
inline auto
format_time_us(double microseconds) -> std::string
{
  if (microseconds < 1000.0) {
    return std::to_string(microseconds) + " μs";
  }
  if (microseconds < 1000000.0) {
    return std::to_string(microseconds / 1000.0) + " ms";
  }
  return std::to_string(microseconds / 1000000.0) + " s";
}

#endif /* SFCGALOP_TIMEUS_HPP */
