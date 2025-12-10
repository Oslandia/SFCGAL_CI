// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_NUMERIC_H_
#define SFCGAL_NUMERIC_H_

#include <cmath>
#include <limits>

#include "SFCGAL/export.h"

#include "SFCGAL/Kernel.h"

namespace SFCGAL {

/// @brief Default epsilon value for floating point comparisons
constexpr double EPSILON = 1e-8;

/// @brief Squared epsilon value for squared distance comparisons
constexpr double EPSILON_SQ = 1e-16;

#ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wfloat-equal"
#elif defined(__GCC__)
  #pragma gcc diagnostic push
  #pragma gcc diagnostic ignored "-Wfloat-equal"
#endif
/**
 * @brief Check if two double values are almost equal within epsilon
 * @param first First value to compare
 * @param second Second value to compare
 * @param epsilon Tolerance for comparison
 * @return true if values are almost equal, false otherwise
 */
inline auto
almostEqual(const double first, const double second, const double epsilon)
    -> bool
{
  // shortcut and handles inf values
  if (first == second) {
    return true;
  }

  if (std::isnan(first) || std::isnan(second)) {
    return std::isnan(first) && std::isnan(second);
  }

  const double absFirst  = std::fabs(first);
  const double absSecond = std::fabs(second);
  const double diff      = std::fabs(first - second);
  // fixed epsilon
  if (diff <= epsilon) {
    return true;
  }

  return diff <= epsilon * std::max(absFirst, absSecond); // adaptative epsilon
}

/**
 * @brief Check if two Kernel::FT values are almost equal within epsilon
 * @param first First value to compare
 * @param second Second value to compare
 * @param epsilon Tolerance for comparison
 * @return true if values are almost equal, false otherwise
 */
inline auto
almostEqual(const Kernel::FT &first, const Kernel::FT &second,
            const Kernel::FT &epsilon) -> bool
{
  // shortcut and handles inf values
  if (first == second) {
    return true;
  }

  const Kernel::FT absFirst  = abs(first);
  const Kernel::FT absSecond = abs(second);
  const Kernel::FT diff      = abs(first - second);
  // fixed epsilon
  if (diff <= epsilon) {
    return true;
  }

  return diff <= epsilon * std::max(absFirst, absSecond); // adaptative epsilon
}

#ifdef __clang__
  #pragma clang diagnostic pop
#elif defined(__GCC__)
  #pragma gcc diagnostic pop
#endif

/**
 * @brief shortcut to get NaN for double
 * @return NaN (Not a Number) value for double
 */
inline auto
NaN() -> double
{
  return std::numeric_limits<double>::quiet_NaN();
}

/**
 * @brief round a double to the nearest integer
 * @param value Value to round
 * @return Rounded value
 */
inline auto
round(const double &value) -> double
{
  if (value < 0.0) {
    return ::ceil(value - 0.5);
  }
  return ::floor(value + 0.5);
}

#ifdef CGAL_USE_GMPXX
/**
 * @brief floor a rational to an integer
 * @param value Rational value to floor
 * @return Floor of the rational as integer
 */
SFCGAL_API auto
floor(const ::mpq_class &value) -> ::mpz_class;
/**
 * @brief ceil a rational to an integer
 * @param value Rational value to ceil
 * @return Ceiling of the rational as integer
 */
SFCGAL_API auto
ceil(const ::mpq_class &value) -> ::mpz_class;
/**
 * @brief round a rational to an integer
 * @param value Rational value to round
 * @return Rounded rational as integer
 */
SFCGAL_API auto
round(const ::mpq_class &value) -> ::mpz_class;
#endif

/**
 * @brief floor a rational to an integer
 * @param value CGAL rational value to floor
 * @return Floor of the rational as CGAL integer
 */
SFCGAL_API auto
floor(const CGAL::Gmpq &value) -> CGAL::Gmpz;
/**
 * @brief ceil a rational to an integer
 * @param value CGAL rational value to ceil
 * @return Ceiling of the rational as CGAL integer
 */
SFCGAL_API auto
ceil(const CGAL::Gmpq &value) -> CGAL::Gmpz;
/**
 * @brief round a rational to an integer
 * @param value CGAL rational value to round
 * @return Rounded rational as CGAL integer
 */
SFCGAL_API auto
round(const CGAL::Gmpq &value) -> CGAL::Gmpz;

/**
 * @brief Normalizes a vector
 * @param vector The vector to normalize
 * @return The normalized vector
 */
inline auto
normalizeVector(const Kernel::Vector_3 &vector) -> Kernel::Vector_3
{
  Kernel::FT length = CGAL::sqrt(CGAL::to_double(vector.squared_length()));
  // clang-tidy wrongly assumes that vector / length might leak
  return (length > 0) ? vector / length : vector;
} // NOLINT(clang-analyzer-cplusplus.NewDeleteLeaks)
} // namespace SFCGAL

#endif
