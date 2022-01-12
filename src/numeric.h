// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _SFCGAL_NUMERIC_H_
#define _SFCGAL_NUMERIC_H_

#include <cmath>
#include <limits>

#include <SFCGAL/export.h>

#include <SFCGAL/Kernel.h>

namespace SFCGAL {
/**
 * shortcut to get NaN for double
 */
inline double
NaN()
{
  return std::numeric_limits<double>::quiet_NaN();
}

/**
 * @brief round a double to the nearest integer
 */
inline double
round(const double &v)
{
  if (v < 0.0) {
    return ::ceil(v - 0.5);
  } else {
    return ::floor(v + 0.5);
  }
}

#ifdef CGAL_USE_GMPXX
/**
 * @brief floor a rational to an integer
 */
SFCGAL_API ::mpz_class
floor(const ::mpq_class &v);
/**
 * @brief ceil a rational to an integer
 */
SFCGAL_API ::mpz_class
ceil(const ::mpq_class &v);
/**
 * @brief round a rational to an integer
 */
SFCGAL_API ::mpz_class
round(const ::mpq_class &v);
#endif

/**
 * @brief floor a rational to an integer
 */
SFCGAL_API CGAL::Gmpz
           floor(const CGAL::Gmpq &v);
/**
 * @brief ceil a rational to an integer
 */
SFCGAL_API CGAL::Gmpz
           ceil(const CGAL::Gmpq &v);
/**
 * @brief round a rational to an integer
 */
SFCGAL_API CGAL::Gmpz
           round(const CGAL::Gmpq &v);

} // namespace SFCGAL

#endif
