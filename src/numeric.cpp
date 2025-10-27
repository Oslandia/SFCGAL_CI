// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/numeric.h"

namespace SFCGAL {

/// @private
auto
floor(const CGAL::Gmpq &value) -> CGAL::Gmpz
{
  return value.numerator() / value.denominator();
}

/// @private
auto
ceil(const CGAL::Gmpq &value) -> CGAL::Gmpz
{
  CGAL::Gmpz result(0);
  mpz_cdiv_q(result.mpz(), value.numerator().mpz(), value.denominator().mpz());
  return result;
}

/// @private
auto
round(const CGAL::Gmpq &value) -> CGAL::Gmpz
{
  if (value < 0) {
    // ceil( value - 0.5 ) ;
    return ceil(value - CGAL::Gmpq(1, 2));
  }
  if (value == 0) {
    return 0;
  } // floor( value + 0.5 ) ;
  return floor(value + CGAL::Gmpq(1, 2));
}

#ifdef CGAL_USE_GMPXX
auto
floor(const mpq_class &value) -> mpz_class
{
  return value.get_num() / value.get_den();
}

auto
ceil(const mpq_class &value) -> mpz_class
{
  mpz_class result(0);
  mpz_cdiv_q(result.get_mpz_t(), value.get_num().get_mpz_t(),
             value.get_den().get_mpz_t());
  return result;
}

auto
round(const mpq_class &value) -> mpz_class
{
  if (value < 0) {
    // ceil( value - 0.5 ) ;
    mpq_class const tmp = value - mpq_class(1, 2);
    return ceil(tmp);
  }
  if (value == 0) {
    return 0;
  } // floor( value + 0.5 ) ;
  mpq_class const tmp = value + mpq_class(1, 2);
  return floor(tmp);
}
#endif

} // namespace SFCGAL
