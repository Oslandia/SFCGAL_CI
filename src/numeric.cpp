// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/numeric.h>

namespace SFCGAL {

///
///
///
auto
floor(const CGAL::Gmpq &v) -> CGAL::Gmpz
{
  return v.numerator() / v.denominator();
}

///
///
///
auto
ceil(const CGAL::Gmpq &v) -> CGAL::Gmpz
{
  CGAL::Gmpz result(0);
  mpz_cdiv_q(result.mpz(), v.numerator().mpz(), v.denominator().mpz());
  return result;
}

///
///
///
auto
round(const CGAL::Gmpq &v) -> CGAL::Gmpz
{
  if (v < 0) {
    // ceil( v - 0.5 ) ;
    return ceil(v - CGAL::Gmpq(1, 2));
  } if (v == 0) {
    return 0;
  }     // floor( v + 0.5 ) ;
    return floor(v + CGAL::Gmpq(1, 2));
 
}

#ifdef CGAL_USE_GMPXX
///
///
///
auto
floor(const mpq_class &v) -> mpz_class
{
  return v.get_num() / v.get_den();
}

///
///
///
auto
ceil(const mpq_class &v) -> mpz_class
{
  mpz_class result(0);
  mpz_cdiv_q(result.get_mpz_t(), v.get_num().get_mpz_t(),
             v.get_den().get_mpz_t());
  return result;
}

///
///
///
auto
round(const mpq_class &v) -> mpz_class
{
  if (v < 0) {
    // ceil( v - 0.5 ) ;
    mpq_class const tmp = v - mpq_class(1, 2);
    return ceil(tmp);
  } if (v == 0) {
    return 0;
  }     // floor( v + 0.5 ) ;
    mpq_class const tmp = v + mpq_class(1, 2);
    return floor(tmp);
 
}
#endif

} // namespace SFCGAL
