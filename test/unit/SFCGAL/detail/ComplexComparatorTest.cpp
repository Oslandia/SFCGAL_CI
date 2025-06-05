// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "SFCGAL/detail/ComplexComparator.h"
#include "SFCGAL/numeric.h"

using namespace SFCGAL;
using namespace SFCGAL::detail;

// always after CGAL
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_detail_ComplexComparatorTest)

BOOST_AUTO_TEST_CASE(testComparator)
{
  ComplexComparator const less;

  BOOST_CHECK(
      !less(std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0)));
  BOOST_CHECK(
      !less(std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0)));
  BOOST_CHECK(
      less(std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0)));
  BOOST_CHECK(
      less(std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 1.0)));
}

BOOST_AUTO_TEST_SUITE_END()
