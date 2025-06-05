// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include <exception>

#include "SFCGAL/numeric.h"

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_NumericTest)

BOOST_AUTO_TEST_CASE(testFloorRational)
{
  BOOST_CHECK_EQUAL(SFCGAL::floor(CGAL::Gmpq(0)), 0);
  BOOST_CHECK_EQUAL(SFCGAL::floor(CGAL::Gmpq(1, 2)), 0);
  BOOST_CHECK_EQUAL(SFCGAL::floor(CGAL::Gmpq(1, 3)), 0);
  BOOST_CHECK_EQUAL(SFCGAL::floor(CGAL::Gmpq(2, 3)), 0);
  BOOST_CHECK_EQUAL(SFCGAL::floor(CGAL::Gmpq(1, 1)), 1);
  BOOST_CHECK_EQUAL(SFCGAL::floor(CGAL::Gmpq(4, 3)), 1);
}

BOOST_AUTO_TEST_CASE(testCeilRational)
{
  BOOST_CHECK_EQUAL(SFCGAL::ceil(CGAL::Gmpq(0)), 0);
  BOOST_CHECK_EQUAL(SFCGAL::ceil(CGAL::Gmpq(1, 2)), 1);
  BOOST_CHECK_EQUAL(SFCGAL::ceil(CGAL::Gmpq(1, 3)), 1);
  BOOST_CHECK_EQUAL(SFCGAL::ceil(CGAL::Gmpq(1, 1)), 1);
  BOOST_CHECK_EQUAL(SFCGAL::ceil(CGAL::Gmpq(4, 3)), 2);
}

BOOST_AUTO_TEST_CASE(testRoundRational)
{
  BOOST_CHECK_EQUAL(SFCGAL::round(CGAL::Gmpq(0)), 0);
  BOOST_CHECK_EQUAL(SFCGAL::round(CGAL::Gmpq(1, 2)), 1);
  BOOST_CHECK_EQUAL(SFCGAL::round(CGAL::Gmpq(1, 3)), 0);
  BOOST_CHECK_EQUAL(SFCGAL::round(CGAL::Gmpq(1, 1)), 1);
  BOOST_CHECK_EQUAL(SFCGAL::round(CGAL::Gmpq(4, 3)), 1);
}

BOOST_AUTO_TEST_SUITE_END()
