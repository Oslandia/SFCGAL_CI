// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/io/wkt.h"

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_GeometryTest)

BOOST_AUTO_TEST_CASE(testAlmostEqual)
{
  {
    std::unique_ptr<Geometry> const gA(io::readWkt(
        "TIN Z (((1/1 1/2 0/1,1/2 1/2 0/1,1/1 1/4 0/1,1/1 1/2 0/1)),"
        "((1/2 0/1 0/1,1/1 0/1 0/1,1/1 1/4 0/1,1/2 0/1 0/1)),"
        "((1/2 0/1 0/1,1/1 1/4 0/1,1/2 1/2 0/1,1/2 0/1 0/1)))"));

    // sub geometry order change
    std::unique_ptr<Geometry> const gB(io::readWkt(
        "TIN Z (((1/2 1/2 0/1,1/1 1/4 0/1,1/1 1/2 0/1,1/2 1/2 0/1)),"
        "((1/1 0/1 0/1,1/1 1/4 0/1,1/2 1/2 0/1,1/1 0/1 0/1)),"
        "((1/2 1/2 0/1,1/2 0/1 0/1,1/1 0/1 0/1,1/2 1/2 0/1)))"));

    BOOST_CHECK(SFCGAL::algorithm::isValid(*gA.get()));
    BOOST_CHECK(SFCGAL::algorithm::isValid(*gB.get()));
    // not strict
    BOOST_CHECK(gA->almostEqual(
        *gB, 0.0, algorithm::EqualityStrictness::coverSubGeomNonOrdered()));
  }

  std::unique_ptr<Geometry> const gA(io::readWkt(
      "MULTIPOLYGON (((0.0 0.0, 5.0 -1.0, 6.0 0.0, 5.0 1.0, 0.0 0.0)), "
      "((6.0 0.0, 15.0 -1.0, 16.0 0.0, 15.0 1.0, 6.0 0.0)))"));

  // sub geometry order change
  std::unique_ptr<Geometry> const gB(io::readWkt(
      "MULTIPOLYGON (((6.0 0.0, 15.0 -1.0, 16.0 0.0, 15.0 1.0, 6.0 0.0)), "
      "((0.0 0.0, 5.0 -1.0, 6.0 0.0, 5.0 1.0, 0.0 0.0)))"));

  BOOST_CHECK(SFCGAL::algorithm::isValid(*gA.get()));
  BOOST_CHECK(SFCGAL::algorithm::isValid(*gB.get()));
  // not strict
  BOOST_CHECK(gA->almostEqual(
      *gB, 0.0, algorithm::EqualityStrictness::pointNonOrdered()));
  // only internal point needed to be ordered
  BOOST_CHECK(gA->almostEqual(
      *gB, 0.0, algorithm::EqualityStrictness::InternalPointOrdered));
  // strict, should fail
  BOOST_CHECK(!gA->almostEqual(
      *gB, 0.0, algorithm::EqualityStrictness::allPointOrdered()));
  // not strict, 0.0 tolerance
  BOOST_CHECK(*gA == *gB);

  // polygon point order change
  std::unique_ptr<Geometry> const gC(io::readWkt(
      "MULTIPOLYGON (((5.0 -1.0, 6.0 0.0, 5.0 1.0, 0.0 0.0, 5.0 -1.0)), "
      "((6.0 0.0, 15.0 -1.0, 16.0 0.0, 15.0 1.0, 6.0 0.0)))"));

  BOOST_CHECK(SFCGAL::algorithm::isValid(*gC.get()));
  // not strict
  BOOST_CHECK(gA->almostEqual(
      *gC, 0.0, algorithm::EqualityStrictness::pointNonOrdered()));
  // shifted point allowed but sub part are ordered
  BOOST_CHECK(
      gA->almostEqual(*gC, 0.0,
                      algorithm::EqualityStrictness() |
                          algorithm::EqualityStrictness::InternalPointShifted |
                          algorithm::EqualityStrictness::SubPartOrdered));
  // strict
  BOOST_CHECK(!gA->almostEqual(
      *gC, 0.0, algorithm::EqualityStrictness::allPointOrdered()));
  // not strict, 0.0 tolerance
  BOOST_CHECK(*gA == *gC);

  // not strict
  BOOST_CHECK(gB->almostEqual(
      *gC, 0.0, algorithm::EqualityStrictness::pointNonOrdered()));
  // strict
  BOOST_CHECK(gB->almostEqual(
      *gC, 0.0, algorithm::EqualityStrictness::InternalPointShifted));
  // strict
  BOOST_CHECK(!gB->almostEqual(
      *gC, 0.0, algorithm::EqualityStrictness::allPointOrdered()));
  // not strict, 0.0 tolerance
  BOOST_CHECK(*gB == *gC);

  // slight change
  std::unique_ptr<Geometry> const gD(io::readWkt(
      "MULTIPOLYGON (((0.1 0.0, 5.0 -1.1, 6.1 0.0, 5.0 1.0, 0.1 0.0)), "
      "((6.1 0.0, 15.0 -1.0, 16.0 0.1, 15.0 1.0, 6.1 0.0)))"));

  BOOST_CHECK(SFCGAL::algorithm::isValid(*gD.get()));
  // not strict, 0.0 tolerance, should fail
  BOOST_CHECK(!gA->almostEqual(
      *gD, 0.0, algorithm::EqualityStrictness::pointNonOrdered()));
  // not strict, 0.1 tolerance
  BOOST_CHECK(gA->almostEqual(
      *gD, 0.1, algorithm::EqualityStrictness::pointNonOrdered()));
  // strict, 0.1 tolerance
  BOOST_CHECK(gA->almostEqual(
      *gD, 0.1, algorithm::EqualityStrictness::allPointOrdered()));
  // not strict, 0.0 tolerance, should fail
  BOOST_CHECK(!(*gA == *gD));

  // not strict, 0.1 tolerance
  BOOST_CHECK(gB->almostEqual(
      *gD, 0.1, algorithm::EqualityStrictness::pointNonOrdered()));
  // 0.1 tolerance, only internal point needed to be ordered
  BOOST_CHECK(gB->almostEqual(
      *gD, 0.1, algorithm::EqualityStrictness::InternalPointOrdered));
  // strict, 0.1 tolerance, should fail
  BOOST_CHECK(!gB->almostEqual(
      *gD, 0.1, algorithm::EqualityStrictness::allPointOrdered()));
  // not strict, 0.0 tolerance, should fail
  BOOST_CHECK(!(*gB == *gD));

  // really not the same polygons
  std::unique_ptr<Geometry> const gE(io::readWkt(
      "MULTIPOLYGON (((-5.0 80.0, 0.0 0.0, 5.0 -1.0, 5.0 1.0, -5.0 80.0)), "
      "((6.0 0.0, 15.0 -1.0, 16.0 0.0, 15.0 1.0, 6.0 0.0)))"));

  BOOST_CHECK(SFCGAL::algorithm::isValid(*gE.get()));
  // not strict, should fail
  BOOST_CHECK(!gA->almostEqual(
      *gE, 0.0, algorithm::EqualityStrictness::pointNonOrdered()));
  // not strict (inverted caller), should fail
  BOOST_CHECK(!gE->almostEqual(
      *gA, 0.0, algorithm::EqualityStrictness::pointNonOrdered()));
  // strict, should fail
  BOOST_CHECK(!gA->almostEqual(
      *gE, 0.0, algorithm::EqualityStrictness::allPointOrdered()));
  // not strict, 0.0 tolerance, should fail
  BOOST_CHECK(!(*gA == *gE));

  // can be the same but not exactly
  std::unique_ptr<Geometry> const gF(
      io::readWkt("MULTIPOLYGON (((0.0 0.0, 5.0 -1.0, 6.0 0.0, 5.5 1.0, 5.0 "
                  "0.5, 4.5 1.0, 0.0 0.0)), "
                  "((6.0 0.0, 15.0 -1.0, 16.0 0.0, 15.0 1.0, 6.0 0.0)))"));

  std::unique_ptr<Geometry> const gG(
      io::readWkt("MULTIPOLYGON (((0.0 0.0, 5.0 -1.0, 6.0 0.0, 5.0 0.5, 5.5 "
                  "1.0, 4.5 1.0, 0.0 0.0)), "
                  "((6.0 0.0, 15.0 -1.0, 16.0 0.0, 15.0 1.0, 6.0 0.0)))"));

  BOOST_CHECK(SFCGAL::algorithm::isValid(*gF.get()));
  BOOST_CHECK(SFCGAL::algorithm::isValid(*gG.get()));
  // not strict
  BOOST_CHECK(gF->almostEqual(
      *gG, 0.0, algorithm::EqualityStrictness::pointNonOrdered()));
  // not strict
  BOOST_CHECK(gG->almostEqual(
      *gF, 0.0, algorithm::EqualityStrictness::pointNonOrdered()));
  // only internal point needed to be ordered, should fail
  BOOST_CHECK(!gF->almostEqual(
      *gG, 0.0, algorithm::EqualityStrictness::InternalPointOrdered));
  // shifted point allowed, should fail
  BOOST_CHECK(!gF->almostEqual(
      *gG, 0.0, algorithm::EqualityStrictness::InternalPointShifted));
  // strict, should fail
  BOOST_CHECK(!gF->almostEqual(
      *gG, 0.0, algorithm::EqualityStrictness::allPointOrdered()));
  // not strict
  BOOST_CHECK(*gF == *gG);
}

BOOST_AUTO_TEST_SUITE_END()
