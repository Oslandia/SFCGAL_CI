// Copyright (c) 2012-2023, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include <SFCGAL/LineString.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/algorithm/partition_2.h>

using namespace SFCGAL;

#include "../../../test_config.h"
// always after CGAL
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_Partition2Test)

BOOST_AUTO_TEST_CASE(testPartition2_NoPolygon)
{
  LineString exteriorRing;
  exteriorRing.addPoint(Point(0.0, 0.0));
  exteriorRing.addPoint(Point(2.0, 0.0));
  exteriorRing.addPoint(Point(2.0, 2.0));
  exteriorRing.addPoint(Point(1.0, 1.0));
  exteriorRing.addPoint(Point(0.0, 2.0));

  std::unique_ptr<Geometry> result(algorithm::partition_2(exteriorRing, SFCGAL::algorithm::y_monotone));
  BOOST_CHECK(result->isEmpty());

  std::string expectedWKT(
      "GEOMETRYCOLLECTION EMPTY");
  BOOST_CHECK_EQUAL(result->asText(1), expectedWKT);
}

BOOST_AUTO_TEST_CASE(testPartition2_Empty)
{
  Polygon g;

  std::unique_ptr<Geometry> result(algorithm::partition_2(g, SFCGAL::algorithm::y_monotone));
  BOOST_CHECK(result->isEmpty());

  std::string expectedWKT(
      "GEOMETRYCOLLECTION EMPTY");
  BOOST_CHECK_EQUAL(result->asText(1), expectedWKT);
}

BOOST_AUTO_TEST_CASE(testPartition2_YMonotonePartition2)
{
  LineString exteriorRing;
  exteriorRing.addPoint(Point(0.0, 0.0));
  exteriorRing.addPoint(Point(2.0, 0.0));
  exteriorRing.addPoint(Point(2.0, 2.0));
  exteriorRing.addPoint(Point(1.0, 1.0));
  exteriorRing.addPoint(Point(0.0, 2.0));

  Polygon g(exteriorRing);

  std::unique_ptr<Geometry> result(algorithm::partition_2(g, SFCGAL::algorithm::y_monotone));
  BOOST_CHECK(!result->isEmpty());

  std::string expectedWKT(
      "GEOMETRYCOLLECTION(POLYGON((0.0 0.0,2.0 0.0,2.0 2.0,1.0 1.0,0.0 0.0)))");
  BOOST_CHECK_EQUAL(result->asText(1), expectedWKT);
}

BOOST_AUTO_TEST_CASE(testPartition2_YMonotonePartition2_gross)
{
  LineString exteriorRing;
  exteriorRing.addPoint(Point(391, 374));
  exteriorRing.addPoint(Point(240, 431));
  exteriorRing.addPoint(Point(252, 340));
  exteriorRing.addPoint(Point(374, 320));
  exteriorRing.addPoint(Point(289, 214));
  exteriorRing.addPoint(Point(134, 390));
  exteriorRing.addPoint(Point(68, 186));
  exteriorRing.addPoint(Point(154, 259));
  exteriorRing.addPoint(Point(161, 107));
  exteriorRing.addPoint(Point(435, 108));
  exteriorRing.addPoint(Point(208, 148));
  exteriorRing.addPoint(Point(295, 160));
  exteriorRing.addPoint(Point(421, 212));
  exteriorRing.addPoint(Point(441, 303));
  exteriorRing.addPoint(Point(391, 374));

  Polygon g(exteriorRing);

  std::unique_ptr<Geometry> result(algorithm::partition_2(g));
  BOOST_CHECK(!result->isEmpty());

  std::string expectedWKT(
      "GEOMETRYCOLLECTION(POLYGON((134.0 390.0,68.0 186.0,154.0 259.0,134.0 "
      "390.0)),POLYGON((289.0 214.0,134.0 390.0,154.0 259.0,161.0 107.0,435.0 "
      "108.0,208.0 148.0,295.0 160.0,421.0 212.0,289.0 214.0)),POLYGON((391.0 "
      "374.0,240.0 431.0,252.0 340.0,374.0 320.0,289.0 214.0,421.0 212.0,441.0 "
      "303.0,391.0 374.0)))");
  BOOST_CHECK_EQUAL(result->asText(1), expectedWKT);
}

BOOST_AUTO_TEST_CASE(testPartition2_ApproxConvexPartition2_gross)
{
  LineString exteriorRing;
  exteriorRing.addPoint(Point(391, 374));
  exteriorRing.addPoint(Point(240, 431));
  exteriorRing.addPoint(Point(252, 340));
  exteriorRing.addPoint(Point(374, 320));
  exteriorRing.addPoint(Point(289, 214));
  exteriorRing.addPoint(Point(134, 390));
  exteriorRing.addPoint(Point(68, 186));
  exteriorRing.addPoint(Point(154, 259));
  exteriorRing.addPoint(Point(161, 107));
  exteriorRing.addPoint(Point(435, 108));
  exteriorRing.addPoint(Point(208, 148));
  exteriorRing.addPoint(Point(295, 160));
  exteriorRing.addPoint(Point(421, 212));
  exteriorRing.addPoint(Point(441, 303));
  exteriorRing.addPoint(Point(391, 374));

  Polygon g(exteriorRing);

  std::unique_ptr<Geometry> result(algorithm::partition_2(g, SFCGAL::algorithm::approx_convex));
  BOOST_CHECK(!result->isEmpty());

  std::string expectedWKT( "GEOMETRYCOLLECTION(POLYGON((391.0 374.0,240.0 431.0,252.0 340.0,374.0 320.0,391.0 374.0)),POLYGON((134.0 390.0,68.0 186.0,154.0 259.0,134.0 390.0)),POLYGON((289.0 214.0,134.0 390.0,154.0 259.0,289.0 214.0)),POLYGON((161.0 107.0,435.0 108.0,208.0 148.0,161.0 107.0)),POLYGON((154.0 259.0,161.0 107.0,208.0 148.0,154.0 259.0)),POLYGON((289.0 214.0,154.0 259.0,208.0 148.0,295.0 160.0,289.0 214.0)),POLYGON((374.0 320.0,289.0 214.0,295.0 160.0,421.0 212.0,374.0 320.0)),POLYGON((391.0 374.0,374.0 320.0,421.0 212.0,441.0 303.0,391.0 374.0)))"
      );
  BOOST_CHECK_EQUAL(result->asText(1), expectedWKT);
}

BOOST_AUTO_TEST_CASE(testPartition2_GreeneApproxConvexPartition2_gross)
{
  LineString exteriorRing;
  exteriorRing.addPoint(Point(391, 374));
  exteriorRing.addPoint(Point(240, 431));
  exteriorRing.addPoint(Point(252, 340));
  exteriorRing.addPoint(Point(374, 320));
  exteriorRing.addPoint(Point(289, 214));
  exteriorRing.addPoint(Point(134, 390));
  exteriorRing.addPoint(Point(68, 186));
  exteriorRing.addPoint(Point(154, 259));
  exteriorRing.addPoint(Point(161, 107));
  exteriorRing.addPoint(Point(435, 108));
  exteriorRing.addPoint(Point(208, 148));
  exteriorRing.addPoint(Point(295, 160));
  exteriorRing.addPoint(Point(421, 212));
  exteriorRing.addPoint(Point(441, 303));
  exteriorRing.addPoint(Point(391, 374));

  Polygon g(exteriorRing);

  std::unique_ptr<Geometry> result(algorithm::partition_2(g, SFCGAL::algorithm::greene_approx_convex));
  BOOST_CHECK(!result->isEmpty());

  std::string expectedWKT( "GEOMETRYCOLLECTION(POLYGON((134.0 390.0,68.0 186.0,154.0 259.0,134.0 390.0)),POLYGON((161.0 107.0,435.0 108.0,208.0 148.0,161.0 107.0)),POLYGON((208.0 148.0,295.0 160.0,421.0 212.0,289.0 214.0,208.0 148.0)),POLYGON((154.0 259.0,161.0 107.0,208.0 148.0,154.0 259.0)),POLYGON((289.0 214.0,134.0 390.0,154.0 259.0,208.0 148.0,289.0 214.0)),POLYGON((374.0 320.0,289.0 214.0,421.0 212.0,374.0 320.0)),POLYGON((374.0 320.0,421.0 212.0,441.0 303.0,391.0 374.0,374.0 320.0)),POLYGON((391.0 374.0,240.0 431.0,252.0 340.0,374.0 320.0,391.0 374.0)))"
      );
  BOOST_CHECK_EQUAL(result->asText(1), expectedWKT);
}

BOOST_AUTO_TEST_CASE(testPartition2_OptimalConvexPartition2_gross)
{
  LineString exteriorRing;
  exteriorRing.addPoint(Point(391, 374));
  exteriorRing.addPoint(Point(240, 431));
  exteriorRing.addPoint(Point(252, 340));
  exteriorRing.addPoint(Point(374, 320));
  exteriorRing.addPoint(Point(289, 214));
  exteriorRing.addPoint(Point(134, 390));
  exteriorRing.addPoint(Point(68, 186));
  exteriorRing.addPoint(Point(154, 259));
  exteriorRing.addPoint(Point(161, 107));
  exteriorRing.addPoint(Point(435, 108));
  exteriorRing.addPoint(Point(208, 148));
  exteriorRing.addPoint(Point(295, 160));
  exteriorRing.addPoint(Point(421, 212));
  exteriorRing.addPoint(Point(441, 303));
  exteriorRing.addPoint(Point(391, 374));

  Polygon g(exteriorRing);

  std::unique_ptr<Geometry> result(algorithm::partition_2(g, SFCGAL::algorithm::optimal_convex));
  BOOST_CHECK(!result->isEmpty());

  std::string expectedWKT( "GEOMETRYCOLLECTION(POLYGON((391.0 374.0,240.0 431.0,252.0 340.0,374.0 320.0,391.0 374.0)),POLYGON((134.0 390.0,68.0 186.0,154.0 259.0,134.0 390.0)),POLYGON((161.0 107.0,435.0 108.0,208.0 148.0,161.0 107.0)),POLYGON((154.0 259.0,161.0 107.0,208.0 148.0,154.0 259.0)),POLYGON((289.0 214.0,134.0 390.0,154.0 259.0,208.0 148.0,295.0 160.0,289.0 214.0)),POLYGON((374.0 320.0,289.0 214.0,295.0 160.0,421.0 212.0,441.0 303.0,374.0 320.0)),POLYGON((391.0 374.0,374.0 320.0,441.0 303.0,391.0 374.0)))"
      );
  BOOST_CHECK_EQUAL(result->asText(1), expectedWKT);
}
BOOST_AUTO_TEST_SUITE_END()
