// Copyright (c) 2024-2025, SFCGAL team.
//
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/config.h"

#ifdef SFCGAL_CGAL_VERSION_MAJOR >= 6

  #include "SFCGAL/version.h"
  #include <memory>

  #include <SFCGAL/Exception.h>
  #include <SFCGAL/LineString.h>
  #include <SFCGAL/MultiPolygon.h>
  #include <SFCGAL/Point.h>
  #include <SFCGAL/Polygon.h>
  #include <SFCGAL/Triangle.h>
  #include <SFCGAL/algorithm/polygonRepair.h>
  #include <SFCGAL/io/wkt.h>

  #include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

using namespace SFCGAL;
using namespace SFCGAL::algorithm;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_PolygonRepairTest)

/**
 * Test basic functionality with simple valid polygon
 */
BOOST_AUTO_TEST_CASE(testValidPolygon)
{
  auto polygon = io::readWkt("POLYGON((0 0, 1 0, 1 1, 0 1, 0 0))");
  auto result  = polygonRepair(*polygon);

  BOOST_CHECK(!result->isEmpty());

  // Test with explicit rule
  result = polygonRepair(*polygon, PolygonRepairRule::EVEN_ODD_RULE);
  BOOST_CHECK_EQUAL(
      result->asText(2),
      "MULTIPOLYGON (((0.00 0.00,1.00 0.00,1.00 1.00,0.00 1.00,0.00 0.00)))");
}

/**
 * Test basic functionality with simple valid polygon XYZM
 */
BOOST_AUTO_TEST_CASE(testValidPolygonXYZM)
{
  auto polygon =
      io::readWkt("POLYGON ZM((0 0 1 2, 1 0 2 3, 1 1 3 4, 0 1 4 5, 0 0 1 2))");
  auto result = polygonRepair(*polygon);

  BOOST_CHECK(!result->isEmpty());

  // Test with explicit rule
  result = polygonRepair(*polygon, PolygonRepairRule::EVEN_ODD_RULE);
  BOOST_CHECK_EQUAL(
      result->asText(2),
      "MULTIPOLYGON (((0.00 0.00,1.00 0.00,1.00 1.00,0.00 1.00,0.00 0.00)))");
}

/**
 * Test with self-intersecting polygon
 */
BOOST_AUTO_TEST_CASE(testSelfIntersectingPolygon)
{
  // Figure-8 polygon
  auto polygon = io::readWkt("POLYGON((0 0, 2 2, 2 0, 0 2, 0 0))");
  auto result  = polygonRepair(*polygon);

  // Self-intersecting polygon typically results in multipolygon with even-odd
  // rule
  BOOST_CHECK_EQUAL(result->asText(2),
                    "MULTIPOLYGON (((0.00 0.00,1.00 1.00,0.00 2.00,0.00 "
                    "0.00)),((1.00 1.00,2.00 0.00,2.00 2.00,1.00 1.00)))");
}

/**
 * Test with polygon with holes
 */
BOOST_AUTO_TEST_CASE(testPolygonWithHoles)
{
  auto polygon = io::readWkt(
      "POLYGON((0 0, 4 0, 4 4, 0 4, 0 0), (1 1, 3 1, 3 3, 1 3, 1 1))");
  auto result = polygonRepair(*polygon);

  BOOST_CHECK_EQUAL(result->asText(2),
                    "MULTIPOLYGON (((0.00 0.00,4.00 0.00,4.00 4.00,0.00 "
                    "4.00,0.00 0.00),(1.00 "
                    "1.00,1.00 3.00,3.00 3.00,3.00 1.00,1.00 1.00)))");
}

/**
 * Test with MultiPolygon
 */
BOOST_AUTO_TEST_CASE(testMultiPolygon)
{
  auto multipolygon = io::readWkt(
      "MULTIPOLYGON(((0 0, 1 0, 1 1, 0 1, 0 0)), ((2 2, 3 2, 3 3, 2 3, 2 2)))");
  auto result = polygonRepair(*multipolygon);

  BOOST_CHECK_EQUAL(
      result->asText(2),
      "MULTIPOLYGON (((0.00 0.00,1.00 0.00,1.00 1.00,0.00 1.00,0.00 "
      "0.00)),((2.00 2.00,3.00 2.00,3.00 3.00,2.00 3.00,2.00 2.00)))");
}

/**
 * Test empty geometry handling
 */
BOOST_AUTO_TEST_CASE(testEmptyGeometry)
{
  auto empty  = io::readWkt("POLYGON EMPTY");
  auto result = polygonRepair(*empty);

  BOOST_CHECK_EQUAL(result->asText(2), "MULTIPOLYGON EMPTY");
}

/**
 * Test different repair rules
 */
BOOST_AUTO_TEST_CASE(testRepairRules)
{
  // Self-intersecting polygon for testing different rules
  auto polygon = io::readWkt("POLYGON((0 0, 2 2, 2 0, 0 2, 0 0))");

  // Test even-odd rule
  auto even_odd = polygonRepair(*polygon, PolygonRepairRule::EVEN_ODD_RULE);
  BOOST_CHECK_EQUAL(even_odd->asText(2),
                    "MULTIPOLYGON (((0.00 0.00,1.00 1.00,0.00 2.00,0.00 "
                    "0.00)),((1.00 1.00,2.00 0.00,2.00 2.00,1.00 1.00)))");

  #if SFCGAL_CGAL_VERSION_MAJOR == 6 && SFCGAL_CGAL_VERSION_MINOR >= 1
  // Test non-zero rule
  auto non_zero = polygonRepair(*polygon, PolygonRepairRule::NON_ZERO_RULE);
  BOOST_CHECK_EQUAL(non_zero->asText(2),
                    "MULTIPOLYGON (((0.00 0.00,1.00 1.00,0.00 2.00,0.00 "
                    "0.00)),((1.00 1.00,2.00 0.00,2.00 2.00,1.00 1.00)))");

  // Test union rule
  auto union_result = polygonRepair(*polygon, PolygonRepairRule::UNION_RULE);
  BOOST_CHECK_EQUAL(
      union_result->asText(2),
      "MULTIPOLYGON (((0.00 0.00,1.00 1.00,0.00 2.00,0.00 0.00)))");

  // Test intersection rule
  auto intersection =
      polygonRepair(*polygon, PolygonRepairRule::INTERSECTION_RULE);
  BOOST_CHECK_EQUAL(
      intersection->asText(2),
      "MULTIPOLYGON (((0.00 0.00,1.00 1.00,0.00 2.00,0.00 0.00)))");
  #endif
}

/**
 * Test unsupported geometry type
 */
BOOST_AUTO_TEST_CASE(testUnsupportedGeometry)
{
  auto point = io::readWkt("POINT(1 2)");

  BOOST_CHECK_THROW(polygonRepair(*point), SFCGAL::Exception);
}

  #if SFCGAL_CGAL_VERSION_MAJOR == 6 && SFCGAL_CGAL_VERSION_MINOR >= 1
/**
 * Test overlapping polygons in multipolygon
 */
BOOST_AUTO_TEST_CASE(testOverlappingMultiPolygon)
{
  // Two overlapping squares
  auto multipolygon = io::readWkt(
      "MULTIPOLYGON(((0 0, 2 0, 2 2, 0 2, 0 0)), ((1 1, 3 1, 3 3, 1 3, 1 1)))");
  auto result = polygonRepair(*multipolygon, PolygonRepairRule::UNION_RULE);

  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK_EQUAL(
      result->asText(2),
      "MULTIPOLYGON (((0.00 0.00,2.00 0.00,2.00 1.00,3.00 1.00,3.00 "
      "3.00,1.00 3.00,1.00 2.00,0.00 2.00,0.00 0.00)))");
}
  #endif

/**
 * Test with invalid orientation
 */
BOOST_AUTO_TEST_CASE(testInvalidOrientation)
{
  // Clockwise polygon (invalid in OGC)
  auto polygon = io::readWkt("POLYGON((0 0, 0 1, 1 1, 1 0, 0 0))");
  auto result  = polygonRepair(*polygon);

  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK_EQUAL(
      result->asText(2),
      "MULTIPOLYGON (((0.00 0.00,1.00 0.00,1.00 1.00,0.00 1.00,0.00 0.00)))");
}

/**
 * Test with degenerate polygon
 */
BOOST_AUTO_TEST_CASE(testDegeneratePolygon)
{
  // Polygon collapsing to a line
  auto polygon = io::readWkt("POLYGON((0 0, 1 0, 1 0, 0 0, 0 0))");
  auto result  = polygonRepair(*polygon);
  BOOST_CHECK_EQUAL(result->asText(2), "MULTIPOLYGON EMPTY");
}

/**
 * Test with polygon having duplicate points
 */
BOOST_AUTO_TEST_CASE(testPolygonWithDuplicates)
{
  auto polygon = io::readWkt("POLYGON((0 0, 1 0, 1 0, 1 1, 0 1, 0 1, 0 0))");
  auto result  = polygonRepair(*polygon);

  BOOST_CHECK_EQUAL(
      result->asText(2),
      "MULTIPOLYGON (((0.00 0.00,1.00 0.00,1.00 1.00,0.00 1.00,0.00 0.00)))");
}

BOOST_AUTO_TEST_SUITE_END()

#endif // SFCGAL_CGAL_VERSION_MAJOR >= 6
