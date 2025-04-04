// Copyright (c) 2025-2025, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/simplification.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/io/wkt.h"
#include <boost/test/unit_test.hpp>
using namespace SFCGAL;
#include "../../../test_config.h"
// always after CGAL
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_SimplificationTest)

// Harmonized helper function to test simplification and dimension preservation
void
testSimplification(const std::string &sourceWkt, const std::string &expectedWkt,
                   double threshold, bool preserveTopology)
{
  std::unique_ptr<Geometry> source(io::readWkt(sourceWkt));
  bool                      sourceIs3D       = source->is3D();
  bool                      sourceIsMeasured = source->isMeasured();

  std::unique_ptr<Geometry> simplified =
      algorithm::simplify(*source, threshold, preserveTopology);
  bool resultIs3D       = simplified->is3D();
  bool resultIsMeasured = simplified->isMeasured();

  // Parse expected result
  std::unique_ptr<Geometry> expected(io::readWkt(expectedWkt));

  // Display results
  std::cout << "Source:\t\t" << source->asText(2) << std::endl;
  std::cout << "Source Dimensions: is3D=" << (sourceIs3D ? "true" : "false")
            << ", isMeasured=" << (sourceIsMeasured ? "true" : "false")
            << std::endl;
  std::cout << "Result:\t\t" << simplified->asText(2) << std::endl;
  std::cout << "Expected:\t" << expected->asText(2) << std::endl;
  std::cout << "Result Dimensions: is3D=" << (resultIs3D ? "true" : "false")
            << ", isMeasured=" << (resultIsMeasured ? "true" : "false")
            << std::endl;
  std::cout << "Threshold: " << threshold << std::endl;
  std::cout << "Preserve Topology: " << (preserveTopology ? "true" : "false")
            << std::endl;

  // Verify dimension preservation
  BOOST_CHECK_EQUAL(sourceIs3D, resultIs3D);
  BOOST_CHECK_EQUAL(sourceIsMeasured, resultIsMeasured);

  // Verify the result matches the expected geometry
  BOOST_CHECK_EQUAL(simplified->asText(2), expected->asText(2));

  std::cout << "Dimension preservation: "
            << (sourceIs3D == resultIs3D && sourceIsMeasured == resultIsMeasured
                    ? "PASSED"
                    : "FAILED")
            << std::endl;
  std::cout << "Result matches expected: "
            << ((simplified->asText(2) == expected->asText(2)) ? "PASSED"
                                                               : "FAILED")
            << std::endl;
  std::cout << "-------------------" << std::endl;
}

// Test Point geometries
BOOST_AUTO_TEST_CASE(testSimplify_Point)
{
  // XY Point
  testSimplification("POINT(1 2)", "POINT(1 2)", 5, false);

  // XYZ Point
  testSimplification("POINT Z(1 2 3)", "POINT Z(1 2 3)", 5, false);

  // XYM Point
  testSimplification("POINT M(1 2 100)", "POINT M(1 2 100)", 5, false);

  // XYZM Point
  testSimplification("POINT ZM(1 2 3 100)", "POINT ZM(1 2 3 100)", 5, false);

  // Empty Point
  testSimplification("POINT EMPTY", "POINT EMPTY", 5, false);
}

// Test MultiPoint geometries
BOOST_AUTO_TEST_CASE(testSimplify_MultiPoint)
{
  // XY MultiPoint
  testSimplification("MULTIPOINT((1 2), (3 4), (5 6))",
                     "MULTIPOINT((1 2), (3 4), (5 6))", 5, false);

  // XYZ MultiPoint
  testSimplification("MULTIPOINT Z((1 2 3), (4 5 6), (7 8 9))",
                     "MULTIPOINT Z((1 2 3), (4 5 6), (7 8 9))", 5, false);

  // XYM MultiPoint
  testSimplification("MULTIPOINT M((1 2 100), (3 4 200), (5 6 300))",
                     "MULTIPOINT M((1 2 100), (3 4 200), (5 6 300))", 5, false);

  // XYZM MultiPoint
  testSimplification("MULTIPOINT ZM((1 2 3 100), (4 5 6 200), (7 8 9 300))",
                     "MULTIPOINT ZM((1 2 3 100), (4 5 6 200), (7 8 9 300))", 5,
                     false);

  // Empty MultiPoint
  testSimplification("MULTIPOINT EMPTY", "MULTIPOINT EMPTY", 5, false);
}

BOOST_AUTO_TEST_CASE(testSimplify_LineString)
{
  // XY LineString - small threshold
  testSimplification("LINESTRING(1 4, 4 9, 4 12, 4 16, 2 19, -4 20)",
                     "LINESTRING(1 4, 4 9, 4 16, 2 19, -4 20)", 1, false);

  // XY LineString - larger threshold
  testSimplification("LINESTRING(1 4, 4 9, 4 12, 4 16, 2 19, -4 20)",
                     "LINESTRING(1 4, 2 19, -4 20)", 5, false);

  // Empty LineString
  testSimplification("LINESTRING EMPTY", "LINESTRING EMPTY", 5, false);

  // Create LineString with 3D coordinates
  testSimplification(
      "LINESTRING Z(1 4 10, 4 9 20, 4 12 15, 4 16 30, 2 19 25, -4 20 15)",
      "LINESTRING Z(1 4 10, 2 19 25, -4 20 15)", 5, false);

  // Create LineString with measure
  testSimplification(
      "LINESTRING M(1 4 100, 4 9 200, 4 12 150, 4 16 300, 2 19 250, -4 20 150)",
      "LINESTRING M(1 4 100, 2 19 250, -4 20 150)", 5, false);

  // Create LineString with 3D coordinates and measure
  testSimplification("LINESTRING ZM(1 4 10 100, 4 9 20 200, 4 12 15 150, 4 16 "
                     "30 300, 2 19 25 250, -4 20 15 150)",
                     "LINESTRING ZM(1 4 10 100, 2 19 25 250, -4 20 15 150)", 5,
                     false);
}

BOOST_AUTO_TEST_CASE(testSimplify_MultiLineString)
{
  // XY MultiLineString
  testSimplification("MULTILINESTRING((1 4, 4 9, 4 12, 4 16, 2 19, -4 20), (6 "
                     "4, 4 9, 4 12, 4 16, 5 19, 10 20))",
                     "MULTILINESTRING((1 4, 2 19, -4 20), (6 4, 10 20))", 5,
                     false);

  // With preserveTopology
  testSimplification(
      "MULTILINESTRING((1 4, 4 9, 4 12, 4 16, 2 19, -4 20), (6 4, 4 9, 4 12, 4 "
      "16, 5 19, 10 20))",
      "MULTILINESTRING((1 4, 4 9, 4 16, -4 20), (6 4, 4 9, 4 16, 10 20))", 5,
      true);

  // Empty MultiLineString
  testSimplification("MULTILINESTRING EMPTY", "MULTILINESTRING EMPTY", 5,
                     false);

  // Test XYZ
  std::string sourceXYZ =
      "MULTILINESTRING Z((1 4 10, 4 9 20, 4 12 15, 4 16 30, 2 19 25, -4 20 "
      "15), "
      "(6 4 5, 4 9 15, 4 12 25, 4 16 20, 5 19 15, 10 20 10))";
  std::string expectedXYZ_noPreserve =
      "MULTILINESTRING Z((1 4 10, 2 19 25, -4 20 15), (6 4 5, 10 20 10))";
  std::string expectedXYZ_preserve =
      "MULTILINESTRING Z((1 4 10, 4 9 20, 4 16 30, -4 20 15), (6 4 5, 4 9 20, "
      "4 16 30, 10 20 10))";
  // CGAL keep first Z when XY coincide
  testSimplification(sourceXYZ, expectedXYZ_noPreserve, 5, false);
  testSimplification(sourceXYZ, expectedXYZ_preserve, 5, true);

  // Test XYM
  std::string sourceXYM =
      "MULTILINESTRING M((1 4 100, 4 9 200, 4 12 150, 4 16 300, 2 19 250, -4 "
      "20 150), "
      "(6 4 50, 4 9 150, 4 12 250, 4 16 200, 5 19 150, 10 20 100))";
  std::string expectedXYM_noPreserve =
      "MULTILINESTRING M((1 4 100, 2 19 250, -4 20 150), (6 4 50, 10 20 100))";
  std::string expectedXYM_preserve =
      "MULTILINESTRING M((1 4 100, 4 9 200, 4 16 300, -4 20 150), (6 4 50, 4 9 "
      "200, 4 16 300, 10 20 100))";

  testSimplification(sourceXYM, expectedXYM_noPreserve, 5, false);
  testSimplification(sourceXYM, expectedXYM_preserve, 5, true);

  // Test XYZM
  std::string sourceXYZM = "MULTILINESTRING ZM((1 4 10 100, 4 9 20 200, 4 12 "
                           "15 150, 4 16 30 300, 2 19 25 250, -4 20 15 150), "
                           "(6 4 5 50, 4 9 15 150, 4 12 25 250, 4 16 20 200, 5 "
                           "19 15 150, 10 20 10 100))";
  std::string expectedXYZM_noPreserve =
      "MULTILINESTRING ZM((1 4 10 100, 2 19 25 250, -4 20 15 150), (6 4 5 50, "
      "10 20 10 100))";
  std::string expectedXYZM_preserve =
      "MULTILINESTRING ZM((1 4 10 100, 4 9 20 200, 4 16 30 300, -4 20 15 150), "
      "(6 4 5 50, 4 9 20 200, 4 16 30 300, 10 20 10 100))";

  testSimplification(sourceXYZM, expectedXYZM_noPreserve, 5, false);
  testSimplification(sourceXYZM, expectedXYZM_preserve, 5, true);
}

BOOST_AUTO_TEST_CASE(testSimplify_Triangle)
{
  // XY Triangle
  testSimplification("TRIANGLE((0 0, 10 0, 0 10, 0 0))",
                     "TRIANGLE((0 0, 10 0, 0 10, 0 0))", 1, false);

  // XYZ Triangle
  testSimplification("TRIANGLE Z((0 0 0, 10 0 5, 0 10 10, 0 0 0))",
                     "TRIANGLE Z((0 0 0, 10 0 5, 0 10 10, 0 0 0))", 1, false);

  // XYM Triangle
  testSimplification("TRIANGLE M((0 0 0, 10 0 50, 0 10 100, 0 0 0))",
                     "TRIANGLE M((0 0 0, 10 0 50, 0 10 100, 0 0 0))", 1, false);

  // XYZM Triangle
  testSimplification("TRIANGLE ZM((0 0 0 0, 10 0 5 50, 0 10 10 100, 0 0 0 0))",
                     "TRIANGLE ZM((0 0 0 0, 10 0 5 50, 0 10 10 100, 0 0 0 0))",
                     1, false);

  // Empty Triangle
  testSimplification("TRIANGLE EMPTY", "TRIANGLE EMPTY", 1, false);
}

BOOST_AUTO_TEST_CASE(testSimplify_Polygon_Dimensions)
{
  // XY Polygon
  testSimplification(
      "POLYGON((0 0, 10 0, 10 10, 0 10, 0 0), (2 2, 8 2, 8 8, 2 8, 2 2))",
      "POLYGON((0 0, 10 0, 10 10, 0 10, 0 0), (2 2, 8 2, 8 8, 2 8, 2 2))", 1,
      false);

  // Test XYZ
  std::string sourceXYZ = "POLYGON Z((0 0 0, 10 0 5, 10 10 10, 0 10 5, 0 0 0), "
                          "(2 2 1, 8 2 3, 8 8 5, 2 8 3, 2 2 1))";
  std::string expectedXYZ =
      "POLYGON Z((0 0 0, 10 0 5, 10 10 10, 0 10 5, 0 0 0), "
      "(2 2 1, 8 2 3, 8 8 5, 2 8 3, 2 2 1))";
  testSimplification(sourceXYZ, expectedXYZ, 1, false);

  // Test XYM
  std::string sourceXYM =
      "POLYGON M((0 0 0, 10 0 50, 10 10 100, 0 10 50, 0 0 0), "
      "(2 2 10, 8 2 30, 8 8 50, 2 8 30, 2 2 10))";
  std::string expectedXYM =
      "POLYGON M((0 0 0, 10 0 50, 10 10 100, 0 10 50, 0 0 0), "
      "(2 2 10, 8 2 30, 8 8 50, 2 8 30, 2 2 10))";
  testSimplification(sourceXYM, expectedXYM, 1, false);

  // Test XYZM
  std::string sourceXYZM =
      "POLYGON ZM((0 0 0 0, 10 0 5 50, 10 10 10 100, 0 10 5 50, 0 0 0 0), "
      "(2 2 1 10, 8 2 3 30, 8 8 5 50, 2 8 3 30, 2 2 1 10))";
  std::string expectedXYZM =
      "POLYGON ZM((0 0 0 0, 10 0 5 50, 10 10 10 100, 0 10 5 50, 0 0 0 0), "
      "(2 2 1 10, 8 2 3 30, 8 8 5 50, 2 8 3 30, 2 2 1 10))";
  testSimplification(sourceXYZM, expectedXYZM, 1, false);

  // Empty Polygon
  testSimplification("POLYGON EMPTY", "POLYGON EMPTY", 1, false);
}

BOOST_AUTO_TEST_CASE(testSimplify_MultiPolygon_Dimensions)
{
  // XY MultiPolygon
  testSimplification("MULTIPOLYGON(((0 0, 5 0, 5 5, 0 5, 0 0)), ((10 0, 15 0, "
                     "15 5, 10 5, 10 0)))",
                     "MULTIPOLYGON(((0 0, 5 0, 5 5, 0 5, 0 0)), ((10 0, 15 0, "
                     "15 5, 10 5, 10 0)))",
                     1, false);

  // Test XYZ
  std::string sourceXYZ =
      "MULTIPOLYGON Z(((0 0 0, 5 0 2, 5 5 5, 0 5 2, 0 0 0)), "
      "((10 0 0, 15 0 2, 15 5 5, 10 5 2, 10 0 0)))";
  std::string expectedXYZ =
      "MULTIPOLYGON Z(((0 0 0, 5 0 2, 5 5 5, 0 5 2, 0 0 0)), "
      "((10 0 0, 15 0 2, 15 5 5, 10 5 2, 10 0 0)))";
  testSimplification(sourceXYZ, expectedXYZ, 1, false);
  testSimplification(sourceXYZ, expectedXYZ, 1, true);

  // Test XYM
  std::string sourceXYM =
      "MULTIPOLYGON M(((0 0 0, 5 0 20, 5 5 50, 0 5 20, 0 0 0)), "
      "((10 0 0, 15 0 20, 15 5 50, 10 5 20, 10 0 0)))";
  std::string expectedXYM =
      "MULTIPOLYGON M(((0 0 0, 5 0 20, 5 5 50, 0 5 20, 0 0 0)), "
      "((10 0 0, 15 0 20, 15 5 50, 10 5 20, 10 0 0)))";
  std::string expectedXYM_preserve =
      "MULTIPOLYGON M(((0 0 0, 5 0 20, 5 5 50, 0 5 20, 0 0 0)), "
      "((10 0 0, 15 0 20, 15 5 50, 10 5 20, 10 0 0)))";
  testSimplification(sourceXYM, expectedXYM, 1, false);
  testSimplification(sourceXYM, expectedXYM_preserve, 1, true);

  // Test XYZM
  std::string sourceXYZM =
      "MULTIPOLYGON ZM(((0 0 0 0, 5 0 2 20, 5 5 5 50, 0 5 2 20, 0 0 0 0)), "
      "((10 0 0 0, 15 0 2 20, 15 5 5 50, 10 5 2 20, 10 0 0 0)))";
  std::string expectedXYZM =
      "MULTIPOLYGON ZM(((0 0 0 0, 5 0 2 20, 5 5 5 50, 0 5 2 20, 0 0 0 0)), "
      "((10 0 0 0, 15 0 2 20, 15 5 5 50, 10 5 2 20, 10 0 0 0)))";
  testSimplification(sourceXYZM, expectedXYZM, 1, false);
  testSimplification(sourceXYZM, expectedXYZM, 1, true);

  // Empty MultiPolygon
  testSimplification("MULTIPOLYGON EMPTY", "MULTIPOLYGON EMPTY", 1, false);
}

BOOST_AUTO_TEST_CASE(testSimplify_TriangulatedSurface)
{
  // XYZ TriangulatedSurface (TIN)
  testSimplification(
      "TIN Z(((0 0 0, 1 0 0, 0 1 0, 0 0 0)), ((1 0 0, 1 1 0, 0 1 0, 1 0 0)))",
      "TIN Z(((0 0 0, 1 0 0, 0 1 0, 0 0 0)), ((1 0 0, 1 1 0, 0 1 0, 1 0 0)))",
      1, false);

  // XYZM TriangulatedSurface
  testSimplification("TIN ZM(((0 0 0 0, 1 0 0 1, 0 1 0 2, 0 0 0 0)), ((1 0 0 "
                     "1, 1 1 0 3, 0 1 0 2, 1 0 0 1)))",
                     "TIN ZM(((0 0 0 0, 1 0 0 1, 0 1 0 2, 0 0 0 0)), ((1 0 0 "
                     "1, 1 1 0 3, 0 1 0 2, 1 0 0 1)))",
                     1, false);

  // Empty TriangulatedSurface
  testSimplification("TIN Z EMPTY", "TIN Z EMPTY", 1, false);
}

BOOST_AUTO_TEST_CASE(testSimplify_PolyhedralSurface)
{
  // XY PolyhedralSurface
  testSimplification("POLYHEDRALSURFACE (((0 0, 0 5, 5 5, 4.5 2.5, 5 0, 0 0)), "
                     "((5 0, 4.5 2.5, 5 5, 10 5, 10 0, 5 0)))",
                     "POLYHEDRALSURFACE (((0 0, 0 5, 5 5, 5 0, 0 0)), ((5 0, 5 "
                     "5, 10 5, 10 0, 5 0)))",
                     2, false);
  testSimplification("POLYHEDRALSURFACE (((0 0, 0 5, 5 5, 4.5 2.5, 5 0, 0 0)), "
                     "((5 0, 4.5 2.5, 5 5, 10 5, 10 0, 5 0)))",
                     "POLYHEDRALSURFACE (((0 0, 0 5, 5 5, 5 0, 0 0)), ((5 0, 5 "
                     "5, 10 5, 10 0, 5 0)))",
                     2, true);
  // XYZ PolyhedralSurface
  testSimplification(
      "POLYHEDRALSURFACE Z(((0 0 10, 0 5 15, 5 5 20, 4.5 2.5 25, 5 0 30, 0 0 "
      "10)), ((5 0 30, 4.5 2.5 25, 5 5 20, 10 5 35, 10 0 40, 5 0 30)))",
      "POLYHEDRALSURFACE Z(((0 0 10, 0 5 15, 5 5 20, 5 0 30, 0 0 10)), ((5 0 "
      "30, 5 5 20, 10 5 35, 10 0 40, 5 0 30)))",
      2, false);
  testSimplification(
      "POLYHEDRALSURFACE Z(((0 0 10, 0 5 15, 5 5 20, 4.5 2.5 25, 5 0 30, 0 0 "
      "10)), ((5 0 30, 4.5 2.5 25, 5 5 20, 10 5 35, 10 0 40, 5 0 30)))",
      "POLYHEDRALSURFACE Z(((0 0 10, 0 5 15, 5 5 20, 5 0 30, 0 0 10)), ((5 0 "
      "30, 5 5 20, 10 5 35, 10 0 40, 5 0 30)))",
      2, true);

  // XYM PolyhedralSurface
  testSimplification(
      "POLYHEDRALSURFACE M(((0 0 100, 0 5 150, 5 5 200, 4.5 2.5 250, 5 0 300, "
      "0 0 100)), ((5 0 300, 4.5 2.5 250, 5 5 200, 10 5 350, 10 0 400, 5 0 "
      "300)))",
      "POLYHEDRALSURFACE M(((0 0 100, 0 5 150, 5 5 200, 5 0 300, 0 0 100)), "
      "((5 0 300, 5 5 200, 10 5 350, 10 0 400, 5 0 300)))",
      2, false);
  testSimplification(
      "POLYHEDRALSURFACE M(((0 0 100, 0 5 150, 5 5 200, 4.5 2.5 250, 5 0 300, "
      "0 0 100)), ((5 0 300, 4.5 2.5 250, 5 5 200, 10 5 350, 10 0 400, 5 0 "
      "300)))",
      "POLYHEDRALSURFACE M(((0 0 100, 0 5 150, 5 5 200, 5 0 300, 0 0 100)), "
      "((5 0 300, 5 5 200, 10 5 350, 10 0 400, 5 0 300)))",
      2, true);

  // XYZM PolyhedralSurface
  testSimplification(
      "POLYHEDRALSURFACE ZM(((0 0 10 100, 0 5 15 150, 5 5 20 200, 4.5 2.5 25 "
      "250, 5 0 30 300, 0 0 10 100)), ((5 0 30 300, 4.5 2.5 25 250, 5 5 20 "
      "200, 10 5 35 350, 10 0 40 400, 5 0 30 300)))",
      "POLYHEDRALSURFACE ZM(((0 0 10 100, 0 5 15 150, 5 5 20 200, 5 0 30 300, "
      "0 0 10 100)), ((5 0 30 300, 5 5 20 200, 10 5 35 350, 10 0 40 400, 5 0 "
      "30 300)))",
      2, false);
  testSimplification(
      "POLYHEDRALSURFACE ZM(((0 0 10 100, 0 5 15 150, 5 5 20 200, 4.5 2.5 25 "
      "250, 5 0 30 300, 0 0 10 100)), ((5 0 30 300, 4.5 2.5 25 250, 5 5 20 "
      "200, 10 5 35 350, 10 0 40 400, 5 0 30 300)))",
      "POLYHEDRALSURFACE ZM(((0 0 10 100, 0 5 15 150, 5 5 20 200, 5 0 30 300, "
      "0 0 10 100)), ((5 0 30 300, 5 5 20 200, 10 5 35 350, 10 0 40 400, 5 0 "
      "30 300)))",
      2, true);
}

BOOST_AUTO_TEST_CASE(testSimplify_Solid)
{
  // XYZ Solid
  testSimplification("SOLID Z((((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0), (0 0 0, 0 "
                     "0 1, 0 1 1, 0 1 0, 0 0 0))))",
                     "SOLID Z((((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0), (0 0 0, 0 "
                     "0 1, 0 1 1, 0 1 0, 0 0 0))))",
                     1, false);

  // XYZM Solid
  testSimplification("SOLID ZM((((0 0 0 0, 0 1 0 1, 1 1 0 2, 1 0 0 3, 0 0 0 "
                     "0), (0 0 0 0, 0 0 1 4, 0 1 1 5, 0 1 0 1, 0 0 0 0))))",
                     "SOLID ZM((((0 0 0 0, 0 1 0 1, 1 1 0 2, 1 0 0 3, 0 0 0 "
                     "0), (0 0 0 0, 0 0 1 4, 0 1 1 5, 0 1 0 1, 0 0 0 0))))",
                     1, false);

  // Empty Solid
  testSimplification("SOLID Z EMPTY", "SOLID Z EMPTY", 1, false);
}

BOOST_AUTO_TEST_CASE(testSimplify_GeometryCollection_Dimensions)
{
  // XY GeometryCollection
  testSimplification(
      "GEOMETRYCOLLECTION(POINT(1 2), LINESTRING(1 4, 4 9, 4 16, 2 19, -4 20))",
      "GEOMETRYCOLLECTION(POINT(1 2), LINESTRING(1 4, -4 20))", 5, false);

  // Test XY
  std::string sourceXY =
      "GEOMETRYCOLLECTION(LINESTRING(1 4, 4 9, 4 12, 4 16, 2 19, -4 20), "
      "POLYGON((0 0, 10 0, 10 10, 0 10, 0 0)))";
  std::string expectedXY_noPreserve =
      "GEOMETRYCOLLECTION(LINESTRING(1 4, -4 20), "
      "POLYGON((0 0, 10 10, 0 10, 0 0)))";
  std::string expectedXY_preserve =
      "GEOMETRYCOLLECTION(LINESTRING(1 4, 4 10, -4 20), "
      "POLYGON((0 0, 10 0, 10 10, 4 10, 0 10, 0 0)))";
  testSimplification(sourceXY, expectedXY_noPreserve, 5, false);
  testSimplification(sourceXY, expectedXY_preserve, 5, true);

  // Test XYZ
  std::string sourceXYZ =
      "GEOMETRYCOLLECTION Z(LINESTRING Z(1 4 10, 4 9 20, 4 12 15, 4 16 30, 2 "
      "19 25, -4 20 15), "
      "POLYGON Z((0 0 0, 10 0 5, 10 10 10, 0 10 5, 0 0 0)))";
  std::string expectedXYZ_noPreserve =
      "GEOMETRYCOLLECTION Z(LINESTRING Z(1 4 10, -4 20 15), "
      "POLYGON Z((0 0 0, 10 10 10, 0 10 5, 0 0 0)))";
  std::string expectedXYZ_preserve =
      "GEOMETRYCOLLECTION Z (LINESTRING Z (1.00 4.00 10.00,4.00 10.00 "
      "18.33,-4.00 20.00 15.00),"
      "POLYGON Z ((0.00 0.00 0.00,10.00 0.00 5.00,10.00 10.00 10.00,4.00 10.00 "
      "18.33,0.00 10.00 5.00,0.00 0.00 0.00)))";
  testSimplification(sourceXYZ, expectedXYZ_noPreserve, 5, false);
  testSimplification(sourceXYZ, expectedXYZ_preserve, 5, true);

  // Test XYM
  std::string sourceXYM =
      "GEOMETRYCOLLECTION M(LINESTRING M(1 4 100, 4 9 200, 4 12 150, 4 16 300, "
      "2 19 250, -4 20 150), "
      "POLYGON M((0 0 0, 10 0 50, 10 10 100, 0 10 50, 0 0 0)))";
  std::string expectedXYM_noPreserve =
      "GEOMETRYCOLLECTION M(LINESTRING M(1 4 100, -4 20 150), "
      "POLYGON M((0 0 0, 10 10 100, 0 10 50, 0 0 0)))";
  std::string expectedXYM_preserve =
      "GEOMETRYCOLLECTION M (LINESTRING M (1.00 4.00 100.00,4.00 10.00 "
      "183.33,-4.00 20.00 150.00),"
      "POLYGON M ((0.00 0.00 0.00,10.00 0.00 50.00,10.00 10.00 100.00,4.00 "
      "10.00 183.33,0.00 10.00 50.00,0.00 0.00 0.00)))";

  testSimplification(sourceXYM, expectedXYM_noPreserve, 5, false);
  testSimplification(sourceXYM, expectedXYM_preserve, 5, true);

  // Test XYZM
  std::string sourceXYZM =
      "GEOMETRYCOLLECTION ZM(LINESTRING ZM(1 4 10 100, 4 9 20 200, 4 12 15 "
      "150, 4 16 30 300, 2 19 25 250, -4 20 15 150), "
      "POLYGON ZM((0 0 0 0, 10 0 5 50, 10 10 10 100, 0 10 5 50, 0 0 0 0)))";
  std::string expectedXYZM_noPreserve =
      "GEOMETRYCOLLECTION ZM (LINESTRING ZM (1.00 4.00 10.00 100.00,-4.00 "
      "20.00 15.00 150.00),"
      "POLYGON ZM ((0.00 0.00 0.00 0.00,10.00 10.00 10.00 100.00,0.00 10.00 "
      "5.00 50.00,0.00 0.00 0.00 0.00)))";
  std::string expectedXYZM_preserve =
      "GEOMETRYCOLLECTION ZM (LINESTRING ZM (1.00 4.00 10.00 100.00,4.00 10.00 "
      "18.33 183.33,-4.00 20.00 15.00 150.00),"
      "POLYGON ZM ((0.00 0.00 0.00 0.00,10.00 0.00 5.00 50.00,10.00 10.00 "
      "10.00 100.00,4.00 10.00 18.33 183.33,0.00 10.00 5.00 50.00,0.00 0.00 "
      "0.00 0.00)))";

  testSimplification(sourceXYZM, expectedXYZM_noPreserve, 5, false);
  testSimplification(sourceXYZM, expectedXYZM_preserve, 5, true);

  // Empty GeometryCollection
  testSimplification("GEOMETRYCOLLECTION EMPTY", "GEOMETRYCOLLECTION EMPTY", 5,
                     false);

  // Test XYZM Topology
  std::string sourceTopo =
      "GEOMETRYCOLLECTION ZM(LINESTRING ZM(-1 -1 3 4, 0 0 10 100, 1 1 20 200, "
      "0 2 15 150, 0 5 30 300, 2 19 25 250, -4 20 15 150), "
      "POLYGON ZM((0 0 10 100, 1 1 20 200, 0 2 15 150, 0 5 30 300, 2 19 25 "
      "250, -4 20 15 150, 0 0 10 100)))";
  std::string expectedTopo_noPreserve =
      "GEOMETRYCOLLECTION ZM (LINESTRING ZM (-1.00 -1.00 3.00 4.00,2.00 19.00 "
      "25.00 250.00,-4.00 20.00 15.00 150.00),"
      "POLYGON ZM ((0.00 0.00 10.00 100.00,2.00 19.00 25.00 250.00,-4.00 20.00 "
      "15.00 150.00,0.00 0.00 10.00 100.00)))";
  std::string expectedTopo_preserve =
      "GEOMETRYCOLLECTION ZM (LINESTRING ZM (-1.00 -1.00 3.00 4.00,0.00 0.00 "
      "10.00 100.00,2.00 19.00 25.00 250.00,-4.00 20.00 15.00 150.00),"
      "POLYGON ZM ((0.00 0.00 10.00 100.00,2.00 19.00 25.00 250.00,-4.00 20.00 "
      "15.00 150.00,0.00 0.00 10.00 100.00)))";

  testSimplification(sourceTopo, expectedTopo_noPreserve, 2, false);
  testSimplification(sourceTopo, expectedTopo_preserve, 2, true);

  // Test Collection with Empty
  std::string sourceWithEmpty =
      "GEOMETRYCOLLECTION (POINT (0 0), LINESTRING (0 0, 1 1, 1 2, 5 5), "
      "LINESTRING EMPTY, LINESTRING (1 2, 2 3))";
  std::string expectedWithEmpty_noPreserve =
      "GEOMETRYCOLLECTION (POINT (0 0), LINESTRING (0 0, 5 5), LINESTRING "
      "EMPTY, LINESTRING (1 2, 2 3))";
  std::string expectedWithEmpty_preserve =
      "GEOMETRYCOLLECTION (POINT (0 0), LINESTRING (0 0, 1 2, 5 5), LINESTRING "
      "EMPTY, LINESTRING (1 2, 2 3))";
  testSimplification(sourceWithEmpty, expectedWithEmpty_noPreserve, 4, false);
  testSimplification(sourceWithEmpty, expectedWithEmpty_preserve, 4, true);
}

BOOST_AUTO_TEST_SUITE_END()
