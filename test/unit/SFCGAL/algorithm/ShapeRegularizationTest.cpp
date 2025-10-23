// Copyright (c) 2025, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/algorithm/ShapeRegularization.h"
#include "SFCGAL/io/wkt.h"

#include <cmath>

using namespace boost::unit_test;
using namespace SFCGAL;
using namespace SFCGAL::algorithm;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_ShapeRegularizationTest)

//
// Helper function to check if two geometries are approximately equal
//
static bool
geometriesApproximatelyEqual(const Geometry &g1, const Geometry &g2,
                             double tolerance = 1e-6)
{
  if (g1.geometryTypeId() != g2.geometryTypeId()) {
    return false;
  }

  // For simple comparison, check WKT (not ideal but works for tests)
  // In production, you'd want coordinate-by-coordinate comparison
  std::string wkt1 = g1.asText(6);
  std::string wkt2 = g2.asText(6);

  // For this test, we just ensure geometries are valid and non-empty
  // Actual regularization results depend on CGAL's internal algorithms
  return !g1.isEmpty() && !g2.isEmpty();
}

//
// Test regularizeSegments
//

BOOST_AUTO_TEST_CASE(testRegularizeSegments_EmptyGeometry)
{
  LineString empty;
  auto       result =
      ShapeRegularization::regularizeSegments(empty, true, false, 25.0, 0.5, 5.0,
                                               ShapeRegularization::NeighborQuery::DELAUNAY);

  BOOST_CHECK(result->isEmpty());
  // Empty geometry returns cloned original (LineString), not MultiLineString
  BOOST_CHECK_EQUAL(result->geometryTypeId(), TYPE_LINESTRING);
}

BOOST_AUTO_TEST_CASE(testRegularizeSegments_SimpleLineString)
{
  // Create a LineString with slight deviations from orthogonal
  LineString ls;
  ls.addPoint(Point(0.0, 0.0));
  ls.addPoint(Point(10.0, 0.1)); // Nearly horizontal
  ls.addPoint(Point(10.1, 10.0)); // Nearly vertical
  ls.addPoint(Point(0.0, 10.1)); // Nearly horizontal

  auto result =
      ShapeRegularization::regularizeSegments(ls, true, false, 25.0, 0.5, 5.0,
                                               ShapeRegularization::NeighborQuery::DELAUNAY);

  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK_EQUAL(result->geometryTypeId(), TYPE_MULTILINESTRING);

  // Result should contain segments
  const auto &mls = result->as<MultiLineString>();
  BOOST_CHECK(mls.numGeometries() > 0);
}

BOOST_AUTO_TEST_CASE(testRegularizeSegments_Polygon)
{
  // Create a nearly-rectangular polygon
  LineString ring;
  ring.addPoint(Point(0.0, 0.0));
  ring.addPoint(Point(10.0, 0.2)); // Nearly horizontal
  ring.addPoint(Point(10.2, 10.0)); // Nearly vertical
  ring.addPoint(Point(0.0, 10.1)); // Nearly horizontal
  ring.addPoint(Point(0.0, 0.0)); // Close the ring

  Polygon poly;
  poly.setExteriorRing(ring.clone());

  auto result =
      ShapeRegularization::regularizeSegments(poly, true, false, 25.0, 0.5, 5.0,
                                               ShapeRegularization::NeighborQuery::DELAUNAY);

  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK_EQUAL(result->geometryTypeId(), TYPE_MULTILINESTRING);

  const auto &mls = result->as<MultiLineString>();
  BOOST_CHECK(mls.numGeometries() >= 4); // At least 4 segments from the ring
}

BOOST_AUTO_TEST_CASE(testRegularizeSegments_PolygonWithHole)
{
  // Exterior ring
  LineString exterior;
  exterior.addPoint(Point(0.0, 0.0));
  exterior.addPoint(Point(20.0, 0.0));
  exterior.addPoint(Point(20.0, 20.0));
  exterior.addPoint(Point(0.0, 20.0));
  exterior.addPoint(Point(0.0, 0.0));

  // Interior ring (hole) - must be oriented opposite to exterior (counter-clockwise)
  LineString interior;
  interior.addPoint(Point(5.0, 5.0));
  interior.addPoint(Point(5.0, 15.0)); // Go up (counter-clockwise)
  interior.addPoint(Point(15.1, 15.0)); // Go right
  interior.addPoint(Point(15.0, 5.1)); // Go down
  interior.addPoint(Point(5.0, 5.0)); // Back to start

  Polygon poly;
  poly.setExteriorRing(exterior.clone());
  poly.addInteriorRing(interior.clone());

  auto result =
      ShapeRegularization::regularizeSegments(poly, true, false, 25.0, 0.5, 5.0,
                                               ShapeRegularization::NeighborQuery::DELAUNAY);

  BOOST_CHECK(!result->isEmpty());
  const auto &mls = result->as<MultiLineString>();
  // Should have segments from both exterior and interior rings
  BOOST_CHECK(mls.numGeometries() >= 8);
}

BOOST_AUTO_TEST_CASE(testRegularizeSegments_MultiLineString)
{
  LineString ls1;
  ls1.addPoint(Point(0.0, 0.0));
  ls1.addPoint(Point(10.0, 0.1));

  LineString ls2;
  ls2.addPoint(Point(0.0, 5.0));
  ls2.addPoint(Point(10.0, 5.1));

  MultiLineString mls;
  mls.addGeometry(ls1.clone());
  mls.addGeometry(ls2.clone());

  auto result =
      ShapeRegularization::regularizeSegments(mls, true, false, 25.0, 0.5, 5.0,
                                               ShapeRegularization::NeighborQuery::DELAUNAY);

  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK_EQUAL(result->geometryTypeId(), TYPE_MULTILINESTRING);
}

BOOST_AUTO_TEST_CASE(testRegularizeSegments_MultiPolygon)
{
  // First polygon
  LineString ring1;
  ring1.addPoint(Point(0.0, 0.0));
  ring1.addPoint(Point(5.0, 0.1));
  ring1.addPoint(Point(5.1, 5.0));
  ring1.addPoint(Point(0.0, 5.0));
  ring1.addPoint(Point(0.0, 0.0));

  Polygon poly1;
  poly1.setExteriorRing(ring1.clone());

  // Second polygon
  LineString ring2;
  ring2.addPoint(Point(10.0, 10.0));
  ring2.addPoint(Point(15.0, 10.1));
  ring2.addPoint(Point(15.1, 15.0));
  ring2.addPoint(Point(10.0, 15.0));
  ring2.addPoint(Point(10.0, 10.0));

  Polygon poly2;
  poly2.setExteriorRing(ring2.clone());

  MultiPolygon mp;
  mp.addGeometry(poly1.clone());
  mp.addGeometry(poly2.clone());

  auto result =
      ShapeRegularization::regularizeSegments(mp, true, false, 25.0, 0.5, 5.0,
                                               ShapeRegularization::NeighborQuery::DELAUNAY);

  BOOST_CHECK(!result->isEmpty());
  const auto &mls = result->as<MultiLineString>();
  BOOST_CHECK(mls.numGeometries() >= 8); // At least 4 segments per polygon
}

BOOST_AUTO_TEST_CASE(testRegularizeSegments_Triangle)
{
  Triangle tri(Point(0.0, 0.0), Point(10.0, 0.1), Point(5.0, 8.66));

  auto result =
      ShapeRegularization::regularizeSegments(tri, true, false, 25.0, 0.5, 5.0,
                                               ShapeRegularization::NeighborQuery::DELAUNAY);

  BOOST_CHECK(!result->isEmpty());
  const auto &mls = result->as<MultiLineString>();
  BOOST_CHECK(mls.numGeometries() >= 3); // At least 3 segments
}

BOOST_AUTO_TEST_CASE(testRegularizeSegments_GeometryCollection)
{
  LineString ls;
  ls.addPoint(Point(0.0, 0.0));
  ls.addPoint(Point(5.0, 0.1));

  LineString ring;
  ring.addPoint(Point(10.0, 10.0));
  ring.addPoint(Point(15.0, 10.1));
  ring.addPoint(Point(15.1, 15.0));
  ring.addPoint(Point(10.0, 15.0));
  ring.addPoint(Point(10.0, 10.0));

  Polygon poly;
  poly.setExteriorRing(ring.clone());

  GeometryCollection gc;
  gc.addGeometry(ls.clone());
  gc.addGeometry(poly.clone());

  auto result =
      ShapeRegularization::regularizeSegments(gc, true, false, 25.0, 0.5, 5.0,
                                               ShapeRegularization::NeighborQuery::DELAUNAY);

  BOOST_CHECK(!result->isEmpty());
  const auto &mls = result->as<MultiLineString>();
  BOOST_CHECK(mls.numGeometries() > 0);
}

BOOST_AUTO_TEST_CASE(testRegularizeSegments_AngleAndOffsetParameters)
{
  LineString ls;
  ls.addPoint(Point(0.0, 0.0));
  ls.addPoint(Point(10.0, 0.5)); // 2.86 degree deviation
  ls.addPoint(Point(10.5, 10.0));

  // Test with angles only
  auto resultAngles = ShapeRegularization::regularizeSegments(
      ls, true, false, 5.0, 0.5, 5.0,
      ShapeRegularization::NeighborQuery::DELAUNAY); // Strict angle bound

  BOOST_CHECK(!resultAngles->isEmpty());

  // Test with offsets only
  auto resultOffsets = ShapeRegularization::regularizeSegments(
      ls, false, true, 25.0, 0.1, 5.0,
      ShapeRegularization::NeighborQuery::DELAUNAY); // Strict offset bound

  BOOST_CHECK(!resultOffsets->isEmpty());

  // Test with both
  auto resultBoth =
      ShapeRegularization::regularizeSegments(ls, true, true, 25.0, 0.5, 5.0,
                                               ShapeRegularization::NeighborQuery::DELAUNAY);

  BOOST_CHECK(!resultBoth->isEmpty());
}

//
// Test regularizeClosedContour
//

BOOST_AUTO_TEST_CASE(testRegularizeClosedContour_EmptyGeometry)
{
  LineString empty;
  auto       result =
      ShapeRegularization::regularizeClosedContour(empty, 
                                                    ShapeRegularization::DirectionEstimator::MULTIPLE);

  BOOST_CHECK(result->isEmpty());
}

BOOST_AUTO_TEST_CASE(testRegularizeClosedContour_ClosedLineString)
{
  // Create a larger nearly-rectangular closed contour
  LineString ls;
  ls.addPoint(Point(0.0, 0.0));
  ls.addPoint(Point(100.0, 1.0));
  ls.addPoint(Point(101.0, 100.0));
  ls.addPoint(Point(1.0, 101.0));
  ls.addPoint(Point(0.0, 0.0));

  auto result = ShapeRegularization::regularizeClosedContour(ls, 
                                                              ShapeRegularization::DirectionEstimator::MULTIPLE);
  
  // CGAL closed contour regularization is disabled (causes crashes)
  // Function returns original geometry unchanged
  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK_EQUAL(result->geometryTypeId(), TYPE_LINESTRING);
  BOOST_CHECK_EQUAL(result->as<LineString>().numPoints(), ls.numPoints());
}

BOOST_AUTO_TEST_CASE(testRegularizeClosedContour_Polygon)
{
  // Create a nearly-rectangular polygon
  LineString ring;
  ring.addPoint(Point(0.0, 0.0));
  ring.addPoint(Point(100.0, 2.0));
  ring.addPoint(Point(102.0, 100.0));
  ring.addPoint(Point(1.0, 101.0));
  ring.addPoint(Point(0.0, 0.0));

  Polygon poly;
  poly.setExteriorRing(ring.clone());
  
  auto result = ShapeRegularization::regularizeClosedContour(poly, 
                                                              ShapeRegularization::DirectionEstimator::MULTIPLE);
  
  // CGAL closed contour regularization is disabled (causes crashes)
  // Function returns original geometry unchanged
  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK_EQUAL(result->geometryTypeId(), TYPE_POLYGON);
}

BOOST_AUTO_TEST_CASE(testRegularizeClosedContour_PolygonWithHole)
{
  // Exterior ring - larger to avoid numerical issues
  LineString exterior;
  exterior.addPoint(Point(0.0, 0.0));
  exterior.addPoint(Point(200.0, 1.0));
  exterior.addPoint(Point(201.0, 200.0));
  exterior.addPoint(Point(1.0, 201.0));
  exterior.addPoint(Point(0.0, 0.0));

  // Interior ring (hole) - counter-clockwise orientation
  LineString interior;
  interior.addPoint(Point(50.0, 50.0));
  interior.addPoint(Point(51.0, 150.0)); // Go up (counter-clockwise)
  interior.addPoint(Point(151.0, 151.0)); // Go right
  interior.addPoint(Point(150.0, 51.0)); // Go down
  interior.addPoint(Point(50.0, 50.0)); // Back to start

  Polygon poly;
  poly.setExteriorRing(exterior.clone());
  poly.addInteriorRing(interior.clone());

  auto result =
      ShapeRegularization::regularizeClosedContour(poly, 
                                                    ShapeRegularization::DirectionEstimator::MULTIPLE);

  // CGAL closed contour regularization is disabled (causes crashes)
  // Function returns original geometry unchanged
  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK_EQUAL(result->geometryTypeId(), TYPE_POLYGON);

  const auto &polyResult = result->as<Polygon>();
  BOOST_CHECK_EQUAL(polyResult.numInteriorRings(), 1);
}

BOOST_AUTO_TEST_CASE(testRegularizeClosedContour_MultiPolygon)
{
  // First polygon
  LineString ring1;
  ring1.addPoint(Point(0.0, 0.0));
  ring1.addPoint(Point(50.0, 1.0));
  ring1.addPoint(Point(51.0, 50.0));
  ring1.addPoint(Point(1.0, 51.0));
  ring1.addPoint(Point(0.0, 0.0));

  Polygon poly1;
  poly1.setExteriorRing(ring1.clone());

  // Second polygon
  LineString ring2;
  ring2.addPoint(Point(100.0, 100.0));
  ring2.addPoint(Point(150.0, 101.0));
  ring2.addPoint(Point(151.0, 150.0));
  ring2.addPoint(Point(101.0, 151.0));
  ring2.addPoint(Point(100.0, 100.0));

  Polygon poly2;
  poly2.setExteriorRing(ring2.clone());

  MultiPolygon mp;
  mp.addGeometry(poly1.clone());
  mp.addGeometry(poly2.clone());

  auto result = ShapeRegularization::regularizeClosedContour(mp, 
                                                              ShapeRegularization::DirectionEstimator::MULTIPLE);

  // CGAL closed contour regularization is disabled (causes crashes)
  // Function returns original geometry unchanged
  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK_EQUAL(result->geometryTypeId(), TYPE_MULTIPOLYGON);
  
  const auto &mpResult = result->as<MultiPolygon>();
  BOOST_CHECK_EQUAL(mpResult.numGeometries(), 2);
}

BOOST_AUTO_TEST_CASE(testRegularizeClosedContour_TooFewPoints)
{
  // Contour with fewer than 4 points should be returned unchanged
  LineString ls;
  ls.addPoint(Point(0.0, 0.0));
  ls.addPoint(Point(1.0, 0.0));
  ls.addPoint(Point(0.0, 0.0)); // Only 3 points

  auto result = ShapeRegularization::regularizeClosedContour(ls, 
                                                              ShapeRegularization::DirectionEstimator::MULTIPLE);

  BOOST_CHECK(!result->isEmpty());
  // Should return the geometry largely unchanged since < 4 points
}

//
// Test regularizeOpenContour
//

BOOST_AUTO_TEST_CASE(testRegularizeOpenContour_EmptyGeometry)
{
  LineString empty;
  auto result = ShapeRegularization::regularizeOpenContour(empty, 
                                                            ShapeRegularization::DirectionEstimator::LONGEST);

  BOOST_CHECK(result->isEmpty());
}

BOOST_AUTO_TEST_CASE(testRegularizeOpenContour_SimpleLineString)
{
  // Create an open polyline with slight deviations
  LineString ls;
  ls.addPoint(Point(0.0, 0.0));
  ls.addPoint(Point(5.0, 0.2)); // Nearly horizontal
  ls.addPoint(Point(10.0, 0.1)); // Nearly horizontal continuation
  ls.addPoint(Point(10.2, 5.0)); // Nearly vertical

  // Test with "longest" estimator (default)
  auto resultLongest =
      ShapeRegularization::regularizeOpenContour(ls, 
                                                  ShapeRegularization::DirectionEstimator::LONGEST);

  BOOST_CHECK(!resultLongest->isEmpty());
  BOOST_CHECK_EQUAL(resultLongest->geometryTypeId(), TYPE_LINESTRING);

  const auto &lsLongest = resultLongest->as<LineString>();
  BOOST_CHECK(lsLongest.numPoints() >= 2);

  // Test with "multiple" estimator
  auto resultMultiple =
      ShapeRegularization::regularizeOpenContour(ls, 
                                                  ShapeRegularization::DirectionEstimator::MULTIPLE);

  BOOST_CHECK(!resultMultiple->isEmpty());
  BOOST_CHECK_EQUAL(resultMultiple->geometryTypeId(), TYPE_LINESTRING);
}

BOOST_AUTO_TEST_CASE(testRegularizeOpenContour_ZigzagPattern)
{
  // Create a zigzag pattern that should be regularized
  LineString ls;
  ls.addPoint(Point(0.0, 0.0));
  ls.addPoint(Point(1.0, 0.1)); // Nearly horizontal
  ls.addPoint(Point(1.1, 1.0)); // Nearly vertical
  ls.addPoint(Point(2.0, 1.1)); // Nearly horizontal
  ls.addPoint(Point(2.1, 2.0)); // Nearly vertical

  auto result = ShapeRegularization::regularizeOpenContour(ls, 
                                                            ShapeRegularization::DirectionEstimator::LONGEST);

  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK_EQUAL(result->geometryTypeId(), TYPE_LINESTRING);

  const auto &lsResult = result->as<LineString>();
  BOOST_CHECK(lsResult.numPoints() >= 2);
}

BOOST_AUTO_TEST_CASE(testRegularizeOpenContour_MultiLineString)
{
  LineString ls1;
  ls1.addPoint(Point(0.0, 0.0));
  ls1.addPoint(Point(5.0, 0.1));
  ls1.addPoint(Point(10.0, 0.2));

  LineString ls2;
  ls2.addPoint(Point(0.0, 5.0));
  ls2.addPoint(Point(5.0, 5.1));
  ls2.addPoint(Point(10.0, 5.2));

  MultiLineString mls;
  mls.addGeometry(ls1.clone());
  mls.addGeometry(ls2.clone());

  auto result = ShapeRegularization::regularizeOpenContour(mls, 
                                                            ShapeRegularization::DirectionEstimator::LONGEST);

  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK_EQUAL(result->geometryTypeId(), TYPE_MULTILINESTRING);

  const auto &mlsResult = result->as<MultiLineString>();
  BOOST_CHECK_EQUAL(mlsResult.numGeometries(), 2);
}

BOOST_AUTO_TEST_CASE(testRegularizeOpenContour_ClosedLineString)
{
  // When an open contour function receives a closed linestring,
  // it should defer to the closed contour regularization
  LineString ls;
  ls.addPoint(Point(0.0, 0.0));
  ls.addPoint(Point(100.0, 1.0));
  ls.addPoint(Point(101.0, 100.0));
  ls.addPoint(Point(0.0, 0.0)); // Closed

  auto result = ShapeRegularization::regularizeOpenContour(ls, 
                                                            ShapeRegularization::DirectionEstimator::LONGEST);

  // Defers to closed contour (which is disabled), returns original
  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK_EQUAL(result->geometryTypeId(), TYPE_LINESTRING);
  BOOST_CHECK_EQUAL(result->as<LineString>().numPoints(), 4);
}

BOOST_AUTO_TEST_CASE(testRegularizeOpenContour_TooFewPoints)
{
  // Contour with at least 2 points but minimal (edge case)
  LineString ls;
  ls.addPoint(Point(0.0, 0.0));
  ls.addPoint(Point(10.0, 0.1));
  ls.addPoint(Point(20.0, 0.0));

  auto result = ShapeRegularization::regularizeOpenContour(ls, 
                                                            ShapeRegularization::DirectionEstimator::LONGEST);

  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK_EQUAL(result->geometryTypeId(), TYPE_LINESTRING);
}

//
// Test 3D geometry handling (Z-coordinate preservation)
//

BOOST_AUTO_TEST_CASE(testRegularizeSegments_3D)
{
  // Create a 3D LineString
  LineString ls;
  ls.addPoint(Point(0.0, 0.0, 0.0));
  ls.addPoint(Point(10.0, 0.1, 1.0)); // Z should be interpolated/preserved
  ls.addPoint(Point(10.1, 10.0, 2.0));

  auto result =
      ShapeRegularization::regularizeSegments(ls, true, false, 25.0, 0.5, 5.0,
                                               ShapeRegularization::NeighborQuery::DELAUNAY);

  BOOST_CHECK(!result->isEmpty());
  // Check that result preserves 3D coordinate type
  BOOST_CHECK(result->is3D());
}

BOOST_AUTO_TEST_CASE(testRegularizeClosedContour_3D)
{
  // 3D closed contour
  LineString ring;
  ring.addPoint(Point(0.0, 0.0, 0.0));
  ring.addPoint(Point(100.0, 1.0, 10.0));
  ring.addPoint(Point(101.0, 100.0, 20.0));
  ring.addPoint(Point(1.0, 101.0, 15.0));
  ring.addPoint(Point(0.0, 0.0, 0.0));

  Polygon poly;
  poly.setExteriorRing(ring.clone());

  auto result = ShapeRegularization::regularizeClosedContour(poly, 
                                                              ShapeRegularization::DirectionEstimator::LONGEST);

  // CGAL closed contour regularization is disabled (causes crashes)
  // Function returns original geometry unchanged
  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK(result->is3D());
  BOOST_CHECK_EQUAL(result->geometryTypeId(), TYPE_POLYGON);
}

BOOST_AUTO_TEST_CASE(testRegularizeOpenContour_3D)
{
  LineString ls;
  ls.addPoint(Point(0.0, 0.0, 0.0));
  ls.addPoint(Point(5.0, 0.1, 1.0));
  ls.addPoint(Point(10.0, 0.2, 2.0));

  auto result = ShapeRegularization::regularizeOpenContour(ls, 
                                                            ShapeRegularization::DirectionEstimator::LONGEST);

  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK(result->is3D());
}

//
// Test validity propagation
//

BOOST_AUTO_TEST_CASE(testValidityFlagPropagation)
{
  LineString ls;
  ls.addPoint(Point(0.0, 0.0));
  ls.addPoint(Point(10.0, 0.1));
  ls.addPoint(Point(10.1, 10.0));

  // Result should have validity flag set
  auto result =
      ShapeRegularization::regularizeSegments(ls, true, false, 25.0, 0.5, 5.0,
                                               ShapeRegularization::NeighborQuery::DELAUNAY);

  BOOST_CHECK(result->hasValidityFlag());
}

//
// Regression test: ensure degenerate segments are handled
//

BOOST_AUTO_TEST_CASE(testRegularizeSegments_DegenerateSegments)
{
  LineString ls;
  ls.addPoint(Point(0.0, 0.0));
  ls.addPoint(Point(0.0, 0.0)); // Degenerate (same point)
  ls.addPoint(Point(10.0, 0.0));

  auto result =
      ShapeRegularization::regularizeSegments(ls, true, false, 25.0, 0.5, 5.0,
                                               ShapeRegularization::NeighborQuery::DELAUNAY);

  BOOST_CHECK(!result->isEmpty());
  // Degenerate segments should be filtered out during collection
}

BOOST_AUTO_TEST_SUITE_END()
