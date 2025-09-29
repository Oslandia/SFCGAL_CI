// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/NURBSCurve.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/algorithm/isClosed.h"
#include "SFCGAL/io/wkt.h"

using namespace SFCGAL;
using namespace SFCGAL::algorithm;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_IsClosedTest)

// ----------------------------------------------------------------------------------
// Point tests
// ----------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testPointAlwaysClosed)
{
  Point   point(1.0, 2.0, 3.0);
  Closure result = isClosed(point);
  BOOST_CHECK(result);

  // Empty point
  Point emptyPoint;
  BOOST_CHECK(isClosed(emptyPoint));
}

// ----------------------------------------------------------------------------------
// LineString tests
// ----------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testLineStringClosed)
{
  // Closed LineString (ring)
  LineString line;
  line.addPoint(Point(0, 0));
  line.addPoint(Point(1, 0));
  line.addPoint(Point(1, 1));
  line.addPoint(Point(0, 1));
  line.addPoint(Point(0, 0)); // Closing point

  Closure result = isClosed(line);
  BOOST_CHECK(result);
}

BOOST_AUTO_TEST_CASE(testLineStringOpen)
{
  // Open LineString
  LineString line;
  line.addPoint(Point(0, 0));
  line.addPoint(Point(1, 0));
  line.addPoint(Point(1, 1));

  Closure result = isClosed(line);
  BOOST_CHECK(!result);
  BOOST_CHECK(result.reason().find("not identical") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(testLineStringEmpty)
{
  LineString line;
  BOOST_CHECK(isClosed(line));
}

BOOST_AUTO_TEST_CASE(testLineStringSinglePoint)
{
  LineString line;
  line.addPoint(Point(0, 0));

  Closure result = isClosed(line);
  BOOST_CHECK(!result);
  BOOST_CHECK(result.reason().find("less than 2 points") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(testLineString3DClosed)
{
  // 3D closed LineString
  std::string wkt = "LINESTRING Z(0 0 0, 1 0 0, 1 1 1, 0 1 1, 0 0 0)";
  std::unique_ptr<Geometry> geom(io::readWkt(wkt));

  BOOST_CHECK(isClosed(*geom));
}

// ----------------------------------------------------------------------------------
// Polygon tests
// ----------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testPolygonAlwaysClosed)
{
  // Simple polygon
  std::string               wkt = "POLYGON((0 0, 1 0, 1 1, 0 1, 0 0))";
  std::unique_ptr<Geometry> geom(io::readWkt(wkt));

  BOOST_CHECK(isClosed(*geom));
}

BOOST_AUTO_TEST_CASE(testPolygonWithHole)
{
  // Polygon with hole
  std::string wkt =
      "POLYGON((0 0, 4 0, 4 4, 0 4, 0 0), (1 1, 3 1, 3 3, 1 3, 1 1))";
  std::unique_ptr<Geometry> geom(io::readWkt(wkt));

  BOOST_CHECK(isClosed(*geom));
}

BOOST_AUTO_TEST_CASE(testPolygonEmpty)
{
  Polygon polygon;
  BOOST_CHECK(isClosed(polygon));
}

// ----------------------------------------------------------------------------------
// Triangle tests
// ----------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testTriangleAlwaysClosed)
{
  std::string               wkt = "TRIANGLE((0 0, 1 0, 0 1, 0 0))";
  std::unique_ptr<Geometry> geom(io::readWkt(wkt));

  BOOST_CHECK(isClosed(*geom));
}

// ----------------------------------------------------------------------------------
// PolyhedralSurface tests
// ----------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testPolyhedralSurfaceClosed)
{
  // Closed cube
  std::string wkt = "POLYHEDRALSURFACE Z("
                    "((0 0 0, 1 0 0, 1 1 0, 0 1 0, 0 0 0))," // Bottom
                    "((0 0 1, 0 1 1, 1 1 1, 1 0 1, 0 0 1))," // Top
                    "((0 0 0, 0 0 1, 1 0 1, 1 0 0, 0 0 0))," // Front
                    "((0 1 0, 1 1 0, 1 1 1, 0 1 1, 0 1 0))," // Back
                    "((0 0 0, 0 1 0, 0 1 1, 0 0 1, 0 0 0))," // Left
                    "((1 0 0, 1 0 1, 1 1 1, 1 1 0, 1 0 0))"  // Right
                    ")";

  std::unique_ptr<Geometry> geom(io::readWkt(wkt));

  Closure result = isClosed(*geom);
  BOOST_CHECK(result);
}

BOOST_AUTO_TEST_CASE(testPolyhedralSurfaceOpen)
{
  // Open surface (cube missing one face)
  std::string wkt = "POLYHEDRALSURFACE Z("
                    "((0 0 0, 1 0 0, 1 1 0, 0 1 0, 0 0 0))," // Bottom
                    "((0 0 1, 0 1 1, 1 1 1, 1 0 1, 0 0 1))," // Top
                    "((0 0 0, 0 0 1, 1 0 1, 1 0 0, 0 0 0))," // Front
                    "((0 1 0, 1 1 0, 1 1 1, 0 1 1, 0 1 0))," // Back
                    "((0 0 0, 0 1 0, 0 1 1, 0 0 1, 0 0 0))"  // Left
                    // Missing right face
                    ")";

  std::unique_ptr<Geometry> geom(io::readWkt(wkt));

  Closure result = isClosed(*geom);
  BOOST_CHECK(!result);
  BOOST_CHECK(result.reason().find("boundary") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(testPolyhedralSurfaceEmpty)
{
  PolyhedralSurface phs;
  BOOST_CHECK(isClosed(phs));
}

// ----------------------------------------------------------------------------------
// TriangulatedSurface tests
// ----------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testTriangulatedSurfaceClosed)
{
  // Closed tetrahedron
  std::string wkt =
      "TIN Z("
      "((0 0 0, 1 0 0, 0.5 0.866 0, 0 0 0)),"                // Base
      "((0 0 0, 0.5 0.289 0.816, 1 0 0, 0 0 0)),"            // Face 1
      "((1 0 0, 0.5 0.289 0.816, 0.5 0.866 0, 1 0 0)),"      // Face 2
      "((0.5 0.866 0, 0.5 0.289 0.816, 0 0 0, 0.5 0.866 0))" // Face 3
      ")";

  std::unique_ptr<Geometry> geom(io::readWkt(wkt));

  Closure result = isClosed(*geom);
  BOOST_CHECK(result);
}

BOOST_AUTO_TEST_CASE(testTriangulatedSurfaceOpen)
{
  // Open surface (single triangle)
  std::string wkt = "TIN Z(((0 0 0, 1 0 0, 0.5 0.866 0, 0 0 0)))";

  std::unique_ptr<Geometry> geom(io::readWkt(wkt));

  Closure result = isClosed(*geom);
  BOOST_CHECK(!result);
  BOOST_CHECK(result.reason().find("boundary") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(testTriangulatedSurfaceEmpty)
{
  TriangulatedSurface tin;
  BOOST_CHECK(isClosed(tin));
}

// ----------------------------------------------------------------------------------
// Solid tests
// ----------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testSolidAlwaysClosed)
{
  std::string wkt = "SOLID Z ((((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)),((0 0 1,1 0 "
                    "1,1 1 1,0 1 1,0 0 1)),((0 0 0,0 0 1,0 1 1,0 1 0,0 0 "
                    "0)),((0 1 0,0 1 1,1 1 1,1 1 0,0 1 0)),((1 1 0,1 1 1,1 0 "
                    "1,1 0 0,1 1 0)),((1 0 0,1 0 1,0 0 1,0 0 0,1 0 0))))";

  std::unique_ptr<Geometry> geom(io::readWkt(wkt));

  BOOST_CHECK(isClosed(*geom));
}

BOOST_AUTO_TEST_CASE(testSolidEmpty)
{
  Solid solid;
  BOOST_CHECK(isClosed(solid));
}

// ----------------------------------------------------------------------------------
// MultiPoint tests
// ----------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testMultiPointAlwaysClosed)
{
  std::string               wkt = "MULTIPOINT((0 0), (1 1), (2 2))";
  std::unique_ptr<Geometry> geom(io::readWkt(wkt));

  BOOST_CHECK(isClosed(*geom));
}

BOOST_AUTO_TEST_CASE(testMultiPointEmpty)
{
  MultiPoint multipoint;
  BOOST_CHECK(isClosed(multipoint));
}

// ----------------------------------------------------------------------------------
// MultiLineString tests
// ----------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testMultiLineStringAllClosed)
{
  std::string wkt = "MULTILINESTRING("
                    "(0 0, 1 0, 1 1, 0 1, 0 0)," // Closed
                    "(2 2, 3 2, 3 3, 2 3, 2 2)"  // Closed
                    ")";

  std::unique_ptr<Geometry> geom(io::readWkt(wkt));

  BOOST_CHECK(isClosed(*geom));
}

BOOST_AUTO_TEST_CASE(testMultiLineStringSomeOpen)
{
  std::string wkt = "MULTILINESTRING("
                    "(0 0, 1 0, 1 1, 0 1, 0 0)," // Closed
                    "(2 2, 3 2, 3 3)"            // Open
                    ")";

  std::unique_ptr<Geometry> geom(io::readWkt(wkt));

  Closure result = isClosed(*geom);
  BOOST_CHECK(!result);
  BOOST_CHECK(result.reason().find("LineString 1") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(testMultiLineStringEmpty)
{
  MultiLineString multiline;
  BOOST_CHECK(isClosed(multiline));
}

// ----------------------------------------------------------------------------------
// MultiPolygon tests
// ----------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testMultiPolygonAlwaysClosed)
{
  std::string wkt = "MULTIPOLYGON("
                    "((0 0, 1 0, 1 1, 0 1, 0 0)),"
                    "((2 2, 3 2, 3 3, 2 3, 2 2))"
                    ")";

  std::unique_ptr<Geometry> geom(io::readWkt(wkt));

  BOOST_CHECK(isClosed(*geom));
}

BOOST_AUTO_TEST_CASE(testMultiPolygonEmpty)
{
  MultiPolygon multipoint;
  BOOST_CHECK(isClosed(multipoint));
}

// ----------------------------------------------------------------------------------
// MultiSolid tests
// ----------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testMultiSolidAlwaysClosed)
{
  // Create a simple MultiSolid with one cube using WKT
  std::string wkt = "MULTISOLID Z (((((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)),((0 0 "
                    "1,1 0 1,1 1 1,0 1 1,0 0 1)),((0 0 0,0 0 1,0 1 1,0 1 0,0 0 "
                    "0)),((0 1 0,0 1 1,1 1 1,1 1 0,0 1 0)),((1 1 0,1 1 1,1 0 "
                    "1,1 0 0,1 1 0)),((1 0 0,1 0 1,0 0 1,0 0 0,1 0 0)))))";

  std::unique_ptr<Geometry> geom(io::readWkt(wkt));

  BOOST_CHECK(isClosed(*geom));
}

BOOST_AUTO_TEST_CASE(testMultiSolidEmpty)
{
  MultiSolid multisolid;
  BOOST_CHECK(isClosed(multisolid));
}

// ----------------------------------------------------------------------------------
// GeometryCollection tests
// ----------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testGeometryCollectionAllClosed)
{
  GeometryCollection geometrycollection;

  // Add closed LineString
  LineString linestring;
  linestring.addPoint(Point(0, 0));
  linestring.addPoint(Point(1, 0));
  linestring.addPoint(Point(0, 0));
  geometrycollection.addGeometry(linestring);

  // Add Point (always closed)
  geometrycollection.addGeometry(Point(2, 2));

  // Add Polygon (always closed)
  std::unique_ptr<Geometry> polygon(
      io::readWkt("POLYGON((0 0, 1 0, 1 1, 0 1, 0 0))"));
  geometrycollection.addGeometry(*polygon);

  BOOST_CHECK(isClosed(geometrycollection));
}

BOOST_AUTO_TEST_CASE(testGeometryCollectionSomeOpen)
{
  GeometryCollection geometrycollection;

  // Add open LineString
  LineString linestring;
  linestring.addPoint(Point(0, 0));
  linestring.addPoint(Point(1, 0));
  linestring.addPoint(Point(1, 1));
  geometrycollection.addGeometry(linestring);

  // Add Point (always closed)
  geometrycollection.addGeometry(Point(2, 2));

  Closure result = isClosed(geometrycollection);
  BOOST_CHECK(!result);
  BOOST_CHECK(result.reason().find("Geometry 0") != std::string::npos);
  BOOST_CHECK(result.reason().find("LineString") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(testGeometryCollectionEmpty)
{
  GeometryCollection geometrycollection;
  BOOST_CHECK(isClosed(geometrycollection));
}

BOOST_AUTO_TEST_CASE(testGeometryCollectionNested)
{
  GeometryCollection outer;
  GeometryCollection inner;

  // Inner collection with open LineString
  LineString linestring;
  linestring.addPoint(Point(0, 0));
  linestring.addPoint(Point(1, 1));
  inner.addGeometry(linestring);

  outer.addGeometry(inner);

  Closure result = isClosed(outer);
  BOOST_CHECK(!result);
  BOOST_CHECK(result.reason().find("not closed") != std::string::npos);
}

// ----------------------------------------------------------------------------------
// NURBSCurve tests
// ----------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testNURBSCurveClosedLinear)
{
  // Create a closed NURBS curve (triangle)
  std::vector<Point> controlPoints = {
      Point(0, 0), Point(1, 0), Point(0.5, 1), Point(0, 0) // Closed triangle
  };

  auto curve = NURBSCurve::interpolateCurve(controlPoints, 2);
  BOOST_CHECK(isClosed(*curve));
}

BOOST_AUTO_TEST_CASE(testNURBSCurveOpen)
{
  // Create an open NURBS curve
  std::vector<Point> controlPoints = {
      Point(0, 0), Point(1, 0), Point(2, 1), Point(3, 0) // Open curve
  };

  auto    curve  = NURBSCurve::interpolateCurve(controlPoints, 3);
  Closure result = isClosed(*curve);
  BOOST_CHECK(!result);
  BOOST_CHECK(result.reason().find("not coincident") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(testNURBSCurveClosedCircular)
{
  // Create a closed circular arc (full circle)
  auto curve = NURBSCurve::createCircularArc(Point(0, 0), // center
                                             1.0,         // radius
                                             0.0,         // start angle
                                             2 * M_PI // end angle (full circle)
  );

  // Test that the full circle is properly closed
  Closure result = isClosed(*curve);

  // Verify closure by checking endpoints manually
  Point startPoint = curve->startPoint();
  Point endPoint   = curve->endPoint();

  // Calculate distance between endpoints
  auto dx              = endPoint.x() - startPoint.x();
  auto dy              = endPoint.y() - startPoint.y();
  auto squaredDistance = dx * dx + dy * dy;
  auto distance        = CGAL::sqrt(CGAL::to_double(squaredDistance));

  // The full circle should be closed within numerical tolerance
  // SFCGAL's createCircularArc should produce an exact closed circle
  BOOST_CHECK_MESSAGE(
      result, "Full circle should be closed. Distance between endpoints: "
                  << distance << ". Reason: " << result.reason());

  // Additional verification: endpoints should be very close
  BOOST_CHECK_LT(distance, EPSILON); // Very tight tolerance for exact closure

  // Verify we actually have a circle by checking the start point is at (1,0)
  // for center (0,0), radius 1, start angle 0
  BOOST_CHECK_CLOSE(CGAL::to_double(startPoint.x()), 1.0, EPSILON);
  BOOST_CHECK_CLOSE(CGAL::to_double(startPoint.y()), 0.0, EPSILON);
}

BOOST_AUTO_TEST_CASE(testNURBSCurveOpenCircularArc)
{
  // Create an open circular arc (quarter circle)
  auto curve =
      NURBSCurve::createCircularArc(Point(0, 0), // center
                                    1.0,         // radius
                                    0.0,         // start angle
                                    M_PI / 2     // end angle (quarter circle)
      );

  Closure result = isClosed(*curve);
  BOOST_CHECK(!result);
  BOOST_CHECK(result.reason().find("not coincident") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(testNURBSCurveClosed3D)
{
  // Create a closed 3D NURBS curve
  std::vector<Point> controlPoints = {
      Point(0, 0, 0), Point(1, 0, 1), Point(1, 1, 1), Point(0, 1, 0),
      Point(0, 0, 0) // Closed 3D curve
  };

  auto curve = NURBSCurve::interpolateCurve(controlPoints, 2);
  BOOST_CHECK(isClosed(*curve));
}

BOOST_AUTO_TEST_CASE(testNURBSCurveClosedWithMeasure)
{
  // Create a closed NURBS curve with measure coordinates
  std::vector<Point> controlPoints = {
      Point(0, 0, 0, 0), // XYZ + M
      Point(1, 0, 0, 1), Point(1, 1, 0, 2), Point(0, 1, 0, 3),
      Point(0, 0, 0, 4) // Back to start (XYZ closed, M different)
  };

  auto curve = NURBSCurve::interpolateCurve(controlPoints, 2);
  BOOST_CHECK(isClosed(*curve));
}

BOOST_AUTO_TEST_CASE(testNURBSCurveClosedBSpline)
{
  // Create a closed B-spline
  std::vector<Point> controlPoints = {Point(0, 0), Point(2, 0), Point(2, 2),
                                      Point(0, 2)};

  // Create closed periodic curve
  std::vector<Point> periodicPoints = controlPoints;
  periodicPoints.push_back(controlPoints[0]); // Close explicitly

  auto curve = NURBSCurve::createBSpline(periodicPoints, 3);
  BOOST_CHECK(isClosed(*curve));
}

BOOST_AUTO_TEST_CASE(testNURBSCurveEmptyAlwaysClosed)
{
  // Empty NURBS curve should be considered closed
  NURBSCurve emptyCurve;
  BOOST_CHECK(isClosed(emptyCurve));
}

BOOST_AUTO_TEST_CASE(testNURBSCurveSinglePoint)
{
  // Single control point NURBS curve
  std::vector<Point> controlPoints = {Point(1, 1, 1)};

  auto    curve  = NURBSCurve::createBSpline(controlPoints, 0);
  Closure result = isClosed(*curve);
  BOOST_CHECK(!result); // Single point should be considered open according to
                        // implementation
  BOOST_CHECK(result.reason().find("less than 2 points") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(testNURBSCurveClosedHighDegree)
{
  // Create a closed high-degree NURBS curve
  std::vector<Point> controlPoints = {
      Point(0, 0),  Point(1, 0),  Point(2, 0), Point(3, 1),
      Point(3, 2),  Point(2, 3),  Point(1, 3), Point(0, 3),
      Point(-1, 2), Point(-1, 1), Point(0, 0) // Closed loop
  };

  auto curve = NURBSCurve::interpolateCurve(controlPoints, 5);
  BOOST_CHECK(isClosed(*curve));
}

BOOST_AUTO_TEST_CASE(testNURBSCurveNearlyClosedButOpen)
{
  // Create a NURBS curve that's nearly closed but not quite
  std::vector<Point> controlPoints = {
      Point(0, 0), Point(1, 0), Point(1, 1), Point(0, 1),
      Point(0.001, 0.001) // Very close to start but not identical
  };

  auto    curve  = NURBSCurve::interpolateCurve(controlPoints, 3);
  Closure result = isClosed(*curve);
  BOOST_CHECK(!result);
  BOOST_CHECK(result.reason().find("not coincident") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(testNURBSCurveFromBezierClosed)
{
  // Create a closed Bezier curve
  std::vector<Point> bezierPoints = {Point(0, 0), Point(1, 0), Point(1, 1),
                                     Point(0, 1)};

  auto curve = NURBSCurve::fromBezier(bezierPoints);

  // Bezier curves are typically open, but let's test the closure detection
  Closure result = isClosed(*curve);
  BOOST_CHECK(!result); // Bezier should be open unless explicitly closed
}

// ----------------------------------------------------------------------------------
// Edge cases and error handling
// ----------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testClosureReasonAccess)
{
  // Test Closure class methods
  Closure closed = Closure::closed();
  BOOST_CHECK(closed);
  BOOST_CHECK(closed.reason().empty());

  Closure open = Closure::open("Test reason");
  BOOST_CHECK(!open);
  BOOST_CHECK_EQUAL(open.reason(), "Test reason");
}

BOOST_AUTO_TEST_CASE(testComplexPolyhedralSurface)
{
  // Create a more complex closed surface (pyramid)
  std::string wkt = "POLYHEDRALSURFACE Z("
                    "((0 0 0, 2 0 0, 2 2 0, 0 2 0, 0 0 0))," // Base
                    "((0 0 0, 1 1 2, 2 0 0, 0 0 0)),"        // Front
                    "((2 0 0, 1 1 2, 2 2 0, 2 0 0)),"        // Right
                    "((2 2 0, 1 1 2, 0 2 0, 2 2 0)),"        // Back
                    "((0 2 0, 1 1 2, 0 0 0, 0 2 0))"         // Left
                    ")";

  std::unique_ptr<Geometry> geom(io::readWkt(wkt));

  BOOST_CHECK(isClosed(*geom));
}

BOOST_AUTO_TEST_SUITE_END()
