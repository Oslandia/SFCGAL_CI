// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <CGAL/number_utils.h>
#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/LineString.h>
#include <SFCGAL/MultiLineString.h>
#include <SFCGAL/MultiPoint.h>
#include <SFCGAL/MultiPolygon.h>
#include <SFCGAL/MultiSolid.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Solid.h>
#include <SFCGAL/Triangle.h>
#include <SFCGAL/TriangulatedSurface.h>
#include <SFCGAL/algorithm/rotate.h>
#include <boost/test/unit_test.hpp>
#include <cmath>

using namespace SFCGAL;
using namespace SFCGAL::algorithm;

// Utility function to compare points with a tolerance
auto
pointsAreClose(const Point &p1, const Point &p2, double tolerance = 1e-10)
    -> bool
{
  return std::abs(CGAL::to_double(p1.x()) - CGAL::to_double(p2.x())) <
             tolerance &&
         std::abs(CGAL::to_double(p1.y()) - CGAL::to_double(p2.y())) <
             tolerance &&
         std::abs(CGAL::to_double(p1.z()) - CGAL::to_double(p2.z())) <
             tolerance;
}

BOOST_AUTO_TEST_SUITE(SFCGALRotateTest)

// Test 2D rotation around origin
BOOST_AUTO_TEST_CASE(testRotate2DOrigin)
{
  Point point(1.0, 0.0);
  rotate(point, M_PI / 2); // 90 degree rotation
  BOOST_CHECK(pointsAreClose(point, Point(0.0, 1.0)));
}

// Test 2D rotation around a specified point
BOOST_AUTO_TEST_CASE(testRotate2DPoint)
{
  Point point(2.0, 0.0);
  rotate(point, M_PI / 2, Point(1.0, 0.0)); // 90 degree rotation around (1,0)
  BOOST_CHECK(pointsAreClose(point, Point(1.0, 1.0)));
}

// Test 3D rotation around Z-axis
BOOST_AUTO_TEST_CASE(testRotate3DZ)
{
  Point point(1.0, 0.0, 1.0);
  rotate(point, M_PI / 2,
         Kernel::Vector_3(0, 0, 1)); // 90 degree rotation around Z-axis
  BOOST_CHECK(pointsAreClose(point, Point(0.0, 1.0, 1.0)));
}

// Test 3D rotation around arbitrary axis
BOOST_AUTO_TEST_CASE(testRotate3DArbitraryAxis)
{
  Point point(1.0, 0.0, 0.0);
  rotate(point, M_PI,
         Kernel::Vector_3(1, 1, 1)); // 180 degree rotation around (1,1,1)
  BOOST_CHECK(pointsAreClose(point, Point(-0.33, 0.67, 0.67), 1e-2));
}

// Test rotation around X-axis
BOOST_AUTO_TEST_CASE(testRotateX)
{
  Point point(0.0, 1.0, 0.0);
  rotateX(point, M_PI / 2); // 90 degree rotation around X-axis
  BOOST_CHECK(pointsAreClose(point, Point(0.0, 0.0, 1.0)));
}

// Test rotation around Y-axis
BOOST_AUTO_TEST_CASE(testRotateY)
{
  Point point(1.0, 0.0, 0.0);
  rotateY(point, M_PI / 2); // 90 degree rotation around Y-axis
  BOOST_CHECK(pointsAreClose(point, Point(0.0, 0.0, -1.0)));
}

// Test rotation around Z-axis
BOOST_AUTO_TEST_CASE(testRotateZ)
{
  Point point(1.0, 0.0, 0.0);
  rotateZ(point, M_PI / 2); // 90 degree rotation around Z-axis
  BOOST_CHECK(pointsAreClose(point, Point(0.0, 1.0, 0.0)));
}

// Test rotation of a LineString
BOOST_AUTO_TEST_CASE(testRotateLineString)
{
  LineString line;
  line.addPoint(Point(1.0, 0.0));
  line.addPoint(Point(2.0, 0.0));
  rotate(line, M_PI / 2); // 90 degree rotation
  BOOST_CHECK(pointsAreClose(line.pointN(0), Point(0.0, 1.0)));
  BOOST_CHECK(pointsAreClose(line.pointN(1), Point(0.0, 2.0)));
}

BOOST_AUTO_TEST_CASE(testRotatePolygon)
{
  std::vector<Point> points;
  points.emplace_back(0.0, 0.0);
  points.emplace_back(1.0, 0.0);
  points.emplace_back(1.0, 1.0);
  points.emplace_back(0.0, 1.0);
  points.emplace_back(0.0, 0.0);
  LineString exteriorRing(points);
  Polygon    polygon(exteriorRing);
  rotate(polygon, M_PI / 2); // 90 degree rotation

  BOOST_CHECK(
      pointsAreClose(polygon.exteriorRing().pointN(0), Point(0.0, 0.0)));
  BOOST_CHECK(
      pointsAreClose(polygon.exteriorRing().pointN(1), Point(0.0, 1.0)));
  BOOST_CHECK(
      pointsAreClose(polygon.exteriorRing().pointN(2), Point(-1.0, 1.0)));
  BOOST_CHECK(
      pointsAreClose(polygon.exteriorRing().pointN(3), Point(-1.0, 0.0)));
}

// Test rotation of a PolyhedralSurface
BOOST_AUTO_TEST_CASE(testRotatePolyhedralSurface)
{
  PolyhedralSurface  surface;
  std::vector<Point> points1 = {Point(0, 0, 0), Point(1, 0, 0), Point(1, 1, 0),
                                Point(0, 1, 0), Point(0, 0, 0)};
  std::vector<Point> points2 = {Point(0, 0, 0), Point(1, 0, 0), Point(1, 0, 1),
                                Point(0, 0, 1), Point(0, 0, 0)};
  surface.addPatch(Polygon(LineString(points1)));
  surface.addPatch(Polygon(LineString(points2)));

  rotate(surface, M_PI / 2,
         Kernel::Vector_3(1, 0, 0)); // 90 degree rotation around X-axis

  BOOST_CHECK(pointsAreClose(surface.patchN(0).exteriorRing().pointN(0),
                             Point(0, 0, 0)));
  BOOST_CHECK(pointsAreClose(surface.patchN(0).exteriorRing().pointN(1),
                             Point(1, 0, 0)));
  BOOST_CHECK(pointsAreClose(surface.patchN(0).exteriorRing().pointN(2),
                             Point(1, 0, 1)));
  BOOST_CHECK(pointsAreClose(surface.patchN(0).exteriorRing().pointN(3),
                             Point(0, 0, 1)));

  BOOST_CHECK(pointsAreClose(surface.patchN(1).exteriorRing().pointN(0),
                             Point(0, 0, 0)));
  BOOST_CHECK(pointsAreClose(surface.patchN(1).exteriorRing().pointN(1),
                             Point(1, 0, 0)));
  BOOST_CHECK(pointsAreClose(surface.patchN(1).exteriorRing().pointN(2),
                             Point(1, -1, 0)));
  BOOST_CHECK(pointsAreClose(surface.patchN(1).exteriorRing().pointN(3),
                             Point(0, -1, 0)));
}

// Test rotation of a Solid
BOOST_AUTO_TEST_CASE(testRotateSolid)
{
  std::vector<Point> points = {Point(0, 0, 0), Point(1, 0, 0), Point(1, 1, 0),
                               Point(0, 1, 0), Point(0, 0, 1), Point(1, 0, 1),
                               Point(1, 1, 1), Point(0, 1, 1)};
  PolyhedralSurface  shell;
  shell.addPatch(Polygon(
      LineString({points[0], points[1], points[2], points[3], points[0]})));
  shell.addPatch(Polygon(
      LineString({points[4], points[5], points[6], points[7], points[4]})));
  shell.addPatch(Polygon(
      LineString({points[0], points[1], points[5], points[4], points[0]})));
  shell.addPatch(Polygon(
      LineString({points[1], points[2], points[6], points[5], points[1]})));
  shell.addPatch(Polygon(
      LineString({points[2], points[3], points[7], points[6], points[2]})));
  shell.addPatch(Polygon(
      LineString({points[3], points[0], points[4], points[7], points[3]})));

  Solid solid(shell);

  // Rotate 90 degrees around the Z-axis
  rotate(solid, M_PI / 2, Kernel::Vector_3(0, 0, 1));

  BOOST_CHECK(
      pointsAreClose(solid.exteriorShell().patchN(0).exteriorRing().pointN(0),
                     Point(0, 0, 0), 1e-6));
  BOOST_CHECK(
      pointsAreClose(solid.exteriorShell().patchN(0).exteriorRing().pointN(1),
                     Point(0, 1, 0), 1e-6));
  BOOST_CHECK(
      pointsAreClose(solid.exteriorShell().patchN(0).exteriorRing().pointN(2),
                     Point(-1, 1, 0), 1e-6));
  BOOST_CHECK(
      pointsAreClose(solid.exteriorShell().patchN(0).exteriorRing().pointN(3),
                     Point(-1, 0, 0), 1e-6));
}
// Test 2D rotation with a negative angle
BOOST_AUTO_TEST_CASE(testRotate2DNegativeAngle)
{
  Point point(1.0, 0.0);
  rotate(point, -M_PI / 2); // -90 degree rotation
  BOOST_CHECK(pointsAreClose(point, Point(0.0, -1.0)));
}

// Test 3D rotation with a negative angle
BOOST_AUTO_TEST_CASE(testRotate3DNegativeAngle)
{
  Point point(1.0, 0.0, 0.0);
  rotate(point, -M_PI / 2,
         Kernel::Vector_3(0, 0, 1)); // -90 degree rotation around Z-axis
  BOOST_CHECK(pointsAreClose(point, Point(0.0, -1.0, 0.0)));
}

// Test 2D rotation with an angle > 180 degrees
BOOST_AUTO_TEST_CASE(testRotate2DLargeAngle)
{
  Point point(1.0, 0.0);
  rotate(point, 3 * M_PI / 2); // 270 degree rotation
  BOOST_CHECK(pointsAreClose(point, Point(0.0, -1.0)));
}

// Test 3D rotation with an angle > 180 degrees
BOOST_AUTO_TEST_CASE(testRotate3DLargeAngle)
{
  Point point(1.0, 0.0, 0.0);
  rotate(point, 3 * M_PI / 2,
         Kernel::Vector_3(0, 1, 0)); // 270 degree rotation around Y-axis
  BOOST_CHECK(pointsAreClose(point, Point(0.0, 0.0, 1.0)));
}

// Test 2D rotation with an angle > 360 degrees
BOOST_AUTO_TEST_CASE(testRotate2DExtraLargeAngle)
{
  Point point(1.0, 0.0);
  rotate(point, 5 * M_PI / 2); // 450 degree rotation (equivalent to 90 degree)
  BOOST_CHECK(pointsAreClose(point, Point(0.0, 1.0)));
}

// Test 3D rotation with an angle > 360 degrees
BOOST_AUTO_TEST_CASE(testRotate3DExtraLargeAngle)
{
  Point point(1.0, 0.0, 0.0);
  rotate(point, 5 * M_PI / 2,
         Kernel::Vector_3(0, 0, 1)); // 450 degree rotation around Z-axis
  BOOST_CHECK(pointsAreClose(point, Point(0.0, 1.0, 0.0)));
}

// Test rotation of a LineString with a negative angle
BOOST_AUTO_TEST_CASE(testRotateLineStringNegativeAngle)
{
  LineString line;
  line.addPoint(Point(1.0, 0.0));
  line.addPoint(Point(2.0, 0.0));
  rotate(line, -M_PI / 2); // -90 degree rotation
  BOOST_CHECK(pointsAreClose(line.pointN(0), Point(0.0, -1.0)));
  BOOST_CHECK(pointsAreClose(line.pointN(1), Point(0.0, -2.0)));
}

// Test rotation of a Polygon with an angle > 180 degrees
BOOST_AUTO_TEST_CASE(testRotatePolygonLargeAngle)
{
  // Create a square polygon
  std::vector<Point> points;
  points.emplace_back(0.0, 0.0);
  points.emplace_back(1.0, 0.0);
  points.emplace_back(1.0, 1.0);
  points.emplace_back(0.0, 1.0);
  points.emplace_back(0.0, 0.0);
  LineString exteriorRing(points);
  Polygon    polygon(exteriorRing);

  rotate(polygon, 3 * M_PI / 2); // 270 degree rotation

  // Check if the vertices are in the correct position after rotation
  BOOST_CHECK(
      pointsAreClose(polygon.exteriorRing().pointN(0), Point(0.0, 0.0)));
  BOOST_CHECK(
      pointsAreClose(polygon.exteriorRing().pointN(1), Point(0.0, -1.0)));
  BOOST_CHECK(
      pointsAreClose(polygon.exteriorRing().pointN(2), Point(1.0, -1.0)));
  BOOST_CHECK(
      pointsAreClose(polygon.exteriorRing().pointN(3), Point(1.0, 0.0)));
}

// Test rotation of a Solid with an angle > 360 degrees
BOOST_AUTO_TEST_CASE(testRotateSolidExtraLargeAngle)
{
  // Create a cube
  std::vector<Point> points = {Point(0, 0, 0), Point(1, 0, 0), Point(1, 1, 0),
                               Point(0, 1, 0), Point(0, 0, 1), Point(1, 0, 1),
                               Point(1, 1, 1), Point(0, 1, 1)};
  PolyhedralSurface  shell;

  // Define the faces of the cube
  std::vector<Point> face1 = {points[0], points[1], points[2], points[3],
                              points[0]};
  std::vector<Point> face2 = {points[4], points[5], points[6], points[7],
                              points[4]};
  std::vector<Point> face3 = {points[0], points[1], points[5], points[4],
                              points[0]};
  std::vector<Point> face4 = {points[1], points[2], points[6], points[5],
                              points[1]};
  std::vector<Point> face5 = {points[2], points[3], points[7], points[6],
                              points[2]};
  std::vector<Point> face6 = {points[3], points[0], points[4], points[7],
                              points[3]};

  // Add faces to the polyhedral surface
  shell.addPatch(Polygon(LineString(face1)));
  shell.addPatch(Polygon(LineString(face2)));
  shell.addPatch(Polygon(LineString(face3)));
  shell.addPatch(Polygon(LineString(face4)));
  shell.addPatch(Polygon(LineString(face5)));
  shell.addPatch(Polygon(LineString(face6)));

  Solid solid(shell);

  // Rotate 450 degrees around the Z-axis (equivalent to a 90 degree rotation)
  rotate(solid, 5 * M_PI / 2, Kernel::Vector_3(0, 0, 1));

  // Check if the vertices of the first face are in the correct position after
  // rotation
  BOOST_CHECK(
      pointsAreClose(solid.exteriorShell().patchN(0).exteriorRing().pointN(0),
                     Point(0, 0, 0), 1e-6));
  BOOST_CHECK(
      pointsAreClose(solid.exteriorShell().patchN(0).exteriorRing().pointN(1),
                     Point(0, 1, 0), 1e-6));
  BOOST_CHECK(
      pointsAreClose(solid.exteriorShell().patchN(0).exteriorRing().pointN(2),
                     Point(-1, 1, 0), 1e-6));
  BOOST_CHECK(
      pointsAreClose(solid.exteriorShell().patchN(0).exteriorRing().pointN(3),
                     Point(-1, 0, 0), 1e-6));
}
BOOST_AUTO_TEST_SUITE_END()
