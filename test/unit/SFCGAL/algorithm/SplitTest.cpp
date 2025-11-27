// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

/**
 * Unit tests for polygon split algorithm using CGAL Arrangement_2.
 */

#include <boost/test/unit_test.hpp>

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/algorithm/area.h"
#include "SFCGAL/algorithm/split.h"

#include <memory>
#include <vector>

using namespace SFCGAL;
using namespace SFCGAL::algorithm;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_SplitTest)

// Helper function to count polygons in a GeometryCollection
auto
countPolygons(const Geometry &geom) -> size_t
{
  if (geom.geometryTypeId() == TYPE_GEOMETRYCOLLECTION) {
    const auto &gc = geom.as<GeometryCollection>();
    size_t      count = 0;
    for (size_t i = 0; i < gc.numGeometries(); ++i) {
      if (gc.geometryN(i).geometryTypeId() == TYPE_POLYGON) {
        count++;
      }
    }
    return count;
  }
  return 0;
}

// Helper to create a simple rectangle
auto
makeRectangle(double x1, double y1, double x2, double y2)
    -> std::unique_ptr<Polygon>
{
  std::vector<Point> points = {Point(x1, y1), Point(x2, y1), Point(x2, y2),
                               Point(x1, y2), Point(x1, y1)};
  return std::make_unique<Polygon>(LineString(points));
}

// Helper to create a LineString from 2 points
auto
makeLine(double x1, double y1, double x2, double y2) -> LineString
{
  return LineString(Point(x1, y1), Point(x2, y2));
}

//=============================================================================
// BASIC SPLIT TESTS
// =============================================================================

BOOST_AUTO_TEST_CASE(test_simple_vertical_split)
{
  auto rect = makeRectangle(0, 0, 10, 6);
  auto line = makeLine(5, -1, 5, 7);

  auto result = split(*rect, line);

  BOOST_CHECK_EQUAL(result->geometryTypeId(), TYPE_GEOMETRYCOLLECTION);
  BOOST_CHECK_EQUAL(countPolygons(*result), 2);

  // Check area preservation
  double originalArea = area(*rect);
  double totalArea    = 0.0;
  const auto &gc       = result->as<GeometryCollection>();
  for (size_t i = 0; i < gc.numGeometries(); ++i) {
    totalArea += area(gc.geometryN(i));
  }
  BOOST_CHECK_CLOSE(originalArea, totalArea, 1e-6);
}

BOOST_AUTO_TEST_CASE(test_simple_horizontal_split)
{
  auto rect = makeRectangle(0, 0, 10, 6);
  auto line = makeLine(-1, 3, 11, 3);

  auto result = split(*rect, line);

  BOOST_CHECK_EQUAL(countPolygons(*result), 2);
}

BOOST_AUTO_TEST_CASE(test_diagonal_split)
{
  auto rect = makeRectangle(0, 0, 10, 6);
  auto line = makeLine(0, 0, 10, 6);

  auto result = split(*rect, line);

  BOOST_CHECK_EQUAL(countPolygons(*result), 2);

  // Check area preservation
  double totalArea = 0.0;
  const auto &gc    = result->as<GeometryCollection>();
  for (size_t i = 0; i < gc.numGeometries(); ++i) {
    totalArea += area(gc.geometryN(i));
  }
  BOOST_CHECK_CLOSE(totalArea, 60.0, 1e-6);
}

// =============================================================================
// NO-SPLIT CASES
// =============================================================================

BOOST_AUTO_TEST_CASE(test_no_intersection)
{
  auto rect = makeRectangle(0, 0, 10, 6);
  auto line = makeLine(20, 0, 20, 10);

  auto result = split(*rect, line);

  BOOST_CHECK_EQUAL(countPolygons(*result), 1);
}

BOOST_AUTO_TEST_CASE(test_single_touch_point)
{
  auto rect = makeRectangle(0, 0, 10, 6);
  auto line = makeLine(0, 0, -5, -5);

  auto result = split(*rect, line);

  BOOST_CHECK_EQUAL(countPolygons(*result), 1);
}

BOOST_AUTO_TEST_CASE(test_line_along_edge)
{
  auto rect = makeRectangle(0, 0, 10, 6);
  auto line = makeLine(0, 0, 10, 0);

  auto result = split(*rect, line);

  BOOST_CHECK_EQUAL(countPolygons(*result), 1);
}

// =============================================================================
// VERTEX INTERSECTION TESTS
// =============================================================================

BOOST_AUTO_TEST_CASE(test_line_through_vertices)
{
  auto rect = makeRectangle(0, 0, 10, 6);
  auto line = makeLine(0, 0, 10, 6);

  auto result = split(*rect, line);

  BOOST_CHECK_EQUAL(countPolygons(*result), 2);
}

BOOST_AUTO_TEST_CASE(test_line_through_one_vertex)
{
  auto rect = makeRectangle(0, 0, 10, 6);
  auto line = makeLine(0, 0, 10, 3);

  auto result = split(*rect, line);

  BOOST_CHECK_EQUAL(countPolygons(*result), 2);
}

// =============================================================================
// POLYGON WITH HOLES
// =============================================================================

BOOST_AUTO_TEST_CASE(test_split_polygon_with_hole)
{
  std::vector<Point> outer = {Point(0, 0),   Point(10, 0), Point(10, 10),
                              Point(0, 10), Point(0, 0)};
  // Hole must have opposite orientation (clockwise) from outer ring
  std::vector<Point> hole  = {Point(3, 3), Point(3, 7), Point(7, 7),
                              Point(7, 3), Point(3, 3)};

  auto polygon = std::make_unique<Polygon>(LineString(outer));
  polygon->addInteriorRing(LineString(hole));

  auto line = makeLine(5, -1, 5, 11);

  auto result = split(*polygon, line);

  BOOST_CHECK_EQUAL(countPolygons(*result), 2);

  // Check area preservation
  double originalArea = 100.0 - 16.0;
  double totalArea    = 0.0;
  const auto &gc       = result->as<GeometryCollection>();
  for (size_t i = 0; i < gc.numGeometries(); ++i) {
    totalArea += area(gc.geometryN(i));
  }
  BOOST_CHECK_CLOSE(originalArea, totalArea, 1e-6);
}

// =============================================================================
// MULTI-SEGMENT LINESTRING TESTS
// =============================================================================

BOOST_AUTO_TEST_CASE(test_multi_segment_linestring)
{
  auto rect = makeRectangle(0, 0, 10, 6);

  std::vector<Point> linePoints = {Point(2, -1), Point(5, 3), Point(8, -1),
                                   Point(8, 7)};
  LineString         line(linePoints);

  auto result = split(*rect, line);

  BOOST_CHECK_GE(countPolygons(*result), 2);
}

// =============================================================================
// PRECISION TESTS
// =============================================================================

BOOST_AUTO_TEST_CASE(test_very_small_polygon)
{
  auto tiny = makeRectangle(0, 0, 0.01, 0.01);
  auto line = makeLine(0.005, -0.001, 0.005, 0.011);

  auto result = split(*tiny, line);

  BOOST_CHECK_EQUAL(countPolygons(*result), 2);
}

BOOST_AUTO_TEST_CASE(test_large_coordinates)
{
  auto large = makeRectangle(1000000, 2000000, 1000100, 2000060);
  auto line  = makeLine(1000050, 1999999, 1000050, 2000061);

  auto result = split(*large, line);

  BOOST_CHECK_EQUAL(countPolygons(*result), 2);
}

// =============================================================================
// EDGE CASES
// =============================================================================

BOOST_AUTO_TEST_CASE(test_empty_polygon)
{
  Polygon emptyPoly;
  auto    line = makeLine(0, 0, 10, 10);

  auto result = split(emptyPoly, line);

  BOOST_CHECK_EQUAL(result->geometryTypeId(), TYPE_GEOMETRYCOLLECTION);
  BOOST_CHECK_EQUAL(result->as<GeometryCollection>().numGeometries(), 0);
}

BOOST_AUTO_TEST_CASE(test_empty_linestring)
{
  auto       rect = makeRectangle(0, 0, 10, 6);
  LineString emptyLine;

  auto result = split(*rect, emptyLine);

  BOOST_CHECK_EQUAL(countPolygons(*result), 1);
}

BOOST_AUTO_TEST_CASE(test_asymmetric_split)
{
  auto rect = makeRectangle(0, 0, 100, 10);
  auto line = makeLine(5, -1, 5, 11);

  auto result = split(*rect, line);

  BOOST_CHECK_EQUAL(countPolygons(*result), 2);

  const auto &gc = result->as<GeometryCollection>();
  std::vector<double> areas;
  for (size_t i = 0; i < gc.numGeometries(); ++i) {
    areas.push_back(area(gc.geometryN(i)));
  }
  std::sort(areas.begin(), areas.end());

  BOOST_CHECK_CLOSE(areas[0], 50.0, 1e-6);
  BOOST_CHECK_CLOSE(areas[1], 950.0, 1e-6);
}

// =============================================================================
// COMPLEX POLYGON SHAPES
// =============================================================================

BOOST_AUTO_TEST_CASE(test_l_shaped_polygon)
{
  std::vector<Point> points = {Point(0, 0),  Point(6, 0),  Point(6, 3),
                               Point(3, 3), Point(3, 6), Point(0, 6),
                               Point(0, 0)};

  LineString lShapeRing(points);
  Polygon    lShape(lShapeRing);
  auto       line = makeLine(4, -1, 4, 7);

  auto result = split(lShape, line);

  BOOST_CHECK_GE(countPolygons(*result), 2);

  double originalArea = area(lShape);
  double totalArea    = 0.0;
  const auto &gc       = result->as<GeometryCollection>();
  for (size_t i = 0; i < gc.numGeometries(); ++i) {
    totalArea += area(gc.geometryN(i));
  }
  BOOST_CHECK_CLOSE(originalArea, totalArea, 1e-6);
}

// =============================================================================
// Z AND M COORDINATE INTERPOLATION TESTS
// =============================================================================

BOOST_AUTO_TEST_CASE(test_split_with_z_interpolation)
{
  std::vector<Point> points = {Point(0, 0, 0), Point(10, 0, 5), Point(10, 6, 10),
                               Point(0, 6, 5), Point(0, 0, 0)};
  LineString         rectRing(points);
  Polygon            rect(rectRing);

  LineString line(Point(5, -1, 0), Point(5, 7, 0));

  auto result = split(rect, line);

  BOOST_CHECK_EQUAL(countPolygons(*result), 2);

  const auto &gc = result->as<GeometryCollection>();
  for (size_t i = 0; i < gc.numGeometries(); ++i) {
    const auto &geom = gc.geometryN(i);
    if (geom.geometryTypeId() == TYPE_POLYGON) {
      const auto &poly = geom.as<Polygon>();
      BOOST_CHECK(poly.is3D());

      const auto &ring = poly.exteriorRing();
      for (size_t j = 0; j < ring.numPoints(); ++j) {
        const auto &point = ring.pointN(j);
        if (std::abs(CGAL::to_double(point.x()) - 5.0) < 1e-6) {
          BOOST_CHECK(point.is3D());
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(test_split_with_m_interpolation)
{
  std::vector<Point> points = {Point(0, 0, 0, 0), Point(10, 0, 0, 100),
                               Point(10, 6, 0, 200), Point(0, 6, 0, 50),
                               Point(0, 0, 0, 0)};
  LineString         rectRing(points);
  Polygon            rect(rectRing);

  LineString line(Point(5, -1, 0, 50), Point(5, 7, 0, 150));

  auto result = split(rect, line);

  BOOST_CHECK_EQUAL(countPolygons(*result), 2);

  const auto &gc = result->as<GeometryCollection>();
  for (size_t i = 0; i < gc.numGeometries(); ++i) {
    const auto &geom = gc.geometryN(i);
    if (geom.geometryTypeId() == TYPE_POLYGON) {
      const auto &poly = geom.as<Polygon>();
      BOOST_CHECK(poly.isMeasured());

      const auto &ring = poly.exteriorRing();
      for (size_t j = 0; j < ring.numPoints(); ++j) {
        const auto &point = ring.pointN(j);
        if (std::abs(CGAL::to_double(point.x()) - 5.0) < 1e-6) {
          BOOST_CHECK(point.isMeasured());
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(test_split_xyz_polygon_with_xy_line)
{
  std::vector<Point> points = {Point(0, 0, 10), Point(10, 0, 20),
                               Point(10, 6, 30), Point(0, 6, 15),
                               Point(0, 0, 10)};
  LineString         rectRing(points);
  Polygon            rect(rectRing);

  LineString line(Point(5, -1), Point(5, 7));

  auto result = split(rect, line);

  BOOST_CHECK_EQUAL(countPolygons(*result), 2);

  const auto &gc = result->as<GeometryCollection>();
  for (size_t i = 0; i < gc.numGeometries(); ++i) {
    const auto &geom = gc.geometryN(i);
    if (geom.geometryTypeId() == TYPE_POLYGON) {
      const auto &poly = geom.as<Polygon>();
      BOOST_CHECK(poly.is3D());

      const auto &ring = poly.exteriorRing();
      for (size_t j = 0; j < ring.numPoints(); ++j) {
        const auto &point = ring.pointN(j);
        BOOST_CHECK(point.is3D());
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(test_split_xy_polygon_with_xyz_line)
{
  std::vector<Point> points = {Point(0, 0), Point(10, 0), Point(10, 6),
                               Point(0, 6), Point(0, 0)};
  LineString         rectRing(points);
  Polygon            rect(rectRing);

  LineString line(Point(5, -1, 100), Point(5, 7, 200));

  auto result = split(rect, line);

  BOOST_CHECK_EQUAL(countPolygons(*result), 2);

  const auto &gc = result->as<GeometryCollection>();
  for (size_t i = 0; i < gc.numGeometries(); ++i) {
    const auto &geom = gc.geometryN(i);
    if (geom.geometryTypeId() == TYPE_POLYGON) {
      const auto &poly = geom.as<Polygon>();
      BOOST_CHECK(poly.is3D());
    }
  }
}

BOOST_AUTO_TEST_CASE(test_split_xyz_polygon_with_xym_line)
{
  std::vector<Point> points = {Point(0, 0, 10), Point(10, 0, 20),
                               Point(10, 6, 30), Point(0, 6, 15),
                               Point(0, 0, 10)};
  LineString         rectRing(points);
  Polygon            rect(rectRing);

  LineString line(Point(5, -1, 0, 50), Point(5, 7, 0, 150));

  auto result = split(rect, line);

  BOOST_CHECK_EQUAL(countPolygons(*result), 2);

  const auto &gc = result->as<GeometryCollection>();
  for (size_t i = 0; i < gc.numGeometries(); ++i) {
    const auto &geom = gc.geometryN(i);
    if (geom.geometryTypeId() == TYPE_POLYGON) {
      const auto &poly = geom.as<Polygon>();
      BOOST_CHECK(poly.is3D());
      BOOST_CHECK(poly.isMeasured());

      const auto &ring = poly.exteriorRing();
      for (size_t j = 0; j < ring.numPoints(); ++j) {
        const auto &point = ring.pointN(j);
        BOOST_CHECK(point.is3D());
        BOOST_CHECK(point.isMeasured());
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(test_split_xyzm_polygon_with_xy_line)
{
  std::vector<Point> points = {Point(0, 0, 10, 0), Point(10, 0, 20, 100),
                               Point(10, 6, 30, 200), Point(0, 6, 15, 50),
                               Point(0, 0, 10, 0)};
  LineString         rectRing(points);
  Polygon            rect(rectRing);

  LineString line(Point(5, -1), Point(5, 7));

  auto result = split(rect, line);

  BOOST_CHECK_EQUAL(countPolygons(*result), 2);

  const auto &gc = result->as<GeometryCollection>();
  for (size_t i = 0; i < gc.numGeometries(); ++i) {
    const auto &geom = gc.geometryN(i);
    if (geom.geometryTypeId() == TYPE_POLYGON) {
      const auto &poly = geom.as<Polygon>();
      BOOST_CHECK(poly.is3D());
      BOOST_CHECK(poly.isMeasured());
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
