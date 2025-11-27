// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

/**
 * Unit tests for polygon split algorithm using CGAL Arrangement_2.
 */

#include <boost/test/unit_test.hpp>

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
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

// Triangle tests
BOOST_AUTO_TEST_CASE(test_split_triangle)
{
  Triangle   triangle(Point(0, 0), Point(10, 0), Point(5, 8));
  LineString line(Point(5, -1), Point(5, 10));

  auto result = split(triangle, line);

  BOOST_CHECK_EQUAL(countPolygons(*result), 2);

  const auto &gc = result->as<GeometryCollection>();
  BOOST_CHECK_EQUAL(gc.numGeometries(), 2);
}

BOOST_AUTO_TEST_CASE(test_split_triangle_3d)
{
  Triangle   triangle(Point(0, 0, 10), Point(10, 0, 20), Point(5, 8, 30));
  LineString line(Point(5, -1), Point(5, 10));

  auto result = split(triangle, line);

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

BOOST_AUTO_TEST_CASE(test_split_triangle_no_intersection)
{
  Triangle   triangle(Point(0, 0), Point(10, 0), Point(5, 8));
  LineString line(Point(20, -1), Point(20, 10));

  auto result = split(triangle, line);

  BOOST_CHECK_EQUAL(countPolygons(*result), 1);
}

// MultiPolygon tests
BOOST_AUTO_TEST_CASE(test_split_multipolygon)
{
  auto poly1 = makeRectangle(0, 0, 10, 6);
  auto poly2 = makeRectangle(15, 0, 25, 6);

  MultiPolygon multiPoly;
  multiPoly.addGeometry(poly1.release());
  multiPoly.addGeometry(poly2.release());

  LineString line(Point(5, -1), Point(5, 7));

  auto result = split(multiPoly, line);

  const auto &gc = result->as<GeometryCollection>();
  BOOST_CHECK_EQUAL(gc.numGeometries(), 3);
}

BOOST_AUTO_TEST_CASE(test_split_multipolygon_3d)
{
  std::vector<Point> points1 = {Point(0, 0, 5), Point(10, 0, 10),
                                Point(10, 6, 15), Point(0, 6, 8),
                                Point(0, 0, 5)};
  auto               poly1   = std::make_unique<Polygon>(LineString(points1));

  std::vector<Point> points2 = {Point(15, 0, 20), Point(25, 0, 25),
                                Point(25, 6, 30), Point(15, 6, 22),
                                Point(15, 0, 20)};
  auto               poly2   = std::make_unique<Polygon>(LineString(points2));

  MultiPolygon multiPoly;
  multiPoly.addGeometry(poly1.release());
  multiPoly.addGeometry(poly2.release());

  LineString line(Point(5, -1), Point(5, 7));

  auto result = split(multiPoly, line);

  const auto &gc = result->as<GeometryCollection>();
  BOOST_CHECK_EQUAL(gc.numGeometries(), 3);

  for (size_t i = 0; i < gc.numGeometries(); ++i) {
    const auto &geom = gc.geometryN(i);
    if (geom.geometryTypeId() == TYPE_POLYGON) {
      const auto &poly = geom.as<Polygon>();
      BOOST_CHECK(poly.is3D());
    }
  }
}

// PolyhedralSurface tests
BOOST_AUTO_TEST_CASE(test_split_polyhedralsurface)
{
  PolyhedralSurface surface;

  auto poly1 = makeRectangle(0, 0, 10, 6);
  auto poly2 = makeRectangle(10, 0, 20, 6);

  surface.addPolygon(poly1.release());
  surface.addPolygon(poly2.release());

  LineString line(Point(5, -1), Point(5, 7));

  auto result = split(surface, line);

  const auto &gc = result->as<GeometryCollection>();
  BOOST_CHECK_EQUAL(gc.numGeometries(), 3);
}

BOOST_AUTO_TEST_CASE(test_split_polyhedralsurface_3d)
{
  PolyhedralSurface surface;

  std::vector<Point> points1 = {Point(0, 0, 5), Point(10, 0, 10),
                                Point(10, 6, 15), Point(0, 6, 8),
                                Point(0, 0, 5)};
  auto               poly1   = std::make_unique<Polygon>(LineString(points1));

  std::vector<Point> points2 = {Point(10, 0, 10), Point(20, 0, 15),
                                Point(20, 6, 20), Point(10, 6, 15),
                                Point(10, 0, 10)};
  auto               poly2   = std::make_unique<Polygon>(LineString(points2));

  surface.addPolygon(poly1.release());
  surface.addPolygon(poly2.release());

  LineString line(Point(5, -1), Point(5, 7));

  auto result = split(surface, line);

  const auto &gc = result->as<GeometryCollection>();
  BOOST_CHECK_EQUAL(gc.numGeometries(), 3);

  for (size_t i = 0; i < gc.numGeometries(); ++i) {
    const auto &geom = gc.geometryN(i);
    if (geom.geometryTypeId() == TYPE_POLYGON) {
      const auto &poly = geom.as<Polygon>();
      BOOST_CHECK(poly.is3D());
    }
  }
}

// TriangulatedSurface tests
BOOST_AUTO_TEST_CASE(test_split_triangulatedsurface)
{
  TriangulatedSurface surface;

  Triangle triangle1(Point(0, 0), Point(10, 0), Point(5, 8));
  Triangle triangle2(Point(5, 8), Point(10, 0), Point(10, 8));

  surface.addTriangle(triangle1);
  surface.addTriangle(triangle2);

  LineString line(Point(5, -1), Point(5, 10));

  auto result = split(surface, line);

  const auto &gc = result->as<GeometryCollection>();
  BOOST_CHECK_EQUAL(gc.numGeometries(), 3);
}

BOOST_AUTO_TEST_CASE(test_split_triangulatedsurface_3d)
{
  TriangulatedSurface surface;

  Triangle triangle1(Point(0, 0, 10), Point(10, 0, 20), Point(5, 8, 30));
  Triangle triangle2(Point(5, 8, 30), Point(10, 0, 20), Point(10, 8, 35));

  surface.addTriangle(triangle1);
  surface.addTriangle(triangle2);

  LineString line(Point(5, -1), Point(5, 10));

  auto result = split(surface, line);

  const auto &gc = result->as<GeometryCollection>();
  BOOST_CHECK_EQUAL(gc.numGeometries(), 3);

  for (size_t i = 0; i < gc.numGeometries(); ++i) {
    const auto &geom = gc.geometryN(i);
    if (geom.geometryTypeId() == TYPE_POLYGON) {
      const auto &poly = geom.as<Polygon>();
      BOOST_CHECK(poly.is3D());
    }
  }
}

// GeometryCollection tests
BOOST_AUTO_TEST_CASE(test_split_geometrycollection)
{
  GeometryCollection collection;

  auto poly1 = makeRectangle(0, 0, 10, 6);
  auto poly2 = makeRectangle(15, 0, 25, 6);

  Triangle triangle(Point(30, 0), Point(40, 0), Point(35, 8));

  collection.addGeometry(poly1.release());
  collection.addGeometry(poly2.release());
  collection.addGeometry(triangle.clone().release());

  LineString line(Point(5, -1), Point(5, 10));

  auto result = split(collection, line);

  const auto &gc = result->as<GeometryCollection>();
  BOOST_CHECK_EQUAL(gc.numGeometries(), 4);
}

BOOST_AUTO_TEST_CASE(test_split_geometrycollection_3d)
{
  GeometryCollection collection;

  std::vector<Point> points1 = {Point(0, 0, 5), Point(10, 0, 10),
                                Point(10, 6, 15), Point(0, 6, 8),
                                Point(0, 0, 5)};
  auto               poly1   = std::make_unique<Polygon>(LineString(points1));

  Triangle triangle(Point(15, 0, 20), Point(25, 0, 25), Point(20, 8, 30));

  collection.addGeometry(poly1.release());
  collection.addGeometry(triangle.clone().release());

  LineString line(Point(5, -1), Point(5, 10));

  auto result = split(collection, line);

  const auto &gc = result->as<GeometryCollection>();
  BOOST_CHECK_EQUAL(gc.numGeometries(), 3);

  for (size_t i = 0; i < gc.numGeometries(); ++i) {
    const auto &geom = gc.geometryN(i);
    if (geom.geometryTypeId() == TYPE_POLYGON) {
      const auto &poly = geom.as<Polygon>();
      BOOST_CHECK(poly.is3D());
    }
  }
}

BOOST_AUTO_TEST_CASE(test_split_empty_geometrycollection)
{
  GeometryCollection emptyCollection;
  LineString         line(Point(5, -1), Point(5, 7));

  auto result = split(emptyCollection, line);

  const auto &gc = result->as<GeometryCollection>();
  BOOST_CHECK_EQUAL(gc.numGeometries(), 1);
  BOOST_CHECK(gc.geometryN(0).is<GeometryCollection>());
}

// Empty geometry tests
BOOST_AUTO_TEST_CASE(test_split_empty_polygon)
{
  Polygon    emptyPolygon;
  LineString line(Point(5, -1), Point(5, 7));

  auto result = split(emptyPolygon, line);

  BOOST_CHECK(result->is<GeometryCollection>());
  const auto &gc = result->as<GeometryCollection>();
  BOOST_CHECK(gc.isEmpty());
}

BOOST_AUTO_TEST_CASE(test_split_polygon_with_empty_line)
{
  auto       rect = makeRectangle(0, 0, 10, 6);
  LineString emptyLine;

  auto result = split(*rect, emptyLine);

  const auto &gc = result->as<GeometryCollection>();
  BOOST_CHECK_EQUAL(gc.numGeometries(), 1);
  BOOST_CHECK(gc.geometryN(0).is<Polygon>());
}

BOOST_AUTO_TEST_CASE(test_split_empty_multipolygon)
{
  MultiPolygon emptyMultiPoly;
  LineString   line(Point(5, -1), Point(5, 7));

  auto result = split(emptyMultiPoly, line);

  const auto &gc = result->as<GeometryCollection>();
  BOOST_CHECK_EQUAL(gc.numGeometries(), 1);
  BOOST_CHECK(gc.geometryN(0).is<MultiPolygon>());
}

BOOST_AUTO_TEST_CASE(test_split_empty_polyhedralsurface)
{
  PolyhedralSurface emptySurface;
  LineString        line(Point(5, -1), Point(5, 7));

  auto result = split(emptySurface, line);

  const auto &gc = result->as<GeometryCollection>();
  BOOST_CHECK_EQUAL(gc.numGeometries(), 1);
  BOOST_CHECK(gc.geometryN(0).is<PolyhedralSurface>());
}

BOOST_AUTO_TEST_CASE(test_split_empty_triangulatedsurface)
{
  TriangulatedSurface emptySurface;
  LineString          line(Point(5, -1), Point(5, 7));

  auto result = split(emptySurface, line);

  const auto &gc = result->as<GeometryCollection>();
  BOOST_CHECK_EQUAL(gc.numGeometries(), 1);
  BOOST_CHECK(gc.geometryN(0).is<TriangulatedSurface>());
}

BOOST_AUTO_TEST_SUITE_END()
