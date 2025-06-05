// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

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
#include "SFCGAL/algorithm/convexHull.h"
#include "SFCGAL/io/wkt.h"

using namespace SFCGAL;

// always after CGAL
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_ConvexHullTest)

// algorithm::convexHull

BOOST_AUTO_TEST_CASE(testConvexHull2D_Empty)
{
  GeometryCollection collect;
  collect.addGeometry(Polygon());
  collect.addGeometry(Polygon());
  std::unique_ptr<Geometry> hull(algorithm::convexHull(collect));
  BOOST_CHECK(hull->isEmpty());
}

BOOST_AUTO_TEST_CASE(testConvexHull2D_ColinearProduceLineString)
{
  LineString lineString;
  lineString.addPoint(Point(0.0, 0.0));
  lineString.addPoint(Point(1.0, 1.0));
  lineString.addPoint(Point(2.0, 2.0));

  std::unique_ptr<Geometry> hull(algorithm::convexHull(lineString));
  BOOST_REQUIRE(hull->is<LineString>());
  BOOST_CHECK_EQUAL(hull->as<LineString>().numPoints(), 2U);

  std::string const hullWKT = hull->asText(1);
  BOOST_CHECK((hullWKT == "LINESTRING (0.0 0.0,2.0 2.0)") ||
              (hullWKT == "LINESTRING (2.0 2.0,0.0 0.0)"));
}

BOOST_AUTO_TEST_CASE(testConvexHull2D_Triangle)
{
  std::vector<Point> points;
  points.emplace_back(0.0, 0.0);
  points.emplace_back(0.5, 0.5);
  points.emplace_back(1.0, 0.0);
  points.emplace_back(0.0, 1.0);

  LineString const          lineString(points);
  std::unique_ptr<Geometry> hull(algorithm::convexHull(lineString));
  BOOST_CHECK(hull->is<Triangle>());
}

BOOST_AUTO_TEST_CASE(testConvexHull2D_Polygon)
{
  std::vector<Point> points;
  points.emplace_back(0.0, 0.0);
  points.emplace_back(1.0, 0.0);
  points.emplace_back(1.0, 1.0);
  points.emplace_back(0.0, 1.0);

  LineString const          lineString(points);
  std::unique_ptr<Geometry> hull(algorithm::convexHull(lineString));
  BOOST_CHECK(hull->is<Polygon>());
}

// algorithm::convexHull3D

BOOST_AUTO_TEST_CASE(testConvexHull3D_Empty)
{
  GeometryCollection collect;
  collect.addGeometry(Polygon());
  collect.addGeometry(Polygon());
  std::unique_ptr<Geometry> hull(algorithm::convexHull3D(collect));
  BOOST_CHECK(hull->isEmpty());
}

BOOST_AUTO_TEST_CASE(testConvexHull3D_Point)
{
  Point const               p(1.0, 2.0, 3.0);
  std::unique_ptr<Geometry> hull(algorithm::convexHull3D(p));
  BOOST_CHECK(hull->is<Point>());
  BOOST_CHECK_EQUAL(hull->as<Point>().x(), 1.0);
  BOOST_CHECK_EQUAL(hull->as<Point>().y(), 2.0);
  BOOST_CHECK_EQUAL(hull->as<Point>().z(), 3.0);
}

/*
 * @todo Test if points are collinear
 */
BOOST_AUTO_TEST_CASE(testConvexHull3D_LineStringCollinear)
{
  std::vector<Point> points;
  points.emplace_back(0.0, 0.0, 0.0);
  points.emplace_back(1.0, 1.0, 1.0);
  points.emplace_back(2.0, 2.0, 2.0);
  points.emplace_back(3.0, 3.0, 3.0);

  LineString const          lineString(points);
  std::unique_ptr<Geometry> hull(algorithm::convexHull3D(lineString));
  BOOST_CHECK(hull->is<LineString>());
}

BOOST_AUTO_TEST_CASE(testConvexHull3D_LineStringCoplanar)
{
  std::vector<Point> points;
  points.emplace_back(0.0, 0.0, 1.0);
  points.emplace_back(1.0, 0.0, 1.0);
  points.emplace_back(1.0, 1.0, 1.0);
  points.emplace_back(0.0, 1.0, 1.0);

  LineString const          lineString(points);
  std::unique_ptr<Geometry> hull(algorithm::convexHull3D(lineString));
  BOOST_CHECK(hull->is<PolyhedralSurface>());
  BOOST_CHECK_EQUAL(hull->as<PolyhedralSurface>().numPatches(), 2U);
}

BOOST_AUTO_TEST_CASE(testConvexHull3D_Tetrahedron)
{
  std::vector<Point> points;
  points.emplace_back(0.0, 0.0, 0.0);
  points.emplace_back(1.0, 0.0, 0.0);
  points.emplace_back(0.0, 1.0, 0.0);
  points.emplace_back(0.0, 0.0, 1.0);

  LineString const          lineString(points);
  std::unique_ptr<Geometry> hull(algorithm::convexHull3D(lineString));
  BOOST_CHECK(hull->is<PolyhedralSurface>());
  BOOST_CHECK_EQUAL(hull->as<PolyhedralSurface>().numPatches(), 4U);
}

BOOST_AUTO_TEST_SUITE_END()
