// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
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
#include "SFCGAL/algorithm/extrude.h"
#include "SFCGAL/detail/transform/ForceZ.h"
#include "SFCGAL/io/wkt.h"

using namespace SFCGAL;

// always after CGAL
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_ExtrudeTest)

BOOST_AUTO_TEST_CASE(testExtrudePoint)
{
  Point const               g(0.0, 0.0, 0.0);
  std::unique_ptr<Geometry> ext(algorithm::extrude(g, 0.0, 0.0, 1.0));
  BOOST_CHECK(ext->is<LineString>());
  BOOST_CHECK(ext->as<LineString>().is3D());
  BOOST_CHECK_EQUAL(ext->asText(1), "LINESTRING Z (0.0 0.0 0.0,0.0 0.0 1.0)");
}

BOOST_AUTO_TEST_CASE(testExtrudeLineString)
{
  LineString const          g(Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0));
  std::unique_ptr<Geometry> ext(algorithm::extrude(g, 0.0, 0.0, 1.0));
  BOOST_CHECK(ext->is<PolyhedralSurface>());
  BOOST_CHECK(ext->as<PolyhedralSurface>().is3D());
  BOOST_CHECK_EQUAL(ext->asText(1),
                    "POLYHEDRALSURFACE Z (((0.0 0.0 0.0,1.0 0.0 0.0,1.0 0.0 "
                    "1.0,0.0 0.0 1.0,0.0 0.0 0.0)))");
}

BOOST_AUTO_TEST_CASE(testExtrudeSquare)
{
  std::vector<Point> points;
  points.emplace_back(0.0, 0.0, 0.0);
  points.emplace_back(1.0, 0.0, 0.0);
  points.emplace_back(1.0, 1.0, 0.0);
  points.emplace_back(0.0, 1.0, 0.0);
  points.emplace_back(0.0, 0.0, 0.0);

  LineString const          exteriorRing(points);
  Polygon const             g(exteriorRing);
  std::unique_ptr<Geometry> ext(algorithm::extrude(g, 0.0, 0.0, 1.0));
  BOOST_CHECK(ext->is<Solid>());
  BOOST_CHECK_EQUAL(ext->as<Solid>().numShells(), 1U);
  BOOST_CHECK_EQUAL(ext->as<Solid>().exteriorShell().numPatches(), 6U);
}

BOOST_AUTO_TEST_CASE(testExtrudePolyhedral)
{
  std::unique_ptr<Geometry> const g =
      io::readWkt("POLYHEDRALSURFACE (((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)))");

  std::unique_ptr<Geometry> ext = algorithm::extrude(*g, 0.0, 0.0, 1.0);
  BOOST_CHECK(ext->is<Solid>());
  BOOST_CHECK_EQUAL(ext->as<Solid>().numShells(), 1U);
}

BOOST_AUTO_TEST_CASE(testExtrudeMultiPolygon)
{
  std::vector<Point> points;
  points.emplace_back(0.0, 0.0, 0.0);
  points.emplace_back(1.0, 0.0, 0.0);
  points.emplace_back(1.0, 1.0, 0.0);
  points.emplace_back(0.0, 1.0, 0.0);
  points.emplace_back(0.0, 0.0, 0.0);

  std::vector<Point> points2;
  points2.emplace_back(2.0, 0.0, 0.0);
  points2.emplace_back(3.0, 0.0, 0.0);
  points2.emplace_back(3.0, 1.0, 0.0);
  points2.emplace_back(2.0, 1.0, 0.0);
  points2.emplace_back(2.0, 0.0, 0.0);

  LineString const exteriorRing(points);
  LineString const exteriorRing2(points2);
  Polygon const    g1(exteriorRing);
  Polygon const    g2(exteriorRing2);
  MultiPolygon     mp;
  mp.addGeometry(g1);
  mp.addGeometry(g2);

  std::unique_ptr<Geometry> ext(algorithm::extrude(mp, 0.0, 0.0, 1.0));
  BOOST_CHECK(ext->is<MultiSolid>());
  BOOST_CHECK_EQUAL(ext->as<MultiSolid>().numGeometries(), 2U);
}

BOOST_AUTO_TEST_CASE(testExtrudeSquareWithHole)
{
  std::vector<LineString> rings;
  {
    std::vector<Point> points;
    points.emplace_back(0.0, 0.0, 0.0);
    points.emplace_back(1.0, 0.0, 0.0);
    points.emplace_back(1.0, 1.0, 0.0);
    points.emplace_back(0.0, 1.0, 0.0);
    points.emplace_back(0.0, 0.0, 0.0);
    rings.emplace_back(points);
  }
  {
    std::vector<Point> points;
    points.emplace_back(0.2, 0.2, 0.0);
    points.emplace_back(0.8, 0.2, 0.0);
    points.emplace_back(0.8, 0.8, 0.0);
    points.emplace_back(0.2, 0.8, 0.0);
    points.emplace_back(0.2, 0.2, 0.0);

    std::reverse(points.begin(), points.end());

    rings.emplace_back(points);
  }

  Polygon const             g(rings);
  std::unique_ptr<Geometry> ext(algorithm::extrude(g, 0.0, 0.0, 1.0));
  BOOST_CHECK(ext->is<Solid>());
  BOOST_CHECK_EQUAL(ext->as<Solid>().numShells(), 1U);
  BOOST_CHECK_EQUAL(ext->as<Solid>().exteriorShell().numPatches(), 10U);
}

// SELECT ST_AsText(ST_Extrude(ST_Extrude(ST_Extrude('POINT (0 0)', 1, 0, 0), 0,
// 1, 0), 0, 0, 1));
BOOST_AUTO_TEST_CASE(testChainingExtrude)
{
  std::unique_ptr<Geometry> g(new Point(0.0, 0.0));
  g = algorithm::extrude(*g, 1.0, 0.0, 0.0);
  BOOST_CHECK_EQUAL(g->asText(0), "LINESTRING Z (0 0 0,1 0 0)");
  g = algorithm::extrude(*g, 0.0, 1.0, 0.0);
  BOOST_CHECK_EQUAL(g->asText(0),
                    "POLYHEDRALSURFACE Z (((0 0 0,1 0 0,1 1 0,0 1 0,0 0 0)))");
  g = algorithm::extrude(*g, 0.0, 0.0, 1.0);
  BOOST_CHECK_EQUAL(
      g->asText(0),
      "SOLID Z ((((0 1 0,1 1 0,1 0 0,0 1 0)),((0 1 1,1 0 1,1 1 1,0 1 1)),((0 1 "
      "0,1 0 0,0 0 0,0 1 0)),((0 1 1,0 0 1,1 0 1,0 1 1)),((1 0 0,1 1 0,1 1 1,1 "
      "0 1,1 0 0)),((1 1 0,0 1 0,0 1 1,1 1 1,1 1 0)),((0 1 0,0 0 0,0 0 1,0 1 "
      "1,0 1 0)),((0 0 0,1 0 0,1 0 1,0 0 1,0 0 0))))");
}

// Test for empty geometries - Issue #315
BOOST_AUTO_TEST_CASE(testExtrudeEmptyPoint)
{
  Point emptyPoint;
  BOOST_CHECK(emptyPoint.isEmpty());
  std::unique_ptr<Geometry> ext(algorithm::extrude(emptyPoint, 0.0, 0.0, 1.0));
  BOOST_CHECK(ext->isEmpty());
  BOOST_CHECK(ext->is<LineString>());
}

BOOST_AUTO_TEST_CASE(testExtrudeEmptyLineString)
{
  LineString emptyLineString;
  BOOST_CHECK(emptyLineString.isEmpty());
  std::unique_ptr<Geometry> ext(
      algorithm::extrude(emptyLineString, 0.0, 0.0, 1.0));
  BOOST_CHECK(ext->isEmpty());
  BOOST_CHECK(ext->is<PolyhedralSurface>());
}

BOOST_AUTO_TEST_CASE(testExtrudeEmptyPolygon)
{
  // This tests ForceZOrderPoints::visit(Polygon) with empty polygon
  Polygon emptyPolygon;
  BOOST_CHECK(emptyPolygon.isEmpty());
  std::unique_ptr<Geometry> ext(
      algorithm::extrude(emptyPolygon, 0.0, 0.0, 1.0));
  BOOST_CHECK(ext->isEmpty());
  BOOST_CHECK(ext->is<Solid>());

  // Also test with WKT parsing
  std::unique_ptr<Geometry> geom(io::readWkt("POLYGON EMPTY"));
  BOOST_REQUIRE(geom->is<Polygon>());
  std::unique_ptr<Geometry> ext2(algorithm::extrude(*geom, 0.0, 0.0, 1.0));
  BOOST_CHECK(ext2->isEmpty());
  BOOST_CHECK(ext2->is<Solid>());
}

BOOST_AUTO_TEST_SUITE_END()
