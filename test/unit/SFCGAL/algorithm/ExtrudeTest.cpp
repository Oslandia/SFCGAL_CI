// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "SFCGAL/Envelope.h"
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
#include "SFCGAL/algorithm/roofGeneration.h"
#include "SFCGAL/algorithm/volume.h"
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

// Test extrudeUntil with flat roof (planar surface)
BOOST_AUTO_TEST_CASE(testExtrudeUntilFlatRoof)
{
  // Create a square footprint at z=0
  std::vector<Point> points;
  points.emplace_back(0.0, 0.0, 0.0);
  points.emplace_back(10.0, 0.0, 0.0);
  points.emplace_back(10.0, 10.0, 0.0);
  points.emplace_back(0.0, 10.0, 0.0);
  points.emplace_back(0.0, 0.0, 0.0);
  LineString const exteriorRing(points);
  Polygon const    footprint(exteriorRing);

  // Create a flat roof at z=5
  PolyhedralSurface  roof;
  std::vector<Point> roofPoints;
  roofPoints.emplace_back(0.0, 0.0, 5.0);
  roofPoints.emplace_back(10.0, 0.0, 5.0);
  roofPoints.emplace_back(10.0, 10.0, 5.0);
  roofPoints.emplace_back(0.0, 10.0, 5.0);
  roofPoints.emplace_back(0.0, 0.0, 5.0);
  roof.addPolygon(Polygon(LineString(roofPoints)));

  // Extrude until roof
  std::unique_ptr<Solid> result = algorithm::extrudeUntil(footprint, roof);

  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK_EQUAL(result->numShells(), 1U);
  BOOST_CHECK(result->exteriorShell().numPolygons() >=
              6U); // bottom + 4 sides + top

  std::cout << footprint.asText(2) << "\n";
  std::cout << roof.asText(2) << "\n";
  std::cout << result->asText(2) << "\n";
}

// Test extrudeUntil with gable roof (non-planar surface)
BOOST_AUTO_TEST_CASE(testExtrudeUntilGableRoof)
{
  // Create a rectangular footprint at z=0
  std::vector<Point> points;
  points.emplace_back(0.0, 0.0, 0.0);
  points.emplace_back(10.0, 0.0, 0.0);
  points.emplace_back(10.0, 6.0, 0.0);
  points.emplace_back(0.0, 6.0, 0.0);
  points.emplace_back(0.0, 0.0, 0.0);
  LineString const exteriorRing(points);
  Polygon const    footprint(exteriorRing);

  // Create a gable roof (two triangular faces)
  PolyhedralSurface roof;

  // First slope: from bottom edge to ridge
  std::vector<Point> slope1;
  slope1.emplace_back(0.0, 0.0, 3.0);
  slope1.emplace_back(10.0, 0.0, 3.0);
  slope1.emplace_back(10.0, 3.0, 6.0);
  slope1.emplace_back(0.0, 3.0, 6.0);
  slope1.emplace_back(0.0, 0.0, 3.0);
  roof.addPolygon(Polygon(LineString(slope1)));

  // Second slope: from ridge to top edge
  std::vector<Point> slope2;
  slope2.emplace_back(0.0, 3.0, 6.0);
  slope2.emplace_back(10.0, 3.0, 6.0);
  slope2.emplace_back(10.0, 6.0, 3.0);
  slope2.emplace_back(0.0, 6.0, 3.0);
  slope2.emplace_back(0.0, 3.0, 6.0);
  roof.addPolygon(Polygon(LineString(slope2)));

  // Extrude until roof
  std::unique_ptr<Solid> result = algorithm::extrudeUntil(footprint, roof);

  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK_EQUAL(result->numShells(), 1U);
  // Should have triangulated top surface due to non-planarity
  BOOST_CHECK(result->exteriorShell().numPolygons() >= 6U);

  std::cout << "gable\n";
  std::cout << "footprint " << footprint.asText(2) << "\n";
  std::cout << "roof " << roof.asText(2) << "\n";
  std::cout << "solid " << result->asText(2) << "\n";
}

// Test extrudeUntil with footprint with hole
BOOST_AUTO_TEST_CASE(testExtrudeUntilWithHole)
{
  // Create exterior ring
  std::vector<Point> exterior;
  exterior.emplace_back(0.0, 0.0, 0.0);
  exterior.emplace_back(10.0, 0.0, 0.0);
  exterior.emplace_back(10.0, 10.0, 0.0);
  exterior.emplace_back(0.0, 10.0, 0.0);
  exterior.emplace_back(0.0, 0.0, 0.0);

  // Create interior ring (hole)
  std::vector<Point> interior;
  interior.emplace_back(3.0, 3.0, 0.0);
  interior.emplace_back(7.0, 3.0, 0.0);
  interior.emplace_back(7.0, 7.0, 0.0);
  interior.emplace_back(3.0, 7.0, 0.0);
  interior.emplace_back(3.0, 3.0, 0.0);
  std::reverse(interior.begin(), interior.end());

  std::vector<LineString> rings;
  rings.emplace_back(exterior);
  rings.emplace_back(interior);
  Polygon const footprint(rings);

  // Create a flat roof at z=5
  PolyhedralSurface  roof;
  std::vector<Point> roofPoints;
  roofPoints.emplace_back(0.0, 0.0, 5.0);
  roofPoints.emplace_back(10.0, 0.0, 5.0);
  roofPoints.emplace_back(10.0, 10.0, 5.0);
  roofPoints.emplace_back(0.0, 10.0, 5.0);
  roofPoints.emplace_back(0.0, 0.0, 5.0);
  roof.addPolygon(Polygon(LineString(roofPoints)));

  // Extrude until roof
  std::unique_ptr<Solid> result = algorithm::extrudeUntil(footprint, roof);

  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK_EQUAL(result->numShells(), 1U);
  // Should have bottom + exterior walls + interior walls + top
  BOOST_CHECK(result->exteriorShell().numPolygons() >= 10U);
  std::cout << footprint.asText(2) << "\n";
  std::cout << roof.asText(2) << "\n";
  std::cout << result->asText(2) << "\n";
}

// Test extrudeUntil with empty inputs
BOOST_AUTO_TEST_CASE(testExtrudeUntilEmptyInputs)
{
  Polygon const          emptyFootprint;
  PolyhedralSurface      emptyRoof;
  std::unique_ptr<Solid> result =
      algorithm::extrudeUntil(emptyFootprint, emptyRoof);

  BOOST_CHECK(result->isEmpty());
}

BOOST_AUTO_TEST_CASE(testExtrudeUntilRoofOverhang)
{
  // Footprint: 10x10 square at z=0
  std::vector<Point> points;
  points.emplace_back(0.0, 0.0, 0.0);
  points.emplace_back(10.0, 0.0, 0.0);
  points.emplace_back(10.0, 10.0, 0.0);
  points.emplace_back(0.0, 10.0, 0.0);
  points.emplace_back(0.0, 0.0, 0.0);
  LineString const exteriorRing(points);
  Polygon const    footprint(exteriorRing);

  // Roof: 20x20 flat surface at z=5, centered on footprint
  // Extends from -5 to 15 in both X and Y
  PolyhedralSurface  roof;
  std::vector<Point> roofPoints;
  roofPoints.emplace_back(-5.0, -5.0, 5.0);
  roofPoints.emplace_back(15.0, -5.0, 5.0);
  roofPoints.emplace_back(15.0, 15.0, 5.0);
  roofPoints.emplace_back(-5.0, 15.0, 5.0);
  roofPoints.emplace_back(-5.0, -5.0, 5.0);
  roof.addPolygon(Polygon(LineString(roofPoints)));

  // Extrude until roof
  // Expected: Solid with 10x10 base, 10x10 top at z=5, and vertical walls
  std::unique_ptr<Solid> result = algorithm::extrudeUntil(footprint, roof);

  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK_EQUAL(result->numShells(), 1U);

  // Check volume (10 * 10 * 5 = 500)
  BOOST_CHECK_CLOSE(CGAL::to_double(algorithm::volume(*result)), 500.0, 0.001);

  // Check bounds
  const Envelope &env = result->envelope();
  BOOST_CHECK_CLOSE(env.xMin(), 0.0, 0.001);
  BOOST_CHECK_CLOSE(env.xMax(), 10.0, 0.001);
  BOOST_CHECK_CLOSE(env.yMin(), 0.0, 0.001);
  BOOST_CHECK_CLOSE(env.yMax(), 10.0, 0.001);
  BOOST_CHECK_CLOSE(env.zMin(), 0.0, 0.001);
  BOOST_CHECK_CLOSE(env.zMax(), 5.0, 0.001);
}

BOOST_AUTO_TEST_CASE(testExtrudeUntilGableRoofOverhang)
{
  std::string footprintWKT{"POLYGON((0 0,10 0,10 6,0 6,0 0))"};
  // std::string overhangWKT{"Polygon ((-1 -1, -1 7, 11 7, 11 -1, -1 -1))"};
  // Generated via: ./build/sfcgalop/sfcgalop -a "Polygon ((-1 -1, -1 7, 11 7,
  // 11 -1, -1 -1))" generate_gable_roof "slope_angle=30" |
  // ./build/sfcgalop/sfcgalop translate "dz=5"
  std::string gableRoofWKT{
      "POLYHEDRALSURFACE Z (((7.000000 3.000000 5.577350,-1.000000 -1.000000 "
      "5.000000,11.000000 -1.000000 5.000000,7.000000 3.000000 "
      "5.577350)),((-1.000000 7.000000 5.000000,-1.000000 3.000000 "
      "5.577350,3.000000 3.000000 5.577350,-1.000000 7.000000 "
      "5.000000)),((3.000000 3.000000 5.577350,-1.000000 3.000000 "
      "5.577350,-1.000000 -1.000000 5.000000,3.000000 3.000000 "
      "5.577350)),((7.000000 3.000000 5.577350,3.000000 3.000000 "
      "5.577350,-1.000000 -1.000000 5.000000,7.000000 3.000000 "
      "5.577350)),((7.000000 3.000000 5.577350,-1.000000 7.000000 "
      "5.000000,3.000000 3.000000 5.577350,7.000000 3.000000 "
      "5.577350)),((11.000000 3.000000 5.577350,7.000000 3.000000 "
      "5.577350,11.000000 -1.000000 5.000000,11.000000 3.000000 "
      "5.577350)),((11.000000 7.000000 5.000000,-1.000000 7.000000 "
      "5.000000,7.000000 3.000000 5.577350,11.000000 7.000000 "
      "5.000000)),((11.000000 7.000000 5.000000,7.000000 3.000000 "
      "5.577350,11.000000 3.000000 5.577350,11.000000 7.000000 5.000000)))"};
  std::string expectedWKT{
      "SOLID Z ((((7.000000 3.000000 5.577350,1.000000 0.000000 "
      "5.144338,10.000000 0.000000 5.144338,7.000000 3.000000 "
      "5.577350)),((0.000000 6.000000 5.144338,0.000000 3.000000 "
      "5.577350,3.000000 3.000000 5.577350,0.000000 6.000000 "
      "5.144338)),((0.000000 3.000000 5.577350,0.000000 0.000000 "
      "5.144338,3.000000 3.000000 5.577350,0.000000 3.000000 "
      "5.577350)),((0.000000 0.000000 5.144338,1.000000 0.000000 "
      "5.144338,7.000000 3.000000 5.577350,3.000000 3.000000 5.577350,0.000000 "
      "0.000000 5.144338)),((1.000000 6.000000 5.144338,0.000000 6.000000 "
      "5.144338,3.000000 3.000000 5.577350,7.000000 3.000000 5.577350,1.000000 "
      "6.000000 5.144338)),((7.000000 3.000000 5.577350,10.000000 0.000000 "
      "5.144338,10.000000 3.000000 5.577350,7.000000 3.000000 "
      "5.577350)),((1.000000 6.000000 5.144338,7.000000 3.000000 "
      "5.577350,10.000000 6.000000 5.144338,1.000000 6.000000 "
      "5.144338)),((7.000000 3.000000 5.577350,10.000000 3.000000 "
      "5.577350,10.000000 6.000000 5.144338,7.000000 3.000000 "
      "5.577350)),((0.000000 3.000000 0.000000,0.000000 6.000000 "
      "0.000000,1.000000 6.000000 0.000000,10.000000 6.000000 0.000000,10.000000 "
      "3.000000 0.000000,10.000000 0.000000 0.000000,1.000000 0.000000 "
      "0.000000,0.000000 0.000000 0.000000,0.000000 3.000000 "
      "0.000000)),((0.000000 3.000000 5.577350,0.000000 3.000000 "
      "0.000000,0.000000 0.000000 0.000000,0.000000 0.000000 5.144338,0.000000 "
      "3.000000 5.577350)),((0.000000 0.000000 5.144338,0.000000 0.000000 "
      "0.000000,1.000000 0.000000 0.000000,1.000000 0.000000 5.144338,0.000000 "
      "0.000000 5.144338)),((1.000000 0.000000 5.144338,1.000000 0.000000 "
      "0.000000,10.000000 0.000000 0.000000,10.000000 0.000000 "
      "5.144338,1.000000 0.000000 5.144338)),((10.000000 0.000000 "
      "5.144338,10.000000 0.000000 0.000000,10.000000 3.000000 "
      "0.000000,10.000000 3.000000 5.577350,10.000000 0.000000 "
      "5.144338)),((10.000000 3.000000 5.577350,10.000000 3.000000 "
      "0.000000,10.000000 6.000000 0.000000,10.000000 6.000000 "
      "5.144338,10.000000 3.000000 5.577350)),((10.000000 6.000000 "
      "5.144338,10.000000 6.000000 0.000000,1.000000 6.000000 0.000000,1.000000 "
      "6.000000 5.144338,10.000000 6.000000 5.144338)),((1.000000 6.000000 "
      "5.144338,1.000000 6.000000 0.000000,0.000000 6.000000 0.000000,0.000000 "
      "6.000000 5.144338,1.000000 6.000000 5.144338)),((0.000000 6.000000 "
      "5.144338,0.000000 6.000000 0.000000,0.000000 3.000000 0.000000,0.000000 "
      "3.000000 5.577350,0.000000 6.000000 5.144338))))"};

  std::unique_ptr<Geometry> footprint(io::readWkt(footprintWKT));
  std::unique_ptr<Geometry> gableRoof(io::readWkt(gableRoofWKT));
  // Extrude until roof
  std::unique_ptr<Solid> result = algorithm::extrudeUntil(
      footprint->as<Polygon>(), gableRoof->as<PolyhedralSurface>());

  // Check the result matches expected WKT
  BOOST_CHECK_EQUAL(result->asText(6), expectedWKT);
  // Also verify the solid is valid
  BOOST_CHECK(algorithm::isValid(*result));
}

BOOST_AUTO_TEST_SUITE_END()
