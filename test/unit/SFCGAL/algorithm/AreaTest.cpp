// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Kernel.h"
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
#include "SFCGAL/algorithm/area.h"
#include "SFCGAL/io/wkt.h"

#include "SFCGAL/detail/tools/Registry.h"

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_AreaTest)

BOOST_AUTO_TEST_CASE(testEmpty2D3D)
{
  tools::Registry const         &registry = tools::Registry::instance();
  std::vector<std::string> const typeNames =
      tools::Registry::instance().getGeometryTypes();

  for (const auto &typeName : typeNames) {
    BOOST_TEST_MESSAGE(typeName);

    std::unique_ptr<Geometry> const g(registry.newGeometryByTypeName(typeName));
    BOOST_REQUIRE(g.get() != nullptr);
    BOOST_CHECK_EQUAL(algorithm::area(*g), 0.0);
    BOOST_CHECK_EQUAL(algorithm::area3D(*g), 0.0);
  }
}

BOOST_AUTO_TEST_CASE(testSignedArea2D_lineString)
{
  LineString lineString;
  lineString.addPoint(Point(0.0, 0.0));
  lineString.addPoint(Point(1.0, 0.0));
  lineString.addPoint(Point(1.0, 1.0));
  lineString.addPoint(Point(0.0, 1.0));
  lineString.addPoint(lineString.startPoint());

  BOOST_CHECK_EQUAL(algorithm::signedArea(lineString), 1.0);
  lineString.reverse();
  BOOST_CHECK_EQUAL(algorithm::signedArea(lineString), -1.0);
}

BOOST_AUTO_TEST_CASE(testSignedArea2D_triangle)
{
  Triangle triangle(Point(0.0, 0.0), Point(1.0, 0.0), Point(1.0, 1.0));

  BOOST_CHECK_EQUAL(algorithm::signedArea(triangle), 0.5);
  triangle.reverse();
  BOOST_CHECK_EQUAL(algorithm::signedArea(triangle), -0.5);
}

// must return 0.0
BOOST_AUTO_TEST_CASE(testPoint2D3D)
{
  BOOST_CHECK_EQUAL(algorithm::area(Point(3.0, 4.0)), 0.0);
  BOOST_CHECK_EQUAL(algorithm::area3D(Point(3.0, 4.0, 5.0)), 0.0);
}
// must return 0.0
BOOST_AUTO_TEST_CASE(testLineString2D3D)
{
  BOOST_CHECK_EQUAL(
      algorithm::area(LineString(Point(0.0, 0.0), Point(1.0, 1.0))), 0.0);
  BOOST_CHECK_EQUAL(
      algorithm::area3D(LineString(Point(0.0, 0.0, 0.0), Point(1.0, 1.0, 1.0))),
      0.0);
}

// must return 0.0
BOOST_AUTO_TEST_CASE(testArea2D_PolygonWithHoleWithBadOrientation)
{
  Polygon polygon;

  // exterior ring
  {
    LineString ring;
    ring.addPoint(Point(0.0, 0.0));
    ring.addPoint(Point(5.0, 0.0));
    ring.addPoint(Point(5.0, 5.0));
    ring.addPoint(Point(0.0, 5.0));
    ring.addPoint(ring.startPoint());

    polygon.setExteriorRing(ring);
  }

  // hole 1
  {
    LineString ring;
    ring.addPoint(Point(1.0, 1.0));
    ring.addPoint(Point(2.0, 1.0));
    ring.addPoint(Point(2.0, 2.0));
    ring.addPoint(Point(1.0, 2.0));
    ring.addPoint(ring.startPoint());

    polygon.addRing(ring);
  }

  // hole 2
  {
    LineString ring;
    ring.addPoint(Point(3.0, 3.0));
    ring.addPoint(Point(4.0, 3.0));
    ring.addPoint(Point(4.0, 4.0));
    ring.addPoint(Point(3.0, 4.0));
    ring.addPoint(ring.startPoint());

    polygon.addRing(ring);
  }

  // 5x5 - 1 - 1 = 23
  BOOST_CHECK_EQUAL(algorithm::area3D(polygon), 23.0);
}

BOOST_AUTO_TEST_CASE(testArea3D_Polygon)
{
  Polygon          polygon;
  const LineString exteriorRing({Point(0, 0, 1), Point(0, 10, 2),
                                 Point(10, 10, 3), Point(10, 0, 2),
                                 Point(0, 0, 1)});
  polygon.setExteriorRing(exteriorRing);
  BOOST_CHECK_CLOSE(algorithm::area3D(polygon), 100.995, 0.001);

  const LineString interiorRing({Point(1, 1, 1.2), Point(9, 1, 2),
                                 Point(9, 9, 2.8), Point(1, 9, 2),
                                 Point(1, 1, 1.2)});
  polygon.addInteriorRing(interiorRing);
  BOOST_CHECK_CLOSE(algorithm::area3D(polygon), 36.3582, 0.0001);
}

BOOST_AUTO_TEST_CASE(testArea3D_Triangle1)
{
  Triangle const triangle(Point(0.0, 0.0, 0.0), Point(0.0, 0.0, 1.0),
                          Point(0.0, 1.0, 0.0));
  BOOST_CHECK_EQUAL(algorithm::area3D(triangle), 0.5);
}

BOOST_AUTO_TEST_CASE(testArea3D_Triangle2)
{
  Triangle const triangle(Point(0.0, 0.0, 0.0), Point(0.0, 0.0, 4.0),
                          Point(0.0, 4.0, 0.0));
  BOOST_CHECK_EQUAL(algorithm::area3D(triangle), 8.0);
}

BOOST_AUTO_TEST_CASE(testArea2D_Triangle)
{
  Triangle const triangle1(Point(0.0, 0.0), Point(4.0, 0.0), Point(4.0, 4.0));
  // the same, inverted
  Triangle const triangle2(Point(0.0, 0.0), Point(0.0, 4.0), Point(4.0, 4.0));
  BOOST_CHECK_EQUAL(algorithm::area(triangle1), 8.0);
  BOOST_CHECK_EQUAL(algorithm::area(triangle2), 8.0);
}

BOOST_AUTO_TEST_CASE(testArea3D_GeometryCollection)
{
  GeometryCollection collection;
  LineString         line({Point(0, 0, 0), Point(0, 10, 0), Point(10, 10, 0),
                           Point(10, 0, 0), Point(0, 0, 0)});
  collection.addGeometry(line);
  // area3D of a line is 0
  BOOST_CHECK_EQUAL(algorithm::area3D(line), 0.0);
  BOOST_CHECK_EQUAL(algorithm::area3D(collection), 0.0);

  LineString line2(
      {Point(3, 3, 3), Point(13, 3, 8), Point(13, 13, 13), Point(3, 13, 8)});
  line2.closes();
  Polygon polygon;
  polygon.setExteriorRing(line2);
  collection.addGeometry(polygon.clone());
  BOOST_CHECK_CLOSE(algorithm::area3D(polygon), 122.474487, 0.0001);
  BOOST_CHECK_CLOSE(algorithm::area3D(collection), 122.474487, 0.000001);

  std::unique_ptr<Geometry> polyhedral(io::readWkt(
      "PolyhedralSurface Z (((0 0 1, 0 10 1, 10 10 1, 10 0 1, 0 0 1)))"));
  BOOST_CHECK_EQUAL(algorithm::area3D(polyhedral), 100.0);
  collection.addGeometry(std::move(polyhedral));
  BOOST_CHECK_CLOSE(algorithm::area3D(collection), 222.474487, 0.000001);
}

BOOST_AUTO_TEST_CASE(testArea3D_Tin)
{
  std::string const wkt("TIN Z ("
                        "((0 0 0, 10 0 0, 5 5 2, 0 0 0)),"
                        "((10 0 0, 10 10 1, 5 5 2, 10 0 0)),"
                        "((10 10 1, 0 10 0, 5 5 2, 10 10 1)),"
                        "((0 10 0, 0 0 0, 5 5 2, 0 10 0)))");

  std::unique_ptr<Geometry> const tin(io::readWkt(wkt));
  BOOST_CHECK_CLOSE(algorithm::area3D(*tin), 106.29209, 1e-5);
}

BOOST_AUTO_TEST_CASE(testArea3D_PolyhedralSurface)
{
  std::string const wkt("POLYHEDRALSURFACE Z ("
                        "((0 0 0, 10 0 0, 10 10 0, 0 10 0, 0 0 0)),"
                        "((0 0 10, 0 10 10, 10 10 10, 10 0 10, 0 0 10)),"
                        "((0 0 0, 0 10 0, 0 10 10, 0 0 10, 0 0 0)),"
                        "((10 0 0, 10 0 10, 10 10 10, 10 10 0, 10 0 0)),"
                        "((0 0 0, 0 0 10, 10 0 10, 10 0 0, 0 0 0)),"
                        "((0 10 0, 10 10 0, 10 10 10, 0 10 10, 0 10 0)))");

  std::unique_ptr<Geometry> const polySurface(io::readWkt(wkt));
  BOOST_CHECK_EQUAL(algorithm::area3D(*polySurface), 600.0);
}

BOOST_AUTO_TEST_CASE(testArea3D_Square1x1)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("POLYGON ((0.0 0.0 0.0,0.0 0.0 1.0,0.0 1.0 1.0,0.0 1.0 "
                  "0.0,0.0 0.0 0.0))"));
  BOOST_CHECK_EQUAL(g->asText(1), "POLYGON Z ((0.0 0.0 0.0,0.0 0.0 1.0,0.0 1.0 "
                                  "1.0,0.0 1.0 0.0,0.0 0.0 0.0))");
  BOOST_CHECK_CLOSE(algorithm::area3D(*g), 1.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(testArea3D_Square4X4)
{
  std::string const wkt("POLYGON ((0.0 0.0 0.0,0.0 0.0 4.0,0.0 4.0 4.0,0.0 4.0 "
                        "0.0,0.0 0.0 0.0))");
  std::unique_ptr<Geometry> const g(io::readWkt(wkt));
  BOOST_CHECK_CLOSE(algorithm::area3D(*g), 16.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(testArea3D_Square4X4WithHole)
{
  std::string const wkt(
      "POLYGON ((0.0 0.0 0.0,0.0 0.0 4.0,0.0 4.0 4.0,0.0 4.0 0.0,0.0 0.0 "
      "0.0),(0.0 2.0 2.0,0.0 3.0 2.0,0.0 3.0 3.0,0.0 2.0 3.0,0.0 2.0 2.0))");
  std::unique_ptr<Geometry> const g(io::readWkt(wkt));
  BOOST_CHECK_CLOSE(algorithm::area3D(*g), 15.0, 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()
