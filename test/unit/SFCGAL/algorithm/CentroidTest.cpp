/**
 *   SFCGAL
 *
 *   Copyright (C) 2012-2013 Oslandia <infos@oslandia.com>
 *   Copyright (C) 2012-2013 IGN (http://www.ign.fr)
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Library General Public
 *   License as published by the Free Software Foundation; either
 *   version 2 of the License, or (at your option) any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Library General Public License for more details.

 *   You should have received a copy of the GNU Library General Public
 *   License along with this library; if not, see
 <http://www.gnu.org/licenses/>.
 */
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
#include "SFCGAL/io/wkt.h"

#include "SFCGAL/detail/tools/Registry.h"
#include <SFCGAL/algorithm/extrude.h>
#include <SFCGAL/algorithm/rotate.h>
#include <SFCGAL/algorithm/translate.h>

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_CentroidTest)

BOOST_AUTO_TEST_CASE(testCentroid_Empty)
{
  tools::Registry const &registry = tools::Registry::instance();
  std::vector<std::string> const typeNames = tools::Registry::instance().getGeometryTypes();

  for (const auto &typeName : typeNames) {
    // if (typeName == "Point" || typeName == "Triangle")
    //   continue;
    BOOST_TEST_MESSAGE(typeName);

    std::unique_ptr<Geometry> const g(registry.newGeometryByTypeName(typeName));
    BOOST_REQUIRE(g.get() != nullptr);
    try {
      g->centroid();
      BOOST_CHECK_MESSAGE(false, typeName);
    } catch (Exception &) {
      // ok
    }
  }
}

// must return 0.0
BOOST_AUTO_TEST_CASE(testCentroid2D_Point2D3D)
{
  BOOST_CHECK_EQUAL(Point(3.0, 4.0).centroid().asText(2), Point(3.0, 4.0).asText(2));
  BOOST_CHECK_EQUAL(Point(3.0, 4.0, 5.0).centroid().asText(2), Point(3.0, 4.0, 5.0).asText(2));
}

// must return 0.0
BOOST_AUTO_TEST_CASE(testCentroid2D_LineString2D3D)
{
  BOOST_CHECK_EQUAL(LineString(Point(0.0, 0.0), Point(1.0, 1.0)).centroid().asText(2),
                    Point(0.5, 0.5).asText(2));
  std::vector<Point> points = {Point(0.0, -1.0, 1.0), Point(2.0, -1.0, 1.0), Point(12.0, -1.0, 1.0)};
  BOOST_CHECK_EQUAL(LineString(points).centroid().asText(2), Point(6.0, -1.0, 1.0).asText(2));
}

BOOST_AUTO_TEST_CASE(testCentroid2D_Triangle2D)
{
  Triangle const triangle1(Point(0.0, 0.0), Point(4.0, 0.0), Point(4.0, 4.0));
  // the same, inverted
  Triangle const triangle2(Point(0.0, 0.0), Point(0.0, 4.0), Point(4.0, 4.0));
  BOOST_CHECK_EQUAL(triangle1.centroid().asText(2), Point(2.67, 1.33).asText(2));
  BOOST_CHECK_EQUAL(triangle2.centroid().asText(2), Point(1.33, 2.67).asText(2));
}

BOOST_AUTO_TEST_CASE(testCentroid2D_Triangle3D)
{
  Triangle const triangle(Point(0.0, 0.0, 0.0), Point(0.0, 0.0, 1.0), Point(0.0, 1.0, 0.0));
  BOOST_CHECK_EQUAL(triangle.centroid().asText(2), Point(0.0, 0.33, 0.33).asText(2));

  Triangle const triangle3D(Point(0.0, 0.0, 0.0), Point(0.0, 0.0, 4.0), Point(0.0, 4.0, 0.0));
  BOOST_CHECK_EQUAL(triangle3D.centroid().asText(2), Point(0.0, 1.33, 1.33).asText(2));
}

BOOST_AUTO_TEST_CASE(testCentroid2D_Triangle4D)
{
  Triangle const triangle(Point(0.0, 0.0, 0.0, 0.0),
                          Point(0.0, 0.0, 1.0, 1.0),
                          Point(0.0, 1.0, 0.0, 2.0));
  BOOST_CHECK_EQUAL(triangle.centroid().asText(2), Point(0.0, 0.33, 0.33, 1.0).asText(2));
}

BOOST_AUTO_TEST_CASE(testCentroid2D_Square2D1x1)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("POLYGON ((0.0 0.0, 1.0 0.0, 1.0 1.0, 0.0 1.0, 0.0 0.0))"));
  BOOST_CHECK_EQUAL(g->asText(1), "POLYGON ((0.0 0.0,1.0 0.0,1.0 1.0,0.0 1.0,0.0 0.0))");
  BOOST_CHECK_EQUAL(g->centroid().asText(2), Point(0.5, 0.5).asText(2));
}

BOOST_AUTO_TEST_CASE(testCentroid3D_Square2D1x1)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("POLYGON ((0.0 0.0, 1.0 0.0, 1.0 1.0, 0.0 1.0, 0.0 0.0))"));
  BOOST_CHECK_EQUAL(g->asText(1), "POLYGON ((0.0 0.0,1.0 0.0,1.0 1.0,0.0 1.0,0.0 0.0))");
  BOOST_CHECK_EQUAL(g->centroid3D().asText(2), Point(0.5, 0.5).asText(2));
}

BOOST_AUTO_TEST_CASE(testCentroid2D_Square3D1x1)
{
  { // vertical geometry (no 2D surface)
    std::unique_ptr<Geometry> g(io::readWkt("POLYGON ((0.0 0.0 0.0,0.0 0.0 1.0,0.0 1.0 1.0,0.0 1.0 "
                                            "0.0,0.0 0.0 0.0))"));
    BOOST_CHECK_EQUAL(g->asText(1),
                      "POLYGON Z ((0.0 0.0 0.0,0.0 0.0 1.0,0.0 1.0 "
                      "1.0,0.0 1.0 0.0,0.0 0.0 0.0))");
    try {
      g->centroid();
      BOOST_CHECK_MESSAGE(false, "vertical geometry (no 2D surface) should fail");
    } catch (InappropriateGeometryException &) {
      // ok
    }
  }

  {
    std::unique_ptr<Geometry> g(io::readWkt("POLYGON ((0.0 0.0 0.0,"
                                            "1.0 0.0 1.0,"
                                            "1.0 1.0 1.0,"
                                            "0.0 1.0 0.0,"
                                            "0.0 0.0 0.0))"));
    BOOST_CHECK_EQUAL(g->asText(1),
                      "POLYGON Z ((0.0 0.0 0.0,"
                      "1.0 0.0 1.0,"
                      "1.0 1.0 1.0,"
                      "0.0 1.0 0.0,"
                      "0.0 0.0 0.0))");
    BOOST_CHECK_EQUAL(g->centroid().asText(2), Point(0.5, 0.5, 0.5).asText(2));
  }
}

BOOST_AUTO_TEST_CASE(testCentroid2D_Square3D4X4)
{
  std::unique_ptr<Geometry> g(io::readWkt("POLYGON ((0.0 0.0 0.0,"
                                          "4.0 0.0 4.0,"
                                          "4.0 4.0 4.0,"
                                          "0.0 4.0 0.0,"
                                          "0.0 0.0 0.0))"));
  BOOST_CHECK_EQUAL(g->asText(1),
                    "POLYGON Z ((0.0 0.0 0.0,"
                    "4.0 0.0 4.0,"
                    "4.0 4.0 4.0,"
                    "0.0 4.0 0.0,"
                    "0.0 0.0 0.0))");
  BOOST_CHECK_EQUAL(g->centroid().asText(2), Point(2.0, 2.0, 2.0).asText(2));
}

BOOST_AUTO_TEST_CASE(testCentroid2D_Square3D4X4WithHole)
{
  std::unique_ptr<Geometry> g(io::readWkt("POLYGON ((0.0 0.0 0.0,"
                                          "4.0 0.0 4.0,"
                                          "4.0 4.0 4.0,"
                                          "0.0 4.0 0.0,"
                                          "0.0 0.0 0.0),"
                                          "(1.0 1.0 1.0,"
                                          "3.0 1.0 3.0,"
                                          "3.0 3.0 3.0,"
                                          "1.0 3.0 1.0,"
                                          "1.0 1.0 1.0))"));
  BOOST_CHECK_EQUAL(g->asText(1),
                    "POLYGON Z ((0.0 0.0 0.0,"
                    "4.0 0.0 4.0,"
                    "4.0 4.0 4.0,"
                    "0.0 4.0 0.0,"
                    "0.0 0.0 0.0),"
                    "(1.0 1.0 1.0,"
                    "3.0 1.0 3.0,"
                    "3.0 3.0 3.0,"
                    "1.0 3.0 1.0,"
                    "1.0 1.0 1.0))");
  BOOST_CHECK_EQUAL(g->centroid().asText(2), Point(2.0, 2.0, 2.0).asText(2));
}

BOOST_AUTO_TEST_CASE(testCentroid3D_Square3D4X4WithHole)
{
  std::unique_ptr<Geometry> g(io::readWkt("POLYGON ((0.0 0.0 0.0,"
                                          "4.0 0.0 4.0,"
                                          "4.0 4.0 4.0,"
                                          "0.0 4.0 0.0,"
                                          "0.0 0.0 0.0),"
                                          "(1.0 1.0 1.0,"
                                          "3.0 1.0 3.0,"
                                          "3.0 3.0 3.0,"
                                          "1.0 3.0 1.0,"
                                          "1.0 1.0 1.0))"));
  BOOST_CHECK_EQUAL(g->asText(1),
                    "POLYGON Z ((0.0 0.0 0.0,"
                    "4.0 0.0 4.0,"
                    "4.0 4.0 4.0,"
                    "0.0 4.0 0.0,"
                    "0.0 0.0 0.0),"
                    "(1.0 1.0 1.0,"
                    "3.0 1.0 3.0,"
                    "3.0 3.0 3.0,"
                    "1.0 3.0 1.0,"
                    "1.0 1.0 1.0))");
  BOOST_CHECK_EQUAL(g->centroid3D().asText(2), Point(2.0, 2.0, 2.0).asText(2));
}

BOOST_AUTO_TEST_CASE(testCentroid2D_PerpendicularSquares)
{
  std::unique_ptr<Geometry> g(io::readWkt("MULTIPOLYGON (((0.0 0.0 0.0,"
                                          "0.0 0.0 4.0,"
                                          "4.0 0.0 4.0,"
                                          "4.0 0.0 0.0,"
                                          "0.0 0.0 0.0)),"
                                          "((0.0 0.0 4.0,"
                                          "4.0 0.0 4.0,"
                                          "4.0 4.0 4.0,"
                                          "0.0 4.0 4.0,"
                                          "0.0 0.0 4.0)))"));
  SFCGAL::algorithm::translate(*(g.get()),
                               Kernel::Vector_3(-2.0, -2.0, -2.0)); // translate around 0,0,0
  double angle = M_PI / 90.0;
  SFCGAL::algorithm::rotate(*(g.get()),
                            angle,
                            Kernel::Vector_3(1.0, 0.0, 0.0)); // rotation around (1,0,0)
  BOOST_CHECK_EQUAL(
      g->asText(1),
      "MULTIPOLYGON Z (((-2.0 -1.9 -2.1,-2.0 -2.1 1.9,2.0 -2.1 1.9,2.0 -1.9 -2.1,-2.0 -1.9 -2.1)),"
      "((-2.0 -2.1 1.9,2.0 -2.1 1.9,2.0 1.9 2.1,-2.0 1.9 2.1,-2.0 -2.1 1.9)))");
  BOOST_CHECK_EQUAL(g->centroid().asText(2), Point(0.00, -0.13, 1.93).asText(2));
}

BOOST_AUTO_TEST_CASE(testCentroid3D_PerpendicularSquares)
{
  std::unique_ptr<Geometry> g(io::readWkt("MULTIPOLYGON (((0.0 0.0 0.0,"
                                          "0.0 0.0 4.0,"
                                          "4.0 0.0 4.0,"
                                          "4.0 0.0 0.0,"
                                          "0.0 0.0 0.0)),"
                                          "((0.0 0.0 4.0,"
                                          "4.0 0.0 4.0,"
                                          "4.0 4.0 4.0,"
                                          "0.0 4.0 4.0,"
                                          "0.0 0.0 4.0)))"));
  SFCGAL::algorithm::translate(*(g.get()),
                               Kernel::Vector_3(-2.0, -2.0, -2.0)); // translate around 0,0,0
  Point expected(0.0, -1.0, 1.0);
  for (int i = 0; i < 90; ++i) {
    double angle = static_cast<double>(i) * (M_PI / 90.0);
    SFCGAL::algorithm::rotate(*(g.get()),
                              angle,
                              Kernel::Vector_3(1.0, 0.0, 0.0)); // rotation around (1,0,0)
    SFCGAL::algorithm::rotate(expected,
                              angle,
                              Kernel::Vector_3(1.0, 0.0, 0.0)); // rotation around (1,0,0)
    BOOST_CHECK_EQUAL(g->centroid3D().asText(2), expected.asText(2));
  }
}

BOOST_AUTO_TEST_CASE(testCentroid2D_polyhedral)
{
  std::unique_ptr<Geometry> g(io::readWkt("POLYHEDRALSURFACE (((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)),"
                                          "((1 0 0,1 1 0,2 1 0,2 0 0,1 0 0)))"));
  double z = 2.0;
  SFCGAL::algorithm::translate(*(g.get()), Kernel::Vector_3(0.0, 0.0, z)); // translate around 0,0,0
  Point expected(1.0, 0.5, z);
  for (int i = 0; i < 90; ++i) {
    double angle = static_cast<double>(i) * (M_PI / 90.0);
    SFCGAL::algorithm::rotate(*(g.get()),
                              angle,
                              Kernel::Vector_3(1.0, 0.0, 0.0)); // rotation around (1,0,0)
    SFCGAL::algorithm::rotate(expected,
                              angle,
                              Kernel::Vector_3(1.0, 0.0, 0.0)); // rotation around (1,0,0)
    BOOST_CHECK_EQUAL(g->centroid().asText(2), expected.asText(2));
  }
}

BOOST_AUTO_TEST_CASE(testCentroid3D_polyhedral)
{
  std::unique_ptr<Geometry> g(io::readWkt("POLYHEDRALSURFACE (((0 0 0,0 1 0,5 1 0,5 0 0,0 0 0)),"
                                          "((5 0 0,5 1 0,8 1 4,8 0 4,5 0 0)))"));
  double z = 2.0;
  SFCGAL::algorithm::translate(*(g.get()), Kernel::Vector_3(0.0, 0.0, z)); // translate around 0,0,0
  Point expected(4.5, 0.5, 1.0 + z);
  for (int i = 0; i < 90; ++i) {
    double angle = static_cast<double>(i) * (M_PI / 90.0);
    SFCGAL::algorithm::rotate(*(g.get()),
                              angle,
                              Kernel::Vector_3(1.0, 0.0, 0.0)); // rotation around (1,0,0)
    SFCGAL::algorithm::rotate(expected,
                              angle,
                              Kernel::Vector_3(1.0, 0.0, 0.0)); // rotation around (1,0,0)
    BOOST_CHECK_EQUAL(g->centroid3D().asText(2), expected.asText(2));
  }
}

BOOST_AUTO_TEST_CASE(testCentroid3D_solid)
{
  std::vector<Point> points;
  points.emplace_back(0.0, 0.0, 0.0);
  points.emplace_back(1.0, 0.0, 0.0);
  points.emplace_back(1.0, 1.0, 0.0);
  points.emplace_back(0.0, 1.0, 0.0);
  points.emplace_back(0.0, 0.0, 0.0);

  LineString const exteriorRing(points);
  Polygon const exteriorPoly(exteriorRing);
  std::unique_ptr<Geometry> g(algorithm::extrude(exteriorPoly, 0.0, 0.0, 1.0));
  BOOST_CHECK(g->is<Solid>());

  double z = 2.0;
  SFCGAL::algorithm::translate(*(g.get()), Kernel::Vector_3(0.0, 0.0, z)); // translate around 0,0,0
  Point expected(0.5, 0.5, 0.5 + z);
  for (int i = 0; i < 90; ++i) {
    //int i = 0;
    double angle = static_cast<double>(i) * (M_PI / 90.0);
    SFCGAL::algorithm::rotate(*(g.get()),
                              angle,
                              Kernel::Vector_3(1.0, 0.0, 0.0)); // rotation around (1,0,0)
    SFCGAL::algorithm::rotate(expected,
                              angle,
                              Kernel::Vector_3(1.0, 0.0, 0.0)); // rotation around (1,0,0)
    BOOST_CHECK_EQUAL(g->centroid3D().asText(2), expected.asText(2));
  }
}

BOOST_AUTO_TEST_SUITE_END()
