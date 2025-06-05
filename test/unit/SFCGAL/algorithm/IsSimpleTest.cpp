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
#include "SFCGAL/algorithm/isSimple.h"
#include "SFCGAL/io/wkt.h"

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_IsSimple)

BOOST_AUTO_TEST_CASE(pointIsSimple)
{
  std::string const               wkt = "POINT (3.0 4.0)";
  std::unique_ptr<Geometry> const g(io::readWkt(wkt));
  Simplicity const                s = algorithm::isSimple(*g);
  BOOST_CHECK_MESSAGE(s == Simplicity::simple(),
                      (boost::format("%s should be simple: %s") %
                       g->geometryType() % g->asText()));
}

BOOST_AUTO_TEST_CASE(point3DIsSimple)
{
  std::string const               wkt = "POINT (3.0 4.0 2.0)";
  std::unique_ptr<Geometry> const g(io::readWkt(wkt));
  Simplicity const                s = algorithm::isSimple(*g);
  BOOST_CHECK_MESSAGE(s == Simplicity::simple(),
                      (boost::format("%s should be simple: %s") %
                       g->geometryType() % g->asText()));
}

BOOST_AUTO_TEST_CASE(ShortLinestringIsSimple)
{
  std::string const               wkt = "LINESTRING (0.0 0.0, 2.0 0.0)";
  std::unique_ptr<Geometry> const g(io::readWkt(wkt));
  Simplicity const                s = algorithm::isSimple(*g);
  BOOST_CHECK_MESSAGE(s == Simplicity::simple(),
                      (boost::format("%s should be simple: %s") %
                       g->geometryType() % g->asText()));
}

BOOST_AUTO_TEST_CASE(LongLinestringIsSimple)
{
  std::string const wkt = "LINESTRING (0.0 0.0, 2.0 0.0, 1.0 1.0)";
  std::unique_ptr<Geometry> const g(io::readWkt(wkt));
  Simplicity const                s = algorithm::isSimple(*g);
  BOOST_CHECK_MESSAGE(s == Simplicity::simple(),
                      (boost::format("%s should be simple: %s") %
                       g->geometryType() % g->asText()));
}

BOOST_AUTO_TEST_CASE(ComplexLongLinestringIsSimple)
{
  std::string const wkt = "LINESTRING (0.0 0.0, 2.0 0.0, 1.0 1.0, 1.0 -1.0)";
  std::unique_ptr<Geometry> const g(io::readWkt(wkt));
  Simplicity const                s = algorithm::isSimple(*g);
  BOOST_CHECK_MESSAGE(s == Simplicity::complex("linestring self intersects"),
                      (boost::format("%s should be complex: %s") %
                       g->geometryType() % g->asText()));
}

BOOST_AUTO_TEST_CASE(ClosedLongLinestringIsSimple)
{
  std::string const wkt = "LINESTRING (0.0 0.0, 2.0 0.0, 1.0 1.0, 0.0 0.0)";
  std::unique_ptr<Geometry> const g(io::readWkt(wkt));
  Simplicity const                s = algorithm::isSimple(*g);
  BOOST_CHECK_MESSAGE(s == Simplicity::simple(),
                      (boost::format("%s should be simple: %s") %
                       g->geometryType() % g->asText()));
}

BOOST_AUTO_TEST_CASE(PolygonIsSimple)
{
  std::string const wkt = "POLYGON Z ((0.0 0.0 0.0, 1.0 0.0 0.0, 1.0 1.0 0.0, "
                          "0.0 1.0 0.0, 0.0 0.0 0.0))";
  std::unique_ptr<Geometry> const g(io::readWkt(wkt));
  Simplicity const                s = algorithm::isSimple(*g);
  BOOST_CHECK_MESSAGE(s == Simplicity::simple(),
                      (boost::format("%s should be simple: %s") %
                       g->geometryType() % g->asText()));
}

BOOST_AUTO_TEST_CASE(ComplexPolygonIsSimple)
{
  std::string const wkt = "POLYGON Z ((0.0 0.0 1.0, 1.0 0.0 0.0, 1.0 1.0 0.0, "
                          "0.0 1.0 0.0, 0.0 0.0 1.0))";
  std::unique_ptr<Geometry> const g(io::readWkt(wkt));
  Simplicity const                s = algorithm::isSimple(*g);
  BOOST_CHECK_MESSAGE(
      s == Simplicity::complex("Points don't lie in the same plane."),
      (boost::format("%s should be complex: %s") % g->geometryType() %
       g->asText()));
}

BOOST_AUTO_TEST_CASE(TriangleIsSimple)
{
  std::string const wkt =
      "TRIANGLE Z ((0.0 0.0 0.0, 1.0 0.0 0.0, 1.0 1.0 0.0, 0.0 0.0 0.0))";
  std::unique_ptr<Geometry> const g(io::readWkt(wkt));
  Simplicity const                s = algorithm::isSimple(*g);
  BOOST_CHECK_MESSAGE(s == Simplicity::simple(),
                      (boost::format("%s should be simple: %s") %
                       g->geometryType() % g->asText()));
}

BOOST_AUTO_TEST_CASE(PolyhedralSurfaceIsSimple)
{
  std::string const wkt = "POLYGON Z ((0.0 0.0 0.0, 1.0 0.0 0.0, 1.0 1.0 0.0, "
                          "0.0 1.0 0.0, 0.0 0.0 0.0))";
  std::unique_ptr<Geometry> const poly(io::readWkt(wkt));
  std::unique_ptr<Geometry>       g(algorithm::extrude(*poly, 0.0, 0.0, 10.0));
  Simplicity const                s = algorithm::isSimple(*g);
  BOOST_CHECK_MESSAGE(s == Simplicity::simple(),
                      (boost::format("%s should be simple: %s") %
                       g->geometryType() % g->asText()));
}

BOOST_AUTO_TEST_CASE(ComplexPolyhedralSurfaceIsSimple)
{
  std::string const               wkt = "POLYHEDRALSURFACE (\
                 ((0 0 0, 0 -2 0, -1 -1 0, -1 0 0, 0 0 0)),\
                 ((0 0 0, 0 0 -1, 0 -1 -1, 0 -2 0, 0 0 0)),\
                 ((0 0 0, -1 0 0, -1 0 -1, 0 0 -1, 0 0 0)),\
                 ((-1 -1 -1, 0 -1 -1, 0 0 -1, -1 0 -1, -1 -1 -1)),\
                 ((-1 -1 -1, -1 0 -1, -1 0 0, -1 -1 0, -1 -1 -1)),\
                 ((-1 -1 -1, -1 -1 0, 0 -2 0, 0 -1 -1, -1 -1 -1)))";
  std::unique_ptr<Geometry> const g(io::readWkt(wkt));
  Simplicity const                s = algorithm::isSimple(*g);
  BOOST_CHECK_MESSAGE(
      s == Simplicity::complex(
               "Polygon 0 is complex: Points don't lie in the same plane."),
      (boost::format("%s should be complex: %s") % g->geometryType() %
       g->asText()));
}

BOOST_AUTO_TEST_CASE(TriangulatedSurfaceIsSimple)
{
  std::string const wkt =
      "TIN (((0.0 0.0 0.0, 1.0 0.0 0.0, 0.0 1.0 0.0, 0.0 0.0 0.0)), \
                            ((0.0 1.0 0.0, 0.0 0.0 0.0, 0.0 0.0 1.0, 0.0 1.0 0.0)), \
                            ((0.0 0.0 1.0, 0.0 0.0 0.0, 1.0 0.0 0.0, 0.0 0.0 1.0)), \
                            ((1.0 0.0 0.0, 0.0 1.0 0.0, 0.0 0.0 1.0, 1.0 0.0 0.0)))";
  std::unique_ptr<Geometry> const g(io::readWkt(wkt));
  Simplicity const                s = algorithm::isSimple(*g);
  BOOST_CHECK_MESSAGE(s == Simplicity::simple(),
                      (boost::format("%s should be simple: %s") %
                       g->geometryType() % g->asText()));
}

BOOST_AUTO_TEST_CASE(SolidIsSimple)
{
  std::string const wkt = "SOLID((((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)),\
                 ((0 0 0, 0 0 1, 0 1 1, 0 1 0, 0 0 0)),\
                 ((0 0 0, 1 0 0, 1 0 1, 0 0 1, 0 0 0)),\
                 ((1 1 1, 0 1 1, 0 0 1, 1 0 1, 1 1 1)),\
                 ((1 1 1, 1 0 1, 1 0 0, 1 1 0, 1 1 1)),\
                 ((1 1 1, 1 1 0, 0 1 0, 0 1 1, 1 1 1))))";
  std::unique_ptr<Geometry> const g(io::readWkt(wkt));
  Simplicity const                s = algorithm::isSimple(*g);
  BOOST_CHECK_MESSAGE(s == Simplicity::simple(),
                      (boost::format("%s should be simple: %s") %
                       g->geometryType() % g->asText()));
}

BOOST_AUTO_TEST_CASE(ComplexSolidIsSimple)
{
  std::string const               wkt = "SOLID ((\
                 ((0 0 0, 0 -2 0, -1 -1 0, -1 0 0, 0 0 0)),\
                 ((0 0 0, 0 0 -1, 0 -1 -1, 0 -2 0, 0 0 0)),\
                 ((0 0 0, -1 0 0, -1 0 -1, 0 0 -1, 0 0 0)),\
                 ((-1 -1 -1, 0 -1 -1, 0 0 -1, -1 0 -1, -1 -1 -1)),\
                 ((-1 -1 -1, -1 0 -1, -1 0 0, -1 -1 0, -1 -1 -1)),\
                 ((-1 -1 -1, -1 -1 0, 0 -2 0, 0 -1 -1, -1 -1 -1))))";
  std::unique_ptr<Geometry> const g(io::readWkt(wkt));
  Simplicity const                s = algorithm::isSimple(*g);
  BOOST_CHECK_MESSAGE(
      s == Simplicity::complex("Shell 0 is complex: Polygon 0 is complex: "
                               "Points don't lie in the same plane."),
      (boost::format("%s should be complex: %s") % g->geometryType() %
       g->asText()));
}

BOOST_AUTO_TEST_CASE(MultiPointIsSimple)
{
  std::string const wkt =
      "MULTIPOINT ((0.0 0.0), (1.0 0.0), (1.0 1.0), (0.0 1.0))";
  std::unique_ptr<Geometry> const g(io::readWkt(wkt));
  Simplicity const                s = algorithm::isSimple(*g);
  BOOST_CHECK_MESSAGE(s == Simplicity::simple(),
                      (boost::format("%s should be simple: %s") %
                       g->geometryType() % g->asText()));
}

BOOST_AUTO_TEST_CASE(ComplexMultiPointIsSimple)
{
  std::string const wkt =
      "MULTIPOINT ((0.0 0.0), (1.0 0.0), (1.0 1.0), (0.0 1.0), (0.0 0.0))";
  std::unique_ptr<Geometry> const g(io::readWkt(wkt));
  Simplicity const                s = algorithm::isSimple(*g);
  BOOST_CHECK_MESSAGE(
      s == Simplicity::complex(
               "Points 0 and 4 are duplicated in the MultiPoint."),
      (boost::format("%s should be complex: %s") % g->geometryType() %
       g->asText()));
}

BOOST_AUTO_TEST_CASE(MultiLineStringIsSimple)
{
  std::string const wkt =
      "MULTILINESTRING((0.0 0.0, 2.0 0.0), (2.0 0.0, 1.0 1.0))";
  std::unique_ptr<Geometry> const g(io::readWkt(wkt));
  Simplicity const                s = algorithm::isSimple(*g);
  BOOST_CHECK_MESSAGE(s == Simplicity::simple(),
                      (boost::format("%s should be simple: %s") %
                       g->geometryType() % g->asText()));
}

BOOST_AUTO_TEST_CASE(ComplexMultiLineStringIsSimple)
{
  std::string const wkt = "MULTILINESTRING ((-4.0 0.0, -2.0 0.0), \
                             (3.0 0.0, 2.0 1.0), \
                             (0.0 0.0, 2.0 0.0, 1.0 1.0, 1.0 -1.0))";
  std::unique_ptr<Geometry> const g(io::readWkt(wkt));
  Simplicity const                s = algorithm::isSimple(*g);
  BOOST_CHECK_MESSAGE(
      s == Simplicity::complex(
               "LineString 2 is complex: linestring self intersects"),
      (boost::format("%s should be complex: %s") % g->geometryType() %
       g->asText()));
}

BOOST_AUTO_TEST_CASE(MultiPolygonIsSimple)
{
  std::string const wkt = "MULTIPOLYGON (((0.0 0.0 0.0, 1.0 0.0 0.0, 1.0 1.0 "
                          "0.0, 0.0 1.0 0.0, 0.0 0.0 0.0)))";
  std::unique_ptr<Geometry> const g(io::readWkt(wkt));
  Simplicity const                s = algorithm::isSimple(*g);
  BOOST_CHECK_MESSAGE(s == Simplicity::simple(),
                      (boost::format("%s should be simple: %s") %
                       g->geometryType() % g->asText()));
}

BOOST_AUTO_TEST_CASE(ComplexMultiPolygonIsSimple)
{
  std::string const wkt = "MULTIPOLYGON (((0.0 0.0 0.0, 1.0 0.0 0.0, 1.0 1.0 "
                          "0.0, 0.0 1.0 0.0, 0.0 0.0 0.0)), ((0.0 0.0 1.0, 1.0 "
                          "0.0 0.0, 1.0 1.0 0.0, 0.0 1.0 0.0, 0.0 0.0 1.0)))";
  std::unique_ptr<Geometry> const g(io::readWkt(wkt));
  Simplicity const                s = algorithm::isSimple(*g);
  BOOST_CHECK_MESSAGE(
      s == Simplicity::complex(
               "Polygon 1 is complex: Points don't lie in the same plane."),
      (boost::format("%s should be complex: %s") % g->geometryType() %
       g->asText()));
}

BOOST_AUTO_TEST_CASE(MultiSolidIsSimple)
{
  std::string const               wkt = "MULTISOLID((\
                 (((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)),\
                 ((0 0 0, 0 0 1, 0 1 1, 0 1 0, 0 0 0)),\
                 ((0 0 0, 1 0 0, 1 0 1, 0 0 1, 0 0 0)),\
                 ((1 1 1, 0 1 1, 0 0 1, 1 0 1, 1 1 1)),\
                 ((1 1 1, 1 0 1, 1 0 0, 1 1 0, 1 1 1)),\
                 ((1 1 1, 1 1 0, 0 1 0, 0 1 1, 1 1 1))),\
                 (((0 0 0, 0 -1 0, -1 -1 0, -1 0 0, 0 0 0)),\
                 ((0 0 0, 0 0 -1, 0 -1 -1, 0 -1 0, 0 0 0)),\
                 ((0 0 0, -1 0 0, -1 0 -1, 0 0 -1, 0 0 0)),\
                 ((-1 -1 -1, 0 -1 -1, 0 0 -1, -1 0 -1, -1 -1 -1)),\
                 ((-1 -1 -1, -1 0 -1, -1 0 0, -1 -1 0, -1 -1 -1)),\
                 ((-1 -1 -1, -1 -1 0, 0 -1 0, 0 -1 -1, -1 -1 -1)))\
                 ))";
  std::unique_ptr<Geometry> const g(io::readWkt(wkt));
  Simplicity const                s = algorithm::isSimple(*g);
  BOOST_CHECK_MESSAGE(s == Simplicity::simple(),
                      (boost::format("%s should be simple: %s") %
                       g->geometryType() % g->asText()));
}

BOOST_AUTO_TEST_CASE(ComplexMultiSolidIsSimple)
{
  std::string const               wkt = "MULTISOLID(\
                 ((((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)),\
                 ((0 0 0, 0 0 1, 0 1 1, 0 1 0, 0 0 0)),\
                 ((0 0 0, 1 0 0, 1 0 1, 0 0 1, 0 0 0)),\
                 ((1 1 1, 0 1 1, 0 0 1, 1 0 1, 1 1 1)),\
                 ((1 1 1, 1 0 1, 1 0 0, 1 1 0, 1 1 1)),\
                 ((1 1 1, 1 1 0, 0 1 0, 0 1 1, 1 1 1)))),\
                 ((((0 0 0, 0 -2 0, -1 -1 0, -1 0 0, 0 0 0)),\
                 ((0 0 0, 0 0 -1, 0 -1 -1, 0 -2 0, 0 0 0)),\
                 ((0 0 0, -1 0 0, -1 0 -1, 0 0 -1, 0 0 0)),\
                 ((-1 -1 -1, 0 -1 -1, 0 0 -1, -1 0 -1, -1 -1 -1)),\
                 ((-1 -1 -1, -1 0 -1, -1 0 0, -1 -1 0, -1 -1 -1)),\
                 ((-1 -1 -1, -1 -1 0, 0 -2 0, 0 -1 -1, -1 -1 -1))))\
                 )";
  std::unique_ptr<Geometry> const g(io::readWkt(wkt));
  Simplicity const                s = algorithm::isSimple(*g);
  BOOST_CHECK_MESSAGE(
      s == Simplicity::complex(
               "Solid 1 is complex: Points don't lie in the same plane."),
      (boost::format("%s should be complex: %s") % g->geometryType() %
       g->asText()));
}

BOOST_AUTO_TEST_CASE(GeometryCollectionIsSimple)
{
  std::string const wkt = "GEOMETRYCOLLECTION (POINT (2.0 3.0), TRIANGLE ((0.0 "
                          "0.0,1.0 0.0,1.0 1.0,0.0 0.0)))";
  std::unique_ptr<Geometry> const g(io::readWkt(wkt));
  Simplicity const                s = algorithm::isSimple(*g);
  BOOST_CHECK_MESSAGE(s == Simplicity::simple(),
                      (boost::format("%s should be simple: %s") %
                       g->geometryType() % g->asText()));
}

BOOST_AUTO_TEST_CASE(ComplexGeometryCollectionIsSimple)
{
  std::string const wkt =
      "GEOMETRYCOLLECTION (POINT (2.0 3.0), TRIANGLE ((0.0 0.0,1.0 0.0,1.0 1.0,0.0 0.0)), \
                             LINESTRING (0.0 0.0, 2.0 0.0, 1.0 1.0, 1.0 -1.0))";
  std::unique_ptr<Geometry> const g(io::readWkt(wkt));
  Simplicity const                s = algorithm::isSimple(*g);
  BOOST_CHECK_MESSAGE(
      s == Simplicity::complex(
               "LineString at index 2 is complex: linestring self intersects."),
      (boost::format("%s should be complex: %s") % g->geometryType() %
       g->asText()));
}

BOOST_AUTO_TEST_SUITE_END()
