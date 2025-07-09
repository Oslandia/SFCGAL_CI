// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <memory>
#include <string>

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
#include "SFCGAL/io/wkt.h"

#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

using namespace SFCGAL;
using namespace SFCGAL::io;

BOOST_AUTO_TEST_SUITE(SFCGAL_io_WktReaderTest)

//-- WKT POINT

BOOST_AUTO_TEST_CASE(pointEmpty)
{
  std::unique_ptr<Geometry> g(readWkt("POINT EMPTY"));
  BOOST_CHECK(g->is<Point>());
  BOOST_CHECK(g->isEmpty());
}

BOOST_AUTO_TEST_CASE(pointXY)
{
  std::unique_ptr<Geometry> g(readWkt("POINT (4.0 6.0)"));
  BOOST_CHECK(g->is<Point>());
  BOOST_CHECK(!g->isEmpty());

  BOOST_CHECK_EQUAL(g->as<Point>().x(), 4.0);
  BOOST_CHECK_EQUAL(g->as<Point>().y(), 6.0);
}

BOOST_AUTO_TEST_CASE(pointXYZ_implicit)
{
  std::unique_ptr<Geometry> g(readWkt("POINT (4.0 5.0 6.0)"));
  BOOST_CHECK(g->is<Point>());
  BOOST_CHECK(!g->isEmpty());

  BOOST_CHECK(g->is3D());
  BOOST_CHECK(!g->isMeasured());

  BOOST_CHECK_EQUAL(g->as<Point>().x(), 4.0);
  BOOST_CHECK_EQUAL(g->as<Point>().y(), 5.0);
  BOOST_CHECK_EQUAL(g->as<Point>().z(), 6.0);
}

BOOST_AUTO_TEST_CASE(pointXYZ_explicit)
{
  std::unique_ptr<Geometry> g(readWkt("POINT Z (4.0 5.0 6.0)"));
  BOOST_CHECK(g->is<Point>());
  BOOST_CHECK(!g->isEmpty());

  BOOST_CHECK(g->is3D());
  BOOST_CHECK(!g->isMeasured());

  BOOST_CHECK_EQUAL(g->as<Point>().x(), 4.0);
  BOOST_CHECK_EQUAL(g->as<Point>().y(), 5.0);
  BOOST_CHECK_EQUAL(g->as<Point>().z(), 6.0);
}

BOOST_AUTO_TEST_CASE(pointXYM_explicit)
{
  std::unique_ptr<Geometry> g(readWkt("POINT M (4.0 5.0 6.0)"));
  BOOST_CHECK(g->is<Point>());
  BOOST_CHECK(!g->isEmpty());

  BOOST_CHECK(!g->is3D());
  BOOST_CHECK(g->isMeasured());

  BOOST_CHECK_EQUAL(g->as<Point>().x(), 4.0);
  BOOST_CHECK_EQUAL(g->as<Point>().y(), 5.0);
  BOOST_CHECK_EQUAL(g->as<Point>().m(), 6.0);
}

//-- WKT LINESTRING

BOOST_AUTO_TEST_CASE(lineStringEmpty)
{
  std::unique_ptr<Geometry> g(readWkt("LINESTRING EMPTY"));
  BOOST_CHECK(g->is<LineString>());
  BOOST_CHECK(g->isEmpty());
}

BOOST_AUTO_TEST_CASE(lineString_twoPoints)
{
  std::unique_ptr<Geometry> g(readWkt("LINESTRING (0.0 0.0,1.0 1.0)"));
  BOOST_CHECK(g->is<LineString>());
  BOOST_CHECK(!g->isEmpty());
  BOOST_CHECK_EQUAL(g->as<LineString>().numPoints(), 2U);
}

BOOST_AUTO_TEST_CASE(lineString_twoPoints3D)
{
  std::unique_ptr<Geometry> g(readWkt("LINESTRING (0.0 0.0 0.0,1.0 1.0 1.0)"));
  BOOST_CHECK(g->is<LineString>());
  BOOST_CHECK(!g->isEmpty());
  BOOST_REQUIRE_EQUAL(g->as<LineString>().numPoints(), 2U);
  BOOST_CHECK(g->as<LineString>().pointN(0).is3D());
  BOOST_CHECK(g->as<LineString>().pointN(1).is3D());
}

//-- WKT POLYGON

BOOST_AUTO_TEST_CASE(polygonEmpty)
{
  std::unique_ptr<Geometry> g(readWkt("POLYGON EMPTY"));
  BOOST_CHECK(g->is<Polygon>());
  BOOST_CHECK(g->isEmpty());
}

// 4 points polygon (triangle)
BOOST_AUTO_TEST_CASE(polygonWithFourPoints)
{
  std::unique_ptr<Geometry> g(readWkt("POLYGON ((0 0,1 0,1 1,0 0))"));
  BOOST_CHECK(g->is<Polygon>());
  BOOST_CHECK(!g->isEmpty());
  BOOST_CHECK_EQUAL(g->as<Polygon>().exteriorRing().numPoints(), 4U);
}

//-- WKT MULTIPOINT

BOOST_AUTO_TEST_CASE(multiPointEmpty)
{
  std::unique_ptr<Geometry> g(readWkt("MULTIPOINT EMPTY"));
  BOOST_CHECK(g->is<MultiPoint>());
  BOOST_CHECK(g->isEmpty());
}

BOOST_AUTO_TEST_CASE(multiPointEmpty2)
{
  std::unique_ptr<Geometry> g(readWkt("MULTIPOINT (0 0,1 1,EMPTY)"));
  BOOST_CHECK(g->asText() == "MULTIPOINT ((0/1 0/1),(1/1 1/1))");
  BOOST_CHECK(g->is<MultiPoint>());
  BOOST_CHECK(g->numGeometries() == 2);
}

BOOST_AUTO_TEST_CASE(multiPointEmpty3)
{
  std::unique_ptr<Geometry> g(readWkt("MULTIPOINT (EMPTY,EMPTY)"));
  BOOST_CHECK(g->asText() == "MULTIPOINT EMPTY");
  BOOST_CHECK(g->is<MultiPoint>());
  BOOST_CHECK(g->isEmpty());
}

//-- WKT MULTILINESTRING

BOOST_AUTO_TEST_CASE(multiLineStringEmpty)
{
  std::unique_ptr<Geometry> g(readWkt("MULTILINESTRING EMPTY"));
  BOOST_CHECK(g->is<MultiLineString>());
  BOOST_CHECK(g->isEmpty());
}

//-- WKT MULTIPOLYGON

BOOST_AUTO_TEST_CASE(multiPolygonEmpty)
{
  std::unique_ptr<Geometry> g(readWkt("MULTIPOLYGON EMPTY"));
  BOOST_CHECK(g->is<MultiPolygon>());
  BOOST_CHECK(g->isEmpty());
}

//-- WKT GEOMETRYCOLLECTION

BOOST_AUTO_TEST_CASE(geometryCollectionEmpty)
{
  std::unique_ptr<Geometry> g(readWkt("GEOMETRYCOLLECTION EMPTY"));
  BOOST_CHECK(g->is<GeometryCollection>());
  BOOST_CHECK(g->isEmpty());
}

//-- WKT TRIANGULATEDSURFACE

BOOST_AUTO_TEST_CASE(triangulatedSurface_Empty)
{
  std::unique_ptr<Geometry> g(readWkt("TIN EMPTY"));
  BOOST_CHECK(g->is<TriangulatedSurface>());
  BOOST_CHECK(g->isEmpty());
}

BOOST_AUTO_TEST_CASE(triangulatedSurface_fourTriangles)
{
  std::string const         wkt = "TIN ("
                                  "((0 0 0, 0 0 1, 0 1 0, 0 0 0)),"
                                  "((0 0 0, 0 1 0, 1 0 0, 0 0 0)),"
                                  "((0 0 0, 1 0 0, 0 0 1, 0 0 0)),"
                                  "((1 0 0, 0 1 0, 0 0 1, 1 0 0))"
                                  ")";
  std::unique_ptr<Geometry> g(readWkt(wkt));
  BOOST_CHECK(g->is<TriangulatedSurface>());
  BOOST_CHECK(!g->isEmpty());

  BOOST_CHECK_EQUAL(g->as<TriangulatedSurface>().numGeometries(), 1U);
  BOOST_CHECK_EQUAL(g->as<TriangulatedSurface>().numPatches(), 4U);
}

BOOST_AUTO_TEST_CASE(wkt_exactTest)
{
  std::unique_ptr<Geometry> g(readWkt("LINESTRING (2/3 3/2,5/4 2/3)"));
  BOOST_CHECK(g->is<LineString>());
  BOOST_CHECK(!g->isEmpty());
  BOOST_REQUIRE_EQUAL(g->as<LineString>().numPoints(), 2U);
  Kernel::Exact_kernel::FT const x =
      CGAL::exact(g->as<LineString>().pointN(0).x());
  Kernel::Exact_kernel::FT const y =
      CGAL::exact(g->as<LineString>().pointN(0).y());

  CGAL::Fraction_traits<Kernel::Exact_kernel::FT>::Numerator_type xn;
  CGAL::Fraction_traits<Kernel::Exact_kernel::FT>::Numerator_type xd;
  CGAL::Fraction_traits<Kernel::Exact_kernel::FT>::Numerator_type yn;
  CGAL::Fraction_traits<Kernel::Exact_kernel::FT>::Numerator_type yd;
  CGAL::Fraction_traits<Kernel::Exact_kernel::FT>::Decompose      decomp;
  decomp(x, xn, xd);
  decomp(y, yn, yd);

  BOOST_CHECK_EQUAL(xn, 2);
  BOOST_CHECK_EQUAL(xd, 3);
  BOOST_CHECK_EQUAL(yn, 3);
  BOOST_CHECK_EQUAL(yd, 2);
}

BOOST_AUTO_TEST_CASE(charArrayRead)
{
  char                      str[] = "LINESTRING (0.0 0.0,1.0 1.0)";
  std::unique_ptr<Geometry> g(readWkt(str, strlen(str)));
  BOOST_CHECK(g->is<LineString>());
  BOOST_CHECK(!g->isEmpty());
  BOOST_CHECK_EQUAL(g->as<LineString>().numPoints(), 2U);
}

BOOST_AUTO_TEST_CASE(wktExtraCharacters)
{
  bool threw = false;
  try {
    std::unique_ptr<Geometry> const g(readWkt("POINT (0 0)POINT (1 0)"));
  } catch (WktParseException &e) {
    std::string const err(e.what());
    BOOST_CHECK_EQUAL(err, "Extra characters in WKT: POINT (1 0)");
    threw = true;
  }
  BOOST_CHECK(threw);

  threw = false;
  try {
    char                            str[] = "POINT (0 0)POINT (1 0)";
    std::unique_ptr<Geometry> const g(readWkt(str, strlen(str)));
  } catch (WktParseException &e) {
    std::string const err(e.what());
    BOOST_CHECK_EQUAL(err, "Extra characters in WKT: POINT (1 0)");
    threw = true;
  }
  BOOST_CHECK(threw);
}

BOOST_AUTO_TEST_SUITE_END()
