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
#include "SFCGAL/io/wkt.h"

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_BoundaryTest)

//-- Point

BOOST_AUTO_TEST_CASE(emptyPoint)
{
  Point const               g;
  std::unique_ptr<Geometry> boundary(g.boundary());

  BOOST_CHECK(boundary->is<GeometryCollection>());
  BOOST_CHECK(boundary->isEmpty());
}
BOOST_AUTO_TEST_CASE(pointXY)
{
  Point const               g(2.0, 3.0);
  std::unique_ptr<Geometry> boundary(g.boundary());

  BOOST_CHECK(boundary->is<GeometryCollection>());
  BOOST_CHECK(boundary->isEmpty());
}

//-- LineString

BOOST_AUTO_TEST_CASE(emptyLineString)
{
  LineString const          g;
  std::unique_ptr<Geometry> boundary(g.boundary());

  BOOST_CHECK(boundary->is<GeometryCollection>());
  BOOST_CHECK(boundary->isEmpty());
}

BOOST_AUTO_TEST_CASE(lineStringSegment)
{
  BOOST_TEST_MESSAGE("check that a LineSegment has a MultiPoint boundary "
                     "composed of 2 points");

  LineString g;
  g.addPoint(Point(0.0, 0.0));
  g.addPoint(Point(1.0, 1.0));
  std::unique_ptr<Geometry> boundary(g.boundary());

  BOOST_CHECK(boundary->is<MultiPoint>());
  BOOST_CHECK_EQUAL(boundary->as<MultiPoint>().numGeometries(), 2U);
}

BOOST_AUTO_TEST_CASE(lineStringWithThreePoints)
{
  BOOST_TEST_MESSAGE("check that a LineString with 3 points has a MultiPoint "
                     "boundary composed of 2 points");

  LineString g;
  g.addPoint(Point(0.0, 0.0));
  g.addPoint(Point(1.0, 1.0));
  g.addPoint(Point(2.0, 2.0));
  std::unique_ptr<Geometry> boundary(g.boundary());

  BOOST_CHECK(boundary->is<MultiPoint>());
  BOOST_CHECK_EQUAL(boundary->as<MultiPoint>().numGeometries(), 2U);
}

BOOST_AUTO_TEST_CASE(lineStringClosed)
{
  BOOST_TEST_MESSAGE("check that a closed LineString has empty boundary");

  LineString g;
  g.addPoint(Point(0.0, 0.0));
  g.addPoint(Point(1.0, 1.0));
  g.addPoint(Point(2.0, 2.0));
  g.addPoint(Point(0.0, 0.0));

  std::unique_ptr<Geometry> boundary(g.boundary());

  BOOST_CHECK(boundary->isEmpty());
  BOOST_CHECK(boundary->is<GeometryCollection>());
}

//-- MultiLineString

BOOST_AUTO_TEST_CASE(multiLineStringSimple)
{
  BOOST_TEST_MESSAGE("check common point between two LineStrings");

  MultiLineString g;
  g.addGeometry(LineString(Point(0.0, 0.0), Point(1.0, 1.0)));
  g.addGeometry(LineString(Point(0.0, 0.0), Point(1.0, 2.0)));

  g.addGeometry(LineString(Point(1.0, 1.0), Point(1.0, 2.0)));

  std::unique_ptr<Geometry> boundary(g.boundary());

  BOOST_CHECK(boundary->isEmpty());
  BOOST_CHECK(boundary->is<GeometryCollection>());
}

//-- MultiPolygon

// TODO TDD

BOOST_AUTO_TEST_CASE(multiPolygonSimple)
{
  BOOST_TEST_MESSAGE("check common point between two LineStrings");

  MultiPolygon g;
  {
    LineString ring;
    ring.addPoint(Point(0.0, 0.0));
    ring.addPoint(Point(1.0, 0.0));
    ring.addPoint(Point(1.0, 1.0));
    ring.addPoint(Point(0.0, 1.0));
    ring.addPoint(Point(0.0, 0.0));

    g.addGeometry(Polygon(ring));
  }
  {
    LineString ring;
    ring.addPoint(Point(1.0, 0.0));
    ring.addPoint(Point(2.0, 0.0));
    ring.addPoint(Point(2.0, 1.0));
    ring.addPoint(Point(1.0, 1.0));
    ring.addPoint(Point(1.0, 0.0));

    g.addGeometry(Polygon(ring));
  }

  std::unique_ptr<Geometry> boundary(g.boundary());

  BOOST_CHECK(!boundary->isEmpty());
  BOOST_CHECK(boundary->is<MultiLineString>());
}

BOOST_AUTO_TEST_SUITE_END()
