// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <fstream>

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

#include "../test_config.h"
#include "Bench.h"

#include <boost/format.hpp>
#include <boost/test/unit_test.hpp>

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_BenchWKT)

//
// Test limit case
BOOST_AUTO_TEST_CASE(testReadPoints)
{
  const int N = 100000;

  bench().start(boost::format("READ WKT POINTS"));

  for (int i = 0; i < N; i++) {
    io::readWkt("POINT (50000 50000)");
  }

  bench().stop();
}

//
// Test limit case
BOOST_AUTO_TEST_CASE(testReadLineString)
{
  const int N = 100000;

  bench().start(boost::format("READ WKT LINESTRING"));

  for (int i = 0; i < N; i++) {
    io::readWkt("LINESTRING (0 0,0 1000,1000 1000)");
  }

  bench().stop();
}

//
// Test limit case
BOOST_AUTO_TEST_CASE(testReadPolygon)
{
  const int N = 100000;

  bench().start(boost::format("READ WKT POLYGON"));

  for (int i = 0; i < N; i++) {
    io::readWkt("POLYGON ((0 0,0 1000,1000 1000,1000 0,0 0))");
  }

  bench().stop();
}

//
// Test limit case
BOOST_AUTO_TEST_CASE(testReadTriangle)
{
  const int N = 100000;

  bench().start(boost::format("READ WKT TRIANGLE"));

  for (int i = 0; i < N; i++) {
    io::readWkt("TRIANGLE ((0 0,0 1000,1000 1000,0 0))");
  }

  bench().stop();
}

//
// Test limit case
BOOST_AUTO_TEST_CASE(testReadSolid)
{
  const int N = 10000;

  bench().start(boost::format("READ WKT SOLID"));

  for (int i = 0; i < N; i++) {
    io::readWkt(
        "SOLID ((((0.0 0.0 0.0,0.0 1.0 0.0,1.0 1.0 0.0,1.0 0.0 0.0,0.0 0.0 "
        "0.0),(0.2 0.2 0.0,0.8 0.2 0.0,0.8 0.8 0.0,0.2 0.8 0.0,0.2 0.2 "
        "0.0)),((0.0 0.0 1.0,1.0 0.0 1.0,1.0 1.0 1.0,0.0 1.0 1.0,0.0 0.0 "
        "1.0),(0.2 0.2 1.0,0.2 0.8 1.0,0.8 0.8 1.0,0.8 0.2 1.0,0.2 0.2 "
        "1.0)),((0.0 0.0 0.0,0.0 0.0 1.0,0.0 1.0 1.0,0.0 1.0 0.0,0.0 0.0 "
        "0.0)),((0.0 1.0 0.0,0.0 1.0 1.0,1.0 1.0 1.0,1.0 1.0 0.0,0.0 1.0 "
        "0.0)),((1.0 1.0 0.0,1.0 1.0 1.0,1.0 0.0 1.0,1.0 0.0 0.0,1.0 1.0 "
        "0.0)),((1.0 0.0 0.0,1.0 0.0 1.0,0.0 0.0 1.0,0.0 0.0 0.0,1.0 0.0 "
        "0.0)),((0.2 0.2 0.0,0.2 0.2 1.0,0.8 0.2 1.0,0.8 0.2 0.0,0.2 0.2 "
        "0.0)),((0.8 0.2 0.0,0.8 0.2 1.0,0.8 0.8 1.0,0.8 0.8 0.0,0.8 0.2 "
        "0.0)),((0.8 0.8 0.0,0.8 0.8 1.0,0.2 0.8 1.0,0.2 0.8 0.0,0.8 0.8 "
        "0.0)),((0.2 0.8 0.0,0.2 0.8 1.0,0.2 0.2 1.0,0.2 0.2 0.0,0.2 0.8 "
        "0.0))))");
  }

  bench().stop();
}

BOOST_AUTO_TEST_SUITE_END()
