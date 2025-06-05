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
#include "SFCGAL/io/wkt.h"

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_CoordinateTest)

BOOST_AUTO_TEST_CASE(testRoundPoint)
{
  std::unique_ptr<Geometry> g(io::readWkt("POINT (1.5 2.6 3.4)"));
  g->round();
  BOOST_CHECK_EQUAL(g->asText(), "POINT Z (2/1 3/1 3/1)");
}

BOOST_AUTO_TEST_CASE(testRoundLineString)
{
  std::unique_ptr<Geometry> g(io::readWkt("LINESTRING (0.5 0.5,1.5 1.5)"));
  g->round(10);
  BOOST_CHECK_EQUAL(g->asText(), "LINESTRING (1/2 1/2,3/2 3/2)");
}

BOOST_AUTO_TEST_SUITE_END()
