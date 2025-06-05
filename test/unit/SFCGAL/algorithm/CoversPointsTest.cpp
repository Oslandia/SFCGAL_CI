// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include <cmath>

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
#include "SFCGAL/detail/algorithm/coversPoints.h"
#include "SFCGAL/io/wkt.h"

using namespace SFCGAL;
using namespace SFCGAL::detail;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_CoversPointsTest)

BOOST_AUTO_TEST_CASE(testPointPointCoversPoints)
{
  Point const pta(0.0, 1.0, 0.0);
  Point const ptb(0.0, 1.0, 0.0);
  Point const ptc(0.0, 0.0, 0.0);
  BOOST_CHECK_EQUAL(algorithm::coversPoints(pta, ptb), true);
  BOOST_CHECK_EQUAL(algorithm::coversPoints(pta, ptc), false);
  BOOST_CHECK_EQUAL(algorithm::coversPoints3D(pta, ptb), true);
  BOOST_CHECK_EQUAL(algorithm::coversPoints3D(pta, ptc), false);
}

BOOST_AUTO_TEST_CASE(testPolygonPolygonCoversPoints)
{
  {
    std::unique_ptr<Geometry> const p1 = io::readWkt(
        "POLYGON ((-1.0 -1.0,1.0 -1.0,1.0 1.0,-1.0 1.0,-1.0 -1.0))");
    std::unique_ptr<Geometry> const p2 = io::readWkt(
        "POLYGON ((-0.5 -0.5,-0.5 0.5,0.5 0.5,0.5 -0.5,-0.5 -0.5))");

    BOOST_CHECK_EQUAL(algorithm::coversPoints(*p1, *p2), true);
    BOOST_CHECK_EQUAL(algorithm::coversPoints3D(*p1, *p2), true);
  }

  {
    // a square with a substracted triangle => concave shape
    std::unique_ptr<Geometry> const p1 =
        io::readWkt("POLYGON ((0.4 0,0 0,0 1,1 1,1 0,0.6 0,0.5 0.4,0.4 0))");
    // a smaller square
    std::unique_ptr<Geometry> const p2 =
        io::readWkt("POLYGON ((0.2 0.2,0.8 0.2,0.8 0.8,0.2 0.8,0.2 0.2))");

    // ST_covers would answer false
    BOOST_CHECK_EQUAL(algorithm::coversPoints(*p1, *p2), true);
    BOOST_CHECK_EQUAL(algorithm::coversPoints3D(*p1, *p2), true);
  }
}

BOOST_AUTO_TEST_CASE(testCollectionCoversPoints)
{
#if 0
    {
        std::unique_ptr<Geometry> p1 = io::readWkt( "GEOMETRYCOLLECTION (LINESTRING (0 0 0,0 0 1),LINESTRING (0 0 0,0 1 0),LINESTRING (0 1 1,0 0 1),LINESTRING (0 1 1,0 1 0))" );
        std::unique_ptr<Geometry> p2 = io::readWkt( "TIN (((0 0.5 0.5,0 0 1,0 0 0,0 0.5 0.5)),((0 0 1,0 0.5 0.5,0 1 1,0 0 1)),((0 0.5 0.5,0 0 0,0 1 0,0 0.5 0.5)),((0 0.5 0.5,0 1 0,0 1 1,0 0.5 0.5)))" );
        std::cout << "p1 covers p2 ? " << algorithm::coversPoints3D( *p1, *p2 ) << std::endl;
        std::cout << "p2 covers p1 ? " << algorithm::coversPoints3D( *p2, *p1 ) << std::endl;
    }

    {
        std::unique_ptr<Geometry> p1 = io::readWkt( "GEOMETRYCOLLECTION (TRIANGLE ((1 1,0.5 0.5,0 1,1 1)),TRIANGLE ((1 0,0.5 0.5,1 1,1 0)),TRIANGLE ((0.5 0.5,0 0,0 1,0.5 0.5)))" );
        std::unique_ptr<Geometry> p2 = io::readWkt( "TRIANGLE ((1 0, 0 0,0.5 0.5,1 0))" );
        std::cout << "p1 covers p2 ? " << algorithm::coversPoints3D( *p1, *p2 ) << std::endl;
        std::cout << "p2 covers p1 ? " << algorithm::coversPoints3D( *p2, *p1 ) << std::endl;
    }
#endif
}

BOOST_AUTO_TEST_SUITE_END()
