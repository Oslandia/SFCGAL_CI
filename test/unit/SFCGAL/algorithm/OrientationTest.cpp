// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/tools/old/interface.hpp>
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
#include "SFCGAL/algorithm/orientation.h"
#include "SFCGAL/detail/transform/ForceOrderPoints.h"
#include "SFCGAL/io/wkt.h"

using namespace SFCGAL;

#include "../../../test_config.h"
// always after CGAL
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_OrientationTest)

// TODO COMPLETE TEST

// void makeValidOrientation( CGAL::Polygon_2< Kernel > & polygon ) ;
// void makeValidOrientation( CGAL::Polygon_with_holes_2< Kernel > & polygon ) ;
// void makeValidOrientation( Polygon & polygon ) ;

// bool hasConsistentOrientation3D( const TriangulatedSurface & g ) ;
BOOST_AUTO_TEST_CASE(testHasConsistentOrientation3D_basicTriangles)
{
  TriangulatedSurface triangulatedSurface;
  BOOST_CHECK(algorithm::hasConsistentOrientation3D(triangulatedSurface));
  triangulatedSurface.addPatch(Triangle(
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0)));
  BOOST_CHECK(algorithm::hasConsistentOrientation3D(triangulatedSurface));

  triangulatedSurface.addPatch(Triangle(
      Point(0.0, 0.0, 0.0), Point(0.0, 1.0, 0.0), Point(-1.0, 0.0, 0.0)));
  BOOST_CHECK(algorithm::hasConsistentOrientation3D(triangulatedSurface));

  triangulatedSurface.addPatch(Triangle(
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, -1.0, 0.0)));
  BOOST_CHECK(!algorithm::hasConsistentOrientation3D(triangulatedSurface));
}

// bool hasConsistentOrientation3D( const PolyhedralSurface & g ) ;

BOOST_AUTO_TEST_CASE(testHasConsistentOrientation3D_basicPolygons)
{
  PolyhedralSurface polyhedralSurface;
  BOOST_CHECK(algorithm::hasConsistentOrientation3D(polyhedralSurface));

  // base polygon
  {
    LineString ring;
    ring.addPoint(Point(0.0, 0.0, 0.0));
    ring.addPoint(Point(1.0, 0.0, 0.0));
    ring.addPoint(Point(1.0, 1.0, 0.0));
    ring.addPoint(Point(0.0, 1.0, 0.0));
    ring.addPoint(ring.startPoint());

    polyhedralSurface.addPatch(Polygon(ring));
  }
  BOOST_CHECK(algorithm::hasConsistentOrientation3D(polyhedralSurface));

  // right polygon
  {
    LineString ring;
    ring.addPoint(Point(1.0, 0.0, 0.0));
    ring.addPoint(Point(1.0, 0.0, 1.0));
    ring.addPoint(Point(1.0, 1.0, 1.0));
    ring.addPoint(Point(1.0, 1.0, 0.0));
    ring.addPoint(ring.startPoint());

    polyhedralSurface.addPatch(Polygon(ring));
  }
  BOOST_CHECK(algorithm::hasConsistentOrientation3D(polyhedralSurface));

  // left polygon (bad orientation)
  {
    LineString ring;
    ring.addPoint(Point(0.0, 0.0, 0.0));
    ring.addPoint(Point(0.0, 0.0, 1.0));
    ring.addPoint(Point(0.0, 1.0, 1.0));
    ring.addPoint(Point(0.0, 1.0, 0.0));
    ring.addPoint(ring.startPoint());

    polyhedralSurface.addPatch(Polygon(ring));
  }
  BOOST_CHECK(!algorithm::hasConsistentOrientation3D(polyhedralSurface));
}

// void makeConsistentOrientation3D( TriangulatedSurface & g ) ;

// bool isCounterClockWiseOriented( const Polygon& );
BOOST_AUTO_TEST_CASE(testIsCounterClockWiseOriented_Polygon)
{
  LineString ring;
  ring.addPoint(Point(0.0, 0.0));
  ring.addPoint(Point(1.0, 0.0));
  ring.addPoint(Point(1.0, 1.0));
  ring.addPoint(Point(0.0, 1.0));
  ring.addPoint(ring.startPoint());

  Polygon polygon(ring);
  BOOST_CHECK(algorithm::isCounterClockWiseOriented(polygon));
  polygon.exteriorRing().reverse();
  BOOST_CHECK(!algorithm::isCounterClockWiseOriented(polygon));
}

// Test for empty polygon - Issue #315
BOOST_AUTO_TEST_CASE(testIsCounterClockWiseOriented_EmptyPolygon)
{
  // Empty polygon should not crash and should return false
  Polygon emptyPolygon;
  BOOST_CHECK(emptyPolygon.isEmpty());
  BOOST_CHECK(!algorithm::isCounterClockWiseOriented(emptyPolygon));

  // Also test with WKT parsing
  std::unique_ptr<Geometry> geom(io::readWkt("POLYGON EMPTY"));
  BOOST_REQUIRE(geom->is<Polygon>());
  BOOST_CHECK(geom->isEmpty());
  BOOST_CHECK(!algorithm::isCounterClockWiseOriented(geom->as<Polygon>()));
}

// bool isCounterClockWiseOriented( const Triangle& );
BOOST_AUTO_TEST_CASE(testIsCounterClockWiseOriented_Triangle)
{
  Triangle triangle(Point(0.0, 0.0), Point(1.0, 0.0), Point(1.0, 1.0));
  BOOST_CHECK(algorithm::isCounterClockWiseOriented(triangle));
  triangle.reverse();
  BOOST_CHECK(!algorithm::isCounterClockWiseOriented(triangle));
}

// bool isCounterClockWiseOriented( const LineString& );
BOOST_AUTO_TEST_CASE(testIsCounterClockWiseOriented_LineString)
{
  LineString ring;
  ring.addPoint(Point(0.0, 0.0));
  ring.addPoint(Point(1.0, 0.0));
  ring.addPoint(Point(1.0, 1.0));
  ring.addPoint(Point(0.0, 1.0));
  ring.addPoint(ring.startPoint());

  BOOST_CHECK(algorithm::isCounterClockWiseOriented(ring));
  ring.reverse();
  BOOST_CHECK(!algorithm::isCounterClockWiseOriented(ring));
}

BOOST_AUTO_TEST_CASE(testPostGISIssue3651)
{
  std::string inputData(SFCGAL_TEST_DIRECTORY);
  inputData += "/data/postgis_crash_orientation_3651.wkt";
  std::ifstream ifs(inputData.c_str());
  BOOST_REQUIRE(ifs.good());

  std::string inputWkt;
  std::getline(ifs, inputWkt);
  std::unique_ptr<Geometry> g(io::readWkt(inputWkt));
  algorithm::isCounterClockWiseOriented(g->as<Polygon>());
}

BOOST_AUTO_TEST_CASE(testPostGISIssue6029)
{
  std::unique_ptr<Geometry> geom(
      io::readWkt("MULTIPOLYGON (((9 9, 9 1, 1 1, 2 4, 7 7, 9 9)), EMPTY)"));
  SFCGAL::transform::ForceOrderPoints force(/* ccw */ true);
  geom->accept(force);
  BOOST_CHECK_EQUAL(geom->asText(0),
                    "MULTIPOLYGON (((9 9,7 7,2 4,1 1,9 1,9 9)))");
}

BOOST_AUTO_TEST_SUITE_END()
