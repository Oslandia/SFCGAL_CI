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
#include "SFCGAL/algorithm/distance.h"
#include "SFCGAL/io/wkt.h"

#include "SFCGAL/detail/tools/Log.h"
#include "SFCGAL/detail/tools/Registry.h"

using namespace SFCGAL;

// always after CGAL
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_DistanceTest)

/*
 * check that distance between empty points is infinity
 */
BOOST_AUTO_TEST_CASE(testDistanceBetweenEmptyPointsIsInfinity)
{
  BOOST_CHECK_EQUAL(Point().distance(Point()),
                    std::numeric_limits<double>::infinity());
}

// TODO enable when implement is complete
#if 0

/*
 * check that distance between all kinds of empty geometry is infinity
 */
BOOST_AUTO_TEST_CASE( testDistanceBetweenEmptyGeometriesIsDefined )
{
    tools::Registry& registry = tools::Registry::instance() ;

    std::vector< std::string > geometryTypes = tools::Registry::instance().getGeometryTypes() ;

    for ( size_t i = 0; i < geometryTypes.size(); i++ ) {
        for ( size_t j = 0; j < geometryTypes.size(); j++ ) {
            BOOST_TEST_MESSAGE( boost::format( "distance(%s,%s)" ) % geometryTypes[i] % geometryTypes[j] );

            std::unique_ptr< Geometry > gA( registry.newGeometryByTypeName( geometryTypes[i] ) );
            std::unique_ptr< Geometry > gB( registry.newGeometryByTypeName( geometryTypes[j] ) );

            double dAB ;
            BOOST_CHECK_NO_THROW( dAB = gA->distance( *gB ) ) ;
            BOOST_CHECK_EQUAL( dAB, std::numeric_limits< double >::infinity() );
        }
    }
}
/*
 * check that distance3D between all kinds of empty geometry is infinity
 */
BOOST_AUTO_TEST_CASE( testDistance3DBetweenEmptyGeometriesIsDefined )
{
    tools::Registry& registry = tools::Registry::instance() ;

    std::vector< std::string > geometryTypes = tools::Registry::instance().getGeometryTypes() ;

    for ( size_t i = 0; i < geometryTypes.size(); i++ ) {
        for ( size_t j = 0; j < geometryTypes.size(); j++ ) {
            BOOST_TEST_MESSAGE( boost::format( "distance3D(%s,%s)" ) % geometryTypes[i] % geometryTypes[j] );

            std::unique_ptr< Geometry > gA( registry.newGeometryByTypeName( geometryTypes[i] ) );
            std::unique_ptr< Geometry > gB( registry.newGeometryByTypeName( geometryTypes[j] ) );

            double dAB ;
            BOOST_CHECK_NO_THROW( dAB = gA->distance3D( *gB ) ) ;
            BOOST_CHECK_EQUAL( dAB, std::numeric_limits< double >::infinity() );
        }
    }
}

#endif

BOOST_AUTO_TEST_CASE(testDistancePointPoint)
{
  BOOST_CHECK_EQUAL(Point(0.0, 0.0).distance(Point(0.0, 0.0)), 0.0);
  BOOST_CHECK_EQUAL(Point(1.0, 1.0).distance(Point(4.0, 5.0)), 5.0);
}
BOOST_AUTO_TEST_CASE(testDistancePointPoint3D)
{
  BOOST_CHECK_EQUAL(Point(0.0, 0.0, 0.0).distance3D(Point(0.0, 0.0, 0.0)), 0.0);
  BOOST_CHECK_EQUAL(Point(1.0, 1.0, 1.0).distance3D(Point(4.0, 1.0, 5.0)), 5.0);
}

// testPointLineString
BOOST_AUTO_TEST_CASE(testDistancePointLineString_pointOnLineString)
{
  Point const      point(1.0, 1.0);
  LineString const lineString(Point(0.0, 0.0), Point(2.0, 2.0));
  BOOST_CHECK_EQUAL(point.distance(lineString), 0.0);
}
BOOST_AUTO_TEST_CASE(
    testDistancePointLineString_pointOnLineString_badLineStringDefinition)
{
  Point const point(3.0, 4.0);
  LineString  lineString;
  lineString.addPoint(Point(0.0, 0.0));
  BOOST_CHECK_THROW(point.distance(lineString), GeometryInvalidityException);
}
BOOST_AUTO_TEST_CASE(
    testDistancePointLineString_pointOnLineString_collapsedSegments)
{
  Point const point(3.0, 4.0);
  LineString  lineString;
  lineString.addPoint(Point(0.0, 0.0));
  lineString.addPoint(Point(0.0, 0.0));
  BOOST_CHECK_THROW(point.distance(lineString), GeometryInvalidityException);
}
BOOST_AUTO_TEST_CASE(
    testDistancePointLineString3D_pointOnLineString_collapsedSegments)
{
  Point const point(0.0, 3.0, 4.0);
  LineString  lineString;
  lineString.addPoint(Point(0.0, 0.0, 0.0));
  lineString.addPoint(Point(0.0, -1.0, -1.0));
  BOOST_CHECK_EQUAL(point.distance3D(lineString), 5.0);
}

BOOST_AUTO_TEST_CASE(testDistancePointLineString_pointOutOfLineString)
{
  Point const      point(0.0, 1.0);
  LineString const lineString(Point(0.0, 0.0), Point(2.0, 2.0));
  BOOST_CHECK_EQUAL(point.distance(lineString), sqrt(2.0) / 2.0);
}

// testPointPolygon
BOOST_AUTO_TEST_CASE(testDistancePointPolygon_pointInPolygon)
{
  std::unique_ptr<Geometry>       gA(io::readWkt("POINT (0.5 0.5)"));
  std::unique_ptr<Geometry> const gB(
      io::readWkt("POLYGON ((0.0 0.0,1.0 0.0,1.0 1.0,0.0 1.0,0.0 0.0))"));
  BOOST_CHECK_EQUAL(gA->distance(*gB), 0.0);
}
BOOST_AUTO_TEST_CASE(testDistancePointPolygon_pointOutOfPolygon)
{
  std::unique_ptr<Geometry>       gA(io::readWkt("POINT (0.0 1.0)"));
  std::unique_ptr<Geometry> const gB(
      io::readWkt("POLYGON ((0.0 0.0,2.0 2.0,2.0 0.0,0.0 0.0))"));
  BOOST_CHECK_EQUAL(gA->distance(*gB), sqrt(2.0) / 2.0);
}

// LineString / LineString 2D
BOOST_AUTO_TEST_CASE(testDistanceLineStringLineString_zeroLengthSegments)
{
  std::unique_ptr<Geometry> gA(io::readWkt("LINESTRING (0.0 0.0,-1.0 -1.0)"));
  std::unique_ptr<Geometry> const gB(
      io::readWkt("LINESTRING (3.0 4.0,4.0 5.0)"));
  BOOST_CHECK_EQUAL(gA->distance(*gB), 5.0);
}
// LineString / LineString 3D
BOOST_AUTO_TEST_CASE(testDistanceLineStringLineString3D_zeroLengthSegments)
{
  std::unique_ptr<Geometry> gA(
      io::readWkt("LINESTRING (0.0 0.0 0.0,-1.0 -1.0 -1.0)"));
  std::unique_ptr<Geometry> const gB(
      io::readWkt("LINESTRING (0.0 3.0 4.0,0.0 4.0 5.0)"));
  BOOST_CHECK_EQUAL(gA->distance3D(*gB), 5.0);
}

// LineString / Triangle
BOOST_AUTO_TEST_CASE(testDistance3DLineStringTriangle_lineStringInTriangle)
{
  std::unique_ptr<Geometry> gA(
      io::readWkt("LINESTRING (-1.0 0.0 1.0,1.0 0.0 1.0)"));
  std::unique_ptr<Geometry> const gB(io::readWkt(
      "TRIANGLE ((-4.0 0.0 1.0,4.0 0.0 1.0,0.0 4.0 1.0,-4.0 0.0 1.0))"));
  BOOST_CHECK_EQUAL(gA->distance3D(*gB), 0.0);
}
BOOST_AUTO_TEST_CASE(
    testDistance3DLineStringTriangle_lineStringStartPointIsNearest)
{
  std::unique_ptr<Geometry> gA(
      io::readWkt("LINESTRING (-1.0 0.0 2.0,1.0 0.0 3.0)"));
  std::unique_ptr<Geometry> const gB(io::readWkt(
      "TRIANGLE ((-4.0 0.0 1.0,4.0 0.0 1.0,0.0 4.0 1.0,-4.0 0.0 1.0))"));
  BOOST_CHECK_EQUAL(gA->distance3D(*gB), 1.0);
}

// Triangle / Triangle
BOOST_AUTO_TEST_CASE(testDistance3DTriangleTriangle_contained)
{
  std::unique_ptr<Geometry>       gA(io::readWkt(
      "TRIANGLE ((-3.0 0.0 1.0,3.0 0.0 1.0,0.0 3.0 1.0,-3.0 0.0 1.0))"));
  std::unique_ptr<Geometry> const gB(io::readWkt(
      "TRIANGLE ((-4.0 0.0 1.0,4.0 0.0 1.0,0.0 4.0 1.0,-4.0 0.0 1.0))"));
  BOOST_CHECK_EQUAL(gA->distance3D(*gB), 0.0);
}
BOOST_AUTO_TEST_CASE(testDistance3DTriangleTriangle_parallel)
{
  std::unique_ptr<Geometry>       gA(io::readWkt(
      "TRIANGLE ((-3.0 0.0 1.0,3.0 0.0 1.0,0.0 3.0 1.0,-3.0 0.0 1.0))"));
  std::unique_ptr<Geometry> const gB(io::readWkt(
      "TRIANGLE ((-4.0 0.0 2.0,4.0 0.0 2.0,0.0 4.0 2.0,-4.0 0.0 2.0))"));
  BOOST_CHECK_EQUAL(gA->distance3D(*gB), 1.0);
}

// Polygon / Polygon

BOOST_AUTO_TEST_CASE(testDistancePolygonPolygon_disjoint)
{
  std::unique_ptr<Geometry> gA(
      io::readWkt("POLYGON ((0.0 0.0,1.0 0.0,1.0 1.0,0.0 1.0,0.0 0.0))"));
  std::unique_ptr<Geometry> const gB(
      io::readWkt("POLYGON ((2.0 0.0,3.0 0.0,3.0 1.0,2.0 1.0,2.0 0.0))"));
  BOOST_CHECK_EQUAL(gA->distance(*gB), 1.0);
}

BOOST_AUTO_TEST_CASE(testDistanceMultiPointMultiPoint_disjoint)
{
  std::unique_ptr<Geometry> gA(
      io::readWkt("MULTIPOINT ((0.0 0.0),(1.0 0.0),(1.0 1.0),(0.0 1.0))"));
  std::unique_ptr<Geometry> const gB(
      io::readWkt("MULTIPOINT ((8.0 8.0),(4.0 5.0))"));
  BOOST_CHECK_EQUAL(gA->distance(*gB), 5.0);
}

// Polygon / Solid
BOOST_AUTO_TEST_CASE(testDistancePolygonSolid)
{
  std::unique_ptr<Geometry> gA(
      io::readWkt("POLYGON ((1 -1 -1,1 1 -1,1 1 1,1 -1 1,1 -1 -1))"));
  std::unique_ptr<Geometry> const gB(
      io::readWkt("SOLID ((((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)),((0 0 0,0 0 1,0 1 "
                  "1,0 1 0,0 0 0)),((0 0 0,1 0 0,1 0 1,0 0 1,0 0 0)),((1 1 1,0 "
                  "1 1,0 0 1,1 0 1,1 1 1)),((1 1 1,1 0 1,1 0 0,1 1 0,1 1 "
                  "1)),((1 1 1,1 1 0,0 1 0,0 1 1,1 1 1))))"));
  BOOST_CHECK_EQUAL(gA->distance3D(*gB), 0);
}

BOOST_AUTO_TEST_CASE(testDistancePolygonSolid_disjoint)
{
  std::unique_ptr<Geometry> polygonA(
      io::readWkt("POLYGON Z ((2 2 2, 3 2 2, 3 3 2, 2 3 2, 2 2 2))"));
  std::unique_ptr<Geometry> const solidB(
      io::readWkt("SOLID ((((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)),((0 0 0,0 0 1,0 1 "
                  "1,0 1 0,0 0 0)),((0 0 0,1 0 0,1 0 1,0 0 1,0 0 0)),((1 1 1,0 "
                  "1 1,0 0 1,1 0 1,1 1 1)),((1 1 1,1 0 1,1 0 0,1 1 0,1 1 "
                  "1)),((1 1 1,1 1 0,0 1 0,0 1 1,1 1 1))))"));
  BOOST_CHECK_CLOSE(polygonA->distance3D(*solidB), 1.7320508, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
