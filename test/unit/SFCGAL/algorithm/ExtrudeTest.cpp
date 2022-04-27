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
 *   License along with this library; if not, see <http://www.gnu.org/licenses/>.
 */
#include <boost/test/unit_test.hpp>

#include <SFCGAL/Point.h>
#include <SFCGAL/LineString.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/Triangle.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/TriangulatedSurface.h>
#include <SFCGAL/Solid.h>
#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/MultiPoint.h>
#include <SFCGAL/MultiLineString.h>
#include <SFCGAL/MultiPolygon.h>
#include <SFCGAL/MultiSolid.h>
#include <SFCGAL/io/wkt.h>
#include <SFCGAL/algorithm/extrude.h>
#include <SFCGAL/detail/transform/ForceZ.h>
#include <SFCGAL/io/wkt.h>

using namespace SFCGAL ;

// always after CGAL
using namespace boost::unit_test ;

BOOST_AUTO_TEST_SUITE( SFCGAL_algorithm_ExtrudeTest )


BOOST_AUTO_TEST_CASE( testExtrudePoint )
{
    Point g( 0.0,0.0,0.0 );
    std::unique_ptr< Geometry > ext( algorithm::extrude( g, 0.0, 0.0, 1.0 ) );
    BOOST_CHECK( ext->is< LineString >() );
    BOOST_CHECK( ext->as< LineString >().is3D() );
    BOOST_CHECK_EQUAL( ext->asText( 1 ), "LINESTRING Z(0.0 0.0 0.0,0.0 0.0 1.0)" );
}


BOOST_AUTO_TEST_CASE( testExtrudeLineString )
{
    LineString g(
        Point( 0.0,0.0,0.0 ),
        Point( 1.0,0.0,0.0 )
    );
    std::unique_ptr< Geometry > ext( algorithm::extrude( g, 0.0, 0.0, 1.0 ) );
    BOOST_CHECK( ext->is< PolyhedralSurface >() );
    BOOST_CHECK( ext->as< PolyhedralSurface >().is3D() );
    BOOST_CHECK_EQUAL( ext->asText( 1 ), "POLYHEDRALSURFACE Z(((0.0 0.0 0.0,1.0 0.0 0.0,1.0 0.0 1.0,0.0 0.0 1.0,0.0 0.0 0.0)))" );
}



BOOST_AUTO_TEST_CASE( testExtrudeSquare )
{
    std::vector< Point > points;
    points.push_back( Point( 0.0,0.0,0.0 ) );
    points.push_back( Point( 1.0,0.0,0.0 ) );
    points.push_back( Point( 1.0,1.0,0.0 ) );
    points.push_back( Point( 0.0,1.0,0.0 ) );
    points.push_back( Point( 0.0,0.0,0.0 ) );

    LineString exteriorRing( points ) ;
    Polygon g( exteriorRing );
    std::unique_ptr< Geometry > ext( algorithm::extrude( g, 0.0, 0.0, 1.0 ) );
    BOOST_CHECK( ext->is< Solid >() );
    BOOST_CHECK_EQUAL( ext->as< Solid >().numShells(), 1U );
    BOOST_CHECK_EQUAL( ext->as< Solid >().exteriorShell().numPolygons(), 6U );
}

BOOST_AUTO_TEST_CASE( testExtrudePolyhedral )
{
    std::unique_ptr<Geometry> g = io::readWkt( "POLYHEDRALSURFACE(((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)))" );

    std::unique_ptr< Geometry > ext = algorithm::extrude( *g, 0.0, 0.0, 1.0 );
    BOOST_CHECK( ext->is< Solid >() );
    BOOST_CHECK_EQUAL( ext->as< Solid >().numShells(), 1U );
}

BOOST_AUTO_TEST_CASE( testExtrudeMultiPolygon )
{
    std::vector< Point > points;
    points.push_back( Point( 0.0,0.0,0.0 ) );
    points.push_back( Point( 1.0,0.0,0.0 ) );
    points.push_back( Point( 1.0,1.0,0.0 ) );
    points.push_back( Point( 0.0,1.0,0.0 ) );
    points.push_back( Point( 0.0,0.0,0.0 ) );

    std::vector< Point > points2;
    points2.push_back( Point( 2.0,0.0,0.0 ) );
    points2.push_back( Point( 3.0,0.0,0.0 ) );
    points2.push_back( Point( 3.0,1.0,0.0 ) );
    points2.push_back( Point( 2.0,1.0,0.0 ) );
    points2.push_back( Point( 2.0,0.0,0.0 ) );

    LineString exteriorRing( points ) ;
    LineString exteriorRing2( points2 ) ;
    Polygon g1( exteriorRing );
    Polygon g2( exteriorRing2 );
    MultiPolygon mp;
    mp.addGeometry( g1 );
    mp.addGeometry( g2 );

    std::unique_ptr< Geometry > ext( algorithm::extrude( mp, 0.0, 0.0, 1.0 ) );
    BOOST_CHECK( ext->is< MultiSolid >() );
    BOOST_CHECK_EQUAL( ext->as<MultiSolid>().numGeometries(), 2U );
}


BOOST_AUTO_TEST_CASE( testExtrudeSquareWithHole )
{
    std::vector< LineString > rings;
    {
        std::vector< Point > points;
        points.push_back( Point( 0.0,0.0,0.0 ) );
        points.push_back( Point( 1.0,0.0,0.0 ) );
        points.push_back( Point( 1.0,1.0,0.0 ) );
        points.push_back( Point( 0.0,1.0,0.0 ) );
        points.push_back( Point( 0.0,0.0,0.0 ) );
        rings.push_back( LineString( points ) );
    }
    {
        std::vector< Point > points;
        points.push_back( Point( 0.2,0.2,0.0 ) );
        points.push_back( Point( 0.8,0.2,0.0 ) );
        points.push_back( Point( 0.8,0.8,0.0 ) );
        points.push_back( Point( 0.2,0.8,0.0 ) );
        points.push_back( Point( 0.2,0.2,0.0 ) );

        std::reverse( points.begin(), points.end() );

        rings.push_back( LineString( points ) );
    }

    Polygon g( rings );
    std::unique_ptr< Geometry > ext( algorithm::extrude( g, 0.0, 0.0, 1.0 ) );
    BOOST_CHECK( ext->is< Solid >() );
    BOOST_CHECK_EQUAL( ext->as< Solid >().numShells(), 1U );
    BOOST_CHECK_EQUAL( ext->as< Solid >().exteriorShell().numPolygons(), 10U );
}


//SELECT ST_AsText(ST_Extrude(ST_Extrude(ST_Extrude('POINT(0 0)', 1, 0, 0), 0, 1, 0), 0, 0, 1));
BOOST_AUTO_TEST_CASE( testChainingExtrude )
{
    std::unique_ptr< Geometry > g( new Point( 0.0,0.0 ) );
    g = algorithm::extrude( *g, 1.0, 0.0, 0.0 ) ;
    BOOST_CHECK_EQUAL( g->asText( 0 ), "LINESTRING Z(0 0 0,1 0 0)" ) ;
    g =  algorithm::extrude( *g, 0.0, 1.0, 0.0 ) ;
    BOOST_CHECK_EQUAL( g->asText( 0 ), "POLYHEDRALSURFACE Z(((0 0 0,1 0 0,1 1 0,0 1 0,0 0 0)))" ) ;
    g =  algorithm::extrude( *g, 0.0, 0.0, 1.0 ) ;
    BOOST_CHECK_EQUAL( g->asText( 0 ), "SOLID Z((((0 1 0,1 1 0,1 0 0,0 1 0)),((0 1 1,1 0 1,1 1 1,0 1 1)),((0 1 0,1 0 0,0 0 0,0 1 0)),((0 1 1,0 0 1,1 0 1,0 1 1)),((1 0 0,1 1 0,1 1 1,1 0 1,1 0 0)),((1 1 0,0 1 0,0 1 1,1 1 1,1 1 0)),((0 1 0,0 0 0,0 0 1,0 1 1,0 1 0)),((0 0 0,1 0 0,1 0 1,0 0 1,0 0 0))))" ) ;
}


BOOST_AUTO_TEST_SUITE_END()

