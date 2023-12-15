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
#include <SFCGAL/algorithm/tesselate.h>

#include <SFCGAL/detail/tools/Registry.h>

using namespace SFCGAL ;

// always after CGAL
using namespace boost::unit_test ;

BOOST_AUTO_TEST_SUITE( SFCGAL_algorithm_TesselateTest )


BOOST_AUTO_TEST_CASE( testEmpty )
{
    tools::Registry const& registry = tools::Registry::instance() ;

    std::vector< std::string > const geometryTypes = registry.getGeometryTypes() ;

    for (auto & geometryType : geometryTypes) {
        std::unique_ptr< Geometry > g( registry.newGeometryByTypeName( geometryType ) ) ;
        BOOST_TEST_MESSAGE( boost::format( "tesselate(%s)" ) % g->asText() );
        std::unique_ptr< Geometry > result = algorithm::tesselate( *g ) ;
        BOOST_CHECK( result->isEmpty() );
    }
}

/*
 * test invariants (Point,LineString & co)
 */

BOOST_AUTO_TEST_CASE( testPoint )
{
    std::string const wkt = "POINT(3.0 4.0)" ;
    std::unique_ptr< Geometry > const g( io::readWkt( wkt ) );
    BOOST_CHECK_EQUAL( algorithm::tesselate( *g )->asText( 1 ), wkt );
}
BOOST_AUTO_TEST_CASE( testLineString )
{
    std::string const wkt = "LINESTRING(0.0 0.0,1.0 1.0)" ;
    std::unique_ptr< Geometry > const g( io::readWkt( wkt ) );
    BOOST_CHECK_EQUAL( algorithm::tesselate( *g )->asText( 1 ), wkt );
}
BOOST_AUTO_TEST_CASE( testMultiPoint )
{
    std::string const wkt = "MULTIPOINT((3.0 4.0),(5.0 6.0))" ;
    std::unique_ptr< Geometry > const g( io::readWkt( wkt ) );
    BOOST_CHECK_EQUAL( algorithm::tesselate( *g )->asText( 1 ), wkt );
}
BOOST_AUTO_TEST_CASE( testMultiLineString )
{
    std::string const wkt = "MULTILINESTRING((0.0 0.0,1.0 1.0),(1.0 1.0,2.0 2.0))" ;
    std::unique_ptr< Geometry > const g( io::readWkt( wkt ) );
    BOOST_CHECK_EQUAL( algorithm::tesselate( *g )->asText( 1 ), wkt );
}


/*
 * test with polygon, MultiPolygon & co
 */
BOOST_AUTO_TEST_CASE( testPolygon )
{
    std::string const wkt = "POLYGON((0.0 0.0,1.0 0.0,1.0 1.0,0.0 1.0,0.0 0.0))" ;
    std::string const wktOut = "TIN(((0.0 1.0,1.0 0.0,1.0 1.0,0.0 1.0)),((0.0 1.0,0.0 0.0,1.0 0.0,0.0 1.0)))" ;
    std::unique_ptr< Geometry > const g( io::readWkt( wkt ) );
    std::unique_ptr< Geometry > result( algorithm::tesselate( *g ) );
    BOOST_CHECK_EQUAL( result->asText( 1 ), wktOut );
}
BOOST_AUTO_TEST_CASE( testMultiPolygon )
{
    std::string const wkt = "MULTIPOLYGON(((0.0 0.0,1.0 0.0,1.0 1.0,0.0 1.0,0.0 0.0)),((2.0 0.0,3.0 0.0,3.0 1.0,2.0 1.0,2.0 0.0)))" ;
    std::string const wktOut = "GEOMETRYCOLLECTION(TIN(((0.0 1.0,1.0 0.0,1.0 1.0,0.0 1.0)),((0.0 1.0,0.0 0.0,1.0 0.0,0.0 1.0))),TIN(((2.0 1.0,3.0 0.0,3.0 1.0,2.0 1.0)),((2.0 1.0,2.0 0.0,3.0 0.0,2.0 1.0))))" ;
    std::unique_ptr< Geometry > const g( io::readWkt( wkt ) );
    std::unique_ptr< Geometry > result( algorithm::tesselate( *g ) );
    BOOST_CHECK_EQUAL( result->asText( 1 ), wktOut );
}


BOOST_AUTO_TEST_SUITE_END()

