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

#include <SFCGAL/Kernel.h>
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
#include <SFCGAL/algorithm/force3D.h>

#include <SFCGAL/detail/tools/Registry.h>

using namespace boost::unit_test ;
using namespace SFCGAL ;

// note that it relies on a visitor

BOOST_AUTO_TEST_SUITE( SFCGAL_algorithm_Force3DTest )

BOOST_AUTO_TEST_CASE( testIgnoreEmpty )
{
    tools::Registry& registry = tools::Registry::instance() ;
    std::vector< std::string > typeNames = tools::Registry::instance().getGeometryTypes();

    for ( size_t i = 0; i < typeNames.size(); i++ ) {
        BOOST_TEST_MESSAGE( typeNames[i] ) ;

        std::unique_ptr< Geometry > g( registry.newGeometryByTypeName( typeNames[i] ) ) ;
        BOOST_REQUIRE( g.get() != NULL ) ;
        algorithm::force3D( *g ) ;
        BOOST_CHECK( g->isEmpty() ) ;
    }
}

BOOST_AUTO_TEST_CASE( testPointForceZ )
{
    Point p( 3.0,4.0 );
    algorithm::force3D( p );
    BOOST_CHECK_EQUAL( p.asText( 1 ), "POINT Z(3.0 4.0 0.0)" );
}
BOOST_AUTO_TEST_CASE( testPointForceZWithValue )
{
    Point p( 3.0,4.0 );
    algorithm::force3D( p,-9999.0 );
    BOOST_CHECK_EQUAL( p.asText( 1 ), "POINT Z(3.0 4.0 -9999.0)" );
}


BOOST_AUTO_TEST_CASE( test_MixedLineString2D3D )
{
    LineString lineString ;
    lineString.addPoint( Point( 1.0,1.0 ) );
    lineString.addPoint( Point( 2.0,2.0,1.0 ) );
    lineString.addPoint( Point( 3.0,3.0 ) );
    algorithm::force3D( lineString );
    // should keep 1.0 for the second point
    BOOST_CHECK_EQUAL( lineString.asText( 1 ), "LINESTRING Z(1.0 1.0 0.0,2.0 2.0 1.0,3.0 3.0 0.0)" );
}




BOOST_AUTO_TEST_SUITE_END()

