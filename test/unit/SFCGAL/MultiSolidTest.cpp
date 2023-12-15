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

#include <exception>

#include <SFCGAL/Envelope.h>
#include <SFCGAL/LineString.h>
#include <SFCGAL/Solid.h>
#include <SFCGAL/MultiSolid.h>

using namespace boost::unit_test ;
using namespace SFCGAL ;

BOOST_AUTO_TEST_SUITE( SFCGAL_MultiSolidTest )

BOOST_AUTO_TEST_CASE( defaultConstructor )
{
    MultiSolid const g;
    BOOST_CHECK( g.isEmpty() ) ;
    BOOST_CHECK( ! g.is3D() );
    BOOST_CHECK_EQUAL( g.numGeometries(), 0U );
}

BOOST_AUTO_TEST_CASE( testGeometryTypeId )
{
    MultiSolid const g;
    BOOST_CHECK_EQUAL( g.geometryTypeId(), TYPE_MULTISOLID );
}

//-- addAllowedGeometry
BOOST_AUTO_TEST_CASE( addSolid )
{
    MultiSolid g;
    g.addGeometry( new Solid() );
    BOOST_CHECK_EQUAL( g.numGeometries(), 1U );
}
//-- addForbidenGeometry
BOOST_AUTO_TEST_CASE( addLineStringThrow )
{
    MultiSolid g;
    BOOST_CHECK_THROW( g.addGeometry( LineString() ), std::exception );
}

//-- asText

BOOST_AUTO_TEST_CASE( asTextEmpty )
{
    MultiSolid const g;
    BOOST_CHECK_EQUAL( g.asText( 1 ), "MULTISOLID EMPTY" );
}
BOOST_AUTO_TEST_CASE( asText2d )
{
    MultiSolid g;
    g.addGeometry( Envelope( 0.0,1.0,0.0,1.0,0.0,1.0 ).toSolid().release() );
    g.addGeometry( Envelope( 2.0,3.0,4.0,5.0,6.0,7.0 ).toSolid().release() );
    BOOST_CHECK_EQUAL( g.asText( 0 ), "MULTISOLID Z(((((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)),((0 0 1,1 0 1,1 1 1,0 1 1,0 0 1)),((0 0 0,1 0 0,1 0 1,0 0 1,0 0 0)),((1 1 0,0 1 0,0 1 1,1 1 1,1 1 0)),((1 0 0,1 1 0,1 1 1,1 0 1,1 0 0)),((0 0 0,0 0 1,0 1 1,0 1 0,0 0 0)))),((((2 4 6,2 5 6,3 5 6,3 4 6,2 4 6)),((2 4 7,3 4 7,3 5 7,2 5 7,2 4 7)),((2 4 6,3 4 6,3 4 7,2 4 7,2 4 6)),((3 5 6,2 5 6,2 5 7,3 5 7,3 5 6)),((3 4 6,3 5 6,3 5 7,3 4 7,3 4 6)),((2 4 6,2 4 7,2 5 7,2 5 6,2 4 6)))))" );
}

//-- is< T >

BOOST_AUTO_TEST_CASE( isGeometryCollection )
{
    MultiSolid const g;
    BOOST_CHECK( g.is< GeometryCollection >() );
}

BOOST_AUTO_TEST_CASE( isMultiSolid )
{
    MultiSolid const g;
    BOOST_CHECK( g.is< MultiSolid >() );
}

BOOST_AUTO_TEST_SUITE_END()




