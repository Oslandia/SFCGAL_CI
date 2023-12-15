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
#include <SFCGAL/algorithm/plane.h>

using namespace SFCGAL ;
using namespace boost::unit_test ;

BOOST_AUTO_TEST_SUITE( SFCGAL_algorithm_PlaneTest )

BOOST_AUTO_TEST_CASE( testPlane1 )
{
    std::unique_ptr<Geometry> gA( io::readWkt( "POLYGON((0 0,1 0,1 1,0 1,0 0))" ) );

    CGAL::Plane_3<Kernel> const plane = algorithm::plane3D<Kernel>( gA->as<Polygon>() );
    BOOST_CHECK_EQUAL( plane.a(), 0.0 );
    BOOST_CHECK_EQUAL( plane.b(), 0.0 );
    BOOST_CHECK_EQUAL( plane.c(), 2.0 );
}

BOOST_AUTO_TEST_CASE( testPlane )
{
    struct TestCase {
        const std::string _wkt;
        const bool        _isPlane;
    };
    const TestCase test[] = {
        {"LINESTRING(1 2 3,4 5 6)", true}, // only two points
        {"LINESTRING(1 2 3,1 2 3,1 2 3,1 2 3)", true}, // all points in the same place
        {"LINESTRING(1 2 3,2 4 6,3 6 9,4 8 12)", true}, // all points aliged
        {"LINESTRING(1 2 3,6 5 4,7 8 9)", true}, // triangle must be plane
        {"LINESTRING(0 0 0,1 0 0,1 1 0,0 1 0,0 0 0)", true}, // all point in the plane z=0
        {"LINESTRING(2 1 0,2 0 0,2 1 0,2 1 0,2 0 3)", true}, // all points in the plane x=2
        {"LINESTRING(2 1 0,2 0 0,2 1 1,2 1 0,1 0 3)", false}, // one point out of plane
        {"LINESTRING(0 0 0, 1e-5 0 0, 1e-5 1e-5 0, 0 1e-5 1e-5)", false}, // fix #247
    };
    const size_t numTest = sizeof( test ) / sizeof( TestCase );

    for ( size_t t=0; t != numTest; ++t ) {
        //std::cout << "test = " << t << "\n";
        std::unique_ptr<Geometry> const g( io::readWkt( test[t]._wkt ) );
        const LineString* l = dynamic_cast<LineString*>( g.get() );
        BOOST_CHECK_MESSAGE(
            algorithm::isPlane3D< Kernel >( *l, 1.e-9 ) == test[t]._isPlane,
            ( boost::format( "LineString %d: %s %s" ) % t % test[t]._wkt % ( test[t]._isPlane ? "is plane" : "isn't plane" ) )
        );
    }
}

BOOST_AUTO_TEST_CASE( testPlane3DDivideByZeroCrash )
{
    std::unique_ptr< Geometry > degenerate_polygon = io::readWkt("POLYGON((1 -1 -1,1 0.5 0.5,1 0.5 0.5,1 -1 -1))");
    BOOST_CHECK( degenerate_polygon->geometryTypeId() == TYPE_POLYGON );

    // Should return degenerate plane without throwing
    auto degenerate_plane = algorithm::plane3D< Kernel >( degenerate_polygon->as<Polygon>() );

    // See triangulatePolygon3D for this pattern
    if (algorithm::hasPlane3D< Kernel >(degenerate_polygon->as<Polygon>()))
    {
        // Should not get here, OR plane3D with Plane3DInexactUnsafe should not divide by zero
        auto div_by_zero_check = algorithm::plane3D< Kernel >( degenerate_polygon->as<Polygon>(), algorithm::Plane3DInexactUnsafe() );
    }

    std::unique_ptr< Geometry > ok_polygon = io::readWkt("POLYGON((1 0.5 0.5,1.5 1.5 0.5,1.5 0.5 0.5,1 0.5 0.5))");
    BOOST_CHECK( ok_polygon->geometryTypeId() == TYPE_POLYGON );

    BOOST_CHECK( algorithm::hasPlane3D< Kernel >(ok_polygon->as<Polygon>()) );

    auto valid_plane = algorithm::plane3D< Kernel >( ok_polygon->as<Polygon>(), algorithm::Plane3DInexactUnsafe() );
}


BOOST_AUTO_TEST_SUITE_END()
