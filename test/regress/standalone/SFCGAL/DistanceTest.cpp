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

#include "../../../test_config.h"

#include <fstream>
#include <boost/format.hpp>

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
//#include <SFCGAL/io/vtk.h>

#include <SFCGAL/detail/tools/Log.h>
#include <cmath>

using namespace boost::unit_test ;
using namespace SFCGAL ;


BOOST_AUTO_TEST_SUITE( SFCGAL_DistanceTest )

/**
 * Compute the area for given geometries
 */
BOOST_AUTO_TEST_CASE( testFileDistanceTest )
{
    logger().setLogLevel( Logger::Debug );

    std::string filename( SFCGAL_TEST_DIRECTORY );
    filename += "/data/DistanceTest.txt" ;

    std::ifstream ifs( filename.c_str() );
    BOOST_REQUIRE( ifs.good() ) ;

    int const argc = framework::master_test_suite().argc;
    char** argv = framework::master_test_suite().argv;

    // look for options
    int test_one_line = -1;

    for ( int i = 0; i < argc; ++i ) {
        std::string const argi( argv[i] );

        if ( argi == "--line" ) {
            // only test one line
            if ( argc > i+1 ) {
                sscanf( argv[i+1], "%d", &test_one_line );
                ++i;
                continue;
            }
        }
    }

    std::string line;
    int lineNo = 0;

    while ( std::getline( ifs, line ) ) {
        ++lineNo;

        if ( line[0] == '#' || line.empty() ) {
            continue ;
        }

        if ( -1 != test_one_line && lineNo != test_one_line ) {
            continue;
        }

        BOOST_TEST_MESSAGE( ( boost::format( "%s:%d" ) % filename % lineNo ).str() );

        std::istringstream iss( line );

        std::string distanceDimension ;
        std::string wktGA;
        std::string wktGB ;
        double expectedDistance = NAN ;

        std::getline( iss, distanceDimension, '|' ) ;
        std::getline( iss, wktGA, '|' ) ;
        std::getline( iss, wktGB, '|' ) ;
        iss >> expectedDistance ;

        std::unique_ptr< Geometry > gA( io::readWkt( wktGA ) );
        std::unique_ptr< Geometry > const gB( io::readWkt( wktGB ) );

        //if (43!=lineNo ) continue;
        //if (43==lineNo )
        //{
        //    io::vtk(gA->as<Polygon>(), "/tmp/gA.vtk");
        //    io::vtk(gB->as<MultiPolygon>(), "/tmp/gB.vtk");
        //}


        if ( distanceDimension == "2" ) {
            BOOST_CHECK_CLOSE( gA->distance( *gB ), expectedDistance, 1e-13 );
        }
        else if ( distanceDimension == "3" ) {
            BOOST_CHECK_CLOSE( gA->distance3D( *gB ), expectedDistance, 1e-13 );
        }
        else {
            BOOST_CHECK( false );
        }
    }
}


BOOST_AUTO_TEST_SUITE_END()




