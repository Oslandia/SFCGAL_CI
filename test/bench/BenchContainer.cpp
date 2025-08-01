// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

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

#include "../test_config.h"

#include "Bench.h"

#include <boost/test/unit_test.hpp>

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_BenchContainer)

#if 0

//
// Test limit case
BOOST_AUTO_TEST_CASE( testVectorOfPoint2 )
{
    const int N = 100000 ;

    // create points
    bench().start( "std::vector< Kernel::Point_2 > create" ) ;
    std::vector< Kernel::Point_2 > points ;
    points.reserve( N ) ;

    for ( int i = 0; i < N; i++ ) {
        points.push_back( Kernel::Point_2( i,i )  ) ;
    }

    bench().stop();

    // centroid
    bench().start( "std::vector< Kernel::Point_2 > access" ) ;
    Kernel::FT x = 0 , y = 0 ;

    for ( std::vector< Kernel::Point_2 >::const_iterator it = points.begin(); it != points.end(); ++it ) {
        x += it->x() ;
        y += it->y() ;
    }

    Kernel::Point_2 centroid( x / Kernel::FT( ( int )points.size() ), y / Kernel::FT( ( int )points.size() ) );
    bench().stop() ;

    // reverse
    bench().start( "std::vector< Kernel::Point_2 > reverse" ) ;
    std::reverse( points.begin(), points.end() ) ;
    bench().stop();

    // clear
    bench().start( "std::vector< Kernel::Point_2 > clear" ) ;
    points.clear() ;
    bench().stop() ;
}

//
// Test limit case
BOOST_AUTO_TEST_CASE( testVectorOfPoint )
{
    const int N = 100000 ;

    // create points
    bench().start( "std::vector< Point > create" );
    std::vector< Point > points ;
    points.reserve( N ) ;

    for ( int i = 0; i < N; i++ ) {
        points.push_back( Point( i,i )  ) ;
    }

    bench().stop();

    // centroid
    bench().start( "std::vector< Point > access" );
    Kernel::FT x = 0 , y = 0 ;

    for ( std::vector< Point >::const_iterator it = points.begin(); it != points.end(); ++it ) {
        x += it->x() ;
        y += it->y() ;
    }

    Point centroid( x / Kernel::FT( ( int )points.size() ), y / Kernel::FT( ( int )points.size() ) );
    bench().stop();


    // reverse
    bench().start( "std::vector< Point > reverse" ) ;
    std::reverse( points.begin(), points.end() ) ;
    bench().stop() ;

    // clear
    bench().start( "std::vector< Point > clear" ) ;
    points.clear() ;
    bench().stop() ;
}


//
// Test limit case
BOOST_AUTO_TEST_CASE( testPtrVectorOfPoint )
{
    const int N = 100000 ;

    boost::timer::cpu_timer timer;

    // create points
    bench().start( "boost::ptr_vector< Point > create" ) ;
    boost::ptr_vector< Point > points ;
    points.reserve( N ) ;

    for ( int i = 0; i < N; i++ ) {
        points.push_back( new Point( i,i )  ) ;
    }

    bench().stop() ;

    // centroid
    bench().start( "boost::ptr_vector< Point > access" ) ;
    Kernel::FT x = 0 , y = 0 ;

    for ( boost::ptr_vector< Point >::const_iterator it = points.begin(); it != points.end(); ++it ) {
        x += it->x() ;
        y += it->y() ;
    }

    Point centroid( x / Kernel::FT( ( int )points.size() ), y / Kernel::FT( ( int )points.size() ) );
    bench().stop();


    // reverse
    bench().start( "boost::ptr_vector< Point > reverse" );
    std::reverse( points.begin(), points.end() ) ;
    bench().stop();

    // clear
    bench().start( "boost::ptr_vector< Point > clear" );
    points.clear() ;
    bench().stop();
}

#endif

BOOST_AUTO_TEST_SUITE_END()
