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
#include <fstream>

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
#include <SFCGAL/algorithm/straightSkeleton.h>

#include "../../../test_config.h"

#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>

#include <CGAL/version.h>
#include <CGAL/version_macros.h>

using namespace boost::unit_test ;

using namespace SFCGAL ;

BOOST_AUTO_TEST_SUITE( SFCGAL_StraightSkeletonTest )

namespace {

    void runTest(const boost::filesystem::path::string_type& filename)
    {
        std::ifstream ifs( filename.c_str() );
        BOOST_REQUIRE( ifs.good() ) ;

        std::string line;
        std::string lbl_base = boost::filesystem::basename(filename);

        int lineno = 0;
        while ( std::getline( ifs, line ) ) {
            ++lineno;
            if ( line[0] == '#' || line.empty() ) {
                continue ;
            }
            std::stringstream lblstream;
            lblstream << lbl_base << ':' << lineno << ": ";
            std::string lbl = lblstream.str();

            std::istringstream iss( line );
            std::string inputWkt ;
            std::string outputWkt ;
            std::string obtWkt ;

            if ( ! std::getline( iss, inputWkt, '|' ) ||
                 ! std::getline( iss, outputWkt, '|' ) )
            {
                std::stringstream ss;
                ss << lbl << "missing `|' char in test file";
                BOOST_CHECK_EQUAL("", ss.str());
                continue;
            }

            std::unique_ptr< Geometry > g;
            try {
              g = io::readWkt( inputWkt );
            } catch (const std::exception &e) {
              std::stringstream ss;
              ss << lbl << e.what();
              BOOST_CHECK_EQUAL("", ss.str());
              continue;
            }
            std::unique_ptr< MultiLineString > result;
            try {
              result = algorithm::straightSkeleton( *g ) ;
              obtWkt = result->asText( 6 );
            } catch (const std::exception& e) {
              obtWkt = std::string(e.what());
            }
            std::string obt = lbl + obtWkt;
            std::string exp = lbl + outputWkt;
            if ( (CGAL_VERSION_MAJOR == 5) && (CGAL_VERSION_MINOR >= 2)) {
              BOOST_CHECK_EQUAL( exp, obt );
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( testStraightSkeletonTest )
{
  std::unique_ptr<Geometry>         g;
  std::unique_ptr<MultiLineString>  result;

  std::string wkt   = "POLYGON((256 760,518 760,518 630,674 630,674  239,673 239,127 239,127 240,126 240,126 513,127 513,127 514,126 514,126  630,255 630,256 630,256 760),(128 629,128 423,270 423,270 422,271  422,271 240,672 240,672 629,128 629),(258 759,258 631,516 631,516  759,258 759),(128 421,128 240,269 240,269 421,128 421))";
  std::string skWkt = "MULTILINESTRING((256.0 760.0,257.0 759.0),(256.0 630.0,256.5 629.5),(255.0 630.0,255.0 629.5),(126.0 630.0,127.0 629.0),(126.0 514.0,127.0 515.0),(127.0 514.0,127.5 514.5),(127.0 513.0,127.5 512.5),(126.0 513.0,127.0 512.0),(126.0 240.0,127.0 241.0),(127.0 240.0,127.5 240.5),(127.0 239.0,127.5 239.5),(673.0 239.0,673.0 240.0),(674.0 239.0,673.0 240.0),(674.0 630.0,673.0 629.0),(518.0 630.0,517.5 629.5),(518.0 760.0,517.0 759.0),(128.0 629.0,127.5 629.5),(672.0 629.0,672.5 629.5),(672.0 240.0,672.5 239.5),(271.0 240.0,270.5 239.5),(271.0 422.0,270.0 421.0),(270.0 422.0,269.5 421.5),(270.0 423.0,269.0 422.0),(128.0 423.0,127.0 422.0),(258.0 759.0,257.5 759.5),(516.0 759.0,516.5 759.5),(516.0 631.0,517.0 630.0),(258.0 631.0,257.0 630.0),(128.0 421.0,127.0 422.0),(269.0 421.0,269.5 421.5),(269.0 240.0,269.5 239.5),(128.0 240.0,127.5 239.5),(256.5 629.5,257.0 630.0),(256.5 629.5,255.0 629.5),(127.5 514.5,127.0 515.0),(127.5 514.5,127.5 512.5),(127.5 512.5,127.0 512.0),(269.5 239.5,270.0 240.0),(269.5 239.5,127.5 239.5),(127.5 239.5,127.5 240.5),(255.0 629.5,127.5 629.5),(127.5 240.5,127.0 241.0),(127.5 629.5,127.0 629.0),(517.5 629.5,672.5 629.5),(517.5 629.5,517.0 630.0),(672.5 629.5,673.0 629.0),(672.5 239.5,673.0 240.0),(672.5 239.5,270.5 239.5),(270.5 239.5,270.0 240.0),(269.5 421.5,269.0 422.0),(269.5 421.5,270.0 421.0),(257.5 759.5,257.0 759.0),(257.5 759.5,516.5 759.5),(516.5 759.5,517.0 759.0),(127.0 512.0,127.0 422.0),(127.0 241.0,127.0 422.0),(269.0 422.0,127.0 422.0),(270.0 421.0,270.0 240.0),(673.0 240.0,673.0 629.0),(127.0 515.0,127.0 629.0),(257.0 759.0,257.0 630.0),(257.0 630.0,517.0 630.0),(517.0 759.0,517.0 630.0))";
  std::string resultWkt;

  g = io::readWkt(wkt);
  result = algorithm::straightSkeleton(*g);
  resultWkt = result->asText(1);
  std::cout << resultWkt << "\n";
  BOOST_CHECK_EQUAL( resultWkt, skWkt );
}

BOOST_AUTO_TEST_CASE( testStraightSkeletonTestIssue )
{
    using namespace boost;
    using namespace boost::filesystem;

    std::string testdir( SFCGAL_TEST_DIRECTORY );
    path dirname = testdir + "/data/StraightSkeletonTest" ;
    if ( is_directory(dirname) )
    {
      directory_iterator it = directory_iterator(dirname);
      while ( it != directory_iterator() )
      {
        path f = *it;
        runTest(f.c_str());
        ++it;
      }
    }
}


BOOST_AUTO_TEST_SUITE_END()
