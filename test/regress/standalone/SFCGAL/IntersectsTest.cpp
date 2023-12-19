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
 *   License along with this library; if not, see
 <http://www.gnu.org/licenses/>.
 */
#include <fstream>

#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/LineString.h>
#include <SFCGAL/MultiLineString.h>
#include <SFCGAL/MultiPoint.h>
#include <SFCGAL/MultiPolygon.h>
#include <SFCGAL/MultiSolid.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Solid.h>
#include <SFCGAL/Triangle.h>
#include <SFCGAL/TriangulatedSurface.h>
#include <SFCGAL/algorithm/intersects.h>
#include <SFCGAL/io/wkt.h>

#include "../../../test_config.h"

#include <boost/format.hpp>
#include <boost/test/unit_test.hpp>

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_IntersectsTest)

//
// Load an obj
// Build a polyhedral surface
// Test intersection with some other objects

// FIXME removed wavefront obj loader, transform to wkt
// BOOST_AUTO_TEST_CASE( testIntersectsTeapot )
//{
//    std::string filename( SFCGAL_TEST_DIRECTORY );
//    filename += "/regress/data/teapot.obj" ;
//
//    io::WaveFrontObj obj;
//    obj.load( filename );
//
//    //check points defined in all faces
//    PolyhedralSurface surf;
//    for ( size_t i = 0; i < obj.numFaces(); i++ ){
//	LineString ext_ring;
//	const io::WaveFrontFace & face = obj.faceN(i);
//	for ( size_t j = 0; j < face.size(); j++ ) {
//	    ext_ring.addPoint( obj.vertexN( face[j].v ) );
//	}
//	// close the polygon
//	ext_ring.addPoint( obj.vertexN( face[0].v ) );
//	surf.addPolygon( Polygon(ext_ring) );
//    }
//
//    Triangle tri1( Point(3.0, 1.8, 0.0), Point(-2.99, 1.8, -0.081),
//    Point(-2.99, 1.8, 0.081) ); Triangle tri2( Point(3.0, 1.8, 5.0),
//    Point(-2.99, 1.8, 5.081), Point(-2.99, 1.8, 5.081) );
//
//    bool intersects1 = algorithm::intersects3D( tri1, surf );
//    bool intersects2 = algorithm::intersects3D( tri2, surf );
//    BOOST_TEST_MESSAGE( boost::format("triangle1 intersects the teapot: %1%")
//    % (intersects1 ? "YES" : "NO") ); BOOST_CHECK_EQUAL( intersects1, true );
//    BOOST_TEST_MESSAGE( boost::format("triangle2 intersects the teapot: %1%")
//    % (intersects2 ? "YES" : "NO") ); BOOST_CHECK_EQUAL( intersects2, false );
//}

//
// Test limit case
BOOST_AUTO_TEST_CASE(testLimitsIntersects)
{
  std::string filename(SFCGAL_TEST_DIRECTORY);
  filename += "/data/countries.wkt";

  std::ifstream ifs(filename.c_str());
  BOOST_REQUIRE(ifs.good());

  std::string wkt1;
  std::string wkt2;
  std::getline(ifs, wkt1);
  std::getline(ifs, wkt2);

  std::unique_ptr<Geometry> const g1(io::readWkt(wkt1));
  std::unique_ptr<Geometry> const g2(io::readWkt(wkt2));

  // check that a call to intersects() does not throw
  bool const throws = false;
  algorithm::intersects(*g1, *g2);
  BOOST_CHECK_EQUAL(throws, false);
}

BOOST_AUTO_TEST_SUITE_END()

// https://gitlab.com/Oslandia/SFCGAL/-/issues/253
BOOST_AUTO_TEST_CASE(issue_gitlab_253)
{
  std::unique_ptr<Geometry> const g1(
      io::readWkt("POLYHEDRALSURFACE Z (((2 0.5 0.5, 0 0 1, 0 0 0, 2 0.5 "
                  "0.5)), ((2 0.5 0.5, 0 1 1, 0 0 1, 2 0.5 0.5)), ((2 0.5 0.5, "
                  "0 1 0, 0 1 1, 2 0.5 0.5)), ((2 0.5 0.5, 0 0 0, 0 1 0, 2 0.5 "
                  "0.5)), ((0 0 0, 0 0 1, 0 1 1, 0 1 0, 0 0 0)))"));
  std::unique_ptr<Geometry> const g2(io::readWkt("POINT (0 0.5 0.5)"));

  // check that a call to intersects() does not throw
  bool const throws     = false;
  bool       intersects = algorithm::intersects3D(*g1, *g2);
  BOOST_CHECK_EQUAL(intersects, true);
  intersects = algorithm::intersects3D(*g2, *g1);
  BOOST_CHECK_EQUAL(intersects, true);
  BOOST_CHECK_EQUAL(throws, false);
}
