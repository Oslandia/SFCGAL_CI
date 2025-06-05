// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>
#include <fstream>

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
// #include "SFCGAL/io/WaveFrontObj.h"

#include "../../../test_config.h"

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_io_WaveFrontObjTest)

// FIXME removed wavefront obj loader, transform to wkt
// BOOST_AUTO_TEST_CASE( testTeaPot )
//{
//	std::string filename( SFCGAL_TEST_DIRECTORY );
//	filename += "/regress/data/teapot.obj" ;
//
//	io::WaveFrontObj obj;
//	obj.load( filename );
//
//	BOOST_CHECK_EQUAL( obj.numVertices(), 3644U );
//	BOOST_CHECK_EQUAL( obj.numFaces(),    6320U );
//
//	//check points defined in all faces
//	for ( size_t i = 0; i < obj.numFaces(); i++ ){
//		const io::WaveFrontFace & face = obj.faceN(i);
//		for ( size_t j = 0; j < face.size(); j++ ){
//			BOOST_CHECK_GE( face[j].v, 0 );
//			BOOST_CHECK_LT( face[j].v, (int)obj.numVertices() );
//		}
//	}
//}

BOOST_AUTO_TEST_SUITE_END()
