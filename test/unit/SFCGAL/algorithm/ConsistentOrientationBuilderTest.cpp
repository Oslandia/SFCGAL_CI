// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Kernel.h"
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
#include "SFCGAL/algorithm/ConsistentOrientationBuilder.h"
#include "SFCGAL/algorithm/orientation.h"
#include "SFCGAL/io/wkt.h"

using namespace SFCGAL;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_ConsistentOrientationBuilderTest)

BOOST_AUTO_TEST_CASE(testOppositeTriangle)
{
  algorithm::ConsistentOrientationBuilder builder;
  builder.addTriangle(Triangle(Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
                               Point(0.0, 1.0, 0.0)));
  builder.addTriangle(Triangle(Point(0.0, 0.0, 0.0), Point(-1.0, 0.0, 0.0),
                               Point(0.0, 1.0, 0.0)));
  TriangulatedSurface const triangulatedSurface =
      builder.buildTriangulatedSurface();
  BOOST_CHECK_EQUAL(triangulatedSurface.numGeometries(), 1U);
  BOOST_CHECK_EQUAL(triangulatedSurface.numPatches(), 2U);
  BOOST_CHECK(algorithm::hasConsistentOrientation3D(triangulatedSurface));
}

BOOST_AUTO_TEST_CASE(testFourTriangle)
{
  algorithm::ConsistentOrientationBuilder builder;
  builder.addTriangle(Triangle(Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
                               Point(0.0, 1.0, 0.0)));
  builder.addTriangle(Triangle(Point(0.0, 0.0, 0.0), Point(-1.0, 0.0, 0.0),
                               Point(0.0, 1.0, 0.0)));
  builder.addTriangle(Triangle(Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
                               Point(0.0, -1.0, 0.0)));
  builder.addTriangle(Triangle(Point(0.0, 0.0, 0.0), Point(-1.0, 0.0, 0.0),
                               Point(0.0, -1.0, 0.0)));
  TriangulatedSurface const triangulatedSurface =
      builder.buildTriangulatedSurface();
  BOOST_CHECK_EQUAL(triangulatedSurface.numGeometries(), 1U);
  BOOST_CHECK_EQUAL(triangulatedSurface.numPatches(), 4U);
  BOOST_CHECK(algorithm::hasConsistentOrientation3D(triangulatedSurface));
}

BOOST_AUTO_TEST_SUITE_END()
