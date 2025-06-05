// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

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
#include "SFCGAL/algorithm/convexHull.h"
#include "SFCGAL/algorithm/intersection.h"
#include "SFCGAL/algorithm/intersects.h"
#include "SFCGAL/io/wkt.h"

#include "../test_config.h"
#include "Bench.h"

#include <boost/format.hpp>
#include <boost/test/unit_test.hpp>

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_IntersectionPerfTest)

#define N_POLYGONS 10000
#define N_POINTS 50

//
// Test limit case
BOOST_AUTO_TEST_CASE(testIntersectionPerf)
{

  //
  // generate polygons
  std::vector<Geometry *> polygons;

  for (size_t i = 0; i < N_POLYGONS; ++i) {
    MultiPoint mp;

    for (size_t j = 0; j < N_POINTS; ++j) {
      double x = (rand() + .0) / RAND_MAX * 10.0;
      double y = (rand() + .0) / RAND_MAX * 10.0;
      mp.addGeometry(Point(x, y));
    }

    std::unique_ptr<Geometry> g(algorithm::convexHull(mp));
    polygons.push_back(g.release());
  }

  bench().start("intersects convex hull");

  for (size_t i = 0; i < N_POLYGONS / 2; ++i) {
    algorithm::intersects(*polygons[2 * i], *polygons[2 * i + 1]);
  }

  bench().stop();
}

BOOST_AUTO_TEST_SUITE_END()
