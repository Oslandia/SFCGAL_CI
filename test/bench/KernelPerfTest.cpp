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
#include "SFCGAL/io/wkt.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "../test_config.h"
#include "Bench.h"

#include <boost/format.hpp>
#include <boost/test/unit_test.hpp>

using namespace boost::unit_test;
using namespace SFCGAL;

typedef CGAL::Exact_predicates_exact_constructions_kernel   ExactKernel;
typedef CGAL::Exact_predicates_inexact_constructions_kernel InexactKernel;

BOOST_AUTO_TEST_SUITE(SFCGAL_KernelPerfTest)

#define N_POINTS 10000000
//
// Test limit case
BOOST_AUTO_TEST_CASE(testPointCopyPerf)
{
  bench().start("Inexact kernel copy");
  InexactKernel::Point_2 ip2;

  for (size_t i = 0; i < N_POINTS; ++i) {
    // test copy
    InexactKernel::Point_2 p((double)rand() / RAND_MAX,
                             (double)rand() / RAND_MAX);
    ip2 = p;
  }

  bench().stop();

  bench().start("Exact kernel copy");
  ExactKernel::Point_2 ep2;

  for (size_t i = 0; i < N_POINTS; ++i) {
    // test copy
    ExactKernel::Point_2 p((double)rand() / RAND_MAX,
                           (double)rand() / RAND_MAX);
    ep2 = p;
  }

  bench().stop();
}

BOOST_AUTO_TEST_SUITE_END()
