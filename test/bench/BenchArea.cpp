// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
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

#include "../test_config.h"
#include "Bench.h"

#include <boost/test/unit_test.hpp>

#include "SFCGAL/detail/generator/sierpinski.h"

#include "SFCGAL/algorithm/area.h"

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_BenchArea)

BOOST_AUTO_TEST_CASE(testAreaSierpinski)
{
  std::unique_ptr<MultiPolygon> fractal(generator::sierpinski(9));

  bench().start("area sierpinski");

  for (int i = 0; i < 10; i++) {
    algorithm::area(*fractal);
  }

  bench().stop();
}

BOOST_AUTO_TEST_CASE(testAreaSierpinski3D)
{
  std::unique_ptr<MultiPolygon> fractal(generator::sierpinski(9));

  bench().start("area sierpinski");

  for (int i = 0; i < 10; i++) {
    algorithm::area3D(*fractal);
  }

  bench().stop();
}

BOOST_AUTO_TEST_SUITE_END()
