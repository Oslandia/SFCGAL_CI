// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "SFCGAL/Kernel.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/detail/io/WktWriter.h"
#include "SFCGAL/detail/transform/ForceZOrderPoints.h"
#include "SFCGAL/io/wkt.h"

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_transform_ForceZOrderPointsTest)

BOOST_AUTO_TEST_CASE(simple)
{
  std::unique_ptr<Geometry> g1 = io::readWkt("POLYGON ((0 0,0 1,1 1,1 0,0 0))");

  const Polygon &p = g1->as<Polygon>();
  BOOST_CHECK(!p.isCounterClockWiseOriented());

  transform::ForceZOrderPoints forceZ;
  g1->accept(forceZ);

  BOOST_CHECK(g1->is3D());
  BOOST_CHECK(g1->as<Polygon>().isCounterClockWiseOriented());
}

BOOST_AUTO_TEST_SUITE_END()
