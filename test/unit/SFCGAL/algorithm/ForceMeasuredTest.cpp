// Copyright (c) 2025-2025, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "SFCGAL/LineString.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/algorithm/forceMeasured.h"

#include "SFCGAL/detail/tools/Registry.h"

using namespace boost::unit_test;
using namespace SFCGAL;

// note that it relies on a visitor

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_ForceMTest)

BOOST_AUTO_TEST_CASE(testIgnoreEmpty)
{
  tools::Registry const         &registry = tools::Registry::instance();
  std::vector<std::string> const typeNames =
      tools::Registry::instance().getGeometryTypes();

  for (const auto &typeName : typeNames) {
    BOOST_TEST_MESSAGE(typeName);

    std::unique_ptr<Geometry> geometry(
        registry.newGeometryByTypeName(typeName));
    BOOST_REQUIRE(geometry.get() != nullptr);
    algorithm::forceMeasured(*geometry);
    BOOST_CHECK(geometry->isEmpty());
  }
}

BOOST_AUTO_TEST_CASE(testPointForceM)
{
  Point point2D(3.0, 4.0);
  algorithm::forceMeasured(point2D);
  BOOST_CHECK_EQUAL(point2D.asText(1), "POINT M (3.0 4.0 0.0)");
}
BOOST_AUTO_TEST_CASE(testPointForceMWithValue)
{
  Point point2D(3.0, 4.0);
  algorithm::forceMeasured(point2D, -9999.0);
  BOOST_CHECK_EQUAL(point2D.asText(1), "POINT M (3.0 4.0 -9999.0)");
}

BOOST_AUTO_TEST_CASE(testPointForceMFromZ)
{
  Point point3D(5.0, 7.0, 1.0);
  algorithm::forceMeasured(point3D);
  BOOST_CHECK_EQUAL(point3D.asText(1), "POINT ZM (5.0 7.0 1.0 0.0)");
}

BOOST_AUTO_TEST_CASE(test_MixedLineString3DM)
{
  LineString lineString;
  lineString.addPoint(Point(1.0, 1.0, 3.0));
  lineString.addPoint(Point(2.0, 3.0, 4.0, 5.0));
  lineString.addPoint(Point(3.0, 3.0, 2.0));
  algorithm::forceMeasured(lineString);
  // should keep 5.0 for the second point
  BOOST_CHECK_EQUAL(
      lineString.asText(1),
      "LINESTRING ZM (1.0 1.0 3.0 0.0,2.0 3.0 4.0 5.0,3.0 3.0 2.0 0.0)");
}

BOOST_AUTO_TEST_SUITE_END()
