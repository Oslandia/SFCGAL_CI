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
#include "SFCGAL/algorithm/force3D.h"
#include "SFCGAL/io/wkt.h"

#include "SFCGAL/detail/tools/Registry.h"

using namespace boost::unit_test;
using namespace SFCGAL;

// note that it relies on a visitor

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_Force3DTest)

BOOST_AUTO_TEST_CASE(testIgnoreEmpty)
{
  tools::Registry const         &registry = tools::Registry::instance();
  std::vector<std::string> const typeNames =
      tools::Registry::instance().getGeometryTypes();

  for (const auto &typeName : typeNames) {
    BOOST_TEST_MESSAGE(typeName);

    std::unique_ptr<Geometry> g(registry.newGeometryByTypeName(typeName));
    BOOST_REQUIRE(g.get() != nullptr);
    algorithm::force3D(*g);
    BOOST_CHECK(g->isEmpty());
  }
}

BOOST_AUTO_TEST_CASE(testPointForceZ)
{
  Point p(3.0, 4.0);
  algorithm::force3D(p);
  BOOST_CHECK_EQUAL(p.asText(1), "POINT Z (3.0 4.0 0.0)");
}
BOOST_AUTO_TEST_CASE(testPointForceZWithValue)
{
  Point p(3.0, 4.0);
  algorithm::force3D(p, -9999.0);
  BOOST_CHECK_EQUAL(p.asText(1), "POINT Z (3.0 4.0 -9999.0)");
}

BOOST_AUTO_TEST_CASE(testPointForceZFromM)
{
  std::unique_ptr<Geometry> ptM(io::readWkt("POINT M (2 3 4)"));
  algorithm::force3D(*ptM);
  BOOST_CHECK_EQUAL(ptM->asText(1), "POINT ZM (2.0 3.0 0.0 4.0)");
}

BOOST_AUTO_TEST_CASE(test_MixedLineString2D3D)
{
  LineString lineString;
  lineString.addPoint(Point(1.0, 1.0));
  lineString.addPoint(Point(2.0, 2.0, 1.0));
  lineString.addPoint(Point(3.0, 3.0));
  algorithm::force3D(lineString);
  // should keep 1.0 for the second point
  BOOST_CHECK_EQUAL(lineString.asText(1),
                    "LINESTRING Z (1.0 1.0 0.0,2.0 2.0 1.0,3.0 3.0 0.0)");
}

BOOST_AUTO_TEST_SUITE_END()
