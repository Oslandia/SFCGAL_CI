// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
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
#include "SFCGAL/algorithm/translate.h"
#include "SFCGAL/io/wkt.h"

#include "SFCGAL/detail/tools/Registry.h"

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_TranslateTest)

BOOST_AUTO_TEST_CASE(testEmpty)
{
  tools::Registry const         &registry = tools::Registry::instance();
  std::vector<std::string> const typeNames =
      tools::Registry::instance().getGeometryTypes();

  for (const auto &typeName : typeNames) {
    BOOST_TEST_MESSAGE(typeName);

    std::unique_ptr<Geometry> g(registry.newGeometryByTypeName(typeName));
    BOOST_REQUIRE(g.get() != nullptr);
    algorithm::translate(*g, 1.0, 1.0, 1.0);
    BOOST_CHECK(g->isEmpty());
  }
}

// TODO complete with 2D/3D test after having renamed translate to translate3D

BOOST_AUTO_TEST_SUITE_END()
