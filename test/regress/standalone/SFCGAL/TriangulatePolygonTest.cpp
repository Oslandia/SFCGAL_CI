// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "../../../test_config.h"

#include <boost/format.hpp>
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
#include "SFCGAL/detail/tools/Log.h"
#include "SFCGAL/io/wkt.h"

#include "SFCGAL/algorithm/area.h"
#include "SFCGAL/triangulate/triangulatePolygon.h"

#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_TriangulatePolygonTest)

/**
 * Triangulate polygon and make some checks
 * @todo Check inPolygon.area3D() == outPolygon.area3D();
 */
BOOST_AUTO_TEST_CASE(testTriangulatePolygon)
{
  // logger().setLogLevel( Logger::Debug );

  std::string filename(SFCGAL_TEST_DIRECTORY);
  filename += "/data/TriangulatePolygonTest.txt";

  std::ifstream ifs(filename.c_str());
  BOOST_REQUIRE(ifs.good());

  int         numLine = 0;
  std::string line;

  while (std::getline(ifs, line)) {
    numLine++;

    if (line[0] == '#' || line.empty()) {
      continue;
    }

    BOOST_TEST_MESSAGE(boost::format("[line#%s]%s") % numLine % line);

    std::istringstream iss(line);
    bool               shouldThrowException = false;
    iss >> shouldThrowException;

    std::string inputWkt;
    std::getline(iss, inputWkt);

    /*
     * parse wkt
     */
    std::unique_ptr<Geometry> g(io::readWkt(inputWkt));

    /*
     * check polygon
     */
    BOOST_CHECK(g->is<Polygon>() || g->is<MultiPolygon>());

    /*
     * triangulate polygon
     */
    TriangulatedSurface triangulatedSurface;

    if (shouldThrowException) {
      BOOST_CHECK_THROW(
          triangulate::triangulatePolygon3D(*g, triangulatedSurface),
          Exception);
      continue;
    }

    BOOST_CHECK_NO_THROW(
        triangulate::triangulatePolygon3D(*g, triangulatedSurface));

    BOOST_TEST_MESSAGE(boost::format("#%1% triangle(s)") %
                       triangulatedSurface.numPatches());

    /*
     * make some checks
     */
    if (!g->isEmpty()) {
      BOOST_CHECK(!triangulatedSurface.isEmpty());
    }

    if (g->is3D()) {
      BOOST_CHECK(triangulatedSurface.is3D());
    }

    // check area
    BOOST_CHECK_CLOSE(algorithm::area3D(*g),
                      algorithm::area3D(triangulatedSurface), 0.5);
  }
}

BOOST_AUTO_TEST_SUITE_END()
