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
#include "SFCGAL/io/wkt.h"

#include "../../../test_config.h"

#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_WktTest)

BOOST_AUTO_TEST_CASE(testReadWriter)
{
  std::string filename(SFCGAL_TEST_DIRECTORY);
  filename += "/data/WktTest.txt";

  std::string expectedFilename(SFCGAL_TEST_DIRECTORY);
  expectedFilename += "/data/WktTestExpected.txt";

  std::ifstream ifs(filename.c_str());
  std::ifstream efs(expectedFilename.c_str());

  BOOST_REQUIRE(ifs.good());
  BOOST_REQUIRE(efs.good());

  std::string inputWkt;
  std::string expectedWkt;

  while (std::getline(ifs, inputWkt) && std::getline(efs, expectedWkt)) {
    if (inputWkt[0] == '#' || inputWkt.empty()) {
      continue;
    }

    /*
     * parse wkt and check symmetry
     */
    std::unique_ptr<Geometry> g(io::readWkt(inputWkt));
    std::string const         outputWkt = g->asText(1);
    BOOST_CHECK_EQUAL(expectedWkt, outputWkt);
  }
}

BOOST_AUTO_TEST_SUITE_END()
