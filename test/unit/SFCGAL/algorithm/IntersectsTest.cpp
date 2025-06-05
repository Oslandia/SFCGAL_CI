// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include <cmath>
#include <fstream>

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
#include "SFCGAL/algorithm/intersects.h"
#include "SFCGAL/detail/io/WktWriter.h"
#include "SFCGAL/detail/transform/AffineTransform3.h"
#include "SFCGAL/io/wkt.h"

#include "../../../test_config.h"

using namespace SFCGAL;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_IntersectsTest)

/**
 * Perform tests in test/data/IntersectsTest.txt
 */
BOOST_AUTO_TEST_CASE(testFileIntersectsTest)
{
  // logger().setLogLevel( Logger::Debug );

  std::string filename(SFCGAL_TEST_DIRECTORY);
  filename += "/data/IntersectsTest.txt";

  std::ifstream ifs(filename.c_str());
  BOOST_REQUIRE(ifs.good());

  int         numLine = 0;
  std::string line;

  while (std::getline(ifs, line)) {
    numLine++;

    if (line[0] == '#' || line.empty()) {
      continue;
    }

    BOOST_TEST_MESSAGE(boost::format("line#%s:%s") % numLine % line);

    std::istringstream iss(line);

    std::string distanceDimension;
    std::string wktGA;
    std::string wktGB;
    std::string trueOrFalse;

    std::getline(iss, distanceDimension, '|');
    std::getline(iss, wktGA, '|');
    std::getline(iss, wktGB, '|');
    std::getline(iss, trueOrFalse, '|');

    bool const expected = trueOrFalse == "true";

    std::unique_ptr<Geometry> gA(io::readWkt(wktGA));
    std::unique_ptr<Geometry> gB(io::readWkt(wktGB));

    try {
      if (distanceDimension == "2") {
        BOOST_CHECK_MESSAGE(algorithm::intersects(*gA, *gB) == expected,
                            numLine << ": intersects(" << gA->asText() << ", "
                                    << gB->asText() << ") should be "
                                    << (expected ? "TRUE" : "FALSE"));
        // test symmetric call
        BOOST_CHECK_MESSAGE(algorithm::intersects(*gB, *gA) == expected,
                            numLine << ": intersects(" << gB->asText() << ", "
                                    << gA->asText() << ") should be "
                                    << (expected ? "TRUE" : "FALSE"));
      } else if (distanceDimension == "3") {
        bool got = algorithm::intersects3D(*gA, *gB);
        BOOST_CHECK_MESSAGE(got == expected,
                            numLine << ": intersects3D(" << gA->asText() << ", "
                                    << gB->asText() << ") should be "
                                    << (expected ? "TRUE" : "FALSE"));
        // test symmetric call
        got = algorithm::intersects3D(*gB, *gA);
        BOOST_CHECK_MESSAGE(got == expected,
                            numLine << ": intersects3D(" << gB->asText() << ", "
                                    << gA->asText() << ") should be "
                                    << (expected ? "TRUE" : "FALSE"));
      } else {
        BOOST_CHECK(false);
      }
    } catch (std::exception &e) {
      BOOST_CHECK_MESSAGE(false, numLine << ": " << e.what());
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
