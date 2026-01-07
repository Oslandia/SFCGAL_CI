// Copyright (c) 2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include <cerrno>
#include <climits>
#include <cmath>
#include <fstream>
#include <iostream>

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
#include "SFCGAL/algorithm/insertPointsWithinTolerance.h"
#include "SFCGAL/io/wkt.h"

#include "../../../test_config.h"

using namespace SFCGAL;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_InsertPointsWithinToleranceTest)

/**
 * Perform tests in test/data/InsertPointsWithinToleranceTest.txt
 */
// NOLINTBEGIN(readability-function-cognitive-complexity)
BOOST_AUTO_TEST_CASE(testFileInsertPointsWithinToleranceTest)
{
  int const argc = framework::master_test_suite().argc;
  char    **argv = framework::master_test_suite().argv;

  // look for options
  int test_one_line = -1;

  for (int i = 0; i < argc; ++i) {
    std::string const argi(argv[i]);

    if (argi == "--line") {
      // only test one line
      if (argc >= i + 1) {
        char *endptr;
        errno    = 0;
        long val = strtol(argv[i + 1], &endptr, 10);

        if (errno == 0 && endptr != argv[i + 1] && *endptr == '\0' &&
            val >= INT_MIN && val <= INT_MAX) {
          test_one_line = static_cast<int>(val);
        } else {
          // Handle error case - invalid number format
          std::cerr << "Error: Invalid number format for --line option: "
                    << argv[i + 1] << "\n";
        }
        ++i;
        continue;
      }
    }
  }

  std::string filename(SFCGAL_TEST_DIRECTORY);
  filename += "/data/InsertPointsWithinToleranceTest.txt";

  std::ifstream ifs(filename.c_str());
  BOOST_REQUIRE(ifs.good());

  int         numLine = 0;
  std::string line;

  while (std::getline(ifs, line)) {
    numLine++;

    if (test_one_line != -1 && test_one_line != numLine) {
      continue;
    }

    if (line[0] == '#' || line.empty()) {
      continue;
    }

    BOOST_TEST_MESSAGE(boost::format("line#%s:%s") % numLine % line);

    std::istringstream iss(line);

    std::string wktBase;
    std::string wktSource;
    std::string toleranceStr;
    std::string wktExpected;

    std::getline(iss, wktBase, '|');
    std::getline(iss, wktSource, '|');
    std::getline(iss, toleranceStr, '|');
    std::getline(iss, wktExpected, '|');

    double const tolerance = std::stod(toleranceStr);

    std::unique_ptr<Geometry> baseGeometry(io::readWkt(wktBase));
    std::unique_ptr<Geometry> sourceGeometry(io::readWkt(wktSource));
    std::unique_ptr<Geometry> expectedGeometry(io::readWkt(wktExpected));

    try {
      std::unique_ptr<Geometry> result = algorithm::insertPointsWithinTolerance(
          *baseGeometry, *sourceGeometry, tolerance);

      // Check if the result matches the expected geometry
      BOOST_CHECK_MESSAGE(
          (*result) == (*expectedGeometry),
          numLine << ": insertPointsWithinTolerance(" << baseGeometry->asText()
                  << ", " << sourceGeometry->asText() << ", " << tolerance
                  << ") should return " << expectedGeometry->asText()
                  << ", but got " << result->asText());
    } catch (std::exception &e) {
      BOOST_CHECK_MESSAGE(false, numLine << ": " << e.what());
    }
  }
}
// NOLINTEND(readability-function-cognitive-complexity)

BOOST_AUTO_TEST_SUITE_END()
