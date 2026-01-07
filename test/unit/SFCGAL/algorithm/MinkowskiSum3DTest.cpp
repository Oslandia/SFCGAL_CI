// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
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
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/algorithm/covers.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/algorithm/minkowskiSum3D.h"
#include "SFCGAL/io/wkt.h"

#include "../../../test_config.h"

using namespace SFCGAL;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_MinkowskiSum3DTest)

/**
 * Perform tests in test/data/MinkowskiSum3DTest.txt
 *
 * File format:
 * testName|geometryA|geometryB|checkType|expectedResult|options
 *
 * checkType values:
 *   covers3D  - verifies algorithm::covers3D(result, expectedResult) is true
 *   empty     - verifies result->isEmpty() is true
 *   not_empty - verifies result->isEmpty() is false and geometryTypeId matches
 *               expectedResult
 *   type_only - only verifies geometryTypeId matches expectedResult
 *
 * options values:
 *   nocheck   - use NoValidityCheck variant of minkowskiSum3D
 *   (empty)   - use default minkowskiSum3D
 */
// NOLINTBEGIN(readability-function-cognitive-complexity)
BOOST_AUTO_TEST_CASE(testFileMinkowskiSum3DTest)
{
  int const argc = framework::master_test_suite().argc;
  char    **argv = framework::master_test_suite().argv;

  // look for options
  int         test_one_line = -1;
  std::string test_name_filter;

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
          std::cerr << "Error: Invalid number format for --line option: "
                    << argv[i + 1] << "\n";
        }
        ++i;
        continue;
      }
    }

    if (argi == "--test") {
      // only test by name
      if (argc >= i + 1) {
        test_name_filter = argv[i + 1];
        ++i;
        continue;
      }
    }
  }

  std::string filename(SFCGAL_TEST_DIRECTORY);
  filename += "/data/MinkowskiSum3DTest.txt";

  std::ifstream ifs(filename.c_str());
  BOOST_REQUIRE(ifs.good());

  int         numLine = 0;
  std::string line;

  while (std::getline(ifs, line)) {
    numLine++;

    if (test_one_line != -1 && test_one_line != numLine) {
      continue;
    }

    if (line.empty() || line[0] == '#') {
      continue;
    }

    std::istringstream iss(line);

    std::string testName;
    std::string wktGA;
    std::string wktGB;
    std::string checkType;
    std::string expectedResult;
    std::string options;

    std::getline(iss, testName, '|');
    std::getline(iss, wktGA, '|');
    std::getline(iss, wktGB, '|');
    std::getline(iss, checkType, '|');
    std::getline(iss, expectedResult, '|');
    std::getline(iss, options, '|');

    // Skip if test name filter is set and doesn't match
    if (!test_name_filter.empty() && testName != test_name_filter) {
      continue;
    }

    BOOST_TEST_MESSAGE(boost::format("line#%s [%s]: %s") % numLine % testName %
                       line.substr(0, 100));

    std::unique_ptr<Geometry> gA(io::readWkt(wktGA));
    std::unique_ptr<Geometry> gB(io::readWkt(wktGB));

    BOOST_REQUIRE_MESSAGE(gA != nullptr,
                          numLine << ": Failed to parse geometry A: " << wktGA);
    BOOST_REQUIRE_MESSAGE(gB != nullptr,
                          numLine << ": Failed to parse geometry B: " << wktGB);

    try {
      std::unique_ptr<Geometry> result;

      // Perform Minkowski sum with appropriate variant
      if (options == "nocheck") {
        result =
            algorithm::minkowskiSum3D(*gA, *gB, algorithm::NoValidityCheck());
      } else {
        result = algorithm::minkowskiSum3D(*gA, *gB);
      }

      BOOST_REQUIRE_MESSAGE(result != nullptr, numLine
                                                   << " [" << testName
                                                   << "]: Result is nullptr");

      // Check based on checkType
      if (checkType == "covers3D") {
        std::unique_ptr<Geometry> expected(io::readWkt(expectedResult));
        BOOST_REQUIRE_MESSAGE(expected != nullptr,
                              numLine << " [" << testName
                                      << "]: Failed to parse expected result: "
                                      << expectedResult);

        BOOST_CHECK_MESSAGE(
            !result->isEmpty(),
            numLine << " [" << testName
                    << "]: Result should not be empty for covers3D check");

        BOOST_CHECK_MESSAGE(
            algorithm::covers3D(*result, *expected),
            numLine << " [" << testName << "]: covers3D(" << result->asText()
                    << ", " << expected->asText() << ") should be TRUE");
      } else if (checkType == "empty") {
        BOOST_CHECK_MESSAGE(result->isEmpty(),
                            numLine << " [" << testName
                                    << "]: Result should be empty, got: "
                                    << result->asText());
      } else if (checkType == "not_empty") {
        BOOST_CHECK_MESSAGE(!result->isEmpty(),
                            numLine << " [" << testName
                                    << "]: Result should not be empty");

        // Parse expected type from expectedResult (e.g.,
        // TYPE_POLYHEDRALSURFACE)
        GeometryType expectedType = TYPE_GEOMETRYCOLLECTION;
        if (expectedResult == "TYPE_POLYHEDRALSURFACE") {
          expectedType = TYPE_POLYHEDRALSURFACE;
        } else if (expectedResult == "TYPE_SOLID") {
          expectedType = TYPE_SOLID;
        } else if (expectedResult == "TYPE_TRIANGULATEDSURFACE") {
          expectedType = TYPE_TRIANGULATEDSURFACE;
        } else if (expectedResult == "TYPE_GEOMETRYCOLLECTION") {
          expectedType = TYPE_GEOMETRYCOLLECTION;
        }

        BOOST_CHECK_MESSAGE(result->geometryTypeId() == expectedType,
                            numLine << " [" << testName << "]: Expected type "
                                    << expectedResult << " but got "
                                    << result->geometryType());
      } else if (checkType == "type_only") {
        // Only check the geometry type, not isEmpty
        GeometryType expectedType = TYPE_GEOMETRYCOLLECTION;
        if (expectedResult == "TYPE_POLYHEDRALSURFACE") {
          expectedType = TYPE_POLYHEDRALSURFACE;
        } else if (expectedResult == "TYPE_SOLID") {
          expectedType = TYPE_SOLID;
        } else if (expectedResult == "TYPE_TRIANGULATEDSURFACE") {
          expectedType = TYPE_TRIANGULATEDSURFACE;
        } else if (expectedResult == "TYPE_GEOMETRYCOLLECTION") {
          expectedType = TYPE_GEOMETRYCOLLECTION;
        }

        BOOST_CHECK_MESSAGE(result->geometryTypeId() == expectedType,
                            numLine << " [" << testName << "]: Expected type "
                                    << expectedResult << " but got "
                                    << result->geometryType());
      } else {
        BOOST_CHECK_MESSAGE(false,
                            numLine << " [" << testName
                                    << "]: Unknown checkType: " << checkType);
      }
    } catch (std::exception &e) {
      BOOST_CHECK_MESSAGE(false,
                          numLine << " [" << testName << "]: " << e.what());
    }
  }
}
// NOLINTEND(readability-function-cognitive-complexity)

BOOST_AUTO_TEST_SUITE_END()
