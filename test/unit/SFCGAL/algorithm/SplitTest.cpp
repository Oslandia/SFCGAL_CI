// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

/**
 * Unit tests for polygon split algorithm using CGAL Arrangement_2.
 */

#include <boost/test/unit_test.hpp>

#include <fstream>
#include <string>
#include <vector>

#include "SFCGAL/Geometry.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/algorithm/split.h"
#include "SFCGAL/io/wkt.h"

using namespace SFCGAL;
using namespace SFCGAL::algorithm;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_SplitTest)

// Helper to read non-comment lines from a file
auto
readTestLines(const std::string &filename) -> std::vector<std::string>
{
  std::vector<std::string> lines;
  std::ifstream            file(filename);
  std::string              line;
  while (std::getline(file, line)) {
    if (!line.empty() && line[0] != '#') {
      lines.push_back(line);
    }
  }
  return lines;
}

BOOST_AUTO_TEST_CASE(test_split_from_files)
{
  std::string inputFile = SFCGAL_TEST_DIRECTORY "/data/split/split_tests.wkt";
  std::string expectedFile =
      SFCGAL_TEST_DIRECTORY "/data/split/split_tests.expected";

  auto inputs   = readTestLines(inputFile);
  auto expected = readTestLines(expectedFile);

  BOOST_REQUIRE_EQUAL(inputs.size(), expected.size());

  for (size_t i = 0; i < inputs.size(); ++i) {
    // Parse input: geometry|linestring
    auto        delimPos      = inputs[i].find('|');
    std::string geometryWkt   = inputs[i].substr(0, delimPos);
    std::string linestringWkt = inputs[i].substr(delimPos + 1);

    auto geometry = io::readWkt(geometryWkt);
    auto blade    = io::readWkt(linestringWkt);

    auto result = split(*geometry, blade->as<LineString>());

    BOOST_CHECK_EQUAL(result->asText(6), expected[i]);
  }
}

BOOST_AUTO_TEST_SUITE_END()
