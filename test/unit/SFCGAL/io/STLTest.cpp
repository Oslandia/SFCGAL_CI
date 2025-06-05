// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

// Copyright (c) 2025-2025, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/filesystem.hpp>

#include "SFCGAL/io/STL.h"
#include "SFCGAL/io/wkt.h"
#include <fstream>
#include <iostream>

#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

#include "../../../test_config.h"

BOOST_AUTO_TEST_SUITE(SFCGAL_io_STLWriterTest)

namespace fs = boost::filesystem;

auto
compareFiles(const std::string &file1, const std::string &file2) -> bool
{
  std::ifstream ifstream_file1(file1,
                               std::ifstream::binary | std::ifstream::ate);
  std::ifstream ifstream_file2(file2,
                               std::ifstream::binary | std::ifstream::ate);

  if (ifstream_file1.fail() || ifstream_file2.fail()) {
    return false; // File opening error
  }

  if (ifstream_file1.tellg() != ifstream_file2.tellg()) {
    return false; // Different sizes
  }

  ifstream_file1.seekg(0, std::ifstream::beg);
  ifstream_file2.seekg(0, std::ifstream::beg);
  return std::equal(std::istreambuf_iterator<char>(ifstream_file1.rdbuf()),
                    std::istreambuf_iterator<char>(),
                    std::istreambuf_iterator<char>(ifstream_file2.rdbuf()));
}

BOOST_AUTO_TEST_CASE(test_stl_compatible_geometries)
{
  // STL only supports triangular facets, so we'll test geometries that can be
  // represented as triangles
  std::vector<std::string> wkt_examples = {
      // NOLINTBEGIN(bugprone-suspicious-missing-comma)
      "TRIANGLE Z ((0 0 0, 1 0 0, 0 1 0, 0 0 0))",
      "POLYGON Z ((0 0 0, 1 0 0, 1 1 0, 0 1 0, 0 0 0))",
      "POLYHEDRALSURFACE Z (((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)), ((0 0 "
      "0, 0 1 0, 0 1 1, 0 0 1, 0 0 0)))",
      "SOLID Z ((((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)), ((1 0 0,1 1 0,1 1 1,1 "
      "0 "
      "1,1 0 0)), ((0 1 0,0 1 1,1 1 1,1 1 0,0 1 0)), ((0 0 1,0 1 1,0 1 0,0 "
      "0 "
      "0,0 0 1)), ((1 0 1,1 1 1,0 1 1,0 0 1,1 0 1)), ((1 0 0,1 0 1,0 0 1,0 "
      "0 "
      "0,1 0 0))))",
      "TIN Z (((0 0 0, 0 0 1, 0 1 0, 0 0 0)), ((0 0 0, 0 1 0, 1 0 0, 0 0 "
      "0)))",
      "MULTIPOLYGON Z (((0 0 0, 1 0 0, 1 1 0, 0 1 0, 0 0 0)), ((2 2 2, 3 2 "
      "2, "
      "3 3 2, 2 3 2, 2 2 2)))",
      "MULTISOLID Z (((((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)),((0 0 1,1 0 1,1 1 "
      "1,0 "
      "1 1,0 0 1)),((0 0 0,1 0 0,1 0 1,0 0 1,0 0 0)),((1 1 0,0 1 0,0 1 1,1 "
      "1 "
      "1,1 1 0)),((1 0 0,1 1 0,1 1 1,1 0 1,1 0 0)),((0 0 0,0 0 1,0 1 1,0 1 "
      "0,0 "
      "0 0)))),((((2 4 6,2 5 6,3 5 6,3 4 6,2 4 6)),((2 4 7,3 4 7,3 5 7,2 5 "
      "7,2 "
      "4 7)),((2 4 6,3 4 6,3 4 7,2 4 7,2 4 6)),((3 5 6,2 5 6,2 5 7,3 5 7,3 "
      "5 "
      "6)),((3 4 6,3 5 6,3 5 7,3 4 7,3 4 6)),((2 4 6,2 4 7,2 5 7,2 5 6,2 4 "
      "6)))))",
      "GEOMETRYCOLLECTION Z (TRIANGLE Z ((0.0 0.0 6.0,1.0 0.0 6.0,1.0 1.0 "
      "6.0,0.0 0.0 6.0)))",
      "TRIANGLE ((0 0, 1 0, 0 1, 0 0))",
      "POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))",
      "POLYHEDRALSURFACE (((0 0, 0 1, 1 1, 1 0, 0 0)))",
      "TIN (((0 0, 0 1, 1 0, 0 0)), ((0 0, 1 0, 1 1, 0 0)))",
      "MULTIPOLYGON (((0 0, 1 0, 1 1, 0 1, 0 0)), ((2 2, 3 2, 3 3, 2 3, 2 "
      "2)))",
      "GEOMETRYCOLLECTION (TRIANGLE ((0.0 0.0,1.0 0.0,1.0 1.0,0.0 0.0)))"};
  // NOLINTEND(bugprone-suspicious-missing-comma)

  // Create a temporary directory for generated files
  fs::path temp_dir = fs::temp_directory_path() / fs::unique_path();
  fs::create_directories(temp_dir);

  for (const auto &wkt : wkt_examples) {
    std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));

    // Generate the expected file name
    std::string geomType = geom->geometryType();
    if (geom->is3D()) {
      geomType += "z";
    }
    std::transform(geomType.begin(), geomType.end(), geomType.begin(),
                   ::tolower);

    std::string expectedFile = std::string(SFCGAL_TEST_DIRECTORY) +
                               "/data/stlfiles/" + geomType + ".stl";

    // Check if the expected file exists
    if (!fs::exists(expectedFile)) {
      std::cout << "Expected file does not exist: " << expectedFile << '\n';
      // continue;
    }

    // Generate the file with our function
    fs::path generatedFile = temp_dir / (geomType + ".stl");
    SFCGAL::io::STL::save(*geom, generatedFile.string());

    // Compare the files
    BOOST_CHECK_MESSAGE(compareFiles(expectedFile, generatedFile.string()),
                        "Output for " << geomType
                                      << " does not match the expected file.\n"
                                         "Expected file: "
                                      << expectedFile
                                      << "\n"
                                         "Generated file: "
                                      << generatedFile.string());
  }

  // Clean up the temporary directory
  // fs::remove_all(temp_dir);
}

BOOST_AUTO_TEST_CASE(test_save_to_string)
{
  std::string wkt = "TRIANGLE Z ((0 0 0, 1 0 0, 0 1 0, 0 0 0))";
  std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));

  std::string result = SFCGAL::io::STL::saveToString(*geom);

  std::string expected = "solid SFCGAL_export\n"
                         "  facet normal 0 0 1\n"
                         "    outer loop\n"
                         "      vertex 0 0 0\n"
                         "      vertex 1 0 0\n"
                         "      vertex 0 1 0\n"
                         "    endloop\n"
                         "  endfacet\n"
                         "endsolid SFCGAL_export\n";
  BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_CASE(test_save_to_buffer)
{
  std::string wkt = "TRIANGLE Z ((0 0 0, 1 0 0, 0 1 0, 0 0 0))";
  std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));

  size_t size = 1000;
  char   buffer[1000]; // NOLINT(modernize-avoid-c-arrays)
  SFCGAL::io::STL::saveToBuffer(*geom, buffer, &size);

  std::string result(buffer, size);
  std::string expected = "solid SFCGAL_export\n"
                         "  facet normal 0 0 1\n"
                         "    outer loop\n"
                         "      vertex 0 0 0\n"
                         "      vertex 1 0 0\n"
                         "      vertex 0 1 0\n"
                         "    endloop\n"
                         "  endfacet\n"
                         "endsolid SFCGAL_export\n";
  BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_CASE(test_buffer_size)
{
  std::string wkt = "TRIANGLE Z ((0 0 0, 1 0 0, 0 1 0, 0 0 0))";
  std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));

  size_t size = 0;
  SFCGAL::io::STL::saveToBuffer(*geom, nullptr, &size);

  BOOST_CHECK_GT(size, 0);

  char *buffer = new char[size];
  SFCGAL::io::STL::saveToBuffer(*geom, buffer, &size);

  std::string result(buffer, size);
  std::string expected = "solid SFCGAL_export\n"
                         "  facet normal 0 0 1\n"
                         "    outer loop\n"
                         "      vertex 0 0 0\n"
                         "      vertex 1 0 0\n"
                         "      vertex 0 1 0\n"
                         "    endloop\n"
                         "  endfacet\n"
                         "endsolid SFCGAL_export\n";
  BOOST_CHECK_EQUAL(result, expected);

  delete[] buffer;
}

BOOST_AUTO_TEST_CASE(test_complex_geometry)
{
  std::string wkt =
      "GEOMETRYCOLLECTION Z (TRIANGLE Z ((0 0 0, 1 0 0, 0 1 0, 0 0 0)), "
      "TRIANGLE Z ((1 1 1, 2 1 1, 1 2 1, 1 1 1)))";
  std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));

  std::string result = SFCGAL::io::STL::saveToString(*geom);

  std::string expected = "solid SFCGAL_export\n"
                         "  facet normal 0 0 1\n"
                         "    outer loop\n"
                         "      vertex 0 0 0\n"
                         "      vertex 1 0 0\n"
                         "      vertex 0 1 0\n"
                         "    endloop\n"
                         "  endfacet\n"
                         "  facet normal 0 0 1\n"
                         "    outer loop\n"
                         "      vertex 1 1 1\n"
                         "      vertex 2 1 1\n"
                         "      vertex 1 2 1\n"
                         "    endloop\n"
                         "  endfacet\n"
                         "endsolid SFCGAL_export\n";

  BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_CASE(test_non_stl_geometries)
{
  // Test how the STL writer handles geometries that are not directly
  // representable in STL
  std::vector<std::string> wkt_examples = {
      "POINT Z (1 2 3)", "LINESTRING Z (0 0 0, 1 1 1, 2 2 2)",
      "MULTIPOINT Z ((1 1 1), (2 2 2))",
      "MULTILINESTRING Z ((0 0 0, 1 1 1), (2 2 2, 3 3 3))"};

  for (const auto &wkt : wkt_examples) {
    std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));

    // These geometries should produce an empty STL file
    std::string result   = SFCGAL::io::STL::saveToString(*geom);
    std::string expected = "solid SFCGAL_export\n"
                           "endsolid SFCGAL_export\n";

    BOOST_CHECK_EQUAL(result, expected);
  }
}

BOOST_AUTO_TEST_SUITE_END()
