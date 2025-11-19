// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/io/OBJ.h"
#include "SFCGAL/Exception.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/io/wkt.h"
#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

#include "../../../test_config.h"

BOOST_AUTO_TEST_SUITE(SFCGAL_io_OBJWriterTest)

namespace fs = std::filesystem;

auto
compareFiles(const std::string &file1, const std::string &file2) -> bool
{
  std::ifstream stream1(file1, std::ifstream::binary | std::ifstream::ate);
  std::ifstream stream2(file2, std::ifstream::binary | std::ifstream::ate);

  if (stream1.fail() || stream2.fail()) {
    return false; // File opening error
  }

  if (stream1.tellg() != stream2.tellg()) {
    return false; // Different sizes
  }

  stream1.seekg(0, std::ifstream::beg);
  stream2.seekg(0, std::ifstream::beg);
  return std::equal(std::istreambuf_iterator<char>(stream1.rdbuf()),
                    std::istreambuf_iterator<char>(),
                    std::istreambuf_iterator<char>(stream2.rdbuf()));
}

BOOST_AUTO_TEST_CASE(test_all_geometries)
{
  std::vector<std::string> wkt_examples = {
      "POINT Z (1 2 3)", "LINESTRING Z (0 0 0, 1 1 1, 2 2 2)",
      "POLYGON Z ((0 0 0, 1 0 0, 1 1 0, 0 1 0, 0 0 0))",
      "TRIANGLE Z ((0 0 0, 1 0 0, 0 1 0, 0 0 0))",
      // NOLINTNEXTLINE(bugprone-suspicious-missing-comma)
      "POLYHEDRALSURFACE Z (((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)), ((0 0 0, 0 "
      "1 "
      "0, 0 1 1, 0 0 1, 0 0 0)))",
      "SOLID Z ((((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)), ((1 0 0,1 1 0,1 1 1,1 0 "
      "1,1 "
      "0 0)), ((0 1 0,0 1 1,1 1 1,1 1 0,0 1 0)), ((0 0 1,0 1 1,0 1 0,0 0 0,0 0 "
      "1)), ((1 0 1,1 1 1,0 1 1,0 0 1,1 0 1)), ((1 0 0,1 0 1,0 0 1,0 0 0,1 0 "
      "0))))",
      "TIN Z (((0 0 0, 0 0 1, 0 1 0, 0 0 0)), ((0 0 0, 0 1 0, 1 0 0, 0 0 0)))",
      "MULTIPOINT Z ((1 1 1), (2 2 2))",
      "MULTILINESTRING Z ((0 0 0, 1 1 1), (2 2 2, 3 3 3))",
      "MULTIPOLYGON Z (((0 0 0, 1 0 0, 1 1 0, 0 1 0, 0 0 0)), ((2 2 2, 3 2 2, "
      "3 "
      "3 2, 2 3 2, 2 2 2)))",
      "MULTISOLID Z (((((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)),((0 0 1,1 0 1,1 1 1,0 "
      "1 1,0 0 1)),((0 0 0,1 0 0,1 0 1,0 0 1,0 0 0)),((1 1 0,0 1 0,0 1 1,1 1 "
      "1,1 1 0)),((1 0 0,1 1 0,1 1 1,1 0 1,1 0 0)),((0 0 0,0 0 1,0 1 1,0 1 0,0 "
      "0 0)))),((((2 4 6,2 5 6,3 5 6,3 4 6,2 4 6)),((2 4 7,3 4 7,3 5 7,2 5 7,2 "
      "4 7)),((2 4 6,3 4 6,3 4 7,2 4 7,2 4 6)),((3 5 6,2 5 6,2 5 7,3 5 7,3 5 "
      "6)),((3 4 6,3 5 6,3 5 7,3 4 7,3 4 6)),((2 4 6,2 4 7,2 5 7,2 5 6,2 4 "
      "6)))))",
      "GEOMETRYCOLLECTION Z (POINT Z (2.0 3.0 5.0),TRIANGLE Z ((0.0 0.0 "
      "6.0,1.0 "
      "0.0 6.0,1.0 1.0 6.0,0.0 0.0 6.0)))",
      "POINT (1 2)", "LINESTRING (0 0, 1 1, 2 2)",
      "POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))", "TRIANGLE ((0 0, 1 0, 0 1, 0 0))",
      "POLYHEDRALSURFACE (((0 0, 0 1, 1 1, 1 0, 0 0)), ((0 0, 0 1, 0 0)))",
      "TIN (((0 0, 0 1, 1 0, 0 0)), ((0 0, 1 0, 1 1, 0 0)))",
      "MULTIPOINT ((1 1), (2 2))", "MULTILINESTRING ((0 0, 1 1), (2 2, 3 3))",
      "MULTIPOLYGON (((0 0, 1 0, 1 1, 0 1, 0 0)), ((2 2, 3 2, 3 3, 2 3, 2 2)))",
      "MULTISOLID (((((0 0,0 1,1 1,1 0,0 0)),((0 0,1 0,1 0,0 0,0 0)),((0 0,1 "
      "0,1 0,0 0,0 0)),((1 1,0 1,0 1,1 1,1 1)),((1 0,1 1,1 1,1 0,1 0)),((0 0,0 "
      "0,0 1,0 1,0 0)))),((((2 4,2 5,3 5,3 4,2 4)),((2 4,3 4,3 5,2 5,2 4)),((2 "
      "4,3 4,3 4,2 4,2 4)),((3 5,2 5,2 5,3 5,3 5)),((3 4,3 5,3 5,3 4,3 4)),((2 "
      "4,2 4,2 5,2 5,2 4)))))",
      "GEOMETRYCOLLECTION (POINT (2.0 3.0),TRIANGLE ((0.0 0.0,1.0 0.0,1.0 "
      "1.0,0.0 "
      "0.0)))"};

  // Create a temporary directory for generated files
  fs::path temp_dir = fs::temp_directory_path() / random_string();
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
                               "/data/objfiles/" + geomType + ".obj";

    // Check if the expected file exists
    if (!fs::exists(expectedFile)) {
      std::cout << "Expected file does not exist: " << expectedFile << '\n';
      continue;
    }

    // Generate the file with our function
    fs::path generatedFile = temp_dir / (geomType + ".obj");
    SFCGAL::io::OBJ::save(*geom, generatedFile.string());

    // Compare the files
    bool filesMatch = compareFiles(expectedFile, generatedFile.string());
    if (!filesMatch) {
      // Display the generated content for debugging
      std::ifstream generatedStream(generatedFile.string());
      std::string   generatedContent(
          (std::istreambuf_iterator<char>(generatedStream)),
          std::istreambuf_iterator<char>());
      std::cout << "Generated content for " << geomType << ":\n"
                << generatedContent << "\n";
    }
    BOOST_CHECK_MESSAGE(filesMatch,
                        "Output for " << geomType
                                      << " does not match the expected file.\n"
                                         "Expected file: "
                                      << expectedFile
                                      << "\n"
                                         "Generated file: "
                                      << generatedFile.string());
  }

  // Clean up the temporary directory
  fs::remove_all(temp_dir);
}

BOOST_AUTO_TEST_CASE(test_save_to_string)
{
  std::string                       wkt = "POINT (1 2)";
  std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));

  std::string result = SFCGAL::io::OBJ::saveToString(*geom);

  std::string expected = "v 1 2 0\np 1\n";
  BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_CASE(test_save_to_buffer)
{
  std::string                       wkt = "LINESTRING (0 0, 1 1)";
  std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));

  size_t                size = 100;
  std::array<char, 100> buffer{};
  SFCGAL::io::OBJ::saveToBuffer(*geom, buffer.data(), &size);

  std::string result(buffer.data(), size);
  std::string expected = "v 0 0 0\nv 1 1 0\nl 1 2\n";
  BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_CASE(test_buffer_size)
{
  std::string                       wkt = "POINT (1 2 3)";
  std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));

  size_t size = 0;
  SFCGAL::io::OBJ::saveToBuffer(*geom, nullptr, &size);

  BOOST_CHECK_GT(size, 0);

  char *buffer = new char[size];
  SFCGAL::io::OBJ::saveToBuffer(*geom, buffer, &size);

  std::string result(buffer, size);
  std::string expected = "v 1 2 3\np 1\n";
  BOOST_CHECK_EQUAL(result, expected);

  delete[] buffer;
}

BOOST_AUTO_TEST_CASE(test_complex_geometry)
{
  std::string wkt = "GEOMETRYCOLLECTION (POINT (1 1),LINESTRING (0 0, 1 "
                    "1),POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0)))";
  std::unique_ptr<SFCGAL::Geometry> geom(SFCGAL::io::readWkt(wkt));

  std::string result = SFCGAL::io::OBJ::saveToString(*geom);

  std::string expected =
      "v 1 1 0\nv 0 0 0\nv 1 0 0\nv 0 1 0\np 1\nl 2 1\nf 2 3 1 4\n";
  BOOST_CHECK_EQUAL(result, expected);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SFCGAL_io_OBJReaderTest)

BOOST_AUTO_TEST_CASE(test_load_from_string_point)
{
  std::string obj_content = "v 1 2 3\np 1\n";

  std::unique_ptr<SFCGAL::Geometry> geom = SFCGAL::io::OBJ::load(obj_content);

  BOOST_CHECK_EQUAL(geom->geometryType(), "MultiPoint");
  const auto &multipoint = geom->as<SFCGAL::MultiPoint>();
  BOOST_CHECK_EQUAL(multipoint.numGeometries(), 1);

  const auto &point = multipoint.geometryN(0).as<SFCGAL::Point>();
  BOOST_CHECK_EQUAL(point.x(), 1.0);
  BOOST_CHECK_EQUAL(point.y(), 2.0);
  BOOST_CHECK_EQUAL(point.z(), 3.0);
}

BOOST_AUTO_TEST_CASE(test_load_from_string_line)
{
  std::string obj_content = "v 0 0 0\nv 1 1 1\nl 1 2\n";

  std::unique_ptr<SFCGAL::Geometry> geom = SFCGAL::io::OBJ::load(obj_content);

  BOOST_CHECK_EQUAL(geom->geometryType(), "MultiLineString");
  const auto &multilinestring = geom->as<SFCGAL::MultiLineString>();
  BOOST_CHECK_EQUAL(multilinestring.numGeometries(), 1);

  const auto &linestring =
      multilinestring.geometryN(0).as<SFCGAL::LineString>();
  BOOST_CHECK_EQUAL(linestring.numPoints(), 2);

  const auto &point0 = linestring.pointN(0);
  BOOST_CHECK_EQUAL(point0.x(), 0.0);
  BOOST_CHECK_EQUAL(point0.y(), 0.0);
  BOOST_CHECK_EQUAL(point0.z(), 0.0);

  const auto &point1 = linestring.pointN(1);
  BOOST_CHECK_EQUAL(point1.x(), 1.0);
  BOOST_CHECK_EQUAL(point1.y(), 1.0);
  BOOST_CHECK_EQUAL(point1.z(), 1.0);
}

BOOST_AUTO_TEST_CASE(test_load_from_string_triangle)
{
  std::string obj_content = "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n";

  std::unique_ptr<SFCGAL::Geometry> geom = SFCGAL::io::OBJ::load(obj_content);

  BOOST_CHECK_EQUAL(geom->geometryType(), "TriangulatedSurface");
  const auto &triangulated_surface = geom->as<SFCGAL::TriangulatedSurface>();
  BOOST_CHECK_EQUAL(triangulated_surface.numPatches(), 1);

  const auto &triangle = triangulated_surface.patchN(0).as<SFCGAL::Triangle>();

  const auto &vertex0 = triangle.vertex(0);
  BOOST_CHECK_EQUAL(vertex0.x(), 0.0);
  BOOST_CHECK_EQUAL(vertex0.y(), 0.0);
  BOOST_CHECK_EQUAL(vertex0.z(), 0.0);

  const auto &vertex1 = triangle.vertex(1);
  BOOST_CHECK_EQUAL(vertex1.x(), 1.0);
  BOOST_CHECK_EQUAL(vertex1.y(), 0.0);
  BOOST_CHECK_EQUAL(vertex1.z(), 0.0);

  const auto &vertex2 = triangle.vertex(2);
  BOOST_CHECK_EQUAL(vertex2.x(), 0.0);
  BOOST_CHECK_EQUAL(vertex2.y(), 1.0);
  BOOST_CHECK_EQUAL(vertex2.z(), 0.0);
}

BOOST_AUTO_TEST_CASE(test_load_from_string_quad)
{
  std::string obj_content = "v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nf 1 2 3 4\n";

  std::unique_ptr<SFCGAL::Geometry> geom = SFCGAL::io::OBJ::load(obj_content);

  BOOST_CHECK_EQUAL(geom->geometryType(), "PolyhedralSurface");
  const auto &polyhedral_surface = geom->as<SFCGAL::PolyhedralSurface>();
  BOOST_CHECK_EQUAL(polyhedral_surface.numPatches(), 1);

  const auto &polygon = polyhedral_surface.patchN(0).as<SFCGAL::Polygon>();
  BOOST_CHECK_EQUAL(polygon.exteriorRing().numPoints(),
                    5); // 4 vertices + closing
}

BOOST_AUTO_TEST_CASE(test_load_2d_vertices)
{
  std::string obj_content = "v 1 2\nv 3 4\nl 1 2\n";

  std::unique_ptr<SFCGAL::Geometry> geom = SFCGAL::io::OBJ::load(obj_content);

  BOOST_CHECK_EQUAL(geom->geometryType(), "MultiLineString");
  const auto &multilinestring = geom->as<SFCGAL::MultiLineString>();
  const auto &linestring =
      multilinestring.geometryN(0).as<SFCGAL::LineString>();

  const auto &point0 = linestring.pointN(0);
  BOOST_CHECK_EQUAL(point0.x(), 1.0);
  BOOST_CHECK_EQUAL(point0.y(), 2.0);
  BOOST_CHECK_EQUAL(point0.z(), 0.0); // Default Z

  const auto &point1 = linestring.pointN(1);
  BOOST_CHECK_EQUAL(point1.x(), 3.0);
  BOOST_CHECK_EQUAL(point1.y(), 4.0);
  BOOST_CHECK_EQUAL(point1.z(), 0.0); // Default Z
}

BOOST_AUTO_TEST_CASE(test_load_with_comments)
{
  std::string obj_content = "# This is a comment\nv 0 0 0\n# Another "
                            "comment\nv 1 0 0\nv 0 1 0\nf 1 2 3\n";

  std::unique_ptr<SFCGAL::Geometry> geom = SFCGAL::io::OBJ::load(obj_content);

  BOOST_CHECK_EQUAL(geom->geometryType(), "TriangulatedSurface");
  const auto &triangulated_surface = geom->as<SFCGAL::TriangulatedSurface>();
  BOOST_CHECK_EQUAL(triangulated_surface.numPatches(), 1);
}

BOOST_AUTO_TEST_CASE(test_face_with_texture_coordinates)
{
  std::string obj_content = "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1/1/1 2/2/2 3/3/3\n";

  std::unique_ptr<SFCGAL::Geometry> geom = SFCGAL::io::OBJ::load(obj_content);

  BOOST_CHECK_EQUAL(geom->geometryType(), "TriangulatedSurface");
  const auto &triangulated_surface = geom->as<SFCGAL::TriangulatedSurface>();
  BOOST_CHECK_EQUAL(triangulated_surface.numPatches(), 1);
}

BOOST_AUTO_TEST_CASE(test_roundtrip_triangle)
{
  std::string wkt = "TRIANGLE ((0 0 0, 1 0 0, 0 1 0, 0 0 0))";
  std::unique_ptr<SFCGAL::Geometry> original_geom(SFCGAL::io::readWkt(wkt));

  // Save to OBJ string
  std::string obj_content = SFCGAL::io::OBJ::saveToString(*original_geom);

  // Load back from OBJ string
  std::unique_ptr<SFCGAL::Geometry> loaded_geom =
      SFCGAL::io::OBJ::load(obj_content);

  // Both should be TriangulatedSurface (Triangle gets converted)
  BOOST_CHECK_EQUAL(loaded_geom->geometryType(), "TriangulatedSurface");
  const auto &triangulated_surface =
      loaded_geom->as<SFCGAL::TriangulatedSurface>();
  BOOST_CHECK_EQUAL(triangulated_surface.numPatches(), 1);
}

BOOST_AUTO_TEST_CASE(test_error_invalid_vertex_index)
{
  std::string obj_content =
      "v 0 0 0\nf 1 2 3\n"; // Face references non-existent vertices

  BOOST_CHECK_THROW(SFCGAL::io::OBJ::load(obj_content), SFCGAL::Exception);
}

BOOST_AUTO_TEST_CASE(test_error_zero_vertex_index)
{
  std::string obj_content =
      "v 0 0 0\nf 0 1 2\n"; // OBJ indices start at 1, not 0

  BOOST_CHECK_THROW(SFCGAL::io::OBJ::load(obj_content), SFCGAL::Exception);
}

BOOST_AUTO_TEST_CASE(test_error_no_geometry)
{
  std::string obj_content = "# Just comments\n";

  BOOST_CHECK_THROW(SFCGAL::io::OBJ::load(obj_content), SFCGAL::Exception);
}

BOOST_AUTO_TEST_SUITE_END()
