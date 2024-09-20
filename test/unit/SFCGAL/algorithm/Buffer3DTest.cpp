#include "SFCGAL/algorithm/buffer3D.h"
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
#include "SFCGAL/io/OBJ.h"
#include "SFCGAL/io/wkt.h"
#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <iostream>
#include <memory>

#include "../../../test_config.h"

using namespace SFCGAL;
namespace fs = boost::filesystem;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_Buffer3DTest)

bool
compareFiles(const std::string &file1, const std::string &file2)
{
  std::ifstream f1(file1, std::ifstream::binary | std::ifstream::ate);
  std::ifstream f2(file2, std::ifstream::binary | std::ifstream::ate);

  if (f1.fail() || f2.fail()) {
    return false; // File opening error
  }

  if (f1.tellg() != f2.tellg()) {
    return false; // Different sizes
  }

  f1.seekg(0, std::ifstream::beg);
  f2.seekg(0, std::ifstream::beg);
  return std::equal(std::istreambuf_iterator<char>(f1.rdbuf()),
                    std::istreambuf_iterator<char>(),
                    std::istreambuf_iterator<char>(f2.rdbuf()));
}

std::string
readFileContent(const std::string &filePath)
{
  std::ifstream file(filePath);
  if (!file) {
    return "Unable to read file: " + filePath;
  }
  std::stringstream buffer;
  buffer << file.rdbuf();
  return buffer.str();
}

void
compareAndReportFiles(const std::string &expectedFile,
                      const std::string &generatedFile,
                      const std::string &testName)
{
  bool filesMatch = compareFiles(expectedFile, generatedFile);

  if (!filesMatch) {
    BOOST_TEST_MESSAGE("Warning for test " << testName << ":");
    BOOST_TEST_MESSAGE("Generated file does not match the expected file.");
    BOOST_TEST_MESSAGE("Expected file: " << expectedFile);
    BOOST_TEST_MESSAGE("Generated file: " << generatedFile);
    BOOST_TEST_MESSAGE("Content of the generated file:");
    BOOST_TEST_MESSAGE(readFileContent(generatedFile));
  } else {
    BOOST_TEST_MESSAGE("Test " << testName << " passed: files match.");
  }
}

BOOST_AUTO_TEST_CASE(testBuffer3D_Point)
{
  double              radius   = 10.0;
  int                 segments = 16;
  Point               point(0, 0, 0);
  algorithm::Buffer3D buffer3d(point, radius, segments);

  // Create a temporary directory for generated files
  fs::path temp_dir = fs::temp_directory_path() / fs::unique_path();
  fs::create_directories(temp_dir);

  std::vector<algorithm::Buffer3D::BufferType> bufferTypes = {
      algorithm::Buffer3D::ROUND, algorithm::Buffer3D::CYLSPHERE,
      algorithm::Buffer3D::FLAT};

  for (auto bufferType : bufferTypes) {
    std::unique_ptr<PolyhedralSurface> buffer = buffer3d.compute(bufferType);
    BOOST_CHECK(buffer->is3D());
    BOOST_CHECK(buffer->numGeometries() > 0);

    std::string typeName;
    switch (bufferType) {
    case algorithm::Buffer3D::ROUND:
      typeName = "round";
      break;
    case algorithm::Buffer3D::CYLSPHERE:
      typeName = "cylsphere";
      break;
    case algorithm::Buffer3D::FLAT:
      typeName = "flat";
      break;
    }

    // Generate the file with our function
    fs::path generatedFile =
        temp_dir / ("point_" + typeName + "_buffer_3d.obj");
    SFCGAL::io::OBJ::save(*buffer, generatedFile.string());

    // Generate the expected file name
    std::string expectedFile = std::string(SFCGAL_TEST_DIRECTORY) +
                               "/data/bufferfiles/point_" + typeName +
                               "_buffer_3d.obj";

    // Check if the expected file exists
    if (!fs::exists(expectedFile)) {
      std::cout << "Expected file does not exist: " << expectedFile
                << std::endl;
      continue;
    }

    // Compare the files

    compareAndReportFiles(expectedFile, generatedFile.string(),
                          "point_" + typeName + "_buffer");
  }

  // Clean up the temporary directory
  fs::remove_all(temp_dir);
}

BOOST_AUTO_TEST_CASE(testBuffer3D_LineString)
{
  double              radius   = 10.0;
  int                 segments = 16;
  std::vector<Point>  points   = {Point(-100, 0, 0), Point(40, -70, 0),
                                  Point(40, 50, 40), Point(-90, -60, 60),
                                  Point(0, 0, -100), Point(30, 0, 150)};
  LineString          lineString(points);
  algorithm::Buffer3D buffer3d(lineString, radius, segments);

  // Create a temporary directory for generated files
  fs::path temp_dir = fs::temp_directory_path() / fs::unique_path();
  fs::create_directories(temp_dir);

  std::vector<algorithm::Buffer3D::BufferType> bufferTypes = {
      algorithm::Buffer3D::ROUND, algorithm::Buffer3D::CYLSPHERE,
      algorithm::Buffer3D::FLAT};

  for (auto bufferType : bufferTypes) {
    std::unique_ptr<PolyhedralSurface> buffer = buffer3d.compute(bufferType);
    BOOST_CHECK(buffer->is3D());
    BOOST_CHECK(buffer->numGeometries() > 0);

    std::string typeName;
    switch (bufferType) {
    case algorithm::Buffer3D::ROUND:
      typeName = "round";
      break;
    case algorithm::Buffer3D::CYLSPHERE:
      typeName = "cylsphere";
      break;
    case algorithm::Buffer3D::FLAT:
      typeName = "flat";
      break;
    }

    // Generate the file with our function
    fs::path generatedFile =
        temp_dir / ("linestring_" + typeName + "_buffer_3d.obj");
    SFCGAL::io::OBJ::save(*buffer, generatedFile.string());

    // Generate the expected file name
    std::string expectedFile = std::string(SFCGAL_TEST_DIRECTORY) +
                               "/data/bufferfiles/linestring_" + typeName +
                               "_buffer_3d.obj";

    // Check if the expected file exists
    if (!fs::exists(expectedFile)) {
      std::cout << "Expected file does not exist: " << expectedFile
                << std::endl;
      continue;
    }

    // Compare the files

    compareAndReportFiles(expectedFile, generatedFile.string(),
                          "linestring_" + typeName + "_buffer");
  }

  // Clean up the temporary directory
  fs::remove_all(temp_dir);
}

BOOST_AUTO_TEST_CASE(testBuffer3D_InvalidGeometry)
{
  double   radius   = 10.0;
  int      segments = 16;
  Triangle triangle(Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0));
  BOOST_CHECK_THROW(algorithm::Buffer3D(triangle, radius, segments),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
