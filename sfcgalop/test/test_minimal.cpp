/**
 * @file test_minimal.cpp
 * @brief Minimal unit tests for SFCGALOP - one test per main functionality
 */

/// @brief Test module definition for Boost.Test framework
#define BOOST_TEST_MODULE SFCGALOP_Minimal_Tests
#include <boost/test/unit_test.hpp>

#include "../error_handler.hpp"
#include "../io.hpp"
#include "../operations/operations.hpp"

#include <SFCGAL/LineString.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>

/// @brief Test geometry loading from WKT format
BOOST_AUTO_TEST_CASE(test_load_wkt)
{
  std::string wkt  = "POINT(1.5 2.5)";
  auto        geom = load_geometry(wkt);
  BOOST_CHECK(geom != nullptr);
}

/// @brief Test geometry validation functionality
BOOST_AUTO_TEST_CASE(test_validate)
{
  SFCGAL::Polygon    polygon;
  SFCGAL::LineString ring;
  ring.addPoint(SFCGAL::Point(0, 0));
  ring.addPoint(SFCGAL::Point(0, 1));
  ring.addPoint(SFCGAL::Point(1, 1));
  ring.addPoint(SFCGAL::Point(1, 0));
  ring.addPoint(SFCGAL::Point(0, 0));
  polygon.setExteriorRing(ring);

  auto result = ErrorHandler::validate_geometry(polygon);
  BOOST_CHECK(result.valid);
}

/// @brief Test area calculation operation
BOOST_AUTO_TEST_CASE(test_area)
{
  SFCGAL::Polygon    polygon;
  SFCGAL::LineString ring;
  ring.addPoint(SFCGAL::Point(0, 0));
  ring.addPoint(SFCGAL::Point(0, 10));
  ring.addPoint(SFCGAL::Point(10, 10));
  ring.addPoint(SFCGAL::Point(10, 0));
  ring.addPoint(SFCGAL::Point(0, 0));
  polygon.setExteriorRing(ring);

  auto result = Operations::execute_operation("area", "", &polygon, nullptr);
  BOOST_CHECK(result.has_value());
}

/// @brief Test distance calculation operation
BOOST_AUTO_TEST_CASE(test_distance)
{
  SFCGAL::Point point1(0, 0);
  SFCGAL::Point point2(3, 4);

  auto result = Operations::execute_operation("distance", "", &point1, &point2);
  BOOST_CHECK(result.has_value());
}

/// @brief Test operations list retrieval
BOOST_AUTO_TEST_CASE(test_operations_list)
{
  auto ops = Operations::get_all_operations_info();
  BOOST_CHECK(!ops.empty());
}

/// @brief Test exception handling mechanism
BOOST_AUTO_TEST_CASE(test_exception)
{
  ErrorHandler::SfcgalopException ex("Test",
                                     ErrorHandler::ErrorCode::INVALID_GEOMETRY);
  BOOST_CHECK_EQUAL(std::string(ex.what()), "Test");
}

/// @brief Test WKT output formatting
BOOST_AUTO_TEST_CASE(test_output_wkt)
{
  SFCGAL::Point     point(1, 2);
  std::stringstream result;
  auto             *old_cout = std::cout.rdbuf(result.rdbuf());
  IO::print_result(std::make_optional(1.0), OutputFormat::WKT, 6);
  std::cout.rdbuf(old_cout);
  BOOST_CHECK(!result.str().empty());
}

/// @brief Test geometric intersection operation
BOOST_AUTO_TEST_CASE(test_intersection)
{
  std::string wkt1 = "POLYGON((0 0, 0 4, 4 4, 4 0, 0 0))";
  std::string wkt2 = "POLYGON((2 2, 2 6, 6 6, 6 2, 2 2))";

  auto geom1 = load_geometry(wkt1);
  auto geom2 = load_geometry(wkt2);
  BOOST_CHECK(geom1 != nullptr);
  BOOST_CHECK(geom2 != nullptr);

  auto result = Operations::execute_operation("intersection", "", geom1.get(),
                                              geom2.get());
  BOOST_CHECK(result.has_value());
}

/// @brief Test convex hull calculation
BOOST_AUTO_TEST_CASE(test_convexhull)
{
  std::string wkt  = "MULTIPOINT((0 0),(1 1),(1 0),(0 1))";
  auto        geom = load_geometry(wkt);
  BOOST_CHECK(geom != nullptr);

  auto result =
      Operations::execute_operation("convexhull", "", geom.get(), nullptr);
  BOOST_CHECK(result.has_value());
}

/// @brief Test geometry validity checking
BOOST_AUTO_TEST_CASE(test_is_valid)
{
  SFCGAL::Point point(1, 2);
  auto result = Operations::execute_operation("is_valid", "", &point, nullptr);
  BOOST_CHECK(result.has_value());
}