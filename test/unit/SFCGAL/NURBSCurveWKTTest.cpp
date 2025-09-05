// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "SFCGAL/Exception.h"
#include "SFCGAL/NURBSCurve.h"
#include "SFCGAL/algorithm/distance.h"
#include "SFCGAL/io/wkt.h"
#include <cmath>

using namespace SFCGAL;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_NURBSCurveWKTTest)

//-- Helper functions for WKT tests

std::vector<NURBSCurve::FT>
convertWeights(const std::vector<double> &doubleWeights)
{
  std::vector<NURBSCurve::FT> weights;
  weights.reserve(doubleWeights.size());
  for (const auto &weight : doubleWeights) {
    weights.emplace_back(weight);
  }
  return weights;
}

bool
isNearlyEqual(double valueA, double valueB, double tolerance = 1e-10)
{
  return std::abs(valueA - valueB) < tolerance;
}

bool
isNearlyEqual(const Point &point1, const Point &point2,
              double tolerance = 1e-10)
{
  return algorithm::distance(point1, point2) < NURBSCurve::FT(tolerance);
}

void
checkControlPointsEqual(const NURBSCurve &curve1, const NURBSCurve &curve2,
                        double tolerance = 1e-10)
{
  BOOST_REQUIRE_EQUAL(curve1.numControlPoints(), curve2.numControlPoints());

  for (size_t idx = 0; idx < curve1.numControlPoints(); ++idx) {
    const Point &p1 = curve1.controlPointN(idx);
    const Point &p2 = curve2.controlPointN(idx);

    BOOST_CHECK(isNearlyEqual(p1, p2, tolerance));
  }
}

//-- WKT Writing Tests

BOOST_AUTO_TEST_CASE(testWriteEmptyNURBSCurve)
{
  NURBSCurve emptyCurve;

  std::string wkt = emptyCurve.asText(0);
  BOOST_CHECK_EQUAL(wkt, "NURBSCURVE EMPTY");

  // Test with precision
  std::string wktPrecision = emptyCurve.asText(3);
  BOOST_CHECK_EQUAL(wktPrecision, "NURBSCURVE EMPTY");
}

BOOST_AUTO_TEST_CASE(testWriteBasicNURBSCurve2D)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(5.0, 10.0);
  controlPoints.emplace_back(10.0, 0.0);

  NURBSCurve curve(controlPoints, 2);

  std::string wkt = curve.asText(0);
  BOOST_CHECK_EQUAL(wkt, "NURBSCURVE ((0 0,5 10,10 0),2)");

  // Test with precision
  std::string wktPrecision = curve.asText(1);
  BOOST_CHECK_EQUAL(wktPrecision, "NURBSCURVE ((0.0 0.0,5.0 10.0,10.0 0.0),2)");
}

BOOST_AUTO_TEST_CASE(testWriteWeightedNURBSCurve)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(5.0, 10.0);
  controlPoints.emplace_back(10.0, 0.0);

  auto       weights = convertWeights({1.0, 2.0, 1.0});
  NURBSCurve curve(controlPoints, weights, 2);

  std::string wkt = curve.asText(0);
  BOOST_CHECK_EQUAL(wkt, "NURBSCURVE ((0 0,5 10,10 0),(1,2,1),2)");

  // Test with high precision weights
  auto       preciseWeights = convertWeights({1.0, 2.5, 1.25});
  NURBSCurve preciseCurve(controlPoints, preciseWeights, 2);

  std::string preciseWkt = preciseCurve.asText(2);
  BOOST_CHECK_EQUAL(
      preciseWkt,
      "NURBSCURVE ((0.00 0.00,5.00 10.00,10.00 0.00),(1.00,2.50,1.25),2)");
}

BOOST_AUTO_TEST_CASE(testWriteNURBSCurve3D)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0, 0.0);
  controlPoints.emplace_back(5.0, 5.0, 10.0);
  controlPoints.emplace_back(10.0, 0.0, 0.0);

  NURBSCurve curve(controlPoints, 2);

  std::string wkt = curve.asText(0);
  BOOST_CHECK_EQUAL(wkt, "NURBSCURVE Z ((0 0 0,5 5 10,10 0 0),2)");

  // With weights
  auto       weights = convertWeights({1.0, 1.5, 1.0});
  NURBSCurve weightedCurve(controlPoints, weights, 2);

  std::string weightedWkt = weightedCurve.asText(1);
  BOOST_CHECK_EQUAL(
      weightedWkt,
      "NURBSCURVE Z ((0.0 0.0 0.0,5.0 5.0 10.0,10.0 0.0 0.0),(1.0,1.5,1.0),2)");
}

BOOST_AUTO_TEST_CASE(testWriteNURBSCurveMeasured)
{
  std::vector<Point> controlPoints;
  Point              point1(0.0, 0.0);
  point1.setM(5.0);
  Point point2(5.0, 5.0);
  point2.setM(15.0);
  Point point3(10.0, 0.0);
  point3.setM(25.0);

  controlPoints.push_back(point1);
  controlPoints.push_back(point2);
  controlPoints.push_back(point3);

  NURBSCurve curve(controlPoints, 2);

  std::string wkt = curve.asText(0);
  BOOST_CHECK_EQUAL(wkt, "NURBSCURVE M ((0 0 5,5 5 15,10 0 25),2)");
}

BOOST_AUTO_TEST_CASE(testWriteNURBSCurve3DM)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0, 0.0, 5.0);
  controlPoints.emplace_back(5.0, 5.0, 10.0, 15.0);
  controlPoints.emplace_back(10.0, 0.0, 0.0, 25.0);

  NURBSCurve curve(controlPoints, 2);

  std::string wkt = curve.asText(0);
  BOOST_CHECK_EQUAL(wkt, "NURBSCURVE ZM ((0 0 0 5,5 5 10 15,10 0 0 25),2)");

  // With weights
  auto       weights = convertWeights({1.0, 2.0, 1.0});
  NURBSCurve weightedCurve(controlPoints, weights, 2);

  std::string weightedWkt = weightedCurve.asText(1);
  BOOST_CHECK_EQUAL(weightedWkt, "NURBSCURVE ZM ((0.0 0.0 0.0 5.0,5.0 5.0 10.0 "
                                 "15.0,10.0 0.0 0.0 25.0),(1.0,2.0,1.0),2)");
}

BOOST_AUTO_TEST_CASE(testWriteHighDegreeNURBSCurve)
{
  std::vector<Point> controlPoints;
  for (int idx = 0; idx <= 5; ++idx) {
    controlPoints.emplace_back(static_cast<double>(idx),
                               std::sin(static_cast<double>(idx) * M_PI / 5.0));
  }

  NURBSCurve curve(controlPoints, 4); // Degree 4 with 6 control points

  std::string wkt = curve.asText(3);
  std::cout << wkt << "\n";

  // Should contain all 6 control points and degree 4
  BOOST_CHECK(wkt.find("NURBSCURVE") != std::string::npos);
  BOOST_CHECK(wkt.find(",4)") != std::string::npos); // Should end with degree 4

  // Count control points in WKT (should have 6 coordinate pairs)
  size_t commaCount = std::count(wkt.begin(), wkt.end(), ',');
  BOOST_CHECK(commaCount == 6); // At least 5 commas between points + degree
}

//-- WKT Reading Tests

BOOST_AUTO_TEST_CASE(testReadEmptyNURBSCurve)
{
  std::string wkt = "NURBSCURVE EMPTY";

  auto geometry = io::readWkt(wkt);
  BOOST_REQUIRE(geometry->is<NURBSCurve>());

  const auto &curve = geometry->as<NURBSCurve>();
  BOOST_CHECK(curve.isEmpty());
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 0U);
  BOOST_CHECK_EQUAL(curve.degree(), 0U);
}

BOOST_AUTO_TEST_CASE(testReadBasicNURBSCurve)
{
  std::string wkt = "NURBSCURVE((0 0, 5 10, 10 0), 2)";

  auto geometry = io::readWkt(wkt);
  BOOST_REQUIRE(geometry->is<NURBSCurve>());

  const auto &curve = geometry->as<NURBSCurve>();
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 3U);
  BOOST_CHECK_EQUAL(curve.degree(), 2U);
  BOOST_CHECK(!curve.isRational());
  BOOST_CHECK(!curve.is3D());
  BOOST_CHECK(!curve.isMeasured());

  // Check control points
  BOOST_CHECK_EQUAL(curve.controlPointN(0).x(), NURBSCurve::FT(0.0));
  BOOST_CHECK_EQUAL(curve.controlPointN(0).y(), NURBSCurve::FT(0.0));
  BOOST_CHECK_EQUAL(curve.controlPointN(1).x(), NURBSCurve::FT(5.0));
  BOOST_CHECK_EQUAL(curve.controlPointN(1).y(), NURBSCurve::FT(10.0));
  BOOST_CHECK_EQUAL(curve.controlPointN(2).x(), NURBSCurve::FT(10.0));
  BOOST_CHECK_EQUAL(curve.controlPointN(2).y(), NURBSCurve::FT(0.0));
}

BOOST_AUTO_TEST_CASE(testReadWeightedNURBSCurve)
{
  std::string wkt = "NURBSCURVE((0 0, 5 10, 10 0), (1, 2, 1), 2)";

  auto geometry = io::readWkt(wkt);
  BOOST_REQUIRE(geometry->is<NURBSCurve>());

  const auto &curve = geometry->as<NURBSCurve>();
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 3U);
  BOOST_CHECK_EQUAL(curve.degree(), 2U);
  BOOST_CHECK(curve.isRational());

  // Check weights
  BOOST_CHECK_EQUAL(curve.weight(0), NURBSCurve::FT(1.0));
  BOOST_CHECK_EQUAL(curve.weight(1), NURBSCurve::FT(2.0));
  BOOST_CHECK_EQUAL(curve.weight(2), NURBSCurve::FT(1.0));
}

BOOST_AUTO_TEST_CASE(testReadNURBSCurve3D)
{
  std::string wkt = "NURBSCURVE Z ((0 0 0, 5 5 10, 10 0 0), 2)";

  auto geometry = io::readWkt(wkt);
  BOOST_REQUIRE(geometry->is<NURBSCurve>());

  const auto &curve = geometry->as<NURBSCurve>();
  BOOST_CHECK(curve.is3D());
  BOOST_CHECK(!curve.isMeasured());
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 3U);

  // Check 3D coordinates
  BOOST_CHECK_EQUAL(curve.controlPointN(1).z(), NURBSCurve::FT(10.0));
}

BOOST_AUTO_TEST_CASE(testReadNURBSCurveMeasured)
{
  std::string wkt = "NURBSCURVE M ((0 0 5, 5 5 15, 10 0 25), 2)";

  auto geometry = io::readWkt(wkt);
  BOOST_REQUIRE(geometry->is<NURBSCurve>());

  const auto &curve = geometry->as<NURBSCurve>();
  BOOST_CHECK(!curve.is3D());
  BOOST_CHECK(curve.isMeasured());
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 3U);

  // Check M coordinates
  BOOST_CHECK_EQUAL(curve.controlPointN(0).m(), 5.0);
  BOOST_CHECK_EQUAL(curve.controlPointN(1).m(), 15.0);
  BOOST_CHECK_EQUAL(curve.controlPointN(2).m(), 25.0);
}

BOOST_AUTO_TEST_CASE(testReadNURBSCurve3DM)
{
  std::string wkt =
      "NURBSCURVE ZM ((0 0 0 5, 5 5 10 15, 10 0 0 25), (1, 2, 1), 2)";

  auto geometry = io::readWkt(wkt);
  BOOST_REQUIRE(geometry->is<NURBSCurve>());

  const auto &curve = geometry->as<NURBSCurve>();
  BOOST_CHECK(curve.is3D());
  BOOST_CHECK(curve.isMeasured());
  BOOST_CHECK(curve.isRational());
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 3U);

  // Check all coordinates
  BOOST_CHECK_EQUAL(curve.controlPointN(1).x(), NURBSCurve::FT(5.0));
  BOOST_CHECK_EQUAL(curve.controlPointN(1).y(), NURBSCurve::FT(5.0));
  BOOST_CHECK_EQUAL(curve.controlPointN(1).z(), NURBSCurve::FT(10.0));
  BOOST_CHECK_EQUAL(curve.controlPointN(1).m(), 15.0);
  BOOST_CHECK_EQUAL(curve.weight(1), NURBSCurve::FT(2.0));
}

BOOST_AUTO_TEST_CASE(testReadComplexNURBSCurve)
{
  std::string wkt = "NURBSCURVE((0 0, 3 6, 6 3, 9 0), (1, 1, 1, 1), (0, 0, 0, "
                    "0.5, 1, 1, 1), 2)";

  auto geometry = io::readWkt(wkt);
  BOOST_REQUIRE(geometry->is<NURBSCurve>());

  const auto &curve = geometry->as<NURBSCurve>();
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 4U);
  BOOST_CHECK_EQUAL(curve.degree(), 2U);
  BOOST_CHECK_EQUAL(curve.numKnots(), 7U);

  // Check knot vector
  BOOST_CHECK_EQUAL(curve.knot(0), NURBSCurve::FT(0.0));
  BOOST_CHECK_EQUAL(curve.knot(3), NURBSCurve::FT(0.5));
  BOOST_CHECK_EQUAL(curve.knot(6), NURBSCurve::FT(1.0));
}

//-- WKT Error Handling Tests

BOOST_AUTO_TEST_CASE(testReadInvalidWKTFormat)
{
  // Missing parentheses
  std::string invalidWkt1 = "NURBSCURVE 0 0, 1 1, 2)";
  BOOST_CHECK_THROW(io::readWkt(invalidWkt1), Exception);

  // Wrong geometry type
  std::string invalidWkt3 = "POINT(0 0)";
  auto        pointGeom   = io::readWkt(invalidWkt3);
  BOOST_CHECK(!pointGeom->is<NURBSCurve>());
}

BOOST_AUTO_TEST_CASE(testReadInvalidNURBSData)
{
  // Degree >= number of control points
  std::string invalidDegree = "NURBSCURVE((0 0, 1 1), 3)";
  BOOST_CHECK_THROW(io::readWkt(invalidDegree), Exception);

  // Mismatched weight count
  std::string invalidWeights = "NURBSCURVE((0 0, 1 1, 2 0), (1, 2), 2)";
  BOOST_CHECK_THROW(io::readWkt(invalidWeights), Exception);

  // Negative weight
  std::string negativeWeight = "NURBSCURVE((0 0, 1 1, 2 0), (1, -0.5, 1), 2)";
  BOOST_CHECK_THROW(io::readWkt(negativeWeight), Exception);

  // Zero weight
  std::string zeroWeight = "NURBSCURVE((0 0, 1 1, 2 0), (1, 0, 1), 2)";
  BOOST_CHECK_THROW(io::readWkt(zeroWeight), Exception);

  // Invalid knot vector (wrong size)
  std::string invalidKnots =
      "NURBSCURVE((0 0, 1 1, 2 0), (1, 1, 1), (0, 1), 2)";
  BOOST_CHECK_THROW(io::readWkt(invalidKnots), Exception);

  // Non-decreasing knots
  std::string badKnotOrder =
      "NURBSCURVE((0 0, 1 1, 2 0), (1, 1, 1), (0, 0, 0.5, 0.3, 1, 1), 2)";
  BOOST_CHECK_THROW(io::readWkt(badKnotOrder), Exception);
}

//-- Round-trip tests

BOOST_AUTO_TEST_CASE(testRoundTripBasicNURBS)
{
  std::vector<Point> originalPoints;
  originalPoints.emplace_back(1.5, 2.5);
  originalPoints.emplace_back(3.7, 8.2);
  originalPoints.emplace_back(7.1, 4.9);
  originalPoints.emplace_back(9.3, 1.8);

  NURBSCurve originalCurve(originalPoints, 3);

  // Write to WKT
  std::string wkt = originalCurve.asText(6); // High precision

  // Read back from WKT
  auto readGeometry = io::readWkt(wkt);
  BOOST_REQUIRE(readGeometry->is<NURBSCurve>());

  const auto &readCurve = readGeometry->as<NURBSCurve>();

  // Compare properties
  BOOST_CHECK_EQUAL(originalCurve.numControlPoints(),
                    readCurve.numControlPoints());
  BOOST_CHECK_EQUAL(originalCurve.degree(), readCurve.degree());
  BOOST_CHECK_EQUAL(originalCurve.isRational(), readCurve.isRational());
  BOOST_CHECK_EQUAL(originalCurve.is3D(), readCurve.is3D());
  BOOST_CHECK_EQUAL(originalCurve.isMeasured(), readCurve.isMeasured());

  // Compare control points
  checkControlPointsEqual(originalCurve, readCurve, 1e-6);
}

BOOST_AUTO_TEST_CASE(testRoundTripWeightedNURBS)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.123456, 0.789012);
  controlPoints.emplace_back(2.345678, 4.567890);
  controlPoints.emplace_back(6.789012, 1.234567);

  auto       weights = convertWeights({1.111111, 2.222222, 0.555555});
  NURBSCurve originalCurve(controlPoints, weights, 2);

  std::string wkt          = originalCurve.asText(6);
  auto        readGeometry = io::readWkt(wkt);

  BOOST_REQUIRE(readGeometry->is<NURBSCurve>());
  const auto &readCurve = readGeometry->as<NURBSCurve>();

  BOOST_CHECK(readCurve.isRational());
  checkControlPointsEqual(originalCurve, readCurve, 1e-6);

  // Check weights
  for (size_t idx = 0; idx < originalCurve.numControlPoints(); ++idx) {
    BOOST_CHECK_CLOSE(CGAL::to_double(originalCurve.weight(idx)),
                      CGAL::to_double(readCurve.weight(idx)), 1e-5);
  }
}

BOOST_AUTO_TEST_CASE(testRoundTrip3DNURBS)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(1.1, 2.2, 3.3);
  controlPoints.emplace_back(4.4, 5.5, 6.6);
  controlPoints.emplace_back(7.7, 8.8, 9.9);
  controlPoints.emplace_back(10.1, 11.2, 12.3);

  NURBSCurve originalCurve(controlPoints, 2);

  std::string wkt          = originalCurve.asText(5);
  auto        readGeometry = io::readWkt(wkt);

  BOOST_REQUIRE(readGeometry->is<NURBSCurve>());
  const auto &readCurve = readGeometry->as<NURBSCurve>();

  BOOST_CHECK(readCurve.is3D());
  BOOST_CHECK(!readCurve.isMeasured());
  checkControlPointsEqual(originalCurve, readCurve, 1e-5);
}

BOOST_AUTO_TEST_CASE(testRoundTripMeasuredNURBS)
{
  std::vector<Point> controlPoints;
  Point              point1(1.1, 2.2);
  point1.setM(10.5);
  Point point2(3.3, 4.4);
  point2.setM(20.7);
  Point point3(5.5, 6.6);
  point3.setM(30.9);

  controlPoints.push_back(point1);
  controlPoints.push_back(point2);
  controlPoints.push_back(point3);

  auto       weights = convertWeights({1.2, 1.8, 0.9});
  NURBSCurve originalCurve(controlPoints, weights, 2);

  std::string wkt          = originalCurve.asText(4);
  auto        readGeometry = io::readWkt(wkt);

  BOOST_REQUIRE(readGeometry->is<NURBSCurve>());
  const auto &readCurve = readGeometry->as<NURBSCurve>();

  BOOST_CHECK(!readCurve.is3D());
  BOOST_CHECK(readCurve.isMeasured());
  BOOST_CHECK(readCurve.isRational());
  checkControlPointsEqual(originalCurve, readCurve, 1e-4);
}

BOOST_AUTO_TEST_CASE(testRoundTrip3DMNURBS)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(1.1, 2.2, 3.3, 10.5);
  controlPoints.emplace_back(4.4, 5.5, 6.6, 20.7);
  controlPoints.emplace_back(7.7, 8.8, 9.9, 30.9);

  NURBSCurve originalCurve(controlPoints, 2);

  std::string wkt          = originalCurve.asText(3);
  auto        readGeometry = io::readWkt(wkt);

  BOOST_REQUIRE(readGeometry->is<NURBSCurve>());
  const auto &readCurve = readGeometry->as<NURBSCurve>();

  BOOST_CHECK(readCurve.is3D());
  BOOST_CHECK(readCurve.isMeasured());
  checkControlPointsEqual(originalCurve, readCurve, 1e-3);
}

//-- Compatibility tests with common WKT variations

BOOST_AUTO_TEST_CASE(testWKTVariationsAndWhitespace)
{
  // Test various whitespace and formatting variations
  std::vector<std::string> variations = {
      "NURBSCURVE((0 0,1 1,2 0),2)",
      "NURBSCURVE( ( 0 0 , 1 1 , 2 0 ) , 2 )",
      "nurbscurve((0 0,1 1,2 0),2)",             // lowercase
      "NurbsCurve((0 0,1 1,2 0),2)",             // mixed case
      "NURBSCURVE((0.0 0.0,1.0 1.0,2.0 0.0),2)", // explicit decimals
  };

  for (const auto &wkt : variations) {
    auto geometry = io::readWkt(wkt);
    BOOST_REQUIRE(geometry->is<NURBSCurve>());

    const auto &curve = geometry->as<NURBSCurve>();
    BOOST_CHECK_EQUAL(curve.numControlPoints(), 3U);
    BOOST_CHECK_EQUAL(curve.degree(), 2U);

    // Check first and last control points
    BOOST_CHECK_EQUAL(curve.controlPointN(0).x(), NURBSCurve::FT(0.0));
    BOOST_CHECK_EQUAL(curve.controlPointN(2).x(), NURBSCurve::FT(2.0));
  }
}

BOOST_AUTO_TEST_CASE(testWKTWithScientificNotation)
{
  std::string wkt = "NURBSCURVE((0 0,1e2 1.5e-1,2e0 0),(1e0,2.5e0,1e0),2)";

  auto geometry = io::readWkt(wkt);
  BOOST_REQUIRE(geometry->is<NURBSCurve>());

  const auto &curve = geometry->as<NURBSCurve>();
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 3U);

  // Check scientific notation parsing
  BOOST_CHECK(
      isNearlyEqual(CGAL::to_double(curve.controlPointN(1).x()), 100.0)); // 1e2
  BOOST_CHECK(isNearlyEqual(
      CGAL::to_double(curve.controlPointN(1).y()),
      0.15)); // 1.5e-1, but why "curve.controlPointN(1).y() ==
              // NURBSCurve::FT(0.15) has failed [0.15 != 0.15]"...
  BOOST_CHECK(
      isNearlyEqual(CGAL::to_double(curve.controlPointN(2).x()), 2.0)); // 2e0

  BOOST_CHECK(isNearlyEqual(CGAL::to_double(curve.weight(1)), 2.5)); // 2.5e0
}

//-- PostGIS compatibility tests

BOOST_AUTO_TEST_CASE(testPostGISCompatibleExamples)
{
  // Test examples that should be compatible with PostGIS
  std::vector<std::string> postgisExamples = {
      "NURBSCURVE ((0 0,10 0,20 0),1)",
      "NURBSCURVE ((0 0,5 10,10 0),2)",
      "NURBSCURVE ((0 0,3 7,7 7,10 0),3)",
      "NURBSCURVE ((0 0,5 10,10 0),(1,3,1),2)",
      "NURBSCURVE Z ((0 0 0,5 5 10,10 0 0),(1,2,1),2)",
  };

  for (const auto &wkt : postgisExamples) {
    BOOST_CHECK_NO_THROW(auto geometry = io::readWkt(wkt));

    auto geometry = io::readWkt(wkt);
    BOOST_REQUIRE(geometry->is<NURBSCurve>());

    const auto &curve = geometry->as<NURBSCurve>();
    BOOST_CHECK(!curve.isEmpty());

    // Should be able to write back in compatible format
    std::string writtenWkt = curve.asText(0);
    BOOST_CHECK(!writtenWkt.empty());
    BOOST_CHECK_EQUAL(writtenWkt, wkt);
  }
}

//-- Performance tests for WKT operations

BOOST_AUTO_TEST_CASE(testWKTLargeDataPerformance)
{
  // Create curve with many control points
  std::vector<Point>  manyPoints;
  std::vector<double> manyWeights;

  for (int idx = 0; idx < 100; ++idx) {
    double angle = 2.0 * M_PI * idx / 99.0;
    manyPoints.emplace_back(std::cos(angle), std::sin(angle));
    manyWeights.push_back(
        1.0 + 0.1 * std::sin(3.0 * angle)); // Slight weight variation
  }

  auto       weights = convertWeights(manyWeights);
  NURBSCurve largeCurve(manyPoints, weights, 3);

  // Writing large curve should not fail
  BOOST_CHECK_NO_THROW(std::string wkt = largeCurve.asText(3));

  std::string wkt = largeCurve.asText(3);
  BOOST_CHECK(!wkt.empty());

  // Reading back should work
  BOOST_CHECK_NO_THROW(auto geometry = io::readWkt(wkt));

  auto readGeometry = io::readWkt(wkt);
  BOOST_REQUIRE(readGeometry->is<NURBSCurve>());

  const auto &readCurve = readGeometry->as<NURBSCurve>();
  BOOST_CHECK_EQUAL(readCurve.numControlPoints(), 100U);
  BOOST_CHECK_EQUAL(readCurve.degree(), 3U);
  BOOST_CHECK(readCurve.isRational());
}

BOOST_AUTO_TEST_CASE(testWKTPrecisionHandling)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.123456789123456, 0.987654321987654);
  controlPoints.emplace_back(1.111111111111111, 2.222222222222222);
  controlPoints.emplace_back(3.333333333333333, 0.555555555555555);

  auto       weights = convertWeights({1.123456789, 2.987654321, 0.555555555});
  NURBSCurve preciseCurve(controlPoints, weights, 2);

  // Test different precision levels
  for (int precision = 0; precision <= 10; ++precision) {
    std::string wkt = preciseCurve.asText(precision);
    BOOST_CHECK(!wkt.empty());

    auto readGeometry = io::readWkt(wkt);
    BOOST_REQUIRE(readGeometry->is<NURBSCurve>());

    const auto &readCurve = readGeometry->as<NURBSCurve>();
    BOOST_CHECK_EQUAL(readCurve.numControlPoints(), 3U);

    // Higher precision should give more accurate results
    if (precision >= 6) {
      checkControlPointsEqual(preciseCurve, readCurve, 1e-6);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
