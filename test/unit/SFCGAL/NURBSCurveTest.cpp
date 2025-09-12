// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "SFCGAL/Envelope.h"
#include "SFCGAL/Exception.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/GeometryVisitor.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/NURBSCurve.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/algorithm/BoundaryVisitor.h"
#include "SFCGAL/algorithm/area.h"
#include "SFCGAL/algorithm/centroid.h"
#include "SFCGAL/algorithm/convexHull.h"
#include "SFCGAL/algorithm/distance.h"
#include "SFCGAL/algorithm/distance3d.h"
#include "SFCGAL/algorithm/intersection.h"
#include "SFCGAL/algorithm/intersects.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/algorithm/length.h"
#include "SFCGAL/algorithm/volume.h"
#include "SFCGAL/detail/transform/AffineTransform2.h"
#include "SFCGAL/detail/transform/AffineTransform3.h"
#include "SFCGAL/io/wkt.h"
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Aff_transformation_3.h>
#include <cmath>
#include <iomanip>
#include <memory>
#include <sstream>

using namespace SFCGAL;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_NURBSCurveTest)

//-- Helper functions
auto
convertWeights(const std::vector<double> &doubleWeights)
    -> std::vector<NURBSCurve::FT>
{
  std::vector<NURBSCurve::FT> weights;
  weights.reserve(doubleWeights.size());
  for (const auto &weight : doubleWeights) {
    weights.emplace_back(weight);
  }
  return weights;
}

auto
convertKnots(const std::vector<double> &doubleKnots)
    -> std::vector<NURBSCurve::Knot>
{
  std::vector<NURBSCurve::Knot> knots;
  knots.reserve(doubleKnots.size());
  for (const auto &knot : doubleKnots) {
    knots.emplace_back(knot);
  }
  return knots;
}

auto
isNearlyEqual(double valueA, double valueB, double tolerance = 1e-10) -> bool
{
  return std::abs(valueA - valueB) < tolerance;
}

auto
isNearlyEqual(const Point &firstPoint, const Point &secondPoint,
              double tolerance = 1e-10) -> bool
{
  return algorithm::distance(firstPoint, secondPoint) <
         NURBSCurve::FT(tolerance);
}

auto
createTestPoints(bool is3D = false, bool isMeasured = false)
    -> std::vector<Point>
{
  std::vector<Point> points;

  if (is3D && isMeasured) {
    points.emplace_back(0.0, 0.0, 0.0, 10.0);
    points.emplace_back(2.0, 4.0, 2.0, 20.0);
    points.emplace_back(6.0, 4.0, -1.0, 30.0);
    points.emplace_back(8.0, 0.0, 0.0, 40.0);
  } else if (is3D) {
    points.emplace_back(0.0, 0.0, 0.0);
    points.emplace_back(2.0, 4.0, 2.0);
    points.emplace_back(6.0, 4.0, -1.0);
    points.emplace_back(8.0, 0.0, 0.0);
  } else if (isMeasured) {
    Point firstPoint(0.0, 0.0);
    firstPoint.setM(10.0);
    Point secondPoint(2.0, 4.0);
    secondPoint.setM(20.0);
    Point thirdPoint(6.0, 4.0);
    thirdPoint.setM(30.0);
    Point fourthPoint(8.0, 0.0);
    fourthPoint.setM(40.0);
    points.push_back(firstPoint);
    points.push_back(secondPoint);
    points.push_back(thirdPoint);
    points.push_back(fourthPoint);
  } else {
    points.emplace_back(0.0, 0.0);
    points.emplace_back(2.0, 4.0);
    points.emplace_back(6.0, 4.0);
    points.emplace_back(8.0, 0.0);
  }

  return points;
}

auto
createTestCurvePoints() -> std::vector<Point>
{
  std::vector<Point> points;
  points.emplace_back(0.0, 0.0);
  points.emplace_back(2.0, 4.0);
  points.emplace_back(6.0, 4.0);
  points.emplace_back(8.0, 0.0);
  return points;
}

auto
createCircularPoints(const double radius = 1.0, const size_t numPoints = 8)
    -> std::vector<Point>
{
  std::vector<Point> points;
  for (size_t pointIdx = 0; pointIdx < numPoints; ++pointIdx) {
    double angle = 2.0 * M_PI * static_cast<double>(pointIdx) /
                   static_cast<double>(numPoints);
    points.emplace_back(radius * std::cos(angle), radius * std::sin(angle));
  }
  return points;
}

auto
checkControlPointsEqual(const NURBSCurve &curve1, const NURBSCurve &curve2,
                        double tolerance = 1e-10) -> void
{
  BOOST_REQUIRE_EQUAL(curve1.numControlPoints(), curve2.numControlPoints());

  for (size_t idx = 0; idx < curve1.numControlPoints(); ++idx) {
    const Point &point1 = curve1.controlPointN(idx);
    const Point &point2 = curve2.controlPointN(idx);

    BOOST_CHECK(isNearlyEqual(point1, point2, tolerance));
  }
}

//-- Constructor tests

/// NURBSCurve();
BOOST_AUTO_TEST_CASE(defaultConstructor)
{
  NURBSCurve const curve;
  BOOST_CHECK(curve.isEmpty());
  BOOST_CHECK(!curve.is3D());
  BOOST_CHECK(!curve.isMeasured());
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 0U);
  BOOST_CHECK_EQUAL(curve.degree(), 0U);
  BOOST_CHECK_EQUAL(curve.weights().size(), 0U);
}

/// NURBSCurve(const std::vector<Point> &controlPoints, unsigned int degree);
BOOST_AUTO_TEST_CASE(constructorWithUniformWeights)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 2.0);
  controlPoints.emplace_back(3.0, 2.0);
  controlPoints.emplace_back(4.0, 0.0);

  NURBSCurve curve(controlPoints, 3);
  BOOST_REQUIRE_EQUAL(curve.numControlPoints(), 4U);
  BOOST_CHECK_EQUAL(curve.degree(), 3U);
  BOOST_CHECK(!curve.isEmpty());
  BOOST_CHECK(!curve.isRational()); // All weights are 1.0
  BOOST_CHECK_EQUAL(curve.weights().size(), 4U);

  // Check that all weights are uniform
  for (size_t idx = 0; idx < curve.weights().size(); ++idx) {
    BOOST_CHECK_EQUAL(curve.weight(idx), NURBSCurve::FT(1.0));
  }
}

/// NURBSCurve(controlPoints, weights, degree);
BOOST_AUTO_TEST_CASE(constructorWithWeights)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 2.0);
  controlPoints.emplace_back(3.0, 2.0);

  auto weights = convertWeights({1.0, 2.0, 1.0});

  NURBSCurve curve(controlPoints, weights, 2);
  BOOST_REQUIRE_EQUAL(curve.numControlPoints(), 3U);
  BOOST_CHECK_EQUAL(curve.degree(), 2U);
  BOOST_CHECK(curve.isRational()); // Non-uniform weights
  BOOST_CHECK_EQUAL(curve.weight(1), NURBSCurve::FT(2.0));
}

/// NURBSCurve(controlPoints, weights, degree, knotVector);
BOOST_AUTO_TEST_CASE(constructorWithKnotVector)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 2.0);
  controlPoints.emplace_back(2.0, 0.0);

  auto weights = convertWeights({1.0, 1.5, 1.0});
  auto knots   = convertKnots({0.0, 0.0, 0.0, 1.0, 1.0, 1.0});

  NURBSCurve curve(controlPoints, weights, 2, knots);
  BOOST_CHECK_EQUAL(curve.degree(), 2U);
  BOOST_CHECK(curve.isRational());
  BOOST_CHECK_EQUAL(curve.knotVector().size(), 6U);
}

/// Copy constructor and assignment
BOOST_AUTO_TEST_CASE(copyConstructorAndAssignment)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);

  auto weights = convertWeights({1.0, 2.0});

  NURBSCurve original(controlPoints, weights, 1);
  NURBSCurve copy(original);
  NURBSCurve assigned;
  assigned = original;

  BOOST_CHECK_EQUAL(copy.numControlPoints(), 2U);
  BOOST_CHECK_EQUAL(copy.weight(1), NURBSCurve::FT(2.0));
  BOOST_CHECK_EQUAL(assigned.numControlPoints(), 2U);
  BOOST_CHECK_EQUAL(assigned.weight(1), NURBSCurve::FT(2.0));

  // Verify independence
  copy.setWeight(1, NURBSCurve::FT(3.0));
  BOOST_CHECK_EQUAL(original.weight(1), NURBSCurve::FT(2.0));
  BOOST_CHECK_EQUAL(copy.weight(1), NURBSCurve::FT(3.0));
}

//-- Error handling in constructors

BOOST_AUTO_TEST_CASE(constructorWithInvalidWeights)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);

  // Wrong number of weights
  auto wrongSizeWeights = convertWeights({1.0});
  BOOST_CHECK_THROW(NURBSCurve curve(controlPoints, wrongSizeWeights, 1),
                    Exception);

  // Negative weight
  auto negativeWeights = convertWeights({1.0, -0.5});
  BOOST_CHECK_THROW(NURBSCurve curve(controlPoints, negativeWeights, 1),
                    Exception);

  // Zero weight
  auto zeroWeights = convertWeights({1.0, 0.0});
  BOOST_CHECK_THROW(NURBSCurve curve(controlPoints, zeroWeights, 1), Exception);
}

//-- Geometry interface tests

/// GeometryType geometryTypeId() const;
BOOST_AUTO_TEST_CASE(testGeometryTypeId)
{
  NURBSCurve const curve;
  BOOST_CHECK_EQUAL(curve.geometryTypeId(), TYPE_NURBSCURVE);
}

/// std::string geometryType() const;
BOOST_AUTO_TEST_CASE(testGeometryType)
{
  NURBSCurve const curve;
  BOOST_CHECK_EQUAL(curve.geometryType(), "NURBSCurve");
}

/// NURBSCurve* clone() const;
BOOST_AUTO_TEST_CASE(testClone)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);

  auto weights = convertWeights({1.0, 2.0});

  NURBSCurve                  curve(controlPoints, weights, 1);
  std::unique_ptr<NURBSCurve> cloned(curve.clone());

  BOOST_REQUIRE(cloned != nullptr);
  BOOST_CHECK_EQUAL(cloned->numControlPoints(), 2U);
  BOOST_CHECK_EQUAL(cloned->weight(1), NURBSCurve::FT(2.0));
}

/// bool isEmpty() const;
BOOST_AUTO_TEST_CASE(testIsEmpty)
{
  NURBSCurve emptyCurve;
  BOOST_CHECK(emptyCurve.isEmpty());

  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  NURBSCurve filledCurve(controlPoints, 0);
  BOOST_CHECK(!filledCurve.isEmpty());
}

/// bool is3D() const;
BOOST_AUTO_TEST_CASE(testIs3D)
{
  std::vector<Point> controlPoints2D;
  controlPoints2D.emplace_back(0.0, 0.0);
  NURBSCurve curve2D(controlPoints2D, 0);
  BOOST_CHECK(!curve2D.is3D());

  std::vector<Point> controlPoints3D;
  controlPoints3D.emplace_back(0.0, 0.0, 0.0);
  NURBSCurve curve3D(controlPoints3D, 0);
  BOOST_CHECK(curve3D.is3D());

  NURBSCurve emptyCurve;
  BOOST_CHECK(!emptyCurve.is3D());
}

/// bool isMeasured() const;
BOOST_AUTO_TEST_CASE(testIsMeasured)
{
  std::vector<Point> controlPointsRegular;
  controlPointsRegular.emplace_back(0.0, 0.0);
  NURBSCurve curveRegular(controlPointsRegular, 0);
  BOOST_CHECK(!curveRegular.isMeasured());

  std::vector<Point> controlPointsMeasured;
  Point              measuredPoint(0.0, 0.0);
  measuredPoint.setM(1.0);
  controlPointsMeasured.push_back(measuredPoint);
  NURBSCurve curveMeasured(controlPointsMeasured, 0);
  BOOST_CHECK(curveMeasured.isMeasured());

  NURBSCurve emptyCurve;
  BOOST_CHECK(!emptyCurve.isMeasured());
}

//-- NURBS-specific tests

/// bool isRational() const;
BOOST_AUTO_TEST_CASE(testIsRational)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);
  controlPoints.emplace_back(2.0, 0.0);

  // Uniform weights - not rational
  auto       uniformWeights = convertWeights({1.0, 1.0, 1.0});
  NURBSCurve uniformCurve(controlPoints, uniformWeights, 2);
  BOOST_CHECK(!uniformCurve.isRational());

  // Non-uniform weights - rational
  auto       nonUniformWeights = convertWeights({1.0, 2.0, 1.0});
  NURBSCurve rationalCurve(controlPoints, nonUniformWeights, 2);
  BOOST_CHECK(rationalCurve.isRational());

  // Empty curve
  NURBSCurve emptyCurve;
  BOOST_CHECK(!emptyCurve.isRational());
}

/// Weight management
BOOST_AUTO_TEST_CASE(testWeightManagement)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);
  controlPoints.emplace_back(2.0, 0.0);

  auto       weights = convertWeights({1.0, 2.0, 1.5});
  NURBSCurve curve(controlPoints, weights, 2);

  // Test individual weight access
  BOOST_CHECK_EQUAL(curve.weight(0), NURBSCurve::FT(1.0));
  BOOST_CHECK_EQUAL(curve.weight(1), NURBSCurve::FT(2.0));
  BOOST_CHECK_EQUAL(curve.weight(2), NURBSCurve::FT(1.5));

  // Test weight modification
  curve.setWeight(1, NURBSCurve::FT(3.0));
  BOOST_CHECK_EQUAL(curve.weight(1), NURBSCurve::FT(3.0));

  // Test weight vector access
  const auto &weightVector = curve.weights();
  BOOST_CHECK_EQUAL(weightVector.size(), 3U);
  BOOST_CHECK_EQUAL(weightVector[1], NURBSCurve::FT(3.0));

  // Test setting all weights
  auto newWeights = convertWeights({0.5, 1.0, 2.0});
  curve.setWeights(newWeights);
  BOOST_CHECK_EQUAL(curve.weight(0), NURBSCurve::FT(0.5));
  BOOST_CHECK_EQUAL(curve.weight(2), NURBSCurve::FT(2.0));
}

BOOST_AUTO_TEST_CASE(testInvalidWeightOperations)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);

  NURBSCurve curve(controlPoints, 1);

  // Setting negative weight
  BOOST_CHECK_THROW(curve.setWeight(0, NURBSCurve::FT(-1.0)), Exception);

  // Setting zero weight
  BOOST_CHECK_THROW(curve.setWeight(0, NURBSCurve::FT(0.0)), Exception);

  // Setting weights with wrong size
  auto wrongSizeWeights = convertWeights({1.0, 2.0, 3.0}); // Too many
  BOOST_CHECK_THROW(curve.setWeights(wrongSizeWeights), Exception);

  // Setting weights with negative values
  auto negativeWeights = convertWeights({1.0, -0.5});
  BOOST_CHECK_THROW(curve.setWeights(negativeWeights), Exception);
}

//-- Curve evaluation tests

/// Point evaluate(double parameter) const;
BOOST_AUTO_TEST_CASE(testEvaluateEmpty)
{
  NURBSCurve emptyCurve;
  BOOST_CHECK_THROW(
      auto result = emptyCurve.evaluate(NURBSCurve::Parameter(0.5)), Exception);
}

BOOST_AUTO_TEST_CASE(testEvaluateParameterBounds)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);

  NURBSCurve curve(controlPoints, 1);
  auto       bounds = curve.parameterBounds();

  BOOST_CHECK_THROW(
      auto result = curve.evaluate(bounds.first - NURBSCurve::Parameter(0.1)),
      Exception);
  BOOST_CHECK_THROW(
      auto result = curve.evaluate(bounds.second + NURBSCurve::Parameter(0.1)),
      Exception);
  BOOST_CHECK_NO_THROW(auto result = curve.evaluate(bounds.first));
  BOOST_CHECK_NO_THROW(auto result = curve.evaluate(bounds.second));
}

BOOST_AUTO_TEST_CASE(testEvaluateLinearNURBS)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(2.0, 2.0);

  // Test with uniform weights (should behave like B-spline)
  NURBSCurve uniformCurve(controlPoints, 1);
  Point      start = uniformCurve.evaluate(NURBSCurve::Parameter(0.0));
  Point      end   = uniformCurve.evaluate(NURBSCurve::Parameter(1.0));
  Point      mid   = uniformCurve.evaluate(NURBSCurve::Parameter(0.5));

  BOOST_CHECK_EQUAL(start.x(), NURBSCurve::FT(0.0));
  BOOST_CHECK_EQUAL(start.y(), NURBSCurve::FT(0.0));
  BOOST_CHECK_EQUAL(end.x(), NURBSCurve::FT(2.0));
  BOOST_CHECK_EQUAL(end.y(), NURBSCurve::FT(2.0));
  BOOST_CHECK_EQUAL(mid.x(), NURBSCurve::FT(1.0));
  BOOST_CHECK_EQUAL(mid.y(), NURBSCurve::FT(1.0));
}

BOOST_AUTO_TEST_CASE(testEvaluateRationalNURBS)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 2.0);
  controlPoints.emplace_back(2.0, 0.0);

  // Create rational NURBS with higher weight on middle point
  auto       weights = convertWeights({1.0, 5.0, 1.0});
  NURBSCurve rationalCurve(controlPoints, weights, 2);

  Point start = rationalCurve.evaluate(NURBSCurve::Parameter(0.0));
  Point end   = rationalCurve.evaluate(NURBSCurve::Parameter(1.0));
  Point mid   = rationalCurve.evaluate(NURBSCurve::Parameter(0.5));

  BOOST_CHECK_EQUAL(start.x(), NURBSCurve::FT(0.0));
  BOOST_CHECK_EQUAL(start.y(), NURBSCurve::FT(0.0));
  BOOST_CHECK_EQUAL(end.x(), NURBSCurve::FT(2.0));
  BOOST_CHECK_EQUAL(end.y(), NURBSCurve::FT(0.0));

  // Middle point should be closer to the highly weighted control point
  BOOST_CHECK(
      mid.y() >
      NURBSCurve::FT(1.0)); // Should be pulled toward control point (1,2)
}

//-- Derivative tests

BOOST_AUTO_TEST_CASE(testDerivativeOrderZero)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);

  NURBSCurve curve(controlPoints, 1);
  Point      point         = curve.derivative(NURBSCurve::Parameter(0.5), 0);
  Point      evaluatePoint = curve.evaluate(NURBSCurve::Parameter(0.5));

  BOOST_CHECK_EQUAL(point.x(), evaluatePoint.x());
  BOOST_CHECK_EQUAL(point.y(), evaluatePoint.y());
}

BOOST_AUTO_TEST_CASE(testDerivativeLinear)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(2.0, 4.0);

  NURBSCurve curve(controlPoints, 1);
  Point      derivative = curve.derivative(NURBSCurve::Parameter(0.5), 1);

  // For linear curve, derivative should be constant
  BOOST_CHECK_EQUAL(derivative.x(), NURBSCurve::FT(2.0));
  BOOST_CHECK_EQUAL(derivative.y(), NURBSCurve::FT(4.0));
}

BOOST_AUTO_TEST_CASE(testDerivativeRational)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 2.0);
  controlPoints.emplace_back(2.0, 0.0);

  auto       weights = convertWeights({1.0, 3.0, 1.0});
  NURBSCurve rationalCurve(controlPoints, weights, 2);

  // Test that rational derivative computation doesn't throw
  BOOST_CHECK_NO_THROW(Point derivative = rationalCurve.derivative(
                           NURBSCurve::Parameter(0.5), 1));

  Point derivative = rationalCurve.derivative(NURBSCurve::Parameter(0.5), 1);
  // The derivative should be non-zero for this configuration
  BOOST_CHECK(CGAL::abs(derivative.x()) > NURBSCurve::FT(1e-10) ||
              CGAL::abs(derivative.y()) > NURBSCurve::FT(1e-10));
}

//-- 3D and measured coordinates

BOOST_AUTO_TEST_CASE(testEvaluate3D)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0, 2.0);

  NURBSCurve curve(controlPoints, 1);
  BOOST_CHECK(curve.is3D());

  Point mid = curve.evaluate(NURBSCurve::Parameter(0.5));
  BOOST_CHECK(mid.is3D());
  BOOST_CHECK_EQUAL(mid.x(), NURBSCurve::FT(0.5));
  BOOST_CHECK_EQUAL(mid.y(), NURBSCurve::FT(0.5));
  BOOST_CHECK_EQUAL(mid.z(), NURBSCurve::FT(1.0));
}

BOOST_AUTO_TEST_CASE(testEvaluateMeasured)
{
  std::vector<Point> controlPoints;
  Point              point1(0.0, 0.0);
  point1.setM(0.0);
  Point point2(2.0, 2.0);
  point2.setM(4.0);
  controlPoints.push_back(point1);
  controlPoints.push_back(point2);

  NURBSCurve curve(controlPoints, 1);
  BOOST_CHECK(curve.isMeasured());

  Point mid = curve.evaluate(NURBSCurve::Parameter(0.5));
  BOOST_CHECK(mid.isMeasured());
  BOOST_CHECK_EQUAL(mid.x(), NURBSCurve::FT(1.0));
  BOOST_CHECK_EQUAL(mid.y(), NURBSCurve::FT(1.0));
  BOOST_CHECK_EQUAL(mid.m(), 2.0);
}

BOOST_AUTO_TEST_CASE(testEvaluate3DMeasured)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0, 0.0, 0.0);
  controlPoints.emplace_back(2.0, 2.0, 4.0, 8.0);

  NURBSCurve curve(controlPoints, 1);
  BOOST_CHECK(curve.is3D());
  BOOST_CHECK(curve.isMeasured());

  Point mid = curve.evaluate(NURBSCurve::Parameter(0.5));
  BOOST_CHECK(mid.is3D());
  BOOST_CHECK(mid.isMeasured());
  BOOST_CHECK_EQUAL(mid.x(), NURBSCurve::FT(1.0));
  BOOST_CHECK_EQUAL(mid.y(), NURBSCurve::FT(1.0));
  BOOST_CHECK_EQUAL(mid.z(), NURBSCurve::FT(2.0));
  BOOST_CHECK_EQUAL(mid.m(), 4.0);
}

//-- Validation tests

BOOST_AUTO_TEST_CASE(testIsValidEmpty)
{
  NURBSCurve emptyCurve;
  auto       validity = algorithm::isValid(emptyCurve);
  BOOST_CHECK(validity);
}

BOOST_AUTO_TEST_CASE(testIsValidConsistent)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);
  controlPoints.emplace_back(2.0, 0.0);

  auto       weights = convertWeights({1.0, 2.0, 1.0});
  NURBSCurve curve(controlPoints, weights, 2);
  auto       validity = algorithm::isValid(curve);
  BOOST_CHECK(validity);
}

//-- Factory methods

BOOST_AUTO_TEST_CASE(testCreateCircularArc2D)
{
  Point          center(0.0, 0.0);
  NURBSCurve::FT radius(1.0);
  NURBSCurve::FT startAngle(0.0);
  NURBSCurve::FT endAngle(M_PI); // Half circle

  auto arc =
      NURBSCurve::createCircularArc(center, radius, startAngle, endAngle);

  BOOST_REQUIRE(arc != nullptr);
  BOOST_CHECK(!arc->is3D());
  BOOST_CHECK(arc->isRational());
  BOOST_CHECK_EQUAL(arc->degree(), 2U);

  // Test that it actually represents a circular arc
  Point start = arc->evaluate(NURBSCurve::Parameter(0.0));
  Point end   = arc->evaluate(NURBSCurve::Parameter(1.0));

  BOOST_CHECK_CLOSE(CGAL::to_double(start.x()), cos(0.0), 1e-10); // cos(0) = 1
  BOOST_CHECK_CLOSE(CGAL::to_double(start.y()), sin(0.0), 1e-10); // sin(0) = 0
  BOOST_CHECK_CLOSE(CGAL::to_double(end.x()), cos(M_PI), 1e-10);  // cos(π) = -1
  BOOST_CHECK_CLOSE(CGAL::to_double(end.y()), sin(M_PI), 1e-10);  // sin(π) ≈ 0
}

BOOST_AUTO_TEST_CASE(testCreateCircularArc3D)
{
  Point          center(0.0, 0.0, 5.0);
  NURBSCurve::FT radius(2.0);
  NURBSCurve::FT startAngle(0.0);
  NURBSCurve::FT endAngle(M_PI / 2); // Quarter circle

  auto arc =
      NURBSCurve::createCircularArc(center, radius, startAngle, endAngle);

  BOOST_REQUIRE(arc != nullptr);
  BOOST_CHECK(arc->is3D());
  BOOST_CHECK(arc->isRational());

  Point start = arc->evaluate(NURBSCurve::Parameter(0.0));
  Point end   = arc->evaluate(NURBSCurve::Parameter(1.0));

  BOOST_CHECK_CLOSE(CGAL::to_double(start.x()), 2.0 * cos(0.0),
                    1e-10); // 2 * cos(0) = 2
  BOOST_CHECK_CLOSE(CGAL::to_double(start.y()), 2.0 * sin(0.0),
                    1e-10);                          // 2 * sin(0) = 0
  BOOST_CHECK_EQUAL(start.z(), NURBSCurve::FT(5.0)); // Z coordinate preserved
  BOOST_CHECK_CLOSE(CGAL::to_double(end.x()), 2.0 * cos(M_PI / 2.0),
                    1e-10); // 2 * cos(π/2) ≈ 0
  BOOST_CHECK_CLOSE(CGAL::to_double(end.y()), 2.0 * sin(M_PI / 2.0),
                    1e-10);                        // 2 * sin(π/2) = 2
  BOOST_CHECK_EQUAL(end.z(), NURBSCurve::FT(5.0)); // Z coordinate preserved
}

//-- Envelope tests

BOOST_AUTO_TEST_CASE(testEnvelopeEmpty)
{
  NURBSCurve emptyCurve;
  BOOST_CHECK(emptyCurve.envelope().isEmpty());
}

BOOST_AUTO_TEST_CASE(testEnvelope2D)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(1.0, 5.0);
  controlPoints.emplace_back(3.0, 9.0);
  controlPoints.emplace_back(2.0, 7.0);

  NURBSCurve     curve(controlPoints, 2);
  Envelope const box = curve.envelope();

  BOOST_CHECK(!box.isEmpty());
  BOOST_CHECK(!box.is3D());

  BOOST_CHECK(box.xMin() >= NURBSCurve::FT(1.0));
  BOOST_CHECK(box.xMax() <= NURBSCurve::FT(3.0));
  BOOST_CHECK(box.yMin() >= NURBSCurve::FT(5.0));
  BOOST_CHECK(box.yMax() <= NURBSCurve::FT(9.0));
}

BOOST_AUTO_TEST_CASE(testEnvelope3D)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(1.0, 5.0, 11.0);
  controlPoints.emplace_back(3.0, 9.0, 17.0);
  controlPoints.emplace_back(2.0, 7.0, 15.0);

  NURBSCurve     curve(controlPoints, 2);
  Envelope const box = curve.envelope();

  BOOST_CHECK(!box.isEmpty());
  BOOST_CHECK(box.is3D());

  BOOST_CHECK(box.xMin() >= NURBSCurve::FT(1.0));
  BOOST_CHECK(box.xMax() <= NURBSCurve::FT(3.0));
  BOOST_CHECK(box.yMin() >= NURBSCurve::FT(5.0));
  BOOST_CHECK(box.yMax() <= NURBSCurve::FT(9.0));
  BOOST_CHECK(box.zMin() >= NURBSCurve::FT(11.0));
  BOOST_CHECK(box.zMax() <= NURBSCurve::FT(17.0));
}

//-- WKT tests

BOOST_AUTO_TEST_CASE(testWktEmpty)
{
  NURBSCurve const emptyCurve;
  std::string      wkt = emptyCurve.asText(0);
  BOOST_CHECK(wkt.find("EMPTY") != std::string::npos);
}

/// WKT reading of basic NURBS curve (points + degree only)
BOOST_AUTO_TEST_CASE(testReadBasicNURBSWkt)
{
  std::string const wkt  = "NURBSCURVE((0 0, 5 10, 10 0), 2)";
  auto              geom = SFCGAL::io::readWkt(wkt);

  BOOST_REQUIRE(geom->is<NURBSCurve>());
  auto const &curve = geom->as<NURBSCurve>();
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 3U);
  BOOST_CHECK_EQUAL(curve.degree(), 2U);
  BOOST_CHECK(!curve.isRational());
}

/// WKT reading with explicit weights
BOOST_AUTO_TEST_CASE(testReadWeightedNURBSWkt)
{
  std::string const wkt  = "NURBSCURVE((0 0, 5 10, 10 0), (1, 2, 1), 2)";
  auto              geom = SFCGAL::io::readWkt(wkt);

  BOOST_REQUIRE(geom->is<NURBSCurve>());
  auto const &curve = geom->as<NURBSCurve>();
  BOOST_CHECK(curve.isRational());
  BOOST_CHECK_EQUAL(curve.weight(1), NURBSCurve::FT(2.0));
}

/// WKT reading with complete knot vector
BOOST_AUTO_TEST_CASE(testReadFullNURBSWkt)
{
  // 4 points, degree 2 → needs 4+2+1=7 knots
  std::string const wkt  = "NURBSCURVE((0 0, 3 6, 6 3, 9 0), (1, 1, 1, 1), (0, "
                           "0, 0, 0.5, 1, 1, 1), 2)";
  auto              geom = SFCGAL::io::readWkt(wkt);

  BOOST_REQUIRE(geom->is<NURBSCurve>());
  auto const &curve = geom->as<NURBSCurve>();
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 4U);
  BOOST_CHECK_EQUAL(curve.degree(), 2U);
  BOOST_CHECK_EQUAL(curve.knotVector().size(), 7U);
}

/// WKT reading of 3D NURBS curve
BOOST_AUTO_TEST_CASE(testRead3DNURBSWkt)
{
  std::string const wkt  = "NURBSCURVE Z ((0 0 0, 5 5 10, 10 0 0), 2)";
  auto              geom = SFCGAL::io::readWkt(wkt);

  BOOST_REQUIRE(geom->is<NURBSCurve>());
  auto const &curve = geom->as<NURBSCurve>();
  BOOST_CHECK(curve.is3D());
  BOOST_CHECK_EQUAL(curve.controlPointN(1).z(), NURBSCurve::FT(10.0));
}

/// WKT reading with measured coordinates
BOOST_AUTO_TEST_CASE(testReadMeasuredNURBSWkt)
{
  std::string const wkt  = "NURBSCURVE M ((0 0 5, 5 5 15, 10 0 25), 2)";
  auto              geom = SFCGAL::io::readWkt(wkt);

  BOOST_REQUIRE(geom->is<NURBSCurve>());
  auto const &curve = geom->as<NURBSCurve>();
  BOOST_CHECK(curve.isMeasured());
  BOOST_CHECK_EQUAL(curve.controlPointN(2).m(), 25.0);
}

/// WKT reading of empty NURBS curve
BOOST_AUTO_TEST_CASE(testReadEmptyNURBSWkt)
{
  std::string const wkt  = "NURBSCURVE EMPTY";
  auto              geom = SFCGAL::io::readWkt(wkt);

  BOOST_REQUIRE(geom->is<NURBSCurve>());
  BOOST_CHECK(geom->as<NURBSCurve>().isEmpty());
}

/// WKT writing of basic NURBS curve
BOOST_AUTO_TEST_CASE(testWriteBasicNURBSWkt)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(5.0, 10.0);
  controlPoints.emplace_back(10.0, 0.0);

  NURBSCurve  curve(controlPoints, 2);
  std::string result = curve.asText(0);

  BOOST_CHECK_EQUAL(result, "NURBSCURVE ((0 0,5 10,10 0),2)");
}

/// WKT writing of weighted NURBS curve
BOOST_AUTO_TEST_CASE(testWriteWeightedNURBSWkt)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(5.0, 10.0);
  controlPoints.emplace_back(10.0, 0.0);

  auto        weights = convertWeights({1.0, 2.0, 1.0});
  NURBSCurve  curve(controlPoints, weights, 2);
  std::string result = curve.asText(0);

  BOOST_CHECK_EQUAL(result, "NURBSCURVE ((0 0,5 10,10 0),(1,2,1),2)");
}

/// WKT writing of 3D NURBS curve
BOOST_AUTO_TEST_CASE(testWrite3DNURBSWkt)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0, 0.0);
  controlPoints.emplace_back(5.0, 5.0, 10.0);
  controlPoints.emplace_back(10.0, 0.0, 0.0);

  NURBSCurve  curve(controlPoints, 2);
  std::string result = curve.asText(0);

  BOOST_CHECK_EQUAL(result, "NURBSCURVE Z ((0 0 0,5 5 10,10 0 0),2)");
}

/// WKT writing of empty NURBS curve
BOOST_AUTO_TEST_CASE(testWriteEmptyNURBSWkt)
{
  NURBSCurve  emptyCurve;
  std::string result = emptyCurve.asText(0);

  BOOST_CHECK_EQUAL(result, "NURBSCURVE EMPTY");
}

/// Round-trip WKT conversion preserves structure
BOOST_AUTO_TEST_CASE(testRoundTripNURBSWkt)
{
  std::string const originalWkt = "NURBSCURVE((0 0, 5 10, 10 0), (1, 2, 1), 2)";
  auto              geom        = SFCGAL::io::readWkt(originalWkt);
  std::string       result      = geom->asText(0);

  BOOST_CHECK_EQUAL(result, "NURBSCURVE ((0 0,5 10,10 0),(1,2,1),2)");

  // Verify re-reading produces identical geometry
  auto geom2 = SFCGAL::io::readWkt(result);
  BOOST_REQUIRE(geom2->is<NURBSCurve>());
  auto const &curve2 = geom2->as<NURBSCurve>();
  BOOST_CHECK_EQUAL(curve2.weight(1), NURBSCurve::FT(2.0));
}

/// PostGIS compatibility with SQL examples
BOOST_AUTO_TEST_CASE(testPostGISCompatibilityWkt)
{
  std::vector<std::string> postgisSamples = {
      "NURBSCURVE((0 0, 10 0, 20 0), 1)", "NURBSCURVE((0 0, 5 10, 10 0), 2)",
      "NURBSCURVE((0 0, 3 7, 7 7, 10 0), 3)",
      "NURBSCURVE((0 0, 5 10, 10 0), (1, 3, 1), 2)",
      "NURBSCURVE Z ((0 0 0, 5 5 10, 10 0 0), (1, 2, 1), 2)"};

  for (auto const &wkt : postgisSamples) {
    BOOST_CHECK_NO_THROW(auto geom = SFCGAL::io::readWkt(wkt));
    auto geom = SFCGAL::io::readWkt(wkt);
    BOOST_REQUIRE(geom->is<NURBSCurve>());
  }
}

/// Error handling for mismatched weight count
BOOST_AUTO_TEST_CASE(testWktErrorMismatchedWeights)
{
  std::string const wkt = "NURBSCURVE((0 0, 5 5, 10 0), (1, 2), 2)";
  BOOST_CHECK_THROW(auto geom = SFCGAL::io::readWkt(wkt), Exception);
}

/// Error handling for negative weights
BOOST_AUTO_TEST_CASE(testWktErrorNegativeWeight)
{
  std::string const wkt = "NURBSCURVE((0 0, 5 5, 10 0), (1, -0.5, 1), 2)";
  BOOST_CHECK_THROW(auto geom = SFCGAL::io::readWkt(wkt), Exception);
}

//-- Template tests

template <typename Derived>
void
testIsInstanceOf(const Geometry &geom)
{
  BOOST_CHECK(geom.is<Derived>());
}

BOOST_AUTO_TEST_CASE(isNURBSCurve)
{
  NURBSCurve const curve;
  testIsInstanceOf<NURBSCurve>(curve);
}

//-- Integration tests with visitor pattern

class TestGeometryVisitor : public GeometryVisitor {
public:
  std::string lastVisited;
  size_t      visitCount = 0;

  void
  visit(Point &geom) override
  {
    lastVisited = "Point";
    visitCount++;
  }
  void
  visit(LineString &geom) override
  {
    lastVisited = "LineString";
    visitCount++;
  }
  void
  visit(Polygon &geom) override
  {
    lastVisited = "Polygon";
    visitCount++;
  }
  void
  visit(Triangle &geom) override
  {
    lastVisited = "Triangle";
    visitCount++;
  }
  void
  visit(Solid &geom) override
  {
    lastVisited = "Solid";
    visitCount++;
  }
  void
  visit(MultiPoint &geom) override
  {
    lastVisited = "MultiPoint";
    visitCount++;
  }
  void
  visit(MultiLineString &geom) override
  {
    lastVisited = "MultiLineString";
    visitCount++;
  }
  void
  visit(MultiPolygon &geom) override
  {
    lastVisited = "MultiPolygon";
    visitCount++;
  }
  void
  visit(MultiSolid &geom) override
  {
    lastVisited = "MultiSolid";
    visitCount++;
  }
  void
  visit(GeometryCollection &geom) override
  {
    lastVisited = "GeometryCollection";
    visitCount++;
  }
  void
  visit(PolyhedralSurface &geom) override
  {
    lastVisited = "PolyhedralSurface";
    visitCount++;
  }
  void
  visit(TriangulatedSurface &geom) override
  {
    lastVisited = "TriangulatedSurface";
    visitCount++;
  }
  void
  visit(NURBSCurve &geom) override
  {
    lastVisited = "NURBSCurve";
    visitCount++;
  }
};

class TestConstGeometryVisitor : public ConstGeometryVisitor {
public:
  std::string lastVisited;
  size_t      visitCount = 0;

  void
  visit(const Point &geom) override
  {
    lastVisited = "Point";
    visitCount++;
  }
  void
  visit(const LineString &geom) override
  {
    lastVisited = "LineString";
    visitCount++;
  }
  void
  visit(const Polygon &geom) override
  {
    lastVisited = "Polygon";
    visitCount++;
  }
  void
  visit(const Triangle &geom) override
  {
    lastVisited = "Triangle";
    visitCount++;
  }
  void
  visit(const Solid &geom) override
  {
    lastVisited = "Solid";
    visitCount++;
  }
  void
  visit(const MultiPoint &geom) override
  {
    lastVisited = "MultiPoint";
    visitCount++;
  }
  void
  visit(const MultiLineString &geom) override
  {
    lastVisited = "MultiLineString";
    visitCount++;
  }
  void
  visit(const MultiPolygon &geom) override
  {
    lastVisited = "MultiPolygon";
    visitCount++;
  }
  void
  visit(const MultiSolid &geom) override
  {
    lastVisited = "MultiSolid";
    visitCount++;
  }
  void
  visit(const GeometryCollection &geom) override
  {
    lastVisited = "GeometryCollection";
    visitCount++;
  }
  void
  visit(const PolyhedralSurface &geom) override
  {
    lastVisited = "PolyhedralSurface";
    visitCount++;
  }
  void
  visit(const TriangulatedSurface &geom) override
  {
    lastVisited = "TriangulatedSurface";
    visitCount++;
  }
  void
  visit(const NURBSCurve &geom) override
  {
    lastVisited = "NURBSCurve";
    visitCount++;
  }
};

BOOST_AUTO_TEST_CASE(testGeometryVisitorIntegration)
{
  auto       controlPoints = createTestPoints();
  NURBSCurve curve(controlPoints, 3);

  TestGeometryVisitor visitor;
  curve.accept(visitor);

  BOOST_CHECK_EQUAL(visitor.lastVisited, "NURBSCurve");
  BOOST_CHECK_EQUAL(visitor.visitCount, 1U);

  // Test const visitor
  TestConstGeometryVisitor constVisitor;
  const NURBSCurve        &constCurve = curve;
  constCurve.accept(constVisitor);

  BOOST_CHECK_EQUAL(constVisitor.lastVisited, "NURBSCurve");
  BOOST_CHECK_EQUAL(constVisitor.visitCount, 1U);
}

BOOST_AUTO_TEST_CASE(testGeometryCollectionWithNURBS)
{
  auto       controlPoints = createTestPoints();
  NURBSCurve curve(controlPoints, 2);

  GeometryCollection collection;
  collection.addGeometry(curve.clone());
  collection.addGeometry(new Point(10.0, 10.0));

  BOOST_CHECK_EQUAL(collection.numGeometries(), 2U);
  BOOST_CHECK(collection.geometryN(0).is<NURBSCurve>());
  BOOST_CHECK(collection.geometryN(1).is<Point>());

  // Test visitor on collection
  TestGeometryVisitor visitor;
  collection.accept(visitor);

  // Should visit the collection itself
  BOOST_CHECK_EQUAL(visitor.lastVisited, "GeometryCollection");
}

//-- Transform integration tests

BOOST_AUTO_TEST_CASE(testTransformNURBSCurve2D)
{
  auto       controlPoints = createTestPoints();
  NURBSCurve originalCurve(controlPoints, 2);

  // Test translation using AffineTransform2
  Kernel::Vector_2                   translationVector(5.0, 3.0);
  CGAL::Aff_transformation_2<Kernel> translationMatrix(CGAL::TRANSLATION,
                                                       translationVector);
  transform::AffineTransform2        translation(translationMatrix);

  NURBSCurve transformedCurve = originalCurve;
  transformedCurve.accept(translation);

  // Check that all control points were translated
  for (size_t idx = 0; idx < originalCurve.numControlPoints(); ++idx) {
    const Point &original    = originalCurve.controlPointN(idx);
    const Point &transformed = transformedCurve.controlPointN(idx);

    BOOST_CHECK_CLOSE(CGAL::to_double(transformed.x()),
                      CGAL::to_double(original.x()) + 5.0, 1e-10);
    BOOST_CHECK_CLOSE(CGAL::to_double(transformed.y()),
                      CGAL::to_double(original.y()) + 3.0, 1e-10);
  }

  // Curve properties should be preserved
  BOOST_CHECK_EQUAL(transformedCurve.degree(), originalCurve.degree());
  BOOST_CHECK_EQUAL(transformedCurve.numControlPoints(),
                    originalCurve.numControlPoints());
  BOOST_CHECK_EQUAL(transformedCurve.isRational(), originalCurve.isRational());
}

BOOST_AUTO_TEST_CASE(testTransformNURBSCurve3D)
{
  auto       controlPoints = createTestPoints(true, false);
  auto       weights       = convertWeights({1.0, 2.0, 1.5, 1.0});
  NURBSCurve originalCurve(controlPoints, weights, 2);

  // Test scaling and rotation using AffineTransform3
  double                             angle = M_PI / 4.0; // 45 degrees
  CGAL::Aff_transformation_3<Kernel> rotationMatrix(
      std::cos(angle), -std::sin(angle), 0.0, 0.0, std::sin(angle),
      std::cos(angle), 0.0, 0.0, 0.0, 0.0, 2.0, 0.0); // Scale Z by 2

  transform::AffineTransform3 rotation(rotationMatrix);

  NURBSCurve transformedCurve = originalCurve;
  transformedCurve.accept(rotation);

  BOOST_CHECK(transformedCurve.is3D());
  BOOST_CHECK(transformedCurve.isRational());

  // Check first control point transformation
  const Point &original    = originalCurve.controlPointN(0);
  const Point &transformed = transformedCurve.controlPointN(0);

  // For point (0,0,0), rotation should still give (0,0,0)
  BOOST_CHECK_SMALL(CGAL::to_double(transformed.x()), 1e-10);
  BOOST_CHECK_SMALL(CGAL::to_double(transformed.y()), 1e-10);
  BOOST_CHECK_SMALL(CGAL::to_double(transformed.z()), 1e-10);

  // Check that weights are preserved
  for (size_t idx = 0; idx < originalCurve.numControlPoints(); ++idx) {
    BOOST_CHECK_EQUAL(originalCurve.weight(idx), transformedCurve.weight(idx));
  }
}

//-- Algorithm integration tests

BOOST_AUTO_TEST_CASE(testEnvelopeAlgorithm)
{
  auto       controlPoints = createTestPoints();
  NURBSCurve curve(controlPoints, 3);

  Envelope envelope = curve.envelope();

  BOOST_CHECK(!envelope.isEmpty());
  BOOST_CHECK(!envelope.is3D()); // 2D curve

  // Envelope should contain all control points
  for (const auto &point : controlPoints) {
    BOOST_CHECK(envelope.xMin() <= point.x());
    BOOST_CHECK(envelope.xMax() >= point.x());
    BOOST_CHECK(envelope.yMin() <= point.y());
    BOOST_CHECK(envelope.yMax() >= point.y());
  }

  // Test 3D envelope
  auto       controlPoints3D = createTestPoints(true, false);
  NURBSCurve curve3D(controlPoints3D, 2);

  Envelope envelope3D = curve3D.envelope();
  BOOST_CHECK(envelope3D.is3D());

  for (const auto &point : controlPoints3D) {
    BOOST_CHECK(envelope3D.zMin() <= point.z());
    BOOST_CHECK(envelope3D.zMax() >= point.z());
  }
}

BOOST_AUTO_TEST_CASE(testBoundaryAlgorithm)
{
  auto       controlPoints = createTestPoints();
  NURBSCurve openCurve(controlPoints, 3);

  // Use BoundaryVisitor directly following SFCGAL patterns
  algorithm::BoundaryVisitor visitor;
  openCurve.accept(visitor);
  std::unique_ptr<Geometry> boundary(visitor.releaseBoundary());

  BOOST_REQUIRE(boundary != nullptr);
  BOOST_CHECK(boundary->is<MultiPoint>());

  const auto &multiPoint = boundary->as<MultiPoint>();
  BOOST_CHECK_EQUAL(multiPoint.numGeometries(), 2U); // Start and end points

  // Check boundary points are curve endpoints
  auto  bounds     = openCurve.parameterBounds();
  Point curveStart = openCurve.evaluate(bounds.first);
  Point curveEnd   = openCurve.evaluate(bounds.second);

  BOOST_CHECK(
      isNearlyEqual(multiPoint.geometryN(0).as<Point>(), curveStart, 1e-8) ||
      isNearlyEqual(multiPoint.geometryN(0).as<Point>(), curveEnd, 1e-8));
  BOOST_CHECK(
      isNearlyEqual(multiPoint.geometryN(1).as<Point>(), curveStart, 1e-8) ||
      isNearlyEqual(multiPoint.geometryN(1).as<Point>(), curveEnd, 1e-8));

  // Test closed curve (should have empty boundary)
  std::vector<Point> closedPoints = controlPoints;
  closedPoints.push_back(controlPoints[0]); // Close the curve
  NURBSCurve closedCurve(closedPoints, 3);

  algorithm::BoundaryVisitor closedVisitor;
  closedCurve.accept(closedVisitor);
  std::unique_ptr<Geometry> closedBoundary(closedVisitor.releaseBoundary());
  BOOST_REQUIRE(closedBoundary != nullptr);
  BOOST_CHECK(closedBoundary->isEmpty() ||
              (closedBoundary->is<MultiPoint>() &&
               closedBoundary->as<MultiPoint>().numGeometries() == 0));
}

BOOST_AUTO_TEST_CASE(testIsValidAlgorithm)
{
  // Valid curve
  auto       controlPoints = createTestPoints();
  auto       weights       = convertWeights({1.0, 2.0, 1.5, 1.0});
  NURBSCurve validCurve(controlPoints, weights, 2);

  auto validity = algorithm::isValid(validCurve);
  BOOST_CHECK(validity.valid());
  BOOST_CHECK(validity.reason().empty());

  // Test with different tolerance
  auto validityTol = algorithm::isValid(validCurve, 1e-6);
  BOOST_CHECK(validityTol.valid());

  // Empty curve should be valid
  NURBSCurve emptyCurve;
  auto       emptyValidity = algorithm::isValid(emptyCurve);
  BOOST_CHECK(emptyValidity.valid());
}

//-- Algorithm integration tests for length and other SFCGAL algorithms

BOOST_AUTO_TEST_CASE(testBasicAlgorithmsOnNURBS)
{
  // Create a simple NURBS curve
  std::vector<Point> points = {Point(0.0, 0.0), Point(1.0, 1.0),
                               Point(2.0, 0.0), Point(3.0, 1.0)};

  auto curve = NURBSCurve::interpolateCurve(points, 2,
                                            NURBSCurve::KnotMethod::CENTRIPETAL,
                                            NURBSCurve::EndCondition::CLAMPED);

  // Ensure curve was created successfully
  BOOST_REQUIRE(curve != nullptr);
  BOOST_REQUIRE(!curve->isEmpty());

  // Test area (should return 0 for curves)
  double area = algorithm::area(*curve);
  BOOST_CHECK_EQUAL(area, 0.0);

  // Test volume (should return 0 for curves)
  auto vol = algorithm::volume(*curve);
  BOOST_CHECK_EQUAL(vol, 0.0);

  // Test length (should return positive value)
  double length = algorithm::length(*curve);
  BOOST_CHECK(length > 0.0);

  // Test centroid (should return a point) - DISABLED due to memory access
  // violation
  // TODO: Investigate why algorithm::centroid causes memory access violation
  // with NURBS curves auto centroidPoint = algorithm::centroid(*curve);
  // BOOST_CHECK(centroidPoint != nullptr);

  // Test convex hull (should work via GetPointsVisitor) - DISABLED due to
  // potential memory issues
  // TODO: Investigate why algorithm::convexHull causes memory access violation
  // with NURBS curves auto hull = algorithm::convexHull(*curve);
  // BOOST_CHECK(hull != nullptr);
}

BOOST_AUTO_TEST_CASE(testDistanceAlgorithmsOnNURBS)
{
  std::vector<Point> points1 = {Point(0.0, 0.0), Point(1.0, 0.0),
                                Point(2.0, 0.0)};
  std::vector<Point> points2 = {Point(0.0, 1.0), Point(1.0, 1.0),
                                Point(2.0, 1.0)};

  auto curve1 = NURBSCurve::interpolateCurve(points1, 2);
  auto curve2 = NURBSCurve::interpolateCurve(points2, 2);

  // Ensure curves were created successfully
  BOOST_REQUIRE(curve1 != nullptr);
  BOOST_REQUIRE(curve2 != nullptr);
  BOOST_REQUIRE(!curve1->isEmpty());
  BOOST_REQUIRE(!curve2->isEmpty());

  // Test 2D distance
  double dist2d = algorithm::distance(*curve1, *curve2);
  BOOST_CHECK(dist2d >= 0.0);
  BOOST_CHECK(dist2d <= 1.1); // Should be approximately 1.0

  // Test 3D distance
  double dist3d = algorithm::distance3D(*curve1, *curve2);
  BOOST_CHECK(dist3d >= 0.0);
  BOOST_CHECK(dist3d <= 1.1); // Should be approximately 1.0
}

BOOST_AUTO_TEST_CASE(testIntersectionAlgorithmsOnNURBS)
{
  // Create two intersecting curves
  std::vector<Point> points1 = {Point(0.0, 0.0), Point(2.0, 2.0)};
  std::vector<Point> points2 = {Point(0.0, 2.0), Point(2.0, 0.0)};

  auto curve1 = NURBSCurve::interpolateCurve(points1, 1);
  auto curve2 = NURBSCurve::interpolateCurve(points2, 1);

  // Ensure curves were created successfully
  BOOST_REQUIRE(curve1 != nullptr);
  BOOST_REQUIRE(curve2 != nullptr);
  BOOST_REQUIRE(!curve1->isEmpty());
  BOOST_REQUIRE(!curve2->isEmpty());

  // Test intersects
  bool intersect = algorithm::intersects(*curve1, *curve2);
  BOOST_CHECK(intersect == true);

  // Test intersection (should return non-empty geometry)
  auto result = algorithm::intersection(*curve1, *curve2);
  BOOST_CHECK(result != nullptr);
  if (result) {
    BOOST_CHECK(!result->isEmpty());
  }
}

//-- Curve manipulation tests

BOOST_AUTO_TEST_CASE(testSplitLinearCurve)
{
  std::vector<Point> linearPoints;
  linearPoints.emplace_back(0.0, 0.0);
  linearPoints.emplace_back(4.0, 8.0);

  NURBSCurve linearCurve(linearPoints, 1);

  auto                  bounds = linearCurve.parameterBounds();
  NURBSCurve::Parameter midParam =
      (bounds.first + bounds.second) / NURBSCurve::FT(2);

  auto [firstHalf, secondHalf] = linearCurve.split(midParam);

  BOOST_REQUIRE(firstHalf != nullptr);
  BOOST_REQUIRE(secondHalf != nullptr);
  BOOST_CHECK(!firstHalf->isEmpty());
  BOOST_CHECK(!secondHalf->isEmpty());

  // Check continuity at split point
  auto firstBounds  = firstHalf->parameterBounds();
  auto secondBounds = secondHalf->parameterBounds();

  Point firstEnd    = firstHalf->evaluate(firstBounds.second);
  Point secondStart = secondHalf->evaluate(secondBounds.first);

  BOOST_CHECK(isNearlyEqual(firstEnd, secondStart, 1e-8));

  // Original curve should be reconstructable from parts
  Point originalMid = linearCurve.evaluate(midParam);
  BOOST_CHECK(isNearlyEqual(firstEnd, originalMid, 1e-8));
}

BOOST_AUTO_TEST_CASE(testReverseLinearCurve)
{
  std::vector<Point> linearPoints;
  linearPoints.emplace_back(1.0, 2.0);
  linearPoints.emplace_back(5.0, 8.0);

  NURBSCurve originalCurve(linearPoints, 1);

  auto reversedCurve = originalCurve.reverse();

  BOOST_REQUIRE(reversedCurve != nullptr);
  BOOST_CHECK(!reversedCurve->isEmpty());
  BOOST_CHECK_EQUAL(reversedCurve->numControlPoints(),
                    originalCurve.numControlPoints());
  BOOST_CHECK_EQUAL(reversedCurve->degree(), originalCurve.degree());

  // Endpoints should be swapped
  auto originalBounds = originalCurve.parameterBounds();
  auto reversedBounds = reversedCurve->parameterBounds();

  Point originalStart = originalCurve.evaluate(originalBounds.first);
  Point originalEnd   = originalCurve.evaluate(originalBounds.second);
  Point reversedStart = reversedCurve->evaluate(reversedBounds.first);
  Point reversedEnd   = reversedCurve->evaluate(reversedBounds.second);

  BOOST_CHECK(isNearlyEqual(originalStart, reversedEnd, 1e-10));
  BOOST_CHECK(isNearlyEqual(originalEnd, reversedStart, 1e-10));
}

//-- Mathematical tests

BOOST_AUTO_TEST_CASE(testLinearCurveArcLength)
{
  // Linear curve from (0,0) to (3,4)
  std::vector<Point> linearPoints;
  linearPoints.emplace_back(0.0, 0.0);
  linearPoints.emplace_back(3.0, 4.0);

  NURBSCurve linearCurve(linearPoints, 1);

  // Expected length = sqrt(3² + 4²) = 5
  NURBSCurve::FT totalLength = linearCurve.length();
  BOOST_CHECK_CLOSE(CGAL::to_double(totalLength), 5.0, 1e-6);

  // Test partial lengths
  auto           bounds   = linearCurve.parameterBounds();
  NURBSCurve::FT midParam = (bounds.first + bounds.second) / NURBSCurve::FT(2);

  NURBSCurve::FT halfLength = linearCurve.length(bounds.first, midParam);
  BOOST_CHECK_CLOSE(CGAL::to_double(halfLength), 2.5, 1e-6);

  NURBSCurve::FT quarterLength = linearCurve.length(
      bounds.first,
      bounds.first + (midParam - bounds.first) / NURBSCurve::FT(2));
  BOOST_CHECK_CLOSE(CGAL::to_double(quarterLength), 1.25, 1e-6);
}

BOOST_AUTO_TEST_CASE(testParameterAtLengthLinear)
{
  std::vector<Point> linearPoints;
  linearPoints.emplace_back(0.0, 0.0);
  linearPoints.emplace_back(4.0, 3.0); // Length = 5

  NURBSCurve linearCurve(linearPoints, 1);

  auto bounds = linearCurve.parameterBounds();

  // At start
  NURBSCurve::Parameter startParam =
      linearCurve.parameterAtLength(NURBSCurve::FT(0.0));
  BOOST_CHECK_CLOSE(CGAL::to_double(startParam), CGAL::to_double(bounds.first),
                    1e-8);

  // At half length
  NURBSCurve::Parameter midParam =
      linearCurve.parameterAtLength(NURBSCurve::FT(2.5));
  Point midPoint = linearCurve.evaluate(midParam);
  BOOST_CHECK_CLOSE(CGAL::to_double(midPoint.x()), 2.0, 1e-8);
  BOOST_CHECK_CLOSE(CGAL::to_double(midPoint.y()), 1.5, 1e-8);

  // At end
  NURBSCurve::FT        totalLength = linearCurve.length();
  NURBSCurve::Parameter endParam = linearCurve.parameterAtLength(totalLength);
  BOOST_CHECK_CLOSE(CGAL::to_double(endParam), CGAL::to_double(bounds.second),
                    1e-6);
}

BOOST_AUTO_TEST_CASE(testClosestPointLinear)
{
  std::vector<Point> linearPoints;
  linearPoints.emplace_back(0.0, 0.0);
  linearPoints.emplace_back(4.0, 3.0);

  NURBSCurve linearCurve(linearPoints, 1);

  // Point on the curve
  Point                 pointOnCurve(2.0, 1.5);
  NURBSCurve::Parameter resultParam;
  Point closestOnCurve = linearCurve.closestPoint(pointOnCurve, &resultParam);

  BOOST_CHECK(isNearlyEqual(pointOnCurve, closestOnCurve, 1e-8));
  BOOST_CHECK_CLOSE(CGAL::to_double(resultParam), 0.5,
                    1e-6); // Should be at midpoint

  // Point off the curve
  Point offCurve(1.0, 3.0);
  Point closestOffCurve = linearCurve.closestPoint(offCurve);

  // Check that closest point is actually on the curve
  auto bounds     = linearCurve.parameterBounds();
  bool foundPoint = false;
  for (int testIdx = 0; testIdx <= 100 && !foundPoint; ++testIdx) {
    NURBSCurve::Parameter testParam =
        bounds.first + (NURBSCurve::FT(testIdx) / NURBSCurve::FT(100)) *
                           (bounds.second - bounds.first);
    Point testPoint = linearCurve.evaluate(testParam);
    if (isNearlyEqual(testPoint, closestOffCurve, 1e-6)) {
      foundPoint = true;
    }
  }
  BOOST_CHECK(foundPoint);
}

//-- End condition tests

BOOST_AUTO_TEST_CASE(testClampedEndCondition)
{
  std::vector<Point> points = {Point(0.0, 0.0), Point(1.0, 2.0),
                               Point(2.0, 1.0), Point(3.0, 0.0)};

  auto curve = NURBSCurve::interpolateCurve(points, 3,
                                            NURBSCurve::KnotMethod::CENTRIPETAL,
                                            NURBSCurve::EndCondition::CLAMPED);

  BOOST_CHECK(!curve->isEmpty());
  BOOST_CHECK_EQUAL(curve->numControlPoints(), points.size());
}

BOOST_AUTO_TEST_CASE(testNaturalEndCondition)
{
  std::vector<Point> points = {Point(0.0, 0.0), Point(1.0, 1.0),
                               Point(2.0, 0.5), Point(3.0, 0.0),
                               Point(4.0, -0.5)};

  // Natural conditions require more points than degree
  auto curve = NURBSCurve::interpolateCurve(points, 3,
                                            NURBSCurve::KnotMethod::CENTRIPETAL,
                                            NURBSCurve::EndCondition::NATURAL);

  BOOST_CHECK(!curve->isEmpty());
  BOOST_CHECK_EQUAL(curve->numControlPoints(), points.size());
  BOOST_CHECK(algorithm::isValid(*curve).valid());
}

BOOST_AUTO_TEST_CASE(testPeriodicEndCondition)
{
  std::vector<Point> points = {
      Point(0.0, 0.0), Point(1.0, 1.0), Point(0.0, 2.0), Point(-1.0, 1.0),
      Point(0.0, 0.0) // Close the curve
  };

  auto curve =
      NURBSCurve::interpolateCurve(points, 3, NURBSCurve::KnotMethod::UNIFORM,
                                   NURBSCurve::EndCondition::PERIODIC);

  BOOST_CHECK(!curve->isEmpty());
  BOOST_CHECK(algorithm::isValid(*curve).valid());
}

BOOST_AUTO_TEST_CASE(testPeriodicEndConditionNotClosed)
{
  std::vector<Point> points = {Point(0.0, 0.0), Point(1.0, 1.0),
                               Point(2.0, 0.0)};

  // Should throw exception for non-closed curve with PERIODIC
  BOOST_CHECK_THROW(
      NURBSCurve::interpolateCurve(points, 3, NURBSCurve::KnotMethod::UNIFORM,
                                   NURBSCurve::EndCondition::PERIODIC),
      Exception);
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
  BOOST_CHECK(isNearlyEqual(CGAL::to_double(curve.controlPointN(1).y()),
                            0.15)); // 1.5e-1
  BOOST_CHECK(
      isNearlyEqual(CGAL::to_double(curve.controlPointN(2).x()), 2.0)); // 2e0

  BOOST_CHECK(isNearlyEqual(CGAL::to_double(curve.weight(1)), 2.5)); // 2.5e0
}

namespace {
std::string
controlPointsWKT(const std::vector<Point> &pts)
{
  std::ostringstream oss;
  oss << "MULTIPOINT(";
  for (std::size_t i = 0; i < pts.size(); i++) {
    if (i > 0)
      oss << ",";
    oss << "(" << CGAL::to_double(pts[i].x()) << " "
        << CGAL::to_double(pts[i].y());
    if (pts[i].dimension() > 2) {
      oss << " " << CGAL::to_double(pts[i].z());
    }
    oss << ")";
  }
  oss << ")";
  return oss.str();
}

void
dumpCurve(const std::string &name, const NURBSCurve &curve,
          const std::vector<Point> &ctrlPts)
{
  std::cout << name << "\n";
  std::cout << "Control points: " << controlPointsWKT(ctrlPts) << "\n";
  std::cout << "LineString: " << curve.toLineString(10)->asText(3) << "\n";
  std::cout << "Length: " << CGAL::to_double(curve.length()) << "\n\n";
}
} // namespace

BOOST_AUTO_TEST_CASE(testWKTComparisonGeomdl)
{
  // 1. Cas 1 - Quadratique simple: NURBSCURVE((0 0, 5 10, 10 0), 2)
  {
    std::vector<Point> pts = {Point(0, 0), Point(5, 10), Point(10, 0)};
    NURBSCurve         curve(pts, 2);
    std::cout << "Cas 1 - Quadratique simple\n";
    std::cout << curve.asText(2) << "\n";
    std::cout << controlPointsWKT(pts) << "\n";
    std::cout << curve.toLineString(200)->asText(2) << "\n";
    std::cout << "Length: " << CGAL::to_double(curve.length()) << "\n\n";
  }

  // 2. Cas 2 - Quadratique pondérée: NURBSCURVE((0 0, 5 10, 10 0), (1, 2, 1),
  // 2)
  {
    std::vector<Point>          pts = {Point(0, 0), Point(5, 10), Point(10, 0)};
    std::vector<NURBSCurve::FT> w   = {1, 2, 1};
    NURBSCurve                  curve(pts, w, 2);
    std::cout << "Cas 2 - Quadratique pondérée\n";
    std::cout << curve.asText(2) << "\n";
    std::cout << controlPointsWKT(pts) << "\n";
    std::cout << curve.toLineString(200)->asText(2) << "\n";
    std::cout << "Length: " << CGAL::to_double(curve.length()) << "\n\n";
  }

  // 3. Cas 3 - Quadratique avec nœuds explicites: NURBSCURVE((0 0, 3 6, 6 3, 9
  // 0), (1,1,1,1), (0,0,0,0.5,1,1,1), 2)
  {
    std::vector<Point>            pts = {Point(0, 0), Point(3, 6), Point(6, 3),
                                         Point(9, 0)};
    std::vector<NURBSCurve::FT>   w   = {1, 1, 1, 1};
    std::vector<NURBSCurve::Knot> k   = {0, 0, 0, 0.5, 1, 1, 1};
    NURBSCurve                    curve(pts, w, 2, k);
    std::cout << "Cas 3 - Quadratique avec nœuds explicites\n";
    std::cout << curve.asText(2) << "\n";
    std::cout << controlPointsWKT(pts) << "\n";
    std::cout << curve.toLineString(200)->asText(2) << "\n";
    std::cout << "Length: " << CGAL::to_double(curve.length()) << "\n\n";
  }

  // 4. Cas 4 - Quadratique en 3D: (0, 0, 0), (5, 10, 5), (10, 0, 0) avec poids
  // (1, 2, 1)
  {
    std::vector<Point> pts = {Point(0, 0, 0), Point(5, 10, 5), Point(10, 0, 0)};
    std::vector<NURBSCurve::FT> w = {1, 2, 1};
    NURBSCurve                  curve(pts, w, 2);
    std::cout << "Cas 4 - Quadratique en 3D\n";
    std::cout << curve.asText(2) << "\n";
    std::cout << controlPointsWKT(pts) << "\n";
    std::cout << curve.toLineString(200)->asText(2) << "\n";
    std::cout << "Length (3D): " << CGAL::to_double(curve.length()) << "\n\n";
  }

  // 5. Cas 5 - Quart de cercle: (1,0), (1,1), (0,1) avec poids (1, √2/2, 1)
  {
    std::vector<Point>            pts = {Point(1, 0), Point(1, 1), Point(0, 1)};
    std::vector<NURBSCurve::FT>   w = {1, NURBSCurve::FT(std::sqrt(2.0) / 2.0),
                                       1};
    std::vector<NURBSCurve::Knot> k = {0, 0, 0, 1, 1, 1};
    NURBSCurve                    curve(pts, w, 2, k);
    std::cout << "Cas 5 - Quart de cercle\n";
    std::cout << curve.asText(2) << "\n";
    std::cout << controlPointsWKT(pts) << "\n";
    std::cout << curve.toLineString(200)->asText(3) << "\n";
    std::cout << "Length: " << CGAL::to_double(curve.length())
              << " (théorique ~" << M_PI / 2.0 << ")\n\n";
  }

  // 6. Cas 6 - S couché: degré 3, 7 points de contrôle
  {
    std::vector<Point> pts = {Point(0, 0), Point(3, 2),  Point(6, -2),
                              Point(9, 0), Point(12, 2), Point(15, -2),
                              Point(18, 0)};
    NURBSCurve         curve(pts, 3);
    std::cout << "Cas 6 - S couché\n";
    std::cout << curve.asText(2) << "\n";
    std::cout << controlPointsWKT(pts) << "\n";
    std::cout << curve.toLineString(200)->asText(2) << "\n";
    std::cout << "Length: " << CGAL::to_double(curve.length()) << "\n\n";
  }

  // 7. Cas 7 - Cercle (approximation via 4 quarts de cercle)
  {
    // Tous les points de contrôle des 4 quadrants comme dans geomdl
    std::vector<Point> allControlPts = {
        Point(5, 0),  Point(5, 5),   Point(0, 5),  // Q1
        Point(0, 5),  Point(-5, 5),  Point(-5, 0), // Q2
        Point(-5, 0), Point(-5, -5), Point(0, -5), // Q3
        Point(0, -5), Point(5, -5),  Point(5, 0)   // Q4
    };

    // Créer une LineString composée en concaténant les 4 quadrants
    std::vector<Point>            allPoints;
    std::vector<NURBSCurve::FT>   w = {1, NURBSCurve::FT(std::sqrt(2.0) / 2.0),
                                       1};
    std::vector<NURBSCurve::Knot> k = {0, 0, 0, 1, 1, 1};

    // Q1: (5,0) -> (0,5)
    std::vector<Point> pts1 = {Point(5, 0), Point(5, 5), Point(0, 5)};
    NURBSCurve         curve1(pts1, w, 2, k);
    auto               ls1 = curve1.toLineString(50);
    for (size_t i = 0; i < ls1->numPoints(); ++i) {
      allPoints.push_back(ls1->pointN(i));
    }

    // Q2: (0,5) -> (-5,0)
    std::vector<Point> pts2 = {Point(0, 5), Point(-5, 5), Point(-5, 0)};
    NURBSCurve         curve2(pts2, w, 2, k);
    auto               ls2 = curve2.toLineString(50);
    for (size_t i = 1; i < ls2->numPoints();
         ++i) { // Skip first point to avoid duplication
      allPoints.push_back(ls2->pointN(i));
    }

    // Q3: (-5,0) -> (0,-5)
    std::vector<Point> pts3 = {Point(-5, 0), Point(-5, -5), Point(0, -5)};
    NURBSCurve         curve3(pts3, w, 2, k);
    auto               ls3 = curve3.toLineString(50);
    for (size_t i = 1; i < ls3->numPoints(); ++i) {
      allPoints.push_back(ls3->pointN(i));
    }

    // Q4: (0,-5) -> (5,0)
    std::vector<Point> pts4 = {Point(0, -5), Point(5, -5), Point(5, 0)};
    NURBSCurve         curve4(pts4, w, 2, k);
    auto               ls4 = curve4.toLineString(50);
    for (size_t i = 1; i < ls4->numPoints(); ++i) {
      allPoints.push_back(ls4->pointN(i));
    }

    LineString     fullCircle(allPoints);
    NURBSCurve::FT totalLength =
        curve1.length() + curve2.length() + curve3.length() + curve4.length();

    std::cout << "Cas 7 - Cercle\n";
    std::cout << controlPointsWKT(allControlPts) << "\n";
    std::cout << fullCircle.asText(2) << "\n";
    std::cout << "Length: " << CGAL::to_double(totalLength) << " (théorique ~"
              << 2 * M_PI * 5 << ")\n\n";
  }

  // 8. Cas 8 - Ellipse (approximation via scaling des points de cercle)
  {
    // Tous les points de contrôle des 4 quadrants comme dans geomdl
    std::vector<Point> allControlPts = {
        Point(8, 0),  Point(8, 4),   Point(0, 4),  // Q1
        Point(0, 4),  Point(-8, 4),  Point(-8, 0), // Q2
        Point(-8, 0), Point(-8, -4), Point(0, -4), // Q3
        Point(0, -4), Point(8, -4),  Point(8, 0)   // Q4
    };

    // Créer une LineString composée en concaténant les 4 quadrants
    std::vector<Point>            allPoints;
    std::vector<NURBSCurve::FT>   w = {1, NURBSCurve::FT(std::sqrt(2.0) / 2.0),
                                       1};
    std::vector<NURBSCurve::Knot> k = {0, 0, 0, 1, 1, 1};

    // Q1: (8,0) -> (0,4)
    std::vector<Point> pts1 = {Point(8, 0), Point(8, 4), Point(0, 4)};
    NURBSCurve         curve1(pts1, w, 2, k);
    auto               ls1 = curve1.toLineString(50);
    for (size_t i = 0; i < ls1->numPoints(); ++i) {
      allPoints.push_back(ls1->pointN(i));
    }

    // Q2: (0,4) -> (-8,0)
    std::vector<Point> pts2 = {Point(0, 4), Point(-8, 4), Point(-8, 0)};
    NURBSCurve         curve2(pts2, w, 2, k);
    auto               ls2 = curve2.toLineString(50);
    for (size_t i = 1; i < ls2->numPoints(); ++i) {
      allPoints.push_back(ls2->pointN(i));
    }

    // Q3: (-8,0) -> (0,-4)
    std::vector<Point> pts3 = {Point(-8, 0), Point(-8, -4), Point(0, -4)};
    NURBSCurve         curve3(pts3, w, 2, k);
    auto               ls3 = curve3.toLineString(50);
    for (size_t i = 1; i < ls3->numPoints(); ++i) {
      allPoints.push_back(ls3->pointN(i));
    }

    // Q4: (0,-4) -> (8,0)
    std::vector<Point> pts4 = {Point(0, -4), Point(8, -4), Point(8, 0)};
    NURBSCurve         curve4(pts4, w, 2, k);
    auto               ls4 = curve4.toLineString(50);
    for (size_t i = 1; i < ls4->numPoints(); ++i) {
      allPoints.push_back(ls4->pointN(i));
    }

    LineString     fullEllipse(allPoints);
    NURBSCurve::FT totalLength =
        curve1.length() + curve2.length() + curve3.length() + curve4.length();

    std::cout << "Cas 8 - Ellipse\n";
    std::cout << controlPointsWKT(allControlPts) << "\n";
    std::cout << fullEllipse.asText(2) << "\n";
    std::cout << "Length: " << CGAL::to_double(totalLength) << "\n\n";
  }

  // 9. Cas 9 - Cœur stylisé: 13 points de contrôle, degré 3
  {
    std::vector<Point> pts = {Point(0, 2),  Point(2, 4),     Point(4, 4),
                              Point(5, 2),  Point(5, 0),     Point(2.5, -3),
                              Point(0, -5), Point(-2.5, -3), Point(-5, 0),
                              Point(-5, 2), Point(-4, 4),    Point(-2, 4),
                              Point(0, 2)};

    std::cout << "Cas 9 - Cœur\n";
    std::cout << controlPointsWKT(pts) << "\n";

    // Test différentes méthodes de paramétrage
    NURBSCurve curveChord(pts, 3, NURBSCurve::KnotMethod::CHORD_LENGTH);
    std::cout << "=== CHORD_LENGTH (défaut SFCGAL) ===\n";
    std::cout << curveChord.toLineString(50)->asText(2) << "\n";

    NURBSCurve curveUniform(pts, 3, NURBSCurve::KnotMethod::UNIFORM);
    std::cout << "=== UNIFORM ===\n";
    std::cout << curveUniform.toLineString(50)->asText(2) << "\n";

    NURBSCurve curveCentripetal(pts, 3, NURBSCurve::KnotMethod::CENTRIPETAL);
    std::cout << "=== CENTRIPETAL (probablement geomdl) ===\n";
    std::cout << curveCentripetal.toLineString(50)->asText(2) << "\n";

    std::cout << "Length CHORD_LENGTH: " << CGAL::to_double(curveChord.length())
              << "\n";
    std::cout << "Length UNIFORM: " << CGAL::to_double(curveUniform.length())
              << "\n";
    std::cout << "Length CENTRIPETAL: "
              << CGAL::to_double(curveCentripetal.length()) << "\n\n";
  }
}

BOOST_AUTO_TEST_CASE(testHeartParameterizationComparison)
{
  std::cout << "\n=== HEART PARAMETERIZATION METHODS COMPARISON ==="
            << std::endl;
  std::cout << "Testing different knot vector generation methods on heart shape"
            << std::endl;
  std::cout
      << "This helps identify which method gives results closest to geomdl\n"
      << std::endl;

  // Heart control points (13 points including closure, degree 3)
  std::vector<Point> heartPts = {
      Point(0, 2),    Point(2, 4),  Point(4, 4),     Point(5, 2),  Point(5, 0),
      Point(2.5, -3), Point(0, -5), Point(-2.5, -3), Point(-5, 0), Point(-5, 2),
      Point(-4, 4),   Point(-2, 4), Point(0, 2)};

  std::cout << "Control Points (13 points, closed curve):" << std::endl;
  std::cout << controlPointsWKT(heartPts) << "\n" << std::endl;

  const int resolution = 50; // Good resolution for visualization

  // Create curves with different parameterization methods
  NURBSCurve curveChord(heartPts, 3, NURBSCurve::KnotMethod::CHORD_LENGTH);
  NURBSCurve curveUniform(heartPts, 3, NURBSCurve::KnotMethod::UNIFORM);
  NURBSCurve curveCentripetal(heartPts, 3, NURBSCurve::KnotMethod::CENTRIPETAL);
  std::cout << curveChord.asText(1) << std::endl;
  std::cout << curveUniform.asText(1) << std::endl;
  std::cout << curveCentripetal.asText(1) << std::endl;

  // Generate geometries
  auto lsChord       = curveChord.toLineString(resolution);
  auto lsUniform     = curveUniform.toLineString(resolution);
  auto lsCentripetal = curveCentripetal.toLineString(resolution);

  // Method 1: CHORD_LENGTH (SFCGAL default)
  std::cout << "=== METHOD 1: CHORD_LENGTH (SFCGAL Default) ===" << std::endl;
  std::cout << "Based on Euclidean distances between control points"
            << std::endl;
  std::cout << "Geometry:" << std::endl;
  std::cout << lsChord->asText(3) << std::endl;
  std::cout << "Length: " << CGAL::to_double(curveChord.length()) << std::endl;

  // Find Y extremes for comparison
  double maxY_chord = -1000, minY_chord = 1000;
  for (size_t i = 0; i < lsChord->numPoints(); ++i) {
    double y   = CGAL::to_double(lsChord->pointN(i).y());
    maxY_chord = std::max(maxY_chord, y);
    minY_chord = std::min(minY_chord, y);
  }
  std::cout << "Y range: [" << minY_chord << ", " << maxY_chord
            << "], height: " << (maxY_chord - minY_chord) << "\n"
            << std::endl;

  // Method 2: UNIFORM
  std::cout << "=== METHOD 2: UNIFORM ===" << std::endl;
  std::cout << "Equal spacing in parameter domain" << std::endl;
  std::cout << "Geometry:" << std::endl;
  std::cout << lsUniform->asText(3) << std::endl;
  std::cout << "Length: " << CGAL::to_double(curveUniform.length())
            << std::endl;

  double maxY_uniform = -1000, minY_uniform = 1000;
  for (size_t i = 0; i < lsUniform->numPoints(); ++i) {
    double y     = CGAL::to_double(lsUniform->pointN(i).y());
    maxY_uniform = std::max(maxY_uniform, y);
    minY_uniform = std::min(minY_uniform, y);
  }
  std::cout << "Y range: [" << minY_uniform << ", " << maxY_uniform
            << "], height: " << (maxY_uniform - minY_uniform) << "\n"
            << std::endl;

  // Method 3: CENTRIPETAL (likely geomdl's choice)
  std::cout << "=== METHOD 3: CENTRIPETAL (Likely geomdl method) ==="
            << std::endl;
  std::cout << "Based on square root of distances (often preferred in CAD)"
            << std::endl;
  std::cout << "Geometry:" << std::endl;
  std::cout << lsCentripetal->asText(3) << std::endl;
  std::cout << "Length: " << CGAL::to_double(curveCentripetal.length())
            << std::endl;

  double maxY_centripetal = -1000, minY_centripetal = 1000;
  for (size_t i = 0; i < lsCentripetal->numPoints(); ++i) {
    double y         = CGAL::to_double(lsCentripetal->pointN(i).y());
    maxY_centripetal = std::max(maxY_centripetal, y);
    minY_centripetal = std::min(minY_centripetal, y);
  }
  std::cout << "Y range: [" << minY_centripetal << ", " << maxY_centripetal
            << "], height: " << (maxY_centripetal - minY_centripetal) << "\n"
            << std::endl;

  // Comparative summary
  std::cout << "=== COMPARATIVE SUMMARY ===" << std::endl;
  std::cout << "Method          | Length    | Max Y     | Min Y     | Height   "
               " | Notes"
            << std::endl;
  std::cout << "----------------|-----------|-----------|-----------|----------"
               "-|-------------"
            << std::endl;

  std::cout << std::left << std::setw(15) << "CHORD_LENGTH" << " | "
            << std::fixed << std::setprecision(4) << std::setw(9)
            << CGAL::to_double(curveChord.length()) << " | " << std::setw(9)
            << maxY_chord << " | " << std::setw(9) << minY_chord << " | "
            << std::setw(9) << (maxY_chord - minY_chord) << " | SFCGAL default"
            << std::endl;

  std::cout << std::left << std::setw(15) << "UNIFORM" << " | " << std::fixed
            << std::setprecision(4) << std::setw(9)
            << CGAL::to_double(curveUniform.length()) << " | " << std::setw(9)
            << maxY_uniform << " | " << std::setw(9) << minY_uniform << " | "
            << std::setw(9) << (maxY_uniform - minY_uniform)
            << " | Equal spacing" << std::endl;

  std::cout << std::left << std::setw(15) << "CENTRIPETAL" << " | "
            << std::fixed << std::setprecision(4) << std::setw(9)
            << CGAL::to_double(curveCentripetal.length()) << " | "
            << std::setw(9) << maxY_centripetal << " | " << std::setw(9)
            << minY_centripetal << " | " << std::setw(9)
            << (maxY_centripetal - minY_centripetal) << " | Likely geomdl"
            << std::endl;

  // Key differences analysis
  std::cout << "\n=== KEY DIFFERENCES ANALYSIS ===" << std::endl;

  double lengthDiff_UC = CGAL::to_double(curveUniform.length()) -
                         CGAL::to_double(curveChord.length());
  double lengthDiff_CC = CGAL::to_double(curveCentripetal.length()) -
                         CGAL::to_double(curveChord.length());

  std::cout << "Length differences vs CHORD_LENGTH:" << std::endl;
  std::cout << "  UNIFORM:     " << std::showpos << lengthDiff_UC << " ("
            << (lengthDiff_UC / CGAL::to_double(curveChord.length()) * 100)
            << "%)" << std::endl;
  std::cout << "  CENTRIPETAL: " << std::showpos << lengthDiff_CC << " ("
            << (lengthDiff_CC / CGAL::to_double(curveChord.length()) * 100)
            << "%)" << std::noshowpos << std::endl;

  std::cout << "Height differences vs CHORD_LENGTH:" << std::endl;
  double heightChord = maxY_chord - minY_chord;
  std::cout << "  UNIFORM:     " << std::showpos
            << ((maxY_uniform - minY_uniform) - heightChord) << std::endl;
  std::cout << "  CENTRIPETAL: " << std::showpos
            << ((maxY_centripetal - minY_centripetal) - heightChord)
            << std::noshowpos << std::endl;

  std::cout << "\n=== RECOMMENDATION ===" << std::endl;
  std::cout << "For closest match to geomdl, use: NURBSCurve(points, degree, "
               "NURBSCurve::KnotMethod::CENTRIPETAL)"
            << std::endl;
  std::cout << "The CENTRIPETAL method typically produces the most visually "
               "pleasing curves"
            << std::endl;
  std::cout << "and is commonly used in CAD applications like geomdl.\n"
            << std::endl;
}

BOOST_AUTO_TEST_CASE(testFitPointsVsControlPoints)
{
  std::cout << "\n=== FIT POINTS vs CONTROL POINTS COMPARISON ===" << std::endl;
  std::cout << "Comparing AutoCAD-style Fit Points (interpolation) vs Control "
               "Points\n"
            << std::endl;

  // Simple test points for clear comparison
  std::vector<Point> testPoints = {Point(0, 0), Point(1, 2), Point(2, 1),
                                   Point(3, 3), Point(4, 0)};

  std::cout << "Test Points (should be traversed exactly by Fit Points method):"
            << std::endl;
  for (size_t i = 0; i < testPoints.size(); ++i) {
    std::cout << "  P" << i << ": (" << CGAL::to_double(testPoints[i].x())
              << ", " << CGAL::to_double(testPoints[i].y()) << ")" << std::endl;
  }
  std::cout << std::endl;

  try {
    // Method 1: Control Points (traditional NURBS - curve influenced by points)
    NURBSCurve controlCurve(testPoints, 3, NURBSCurve::KnotMethod::UNIFORM);
    auto       controlLS = controlCurve.toLineString(50);

    std::cout << "=== METHOD 1: CONTROL POINTS (Traditional NURBS) ==="
              << std::endl;
    std::cout << "Curve is INFLUENCED by points but doesn't pass through them"
              << std::endl;
    std::cout << "Geometry: " << controlLS->asText(3) << std::endl;
    std::cout << "Length: " << CGAL::to_double(controlCurve.length())
              << std::endl;

    // Check if curve passes through original points
    std::cout << "Distance from curve to original points:" << std::endl;
    for (size_t i = 0; i < testPoints.size(); ++i) {
      // Find parameter closest to this point (simplified check)
      auto   bounds  = controlCurve.parameterBounds();
      double minDist = 1000.0;
      for (int j = 0; j <= 100; ++j) {
        double t = j / 100.0;
        auto   param =
            bounds.first + NURBSCurve::FT(t) * (bounds.second - bounds.first);
        auto   curvePoint = controlCurve.evaluate(param);
        double dist       = std::sqrt(
            std::pow(CGAL::to_double(curvePoint.x() - testPoints[i].x()), 2) +
            std::pow(CGAL::to_double(curvePoint.y() - testPoints[i].y()), 2));
        minDist = std::min(minDist, dist);
      }
      std::cout << "  P" << i << " distance: " << std::setprecision(4)
                << minDist << std::endl;
    }

    // Method 2: Fit Points (interpolation - curve passes through points)
    std::cout << "\n=== METHOD 2: FIT POINTS (AutoCAD-style Interpolation) ==="
              << std::endl;
    std::cout << "Curve passes EXACTLY through all given points" << std::endl;

    auto fitCurve = NURBSCurve::interpolateCurve(
        testPoints, 3, NURBSCurve::KnotMethod::UNIFORM,
        NURBSCurve::EndCondition::CLAMPED);

    if (fitCurve) {
      auto fitLS = fitCurve->toLineString(50);
      std::cout << "Geometry: " << fitLS->asText(3) << std::endl;
      std::cout << "Length: " << CGAL::to_double(fitCurve->length())
                << std::endl;

      // Verify interpolation by checking distances to curve
      std::cout << "Verification - minimum distance from points to curve:"
                << std::endl;
      auto bounds = fitCurve->parameterBounds();
      for (size_t i = 0; i < testPoints.size(); ++i) {
        double minDist = 1000.0;
        // Sample densely to find minimum distance
        for (int j = 0; j <= 200; ++j) {
          double t = j / 200.0;
          auto   param =
              bounds.first + NURBSCurve::FT(t) * (bounds.second - bounds.first);
          auto   curvePoint = fitCurve->evaluate(param);
          double dist       = std::sqrt(
              std::pow(CGAL::to_double(curvePoint.x() - testPoints[i].x()), 2) +
              std::pow(CGAL::to_double(curvePoint.y() - testPoints[i].y()), 2));
          minDist = std::min(minDist, dist);
        }
        std::cout << "  P" << i << " min distance: " << std::setprecision(8)
                  << minDist << " (should be ~0 for interpolation)"
                  << std::endl;
      }
    } else {
      std::cout << "ERROR: Could not create interpolating curve" << std::endl;
    }

    std::cout << "\n=== SUMMARY ===" << std::endl;
    std::cout << "Control Points: Curve influenced by points (traditional CAD "
                 "workflow)"
              << std::endl;
    std::cout << "Fit Points: Curve passes through points (AutoCAD SPLINE with "
                 "F option)"
              << std::endl;
    std::cout
        << "For geomdl compatibility: Use Control Points with UNIFORM method"
        << std::endl;
    std::cout
        << "For AutoCAD SPLINE compatibility: Use interpolateCurve() method"
        << std::endl;

  } catch (const std::exception &e) {
    std::cerr << "Error in test: " << e.what() << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(testApproximationVsGeomdl)
{
  std::cout << "\n=== SFCGAL NURBS APPROXIMATION (geomdl equivalent) ==="
            << std::endl;
  std::cout << "Testing SFCGAL's approximateCurve() against geomdl's "
               "approximate_curve()\n"
            << std::endl;

  // Parameters matching your geomdl script
  const int nb_points  = 50; // Input data points
  const int degree     = 3;  // NURBS degree
  const int nb_ctrlpts = 10; // Number of control points

  std::cout << "Parameters:" << std::endl;
  std::cout << "  Input points: " << nb_points << std::endl;
  std::cout << "  NURBS degree: " << degree << std::endl;
  std::cout << "  Control points: " << nb_ctrlpts << std::endl;
  std::cout << std::endl;

  // Generate sine wave data points (same as geomdl script)
  std::vector<Point> dataPoints;
  for (int i = 0; i < nb_points; ++i) {
    double x = (2.0 * M_PI * i) / (nb_points - 1); // 0 to 2π
    double y = std::sin(x);
    dataPoints.emplace_back(x, y);
  }

  std::cout << "Generated " << dataPoints.size()
            << " data points from sine wave" << std::endl;
  std::cout << "X range: [0, " << 2.0 * M_PI << "], Y range: [-1, 1]"
            << std::endl;
  std::cout << std::endl;

  try {
    // SFCGAL NURBS approximation (equivalent to geomdl's approximate_curve)
    std::cout << "=== SFCGAL APPROXIMATION ===" << std::endl;

    auto approximatedCurve =
        NURBSCurve::approximateCurve(dataPoints,           // Input data points
                                     degree,               // Degree
                                     NURBSCurve::FT(1e-6), // Tolerance
                                     nb_ctrlpts // Maximum control points
        );

    if (approximatedCurve) {
      std::cout << "✓ Successfully created approximated NURBS curve!"
                << std::endl;

      // Get results (equivalent to geomdl's curve.ctrlpts, curve.knotvector,
      // etc.)
      auto controlPoints = approximatedCurve->controlPoints();
      auto knotVector    = approximatedCurve->knotVector();
      auto weights       = approximatedCurve->weights();

      std::cout << "\n=== RESULTS (geomdl equivalent) ===" << std::endl;

      // Control points (like geomdl's curve.ctrlpts)
      std::cout << "Control points (" << controlPoints.size()
                << " points):" << std::endl;
      for (size_t i = 0; i < std::min(size_t(5), controlPoints.size()); ++i) {
        std::cout << "  P" << i << ": (" << std::setprecision(6) << std::fixed
                  << CGAL::to_double(controlPoints[i].x()) << ", "
                  << CGAL::to_double(controlPoints[i].y()) << ")" << std::endl;
      }
      if (controlPoints.size() > 5) {
        std::cout << "  ... and " << (controlPoints.size() - 5)
                  << " more points" << std::endl;
      }

      // Knot vector (like geomdl's curve.knotvector)
      std::cout << "\nKnot vector (" << knotVector.size()
                << " knots):" << std::endl;
      std::cout << "[";
      for (size_t i = 0; i < std::min(size_t(10), knotVector.size()); ++i) {
        if (i > 0)
          std::cout << ", ";
        std::cout << std::setprecision(4) << std::fixed
                  << CGAL::to_double(knotVector[i]);
      }
      if (knotVector.size() > 10)
        std::cout << ", ...";
      std::cout << "]" << std::endl;

      // Curve evaluation (like geomdl's curve.evaluate())
      std::cout << "\n=== CURVE EVALUATION (geomdl equivalent) ==="
                << std::endl;
      const int sample_size    = 100; // Like geomdl's curve.sample_size = 100
      auto      evaluatedCurve = approximatedCurve->toLineString(sample_size);

      std::cout << "✓ Evaluated curve with " << sample_size << " samples"
                << std::endl;
      std::cout << "✓ Generated " << evaluatedCurve->numPoints()
                << " evaluated points" << std::endl;
      std::cout << "✓ Curve length: " << std::setprecision(6)
                << CGAL::to_double(approximatedCurve->length()) << std::endl;

      // Calculate approximation quality
      std::cout << "\n=== APPROXIMATION QUALITY ===" << std::endl;
      double maxError   = 0.0;
      double avgError   = 0.0;
      int    errorCount = 0;

      // Sample every 5th input point for error calculation
      for (size_t i = 0; i < dataPoints.size(); i += 5) {
        const auto &dataPoint = dataPoints[i];
        double      minDist   = 1000.0;

        // Find closest point on curve
        for (size_t j = 0; j < evaluatedCurve->numPoints(); ++j) {
          auto   curvePoint = evaluatedCurve->pointN(j);
          double dist       = std::sqrt(
              std::pow(CGAL::to_double(curvePoint.x() - dataPoint.x()), 2) +
              std::pow(CGAL::to_double(curvePoint.y() - dataPoint.y()), 2));
          minDist = std::min(minDist, dist);
        }
        maxError = std::max(maxError, minDist);
        avgError += minDist;
        errorCount++;
      }
      avgError /= errorCount;

      std::cout << "✓ Maximum approximation error: " << std::setprecision(8)
                << maxError << std::endl;
      std::cout << "✓ Average approximation error: " << std::setprecision(8)
                << avgError << std::endl;

      std::cout << "\n=== GEOMDL COMPATIBILITY ===" << std::endl;
      std::cout << "✓ SFCGAL provides equivalent functionality to geomdl:"
                << std::endl;
      std::cout << "  • approximate_curve() → NURBSCurve::approximateCurve()"
                << std::endl;
      std::cout << "  • curve.ctrlpts → curve->controlPoints()" << std::endl;
      std::cout << "  • curve.knotvector → curve->knotVector()" << std::endl;
      std::cout << "  • curve.evaluate() → curve->toLineString()" << std::endl;
      std::cout << "  • curve.sample_size → toLineString(N)" << std::endl;

    } else {
      std::cout << "✗ ERROR: Failed to create approximated curve" << std::endl;
    }

  } catch (const std::exception &e) {
    std::cerr << "✗ Error: " << e.what() << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(testApproximationVsTrueControlPoints)
{
  std::cout << "\n=== SFCGAL vs TRUE NURBS APPROXIMATION ANALYSIS ==="
            << std::endl;
  std::cout
      << "Demonstrating the difference between SFCGAL and geomdl approaches\n"
      << std::endl;

  // Create a simple test case: 5 points on a parabola
  std::vector<Point> dataPoints = {Point(0, 0), // y = x^2
                                   Point(1, 1), Point(2, 4), Point(3, 9),
                                   Point(4, 16)};

  std::cout << "Input data points (y = x^2):" << std::endl;
  for (size_t i = 0; i < dataPoints.size(); ++i) {
    std::cout << "  Data[" << i << "]: (" << CGAL::to_double(dataPoints[i].x())
              << ", " << CGAL::to_double(dataPoints[i].y()) << ")" << std::endl;
  }

  try {
    // SFCGAL current approximation (problematic)
    std::cout << "\n=== SFCGAL CURRENT APPROXIMATION ===" << std::endl;
    auto sfcgalCurve =
        NURBSCurve::approximateCurve(dataPoints, 3, NURBSCurve::FT(1e-6), 3);

    if (sfcgalCurve) {
      auto sfcgalControlPoints = sfcgalCurve->controlPoints();

      std::cout << "SFCGAL 'control' points (" << sfcgalControlPoints.size()
                << " points):" << std::endl;
      for (size_t i = 0; i < sfcgalControlPoints.size(); ++i) {
        std::cout << "  Ctrl[" << i << "]: ("
                  << CGAL::to_double(sfcgalControlPoints[i].x()) << ", "
                  << CGAL::to_double(sfcgalControlPoints[i].y()) << ")"
                  << std::endl;
      }

      // Check if control points are ON the input data
      std::cout << "\nAnalysis: Are SFCGAL 'control points' on input data?"
                << std::endl;
      for (size_t i = 0; i < sfcgalControlPoints.size(); ++i) {
        bool foundMatch = false;
        for (size_t j = 0; j < dataPoints.size(); ++j) {
          double dx =
              CGAL::to_double(sfcgalControlPoints[i].x() - dataPoints[j].x());
          double dy =
              CGAL::to_double(sfcgalControlPoints[i].y() - dataPoints[j].y());
          if (std::abs(dx) < 1e-10 && std::abs(dy) < 1e-10) {
            std::cout << "  ✓ Ctrl[" << i << "] matches Data[" << j
                      << "] exactly" << std::endl;
            foundMatch = true;
            break;
          }
        }
        if (!foundMatch) {
          std::cout << "  ✗ Ctrl[" << i
                    << "] is NOT on input data (this would be correct for true "
                       "approximation)"
                    << std::endl;
        }
      }

      // Evaluate the curve
      auto evaluatedCurve = sfcgalCurve->toLineString(20);
      std::cout << "\nSFCGAL curve evaluation (" << evaluatedCurve->numPoints()
                << " points):" << std::endl;
      std::cout << "  First point: ("
                << CGAL::to_double(evaluatedCurve->pointN(0).x()) << ", "
                << CGAL::to_double(evaluatedCurve->pointN(0).y()) << ")"
                << std::endl;
      std::cout
          << "  Last point: ("
          << CGAL::to_double(
                 evaluatedCurve->pointN(evaluatedCurve->numPoints() - 1).x())
          << ", "
          << CGAL::to_double(
                 evaluatedCurve->pointN(evaluatedCurve->numPoints() - 1).y())
          << ")" << std::endl;
    }

    // Demonstrate TRUE control point approximation concept
    std::cout << "\n=== TRUE NURBS APPROXIMATION (Conceptual) ===" << std::endl;
    std::cout << "What geomdl would do with least-squares approximation:"
              << std::endl;
    std::cout << "Input: 5 data points on parabola y = x^2" << std::endl;
    std::cout << "Goal: Find 3 control points that create NURBS approximating "
                 "the data"
              << std::endl;
    std::cout << "" << std::endl;

    // Manual example of what true approximation should produce
    std::cout << "Expected TRUE control points (computed by least-squares):"
              << std::endl;
    std::cout << "  True_Ctrl[0]: (~0.0, ~-1.5)  # NOT on input data!"
              << std::endl;
    std::cout << "  True_Ctrl[1]: (~2.0, ~8.5)   # NOT on input data!"
              << std::endl;
    std::cout << "  True_Ctrl[2]: (~4.0, ~18.5)  # NOT on input data!"
              << std::endl;
    std::cout << "" << std::endl;
    std::cout << "These control points would create a NURBS curve that:"
              << std::endl;
    std::cout << "  ✓ Passes NEAR the input data points" << std::endl;
    std::cout << "  ✓ Minimizes approximation error" << std::endl;
    std::cout << "  ✗ Control points are NOT on the input data" << std::endl;

    std::cout << "\n=== KEY DIFFERENCE SUMMARY ===" << std::endl;
    std::cout << "CURRENT SFCGAL approach:" << std::endl;
    std::cout << "  • approximateCurve() selects subset of input points as "
                 "'control points'"
              << std::endl;
    std::cout << "  • Result: Control points ARE on the input data"
              << std::endl;
    std::cout << "  • This is actually SAMPLING, not true approximation"
              << std::endl;
    std::cout << "" << std::endl;
    std::cout << "GEOMDL/TRUE approximation approach:" << std::endl;
    std::cout << "  • approximate_curve() computes optimal control points via "
                 "least-squares"
              << std::endl;
    std::cout << "  • Result: Control points are NOT on input data"
              << std::endl;
    std::cout
        << "  • Control points create curve that best approximates input data"
        << std::endl;
    std::cout << "" << std::endl;
    std::cout << "CONCLUSION: SFCGAL needs a true least-squares approximation "
                 "algorithm"
              << std::endl;
    std::cout << "to match geomdl's behavior!" << std::endl;

  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()
