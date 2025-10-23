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
#include "SFCGAL/io/wkb.h"
#include "SFCGAL/io/wkt.h"
#include "SFCGAL/numeric.h"
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Aff_transformation_3.h>
#include <cmath>
#include <iomanip>
#include <memory>
#include <sstream>

using namespace SFCGAL;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_NURBSCurveTest)

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
isNearlyEqual(double valueA, double valueB, double tolerance = EPSILON) -> bool
{
  return std::abs(valueA - valueB) < tolerance;
}

auto
isNearlyEqual(const Point &firstPoint, const Point &secondPoint,
              double tolerance = EPSILON) -> bool
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

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
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
// NOLINTEND(bugprone-easily-swappable-parameters)

auto
checkControlPointsEqual(const NURBSCurve &curve1, const NURBSCurve &curve2,
                        double tolerance = EPSILON) -> void
{
  BOOST_REQUIRE_EQUAL(curve1.numControlPoints(), curve2.numControlPoints());

  for (size_t idx = 0; idx < curve1.numControlPoints(); ++idx) {
    const Point &point1 = curve1.controlPointN(idx);
    const Point &point2 = curve2.controlPointN(idx);

    BOOST_CHECK(isNearlyEqual(point1, point2, tolerance));
  }
}

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
  std::unique_ptr<NURBSCurve> cloned = curve.clone();

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
  BOOST_CHECK(CGAL::abs(derivative.x()) > NURBSCurve::FT(EPSILON) ||
              CGAL::abs(derivative.y()) > NURBSCurve::FT(EPSILON));
}

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

  BOOST_CHECK_CLOSE(CGAL::to_double(start.x()), cos(0.0),
                    EPSILON); // cos(0) = 1
  BOOST_CHECK_CLOSE(CGAL::to_double(start.y()), sin(0.0),
                    EPSILON); // sin(0) = 0
  BOOST_CHECK_CLOSE(CGAL::to_double(end.x()), cos(M_PI),
                    EPSILON); // cos(π) = -1
  BOOST_CHECK_CLOSE(CGAL::to_double(end.y()), sin(M_PI), EPSILON); // sin(π) ≈ 0
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
                    EPSILON); // 2 * cos(0) = 2
  BOOST_CHECK_CLOSE(CGAL::to_double(start.y()), 2.0 * sin(0.0),
                    EPSILON);                        // 2 * sin(0) = 0
  BOOST_CHECK_EQUAL(start.z(), NURBSCurve::FT(5.0)); // Z coordinate preserved
  BOOST_CHECK_CLOSE(CGAL::to_double(end.x()), 2.0 * cos(M_PI / 2.0),
                    EPSILON); // 2 * cos(π/2) ≈ 0
  BOOST_CHECK_CLOSE(CGAL::to_double(end.y()), 2.0 * sin(M_PI / 2.0),
                    EPSILON);                      // 2 * sin(π/2) = 2
  BOOST_CHECK_EQUAL(end.z(), NURBSCurve::FT(5.0)); // Z coordinate preserved
}

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

BOOST_AUTO_TEST_CASE(testWktEmpty)
{
  NURBSCurve const emptyCurve;
  std::string      wkt = emptyCurve.asText(0);
  BOOST_CHECK(wkt.find("EMPTY") != std::string::npos);
}

/// WKT reading of basic NURBS curve (points + degree only)
BOOST_AUTO_TEST_CASE(testReadBasicNURBSWkt)
{
  std::string const wkt  = "NURBSCURVE(2, (0 0, 5 10, 10 0))";
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
  std::string const wkt  = "NURBSCURVE(2, (0 0, 5 10, 10 0), (1, 2, 1))";
  auto              geom = SFCGAL::io::readWkt(wkt);

  BOOST_REQUIRE(geom->is<NURBSCurve>());
  auto const &curve = geom->as<NURBSCurve>();
  BOOST_CHECK(curve.isRational());
  BOOST_CHECK_EQUAL(curve.weight(1), NURBSCurve::FT(2.0));
}

/// WKT reading with complete knot vector
BOOST_AUTO_TEST_CASE(testReadFullNURBSWkt)
{
  // 4 points, degree 2 → needs 4+2+1=7 knots (corrected)
  std::string const wkt =
      "NURBSCURVE(2, (0 0, 3 6, 6 3, 9 0), (1, 1, 1, 1), (0, "
      "0, 0, 0.5, 1, 1, 1))";
  auto geom = SFCGAL::io::readWkt(wkt);

  BOOST_REQUIRE(geom->is<NURBSCurve>());
  auto const &curve = geom->as<NURBSCurve>();
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 4U);
  BOOST_CHECK_EQUAL(curve.degree(), 2U);
  BOOST_CHECK_EQUAL(curve.knotVector().size(), 7U);
}

/// WKT reading of 3D NURBS curve
BOOST_AUTO_TEST_CASE(testRead3DNURBSWkt)
{
  std::string const wkt  = "NURBSCURVE Z (2, (0 0 0, 5 5 10, 10 0 0))";
  auto              geom = SFCGAL::io::readWkt(wkt);

  BOOST_REQUIRE(geom->is<NURBSCurve>());
  auto const &curve = geom->as<NURBSCurve>();
  BOOST_CHECK(curve.is3D());
  BOOST_CHECK_EQUAL(curve.controlPointN(1).z(), NURBSCurve::FT(10.0));
}

/// WKT reading with measured coordinates
BOOST_AUTO_TEST_CASE(testReadMeasuredNURBSWkt)
{
  std::string const wkt  = "NURBSCURVE M (2, (0 0 5, 5 5 15, 10 0 25))";
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

  BOOST_CHECK_EQUAL(result,
                    "NURBSCURVE (2,(0 0,5 10,10 0),(1,1,1),(0,0,0,1,1,1))");
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

  BOOST_CHECK_EQUAL(result,
                    "NURBSCURVE (2,(0 0,5 10,10 0),(1,2,1),(0,0,0,1,1,1))");
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

  BOOST_CHECK_EQUAL(
      result, "NURBSCURVE Z (2,(0 0 0,5 5 10,10 0 0),(1,1,1),(0,0,0,1,1,1))");
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
  std::string const originalWkt = "NURBSCURVE(2, (0 0, 5 10, 10 0), (1, 2, 1))";
  auto              geom        = SFCGAL::io::readWkt(originalWkt);
  std::string       result      = geom->asText(0);

  BOOST_CHECK_EQUAL(result,
                    "NURBSCURVE (2,(0 0,5 10,10 0),(1,2,1),(0,0,0,1,1,1))");

  // Verify re-reading produces identical geometry
  auto geom2 = SFCGAL::io::readWkt(result);
  BOOST_REQUIRE(geom2->is<NURBSCurve>());
  auto const &curve2 = geom2->as<NURBSCurve>();
  BOOST_CHECK_EQUAL(curve2.weight(1), NURBSCurve::FT(2.0));
}

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
                      CGAL::to_double(original.x()) + 5.0, EPSILON);
    BOOST_CHECK_CLOSE(CGAL::to_double(transformed.y()),
                      CGAL::to_double(original.y()) + 3.0, EPSILON);
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
  BOOST_CHECK_SMALL(transformed.x(), original.x());
  BOOST_CHECK_SMALL(transformed.y(), original.y());
  BOOST_CHECK_SMALL(transformed.z(), original.z());

  // Check that weights are preserved
  for (size_t idx = 0; idx < originalCurve.numControlPoints(); ++idx) {
    BOOST_CHECK_EQUAL(originalCurve.weight(idx), transformedCurve.weight(idx));
  }
}

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
      isNearlyEqual(multiPoint.geometryN(0).as<Point>(), curveStart, EPSILON) ||
      isNearlyEqual(multiPoint.geometryN(0).as<Point>(), curveEnd, EPSILON));
  BOOST_CHECK(
      isNearlyEqual(multiPoint.geometryN(1).as<Point>(), curveStart, EPSILON) ||
      isNearlyEqual(multiPoint.geometryN(1).as<Point>(), curveEnd, EPSILON));

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
  auto validityTol = algorithm::isValid(validCurve, EPSILON);
  BOOST_CHECK(validityTol.valid());

  // Empty curve should be valid
  NURBSCurve emptyCurve;
  auto       emptyValidity = algorithm::isValid(emptyCurve);
  BOOST_CHECK(emptyValidity.valid());
}

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

  BOOST_CHECK(isNearlyEqual(firstEnd, secondStart, EPSILON));

  // Original curve should be reconstructable from parts
  Point originalMid = linearCurve.evaluate(midParam);
  BOOST_CHECK(isNearlyEqual(firstEnd, originalMid, EPSILON));
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

  BOOST_CHECK(isNearlyEqual(originalStart, reversedEnd, EPSILON));
  BOOST_CHECK(isNearlyEqual(originalEnd, reversedStart, EPSILON));
}

BOOST_AUTO_TEST_CASE(testLinearCurveArcLength)
{
  // Linear curve from (0,0) to (3,4)
  std::vector<Point> linearPoints;
  linearPoints.emplace_back(0.0, 0.0);
  linearPoints.emplace_back(3.0, 4.0);

  NURBSCurve linearCurve(linearPoints, 1);

  // Expected length = sqrt(3² + 4²) = 5
  NURBSCurve::FT totalLength = linearCurve.length();
  BOOST_CHECK_CLOSE(CGAL::to_double(totalLength), 5.0, EPSILON);

  // Test partial lengths
  auto           bounds   = linearCurve.parameterBounds();
  NURBSCurve::FT midParam = (bounds.first + bounds.second) / NURBSCurve::FT(2);

  NURBSCurve::FT halfLength = linearCurve.length(bounds.first, midParam);
  BOOST_CHECK_CLOSE(CGAL::to_double(halfLength), 2.5, EPSILON);

  NURBSCurve::FT quarterLength = linearCurve.length(
      bounds.first,
      bounds.first + (midParam - bounds.first) / NURBSCurve::FT(2));
  BOOST_CHECK_CLOSE(CGAL::to_double(quarterLength), 1.25, EPSILON);
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
                    EPSILON);

  // At half length
  NURBSCurve::Parameter midParam =
      linearCurve.parameterAtLength(NURBSCurve::FT(2.5));
  Point midPoint = linearCurve.evaluate(midParam);
  BOOST_CHECK_CLOSE(CGAL::to_double(midPoint.x()), 2.0, EPSILON);
  BOOST_CHECK_CLOSE(CGAL::to_double(midPoint.y()), 1.5, EPSILON);

  // At end
  NURBSCurve::FT        totalLength = linearCurve.length();
  NURBSCurve::Parameter endParam = linearCurve.parameterAtLength(totalLength);
  BOOST_CHECK_CLOSE(CGAL::to_double(endParam), CGAL::to_double(bounds.second),
                    EPSILON);
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

  BOOST_CHECK(isNearlyEqual(pointOnCurve, closestOnCurve, EPSILON));
  BOOST_CHECK_CLOSE(CGAL::to_double(resultParam), 0.5,
                    EPSILON); // Should be at midpoint

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
    if (isNearlyEqual(testPoint, closestOffCurve, EPSILON)) {
      foundPoint = true;
    }
  }
  BOOST_CHECK(foundPoint);
}

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
  std::string wkt = "NURBSCURVE(2, (0 0,1e2 1.5e-1,2e0 0),(1e0,2.5e0,1e0))";

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

// NOLINTBEGIN(bugprone-suspicious-missing-comma)
BOOST_AUTO_TEST_CASE(testInteroperability)
{
  std::vector<std::string> wkts = {
      "NURBSCURVE(1, (0 0, 10 10))",
      "NURBSCURVE(1, (0 0, 5 5, 10 0))",
      "NURBSCURVE(1, (0 0, 10 10), (1, 1))",
      "NURBSCURVE(1, (0 0, 10 10), (0.5, 2.0))",
      "NURBSCURVE(1, (0 0, 5 5, 10 0), (1, 1.5, 1))",
      "NURBSCURVE(1, (0 0, 10 10), (1, 1), (0, 0, 1, 1))",
      "NURBSCURVE(1, (0 0, 5 5, 10 0), (1, 1, 1), (0, 0, 0.5, 1, 1))",
      "NURBSCURVE(2, (0 0, 5 10, 10 0))",
      "NURBSCURVE(2, (1 0, 1 1, 0 1))",
      "NURBSCURVE(2, (-5 -5, 0 10, 5 -5))",
      "NURBSCURVE(2, (0 0, 2 8, 8 2, 10 10))",
      "NURBSCURVE(2, (0 0, 5 10, 10 0), (1, 1, 1))",
      "NURBSCURVE(2, (1 0, 1 1, 0 1), (1, 0.5, 1))",
      "NURBSCURVE(2, (0 0, 5 10, 10 0), (1, 3, 1))",
      "NURBSCURVE(2, (0 0, 5 10, 10 0), (0.5, 1, 0.5))",
      "NURBSCURVE(2, (0.5 0.25, 2.75 5.5, 5.0 0.25), (1, 1.5, 1))",
      "NURBSCURVE(2, (0 0, 5 10, 10 0), (1, 1, 1), (0, 0, 0, 1, 1, 1))",
      "NURBSCURVE(2, (1 0, 1 1, 0 1), (1, 0.5, 1), (0, 0, 0, 1, 1, 1))",
      "NURBSCURVE(3, (0 0, 3 10, 7 10, 10 0))",
      "NURBSCURVE(3, (0 0, 2 8, 5 12, 8 8, 10 0))",
      "NURBSCURVE(3, (5 0, 10 2.5, 10 7.5, 5 10, 0 7.5, 0 2.5, 5 0))",
      "NURBSCURVE(3, (0 0, 5 5, 10 0, 15 -5, 20 0))",
      "NURBSCURVE(3, (0 0, 3 10, 7 10, 10 0), (1, 1, 1, 1))",
      "NURBSCURVE(3, (0 0, 3 10, 7 10, 10 0), (1, 2, 2, 1))",
      "NURBSCURVE(3, (0 0, 3 10, 7 10, 10 0), (0.5, 1.5, 2.0, 0.8))",
      "NURBSCURVE(3, (0 0, 2 8, 5 12, 8 8, 10 0), (1, 3, 5, 3, 1))",
      "NURBSCURVE(3, (0 0, 3 10, 7 10, 10 0), (1, 1, 1, 1), (0, 0, 0, 0, 1, 1, "
      "1, 1))",
      "NURBSCURVE(3, (0 0, 3 10, 7 10, 10 0), (1, 2, 2, 1), (0, 0, 0, 0, 1, 1, "
      "1, 1))",
      "NURBSCURVE(3, (0 0, 3 10, 7 10, 10 0), (1, 1, 1, 1), (0, 0, 0, 0, 0.3, "
      "0.7, 1, 1))",
      "NURBSCURVE(4, (0 0, 2 8, 5 12, 8 8, 10 0))",
      "NURBSCURVE(4, (0 0, 1 5, 3 8, 7 6, 10 10, 12 0))",
      "NURBSCURVE(4, (0 0, 2 8, 5 12, 8 8, 10 0), (1, 1.5, 2, 1.5, 1))",
      "NURBSCURVE(4, (0 0, 2 8, 5 12, 8 8, 10 0), (1, 1, 1, 1, 1), (0, 0, 0, "
      "0, 0, 1, 1, 1, 1, 1))",
      "NURBSCURVE(5, (0 0, 1 2, 3 4, 5 6, 7 4, 8 2, 10 0), (1, 1, 1, 1, 1, 1, "
      "1))",
      "NURBSCURVE Z(1, (0 0 0, 10 10 5))",
      "NURBSCURVE Z(2, (0 0 0, 5 10 5, 10 0 0))",
      "NURBSCURVE Z(2, (0 0 0, 5 10 5, 10 0 0), (1, 2, 1))",
      "NURBSCURVE Z(3, (0 0 0, 3 10 3, 7 10 7, 10 0 10))",
      "NURBSCURVE Z(3, (0 0 0, 3 10 3, 7 10 7, 10 0 10), (1, 1.5, 1.5, 1), (0, "
      "0, 0, 0, 1, 1, 1, 1))",
      "NURBSCURVE M(1, (0 0 0, 10 10 3600))",
      "NURBSCURVE M(2, (0 0 0, 5 10 1800, 10 0 3600))",
      "NURBSCURVE M(2, (0 0 0, 5 10 1800, 10 0 3600), (1, 2, 1))",
      "NURBSCURVE ZM(1, (0 0 0 0, 10 10 5 3600))",
      "NURBSCURVE ZM(2, (0 0 0 0, 5 10 5 1800, 10 0 0 3600))",
      "NURBSCURVE ZM(2, (0 0 0 0, 5 10 5 1800, 10 0 0 3600), (1, 2, 1))",
      "NURBSCURVE ZM(3, (0 0 0 0, 3 10 3 1200, 7 10 7 2400, 10 0 10 3600), (1, "
      "1.5, 1.5, 1), (0, 0, 0, 0, 1, 1, 1, 1))",
      "NURBSCURVE(2, (1 0, 1 1, 0 1), (1, 0.5, 1), (0, 0, 0, 1, 1, 1))",
      "NURBSCURVE(2, (2 0, 2 1, 0 1), (1, 0.5, 1), (0, 0, 0, 1, 1, 1))",
      "NURBSCURVE(2, (0 0, 1 1, 2 0), (1, 1, 1), (0, 0, 0, 1, 1, 1))",
      "NURBSCURVE(2, (0 0, 5 5, 10 0))",
      "NURBSCURVE(2, (0 0, 1 2, 2 1, 3 3, 4 2, 5 4, 6 3, 7 5, 8 4, 9 6, 10 "
      "5))"};
  for (const auto &wktString : wkts) {
    BOOST_CHECK_NO_THROW(auto geom = SFCGAL::io::readWkt(wktString));
    auto geom         = SFCGAL::io::readWkt(wktString);
    auto wkb          = geom->asWkb(boost::endian::order::native, true);
    auto geom_bak     = SFCGAL::io::readWkb(wkb, true);
    auto wkt_bak      = geom_bak->asText();
    auto geom_wkt_bak = SFCGAL::io::readWkt(wkt_bak);
    auto wkb_bak_bak  = geom_wkt_bak->asWkb(boost::endian::order::native, true);

    BOOST_CHECK_EQUAL(geom->as<NURBSCurve>().degree(),
                      geom_bak->as<NURBSCurve>().degree());
    BOOST_CHECK_EQUAL(geom->as<NURBSCurve>().numControlPoints(),
                      geom_bak->as<NURBSCurve>().numControlPoints());
    BOOST_CHECK_EQUAL(geom->as<NURBSCurve>().isRational(),
                      geom_bak->as<NURBSCurve>().isRational());
    BOOST_CHECK_EQUAL(geom->as<NURBSCurve>().knotVector().size(),
                      geom_bak->as<NURBSCurve>().knotVector().size());
    BOOST_CHECK_EQUAL(geom->as<NURBSCurve>().weights().size(),
                      geom_bak->as<NURBSCurve>().weights().size());

    // Verify round-trip consistency
    BOOST_CHECK_EQUAL(wkb, wkb_bak_bak);
  }
}
// NOLINTEND(bugprone-suspicious-missing-comma)

/// Test new approximation modes
BOOST_AUTO_TEST_CASE(testApproximationWithFixedEndpoints)
{

  // Create test data - simple noisy line
  std::vector<Point> dataPoints;
  dataPoints.emplace_back(0.0, 0.0);
  dataPoints.emplace_back(1.0, 1.1); // slight deviation
  dataPoints.emplace_back(2.0, 1.9); // slight deviation
  dataPoints.emplace_back(3.0, 3.0);

  const unsigned int degree           = 2;
  const size_t       numControlPoints = 3;
  const auto         tolerance        = NURBSCurve::FT(EPSILON);

  // Test approximation (now uses geomdl-like behavior by default)
  auto curve = NURBSCurve::approximateCurve(dataPoints, degree, tolerance,
                                            numControlPoints);

  BOOST_REQUIRE(curve != nullptr);
  BOOST_CHECK_EQUAL(curve->numControlPoints(), numControlPoints);

  // Check endpoints are fixed (geomdl-like behavior)
  auto startPoint = curve->evaluate(0.0);
  auto endPoint   = curve->evaluate(1.0);

  double startDist = std::sqrt(
      std::pow(CGAL::to_double(startPoint.x() - dataPoints[0].x()), 2) +
      std::pow(CGAL::to_double(startPoint.y() - dataPoints[0].y()), 2));
  double endDist = std::sqrt(
      std::pow(CGAL::to_double(endPoint.x() - dataPoints.back().x()), 2) +
      std::pow(CGAL::to_double(endPoint.y() - dataPoints.back().y()), 2));

  // Endpoints should be fixed
  BOOST_CHECK_SMALL(startDist, EPSILON);
  BOOST_CHECK_SMALL(endDist, EPSILON);

  // Test that fitCurve works without mode parameter
  auto fitCurve = NURBSCurve::fitCurve(
      dataPoints, degree, NURBSCurve::FitMethod::APPROXIMATE,
      NURBSCurve::KnotMethod::CHORD_LENGTH, NURBSCurve::EndCondition::CLAMPED,
      tolerance, numControlPoints);

  BOOST_REQUIRE(fitCurve != nullptr);
  BOOST_CHECK_EQUAL(fitCurve->numControlPoints(), numControlPoints);
}

/// Test subcurve extraction
BOOST_AUTO_TEST_CASE(testSubcurve)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(2.0, 4.0);
  controlPoints.emplace_back(4.0, 0.0);
  controlPoints.emplace_back(6.0, 2.0);
  NURBSCurve originalCurve(controlPoints, 3);

  auto                  bounds = originalCurve.parameterBounds();
  NURBSCurve::Parameter quarterParam =
      bounds.first + (bounds.second - bounds.first) / NURBSCurve::FT(4);
  NURBSCurve::Parameter threeQuarterParam =
      bounds.first + 3 * (bounds.second - bounds.first) / NURBSCurve::FT(4);

  // Extract middle half of curve
  auto subcurve = originalCurve.subcurve(quarterParam, threeQuarterParam);

  BOOST_REQUIRE(subcurve != nullptr);
  BOOST_CHECK(!subcurve->isEmpty());

  // Check that subcurve endpoints match original curve at those parameters
  auto  subcurveBounds       = subcurve->parameterBounds();
  Point subcurveStart        = subcurve->evaluate(subcurveBounds.first);
  Point subcurveEnd          = subcurve->evaluate(subcurveBounds.second);
  Point originalQuarter      = originalCurve.evaluate(quarterParam);
  Point originalThreeQuarter = originalCurve.evaluate(threeQuarterParam);

  BOOST_CHECK(isNearlyEqual(subcurveStart, originalQuarter, EPSILON));
  BOOST_CHECK(isNearlyEqual(subcurveEnd, originalThreeQuarter, EPSILON));
}

/// Test subcurve with full range (should be geometrically equivalent)
BOOST_AUTO_TEST_CASE(testSubcurveFullRange)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);
  controlPoints.emplace_back(2.0, 0.0);
  NURBSCurve originalCurve(controlPoints, 2);

  auto bounds   = originalCurve.parameterBounds();
  auto subcurve = originalCurve.subcurve(bounds.first, bounds.second);

  BOOST_REQUIRE(subcurve != nullptr);
  BOOST_CHECK_EQUAL(subcurve->degree(), originalCurve.degree());

  // Note: subcurve may have different internal representation (more control
  // points) but should be geometrically equivalent

  // Check that endpoints are preserved
  auto  subcurveBounds = subcurve->parameterBounds();
  Point originalStart  = originalCurve.evaluate(bounds.first);
  Point originalEnd    = originalCurve.evaluate(bounds.second);
  Point subcurveStart  = subcurve->evaluate(subcurveBounds.first);
  Point subcurveEnd    = subcurve->evaluate(subcurveBounds.second);

  BOOST_CHECK(isNearlyEqual(originalStart, subcurveStart, EPSILON));
  BOOST_CHECK(isNearlyEqual(originalEnd, subcurveEnd, EPSILON));

  // Check that midpoint is preserved (with more tolerance since subcurve may
  // use different parameterization)
  Point originalMid = originalCurve.evaluate((bounds.first + bounds.second) /
                                             NURBSCurve::FT(2));
  Point subcurveMid = subcurve->evaluate(
      (subcurveBounds.first + subcurveBounds.second) / NURBSCurve::FT(2));

  double midDistance =
      CGAL::to_double(algorithm::distance(originalMid, subcurveMid));
  BOOST_CHECK(midDistance <
              0.1); // Allow some tolerance for different parameterizations
}

/// Test reparameterizeByArcLength method
BOOST_AUTO_TEST_CASE(testReparameterizeByArcLength)
{
  // Create a curve where parameter and arc length differ significantly
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(0.5, 10.0); // High point creates non-uniform speed
  controlPoints.emplace_back(1.0, 0.0);

  auto weights = convertWeights({1.0, 5.0, 1.0}); // High weight on middle point
  NURBSCurve originalCurve(controlPoints, weights, 2);

  // Reparameterize by arc length
  auto reparameterizedCurve = originalCurve.reparameterizeByArcLength();

  BOOST_REQUIRE(reparameterizedCurve != nullptr);
  BOOST_CHECK(!reparameterizedCurve->isEmpty());

  // Check that endpoints are preserved
  auto originalBounds = originalCurve.parameterBounds();
  auto reparamBounds  = reparameterizedCurve->parameterBounds();

  Point originalStart = originalCurve.evaluate(originalBounds.first);
  Point originalEnd   = originalCurve.evaluate(originalBounds.second);
  Point reparamStart  = reparameterizedCurve->evaluate(reparamBounds.first);
  Point reparamEnd    = reparameterizedCurve->evaluate(reparamBounds.second);

  BOOST_CHECK(isNearlyEqual(originalStart, reparamStart, EPSILON));
  BOOST_CHECK(isNearlyEqual(originalEnd, reparamEnd, EPSILON));

  // For arc-length parameterized curve, the derivative magnitude should be more
  // uniform We can't easily test this without complex calculations, so we just
  // verify it doesn't crash and produces a valid curve
  BOOST_CHECK(algorithm::isValid(*reparameterizedCurve).valid());
}

BOOST_AUTO_TEST_SUITE_END()
