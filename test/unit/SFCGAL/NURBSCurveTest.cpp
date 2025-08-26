// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "SFCGAL/Envelope.h"
#include "SFCGAL/Exception.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/NURBSCurve.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/io/wkt.h"

using namespace SFCGAL;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_NURBSCurveTest)

//-- Helper function to convert std::vector<double> to std::vector<FT>
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

std::vector<NURBSCurve::Knot>
convertKnots(const std::vector<double> &doubleKnots)
{
  std::vector<NURBSCurve::Knot> knots;
  knots.reserve(doubleKnots.size());
  for (const auto &knot : doubleKnots) {
    knots.emplace_back(knot);
  }
  return knots;
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

BOOST_AUTO_TEST_CASE(testWktEmpty)
{
  NURBSCurve const emptyCurve;
  std::string      wkt = emptyCurve.asText(0);
  BOOST_CHECK(wkt.find("EMPTY") != std::string::npos);
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

BOOST_AUTO_TEST_SUITE_END()
