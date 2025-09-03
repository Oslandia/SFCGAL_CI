// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "SFCGAL/Envelope.h"
#include "SFCGAL/Exception.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/NURBSCurve.h"
#include "SFCGAL/io/wkt.h"

using namespace SFCGAL;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_NURBSCurveTest)

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
  for (size_t i = 0; i < curve.weights().size(); ++i) {
    BOOST_CHECK_EQUAL(curve.weight(i), 1.0);
  }
}

/// NURBSCurve(controlPoints, degree, weights);
BOOST_AUTO_TEST_CASE(constructorWithWeights)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 2.0);
  controlPoints.emplace_back(3.0, 2.0);

  std::vector<double> weights{1.0, 2.0, 1.0};

  NURBSCurve curve(controlPoints, 2, weights);
  BOOST_REQUIRE_EQUAL(curve.numControlPoints(), 3U);
  BOOST_CHECK_EQUAL(curve.degree(), 2U);
  BOOST_CHECK(curve.isRational()); // Non-uniform weights
  BOOST_CHECK_EQUAL(curve.weight(1), 2.0);
}

/// NURBSCurve(controlPoints, degree, weights, knotVector);
BOOST_AUTO_TEST_CASE(constructorWithKnotVector)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 2.0);
  controlPoints.emplace_back(2.0, 0.0);

  std::vector<double> weights{1.0, 1.5, 1.0};
  std::vector<double> knots{0.0, 0.0, 0.0, 1.0, 1.0, 1.0};

  NURBSCurve curve(controlPoints, 2, weights, knots);
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

  std::vector<double> weights{1.0, 2.0};

  NURBSCurve original(controlPoints, 1, weights);
  NURBSCurve copy(original);
  NURBSCurve assigned;
  assigned = original;

  BOOST_CHECK_EQUAL(copy.numControlPoints(), 2U);
  BOOST_CHECK_EQUAL(copy.weight(1), 2.0);
  BOOST_CHECK_EQUAL(assigned.numControlPoints(), 2U);
  BOOST_CHECK_EQUAL(assigned.weight(1), 2.0);

  // Verify independence
  copy.setWeight(1, 3.0);
  BOOST_CHECK_EQUAL(original.weight(1), 2.0);
  BOOST_CHECK_EQUAL(copy.weight(1), 3.0);
}

//-- Error handling in constructors

BOOST_AUTO_TEST_CASE(constructorWithInvalidWeights)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);

  // Wrong number of weights
  std::vector<double> wrongSizeWeights{1.0};
  BOOST_CHECK_THROW(NURBSCurve curve(controlPoints, 1, wrongSizeWeights),
                    Exception);

  // Negative weight
  std::vector<double> negativeWeights{1.0, -0.5};
  BOOST_CHECK_THROW(NURBSCurve curve(controlPoints, 1, negativeWeights),
                    Exception);

  // Zero weight
  std::vector<double> zeroWeights{1.0, 0.0};
  BOOST_CHECK_THROW(NURBSCurve curve(controlPoints, 1, zeroWeights), Exception);
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

  std::vector<double> weights{1.0, 2.0};

  NURBSCurve                  curve(controlPoints, 1, weights);
  std::unique_ptr<NURBSCurve> cloned(curve.clone());

  BOOST_REQUIRE(cloned != nullptr);
  BOOST_CHECK_EQUAL(cloned->numControlPoints(), 2U);
  BOOST_CHECK_EQUAL(cloned->weight(1), 2.0);
}

/// bool isEmpty() const;
BOOST_AUTO_TEST_CASE(testIsEmpty)
{
  NURBSCurve emptyCurve;
  BOOST_CHECK(emptyCurve.isEmpty());

  NURBSCurve         nonEmptyCurve;
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
  controlPointsMeasured.emplace_back(0.0, 0.0, Kernel::FT(NaN()), 1.0);
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
  std::vector<double> uniformWeights{1.0, 1.0, 1.0};
  NURBSCurve          uniformCurve(controlPoints, 2, uniformWeights);
  BOOST_CHECK(!uniformCurve.isRational());

  // Non-uniform weights - rational
  std::vector<double> nonUniformWeights{1.0, 2.0, 1.0};
  NURBSCurve          rationalCurve(controlPoints, 2, nonUniformWeights);
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

  std::vector<double> weights{1.0, 2.0, 1.5};
  NURBSCurve          curve(controlPoints, 2, weights);

  // Test individual weight access
  BOOST_CHECK_EQUAL(curve.weight(0), 1.0);
  BOOST_CHECK_EQUAL(curve.weight(1), 2.0);
  BOOST_CHECK_EQUAL(curve.weight(2), 1.5);

  // Test weight modification
  curve.setWeight(1, 3.0);
  BOOST_CHECK_EQUAL(curve.weight(1), 3.0);

  // Test weight vector access
  const std::vector<double> &weightVector = curve.weights();
  BOOST_CHECK_EQUAL(weightVector.size(), 3U);
  BOOST_CHECK_EQUAL(weightVector[1], 3.0);

  // Test setting all weights
  std::vector<double> newWeights{0.5, 1.0, 2.0};
  curve.setWeights(newWeights);
  BOOST_CHECK_EQUAL(curve.weight(0), 0.5);
  BOOST_CHECK_EQUAL(curve.weight(2), 2.0);
}

BOOST_AUTO_TEST_CASE(testInvalidWeightOperations)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);

  NURBSCurve curve(controlPoints, 1);

  // Setting negative weight
  BOOST_CHECK_THROW(curve.setWeight(0, -1.0), Exception);

  // Setting zero weight
  BOOST_CHECK_THROW(curve.setWeight(0, 0.0), Exception);

  // Setting weights with wrong size
  std::vector<double> wrongSizeWeights{1.0, 2.0, 3.0}; // Too many
  BOOST_CHECK_THROW(curve.setWeights(wrongSizeWeights), Exception);

  // Setting weights with negative values
  std::vector<double> negativeWeights{1.0, -0.5};
  BOOST_CHECK_THROW(curve.setWeights(negativeWeights), Exception);
}

//-- Curve evaluation tests

/// Point evaluate(double parameter) const;
BOOST_AUTO_TEST_CASE(testEvaluateEmpty)
{
  NURBSCurve emptyCurve;
  BOOST_CHECK_THROW(auto result = emptyCurve.evaluate(0.5), Exception);
}

BOOST_AUTO_TEST_CASE(testEvaluateParameterBounds)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);

  NURBSCurve curve(controlPoints, 1);
  auto       bounds = curve.parameterBounds();

  BOOST_CHECK_THROW(auto result = curve.evaluate(bounds.first - 0.1),
                    Exception);
  BOOST_CHECK_THROW(auto result = curve.evaluate(bounds.second + 0.1),
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
  Point      start = uniformCurve.evaluate(0.0);
  Point      end   = uniformCurve.evaluate(1.0);
  Point      mid   = uniformCurve.evaluate(0.5);

  BOOST_CHECK_EQUAL(start.x(), 0.0);
  BOOST_CHECK_EQUAL(start.y(), 0.0);
  BOOST_CHECK_EQUAL(end.x(), 2.0);
  BOOST_CHECK_EQUAL(end.y(), 2.0);
  BOOST_CHECK_EQUAL(mid.x(), 1.0);
  BOOST_CHECK_EQUAL(mid.y(), 1.0);
}

BOOST_AUTO_TEST_CASE(testEvaluateRationalNURBS)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 2.0);
  controlPoints.emplace_back(2.0, 0.0);

  // Create rational NURBS with higher weight on middle point
  std::vector<double> weights{1.0, 5.0, 1.0};
  NURBSCurve          rationalCurve(controlPoints, 2, weights);

  Point start = rationalCurve.evaluate(0.0);
  Point end   = rationalCurve.evaluate(1.0);
  Point mid   = rationalCurve.evaluate(0.5);

  BOOST_CHECK_EQUAL(start.x(), 0.0);
  BOOST_CHECK_EQUAL(start.y(), 0.0);
  BOOST_CHECK_EQUAL(end.x(), 2.0);
  BOOST_CHECK_EQUAL(end.y(), 0.0);

  // Middle point should be closer to the highly weighted control point
  BOOST_CHECK(mid.y() > 1.0); // Should be pulled toward control point (1,2)
}

//-- Derivative tests

BOOST_AUTO_TEST_CASE(testDerivativeOrderZero)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);

  NURBSCurve curve(controlPoints, 1);
  Point      point         = curve.derivative(0.5, 0);
  Point      evaluatePoint = curve.evaluate(0.5);

  BOOST_CHECK_EQUAL(point.x(), evaluatePoint.x());
  BOOST_CHECK_EQUAL(point.y(), evaluatePoint.y());
}

BOOST_AUTO_TEST_CASE(testDerivativeLinear)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(2.0, 4.0);

  NURBSCurve curve(controlPoints, 1);
  Point      derivative = curve.derivative(0.5, 1);

  // For linear curve, derivative should be constant
  BOOST_CHECK_EQUAL(derivative.x(), 2.0);
  BOOST_CHECK_EQUAL(derivative.y(), 4.0);
}

BOOST_AUTO_TEST_CASE(testDerivativeRational)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 2.0);
  controlPoints.emplace_back(2.0, 0.0);

  std::vector<double> weights{1.0, 3.0, 1.0};
  NURBSCurve          rationalCurve(controlPoints, 2, weights);

  // Test that rational derivative computation doesn't throw
  BOOST_CHECK_NO_THROW(Point derivative = rationalCurve.derivative(0.5, 1));

  Point derivative = rationalCurve.derivative(0.5, 1);
  // The derivative should be non-zero for this configuration
  BOOST_CHECK(std::abs(CGAL::to_double(derivative.x())) > 1e-10 ||
              std::abs(CGAL::to_double(derivative.y())) > 1e-10);
}

//-- 3D and measured coordinates

BOOST_AUTO_TEST_CASE(testEvaluate3D)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0, 2.0);

  NURBSCurve curve(controlPoints, 1);
  BOOST_CHECK(curve.is3D());

  Point mid = curve.evaluate(0.5);
  BOOST_CHECK(mid.is3D());
  BOOST_CHECK_EQUAL(mid.x(), 0.5);
  BOOST_CHECK_EQUAL(mid.y(), 0.5);
  BOOST_CHECK_EQUAL(mid.z(), 1.0);
}

BOOST_AUTO_TEST_CASE(testEvaluateMeasured)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0, 0.0, 0.0, COORDINATE_XYM);
  controlPoints.emplace_back(2.0, 2.0, 0.0, 4.0, COORDINATE_XYM);

  NURBSCurve curve(controlPoints, 1);
  BOOST_CHECK(curve.isMeasured());

  Point mid = curve.evaluate(0.5);
  BOOST_CHECK(mid.isMeasured());
  BOOST_CHECK_EQUAL(mid.x(), 1.0);
  BOOST_CHECK_EQUAL(mid.y(), 1.0);
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

  Point mid = curve.evaluate(0.5);
  BOOST_CHECK(mid.is3D());
  BOOST_CHECK(mid.isMeasured());
  BOOST_CHECK_EQUAL(mid.x(), 1.0);
  BOOST_CHECK_EQUAL(mid.y(), 1.0);
  BOOST_CHECK_EQUAL(mid.z(), 2.0);
  BOOST_CHECK_EQUAL(mid.m(), 4.0);
}

//-- Validation tests

BOOST_AUTO_TEST_CASE(testIsValidEmpty)
{
  NURBSCurve emptyCurve;
  BOOST_CHECK(emptyCurve.isValid());
}

BOOST_AUTO_TEST_CASE(testIsValidConsistent)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);
  controlPoints.emplace_back(2.0, 0.0);

  std::vector<double> weights{1.0, 2.0, 1.0};
  NURBSCurve          curve(controlPoints, 2, weights);
  BOOST_CHECK(curve.isValid());
}

//-- Homogeneous coordinates

BOOST_AUTO_TEST_CASE(testHomogeneousControlPoints)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(1.0, 2.0);
  controlPoints.emplace_back(3.0, 4.0);

  std::vector<double> weights{2.0, 0.5};
  NURBSCurve          curve(controlPoints, 1, weights);

  auto homogeneous = curve.homogeneousControlPoints();
  BOOST_CHECK_EQUAL(homogeneous.size(), 2U);

  // First point: (1, 2) with weight 2.0
  BOOST_CHECK_EQUAL(homogeneous[0].size(), 3U); // [wx, wy, w]
  BOOST_CHECK_EQUAL(homogeneous[0][0], 2.0);    // wx = 1.0 * 2.0
  BOOST_CHECK_EQUAL(homogeneous[0][1], 4.0);    // wy = 2.0 * 2.0
  BOOST_CHECK_EQUAL(homogeneous[0][2], 2.0);    // w = 2.0

  // Second point: (3, 4) with weight 0.5
  BOOST_CHECK_EQUAL(homogeneous[1][0], 1.5); // wx = 3.0 * 0.5
  BOOST_CHECK_EQUAL(homogeneous[1][1], 2.0); // wy = 4.0 * 0.5
  BOOST_CHECK_EQUAL(homogeneous[1][2], 0.5); // w = 0.5
}

BOOST_AUTO_TEST_CASE(testHomogeneousControlPoints3D)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(1.0, 2.0, 3.0);

  std::vector<double> weights{2.0};
  NURBSCurve          curve(controlPoints, 0, weights);

  auto homogeneous = curve.homogeneousControlPoints();
  BOOST_CHECK_EQUAL(homogeneous[0].size(), 4U); // [wx, wy, wz, w]
  BOOST_CHECK_EQUAL(homogeneous[0][0], 2.0);    // wx = 1.0 * 2.0
  BOOST_CHECK_EQUAL(homogeneous[0][1], 4.0);    // wy = 2.0 * 2.0
  BOOST_CHECK_EQUAL(homogeneous[0][2], 6.0);    // wz = 3.0 * 2.0
  BOOST_CHECK_EQUAL(homogeneous[0][3], 2.0);    // w = 2.0
}

//-- Knot manipulation tests

BOOST_AUTO_TEST_CASE(testKnotInsertion)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);
  controlPoints.emplace_back(2.0, 0.0);

  NURBSCurve curve(controlPoints, 2);
  size_t     originalNumPoints = curve.numControlPoints();
  auto       originalBounds    = curve.parameterBounds();

  // Insert knot at parameter 0.5
  curve.insertKnot(0.5);

  // Should have one more control point
  BOOST_CHECK_EQUAL(curve.numControlPoints(), originalNumPoints + 1);

  // Parameter bounds should remain the same
  auto newBounds = curve.parameterBounds();
  BOOST_CHECK_EQUAL(newBounds.first, originalBounds.first);
  BOOST_CHECK_EQUAL(newBounds.second, originalBounds.second);
}

BOOST_AUTO_TEST_CASE(testKnotInsertionBounds)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);

  NURBSCurve curve(controlPoints, 1);
  auto       bounds = curve.parameterBounds();

  BOOST_CHECK_THROW(curve.insertKnot(bounds.first - 0.1), Exception);
  BOOST_CHECK_THROW(curve.insertKnot(bounds.second + 0.1), Exception);
  BOOST_CHECK_NO_THROW(curve.insertKnot(bounds.first + 0.1));
}

BOOST_AUTO_TEST_CASE(testKnotRefinement)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);
  controlPoints.emplace_back(2.0, 0.0);

  NURBSCurve curve(controlPoints, 2);
  size_t     originalNumPoints = curve.numControlPoints();

  std::vector<double> newKnots{0.3, 0.7};
  curve.refineKnots(newKnots);

  // Should have added 2 control points
  BOOST_CHECK_EQUAL(curve.numControlPoints(), originalNumPoints + 2);
}

//-- Factory methods

BOOST_AUTO_TEST_CASE(testCreateCircularArc2D)
{
  Point  center(0.0, 0.0);
  double radius     = 1.0;
  double startAngle = 0.0;
  double endAngle   = M_PI; // Half circle

  auto arc = NURBSCurve::createCircularArc(center, radius, startAngle, endAngle,
                                           false);

  BOOST_REQUIRE(arc != nullptr);
  BOOST_CHECK(!arc->is3D());
  BOOST_CHECK(arc->isRational());
  BOOST_CHECK_EQUAL(arc->degree(), 2U);

  // Test that it actually represents a circular arc
  Point start = arc->evaluate(0.0);
  Point end   = arc->evaluate(1.0);

  BOOST_CHECK_CLOSE(start.x(), cos(0.0), 1e-10); // cos(0) = 1
  BOOST_CHECK_CLOSE(start.y(), sin(0.0), 1e-10); // sin(0) = 0
  BOOST_CHECK_CLOSE(end.x(), cos(M_PI), 1e-10);  // cos(π) = -1
  BOOST_CHECK_CLOSE(end.y(), sin(M_PI),
                    1e-10); // sin(π) = 0. But in double: 1.2246467991473532e-16
}

BOOST_AUTO_TEST_CASE(testCreateCircularArc3D)
{
  Point  center(0.0, 0.0, 5.0);
  double radius     = 2.0;
  double startAngle = 0.0;
  double endAngle   = M_PI / 2; // Quarter circle

  auto arc =
      NURBSCurve::createCircularArc(center, radius, startAngle, endAngle, true);

  BOOST_REQUIRE(arc != nullptr);
  BOOST_CHECK(arc->is3D());
  BOOST_CHECK(arc->isRational());

  Point start = arc->evaluate(0.0);
  Point end   = arc->evaluate(1.0);

  BOOST_CHECK_CLOSE(start.x(), 2.0 * cos(0.0), 1e-10); // 2 * cos(0) = 2
  BOOST_CHECK_CLOSE(start.y(), 2.0 * sin(0.0), 1e-10); // 2 * sin(0) = 0
  BOOST_CHECK_EQUAL(start.z(), 5.0);                   // Z coordinate preserved
  BOOST_CHECK_CLOSE(end.x(), 2.0 * cos(M_PI / 2.0), 1e-10); // 2 * cos(π/2) = 0
  BOOST_CHECK_CLOSE(end.y(), 2.0 * sin(M_PI / 2.0), 1e-10); // 2 * sin(π/2) = 2
  BOOST_CHECK_EQUAL(end.z(), 5.0); // Z coordinate preserved
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

  BOOST_CHECK(box.xMin() >= 1.0);
  BOOST_CHECK(box.xMax() <= 3.0);
  BOOST_CHECK(box.yMin() >= 5.0);
  BOOST_CHECK(box.yMax() <= 9.0);
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

  BOOST_CHECK(box.xMin() >= 1.0);
  BOOST_CHECK(box.xMax() <= 3.0);
  BOOST_CHECK(box.yMin() >= 5.0);
  BOOST_CHECK(box.yMax() <= 9.0);
  BOOST_CHECK(box.zMin() >= 11.0);
  BOOST_CHECK(box.zMax() <= 17.0);
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
  BOOST_CHECK_EQUAL(curve.weight(1), 2.0);
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
  BOOST_CHECK_EQUAL(curve.controlPointAt(1).z(), 10.0);
}

/// WKT reading with measured coordinates
BOOST_AUTO_TEST_CASE(testReadMeasuredNURBSWkt)
{
  std::string const wkt  = "NURBSCURVE M ((0 0 5, 5 5 15, 10 0 25), 2)";
  auto              geom = SFCGAL::io::readWkt(wkt);

  BOOST_REQUIRE(geom->is<NURBSCurve>());
  auto const &curve = geom->as<NURBSCurve>();
  BOOST_CHECK(curve.isMeasured());
  BOOST_CHECK_EQUAL(curve.controlPointAt(2).m(), 25.0);
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

  std::vector<double> weights{1.0, 2.0, 1.0};
  NURBSCurve          curve(controlPoints, 2, weights);
  std::string         result = curve.asText(0);

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
  BOOST_CHECK_EQUAL(curve2.weight(1), 2.0);
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
