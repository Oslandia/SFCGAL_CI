// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "SFCGAL/BSplineCurve.h"
#include "SFCGAL/Envelope.h"
#include "SFCGAL/Exception.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/io/wkt.h"

using namespace SFCGAL;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_BSplineCurveTest)

//-- Constructor tests

/// BSplineCurve() ;
BOOST_AUTO_TEST_CASE(defaultConstructor)
{
  BSplineCurve const curve;
  BOOST_CHECK(curve.isEmpty());
  BOOST_CHECK(!curve.is3D());
  BOOST_CHECK(!curve.isMeasured());
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 0U);
  BOOST_CHECK_EQUAL(curve.degree(), 0U);
  BOOST_CHECK(curve.knotVector().empty());
}

/// BSplineCurve(const std::vector<Point> &controlPoints, unsigned int degree) ;
BOOST_AUTO_TEST_CASE(constructorWithControlPointsAndDegree)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 2.0);
  controlPoints.emplace_back(3.0, 2.0);
  controlPoints.emplace_back(4.0, 0.0);

  BSplineCurve curve(controlPoints, 2);
  BOOST_REQUIRE_EQUAL(curve.numControlPoints(), 4U);
  BOOST_CHECK_EQUAL(curve.degree(), 2U);
  BOOST_CHECK(!curve.isEmpty());
  BOOST_CHECK_EQUAL(curve.controlPointAt(0).x(), 0.0);
  BOOST_CHECK_EQUAL(curve.controlPointAt(0).y(), 0.0);
  BOOST_CHECK_EQUAL(curve.controlPointAt(3).x(), 4.0);
  BOOST_CHECK_EQUAL(curve.controlPointAt(3).y(), 0.0);

  // Check knot vector was generated
  BOOST_CHECK(!curve.knotVector().empty());
  BOOST_CHECK_EQUAL(curve.knotVector().size(), 7U); // n + p + 1 = 4 + 2 + 1
}

/// BSplineCurve(const std::vector<Point> &controlPoints, unsigned int degree,
///              const std::vector<double> &knotVector) ;
BOOST_AUTO_TEST_CASE(constructorWithKnotVector)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);
  controlPoints.emplace_back(2.0, 0.0);
  controlPoints.emplace_back(3.0, 1.0);

  std::vector<double> knots = {0, 0, 0, 0.5,
                               1, 1, 1}; // degree 2, 4 control points

  BSplineCurve curve(controlPoints, 2, knots);
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 4U);
  BOOST_CHECK_EQUAL(curve.degree(), 2U);
  BOOST_CHECK_EQUAL(curve.knotVector().size(), 7U);
  BOOST_CHECK_EQUAL(curve.knotVector()[3], 0.5);
}

BOOST_AUTO_TEST_CASE(constructorWithInvalidDegree)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);

  // Degree >= number of control points should throw
  BOOST_CHECK_THROW(BSplineCurve curve(controlPoints, 2), Exception);
  BOOST_CHECK_THROW(BSplineCurve curve(controlPoints, 3), Exception);
}

BOOST_AUTO_TEST_CASE(constructorWithInvalidKnotVector)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);
  controlPoints.emplace_back(2.0, 0.0);

  // Wrong size knot vector
  std::vector<double> wrongSizeKnots = {
      0, 0, 1, 1}; // Should be 5 for degree 1, 3 points
  BOOST_CHECK_THROW(BSplineCurve curve(controlPoints, 1, wrongSizeKnots),
                    Exception);

  // Non-monotonic knot vector
  std::vector<double> nonMonotonicKnots = {0, 0, 0.7, 0.5, 1};
  BOOST_CHECK_THROW(BSplineCurve curve(controlPoints, 1, nonMonotonicKnots),
                    Exception);
}

/// BSplineCurve(const BSplineCurve& other) ;
/// BSplineCurve& operator=(const BSplineCurve& other) ;
BOOST_AUTO_TEST_CASE(copyConstructorAndAssignment)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);
  controlPoints.emplace_back(2.0, 0.0);

  BSplineCurve original(controlPoints, 1);
  BSplineCurve copy(original);
  BSplineCurve assigned;
  assigned = original;

  BOOST_CHECK_EQUAL(copy.numControlPoints(), 3U);
  BOOST_CHECK_EQUAL(copy.degree(), 1U);
  BOOST_CHECK_EQUAL(assigned.numControlPoints(), 3U);
  BOOST_CHECK_EQUAL(assigned.degree(), 1U);

  // Verify independence
  copy.addControlPoint(Point(3.0, 1.0));
  BOOST_CHECK_EQUAL(original.numControlPoints(), 3U);
  BOOST_CHECK_EQUAL(copy.numControlPoints(), 4U);
}

//-- Basic operations tests

/// void clear() ;
BOOST_AUTO_TEST_CASE(testClear)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);

  BSplineCurve curve(controlPoints, 1);
  BOOST_CHECK(!curve.isEmpty());
  BOOST_CHECK(!curve.knotVector().empty());

  curve.clear();
  BOOST_CHECK(curve.isEmpty());
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 0U);
  BOOST_CHECK_EQUAL(curve.degree(), 0U);
  BOOST_CHECK(curve.knotVector().empty());
}

/// void addControlPoint(const Point &point) ;
BOOST_AUTO_TEST_CASE(testAddControlPoint)
{
  BSplineCurve curve;
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 0U);

  curve.addControlPoint(Point(0.0, 0.0));
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 1U);

  curve.addControlPoint(Point(1.0, 1.0));
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 2U);

  curve.addControlPoint(Point(2.0, 0.0));
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 3U);
}

BOOST_AUTO_TEST_CASE(testAddControlPoint_dimensionalConsistency)
{
  BSplineCurve curve;
  curve.addControlPoint(Point(0.0, 0.0));

  // Adding 3D point to 2D curve should throw
  BOOST_CHECK_THROW(curve.addControlPoint(Point(1.0, 1.0, 1.0)), Exception);

  // Adding measured point to non-measured curve should throw
  BOOST_CHECK_THROW(
      curve.addControlPoint(Point(1.0, 1.0, Kernel::FT(NaN()), 1.0)),
      Exception);
}

/// void removeControlPoint(size_t index) ;
BOOST_AUTO_TEST_CASE(testRemoveControlPoint)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);
  controlPoints.emplace_back(2.0, 0.0);
  controlPoints.emplace_back(3.0, 1.0);

  BSplineCurve curve(controlPoints, 2);
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 4U);

  curve.removeControlPoint(2);
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 3U);
  BOOST_CHECK_EQUAL(curve.controlPointAt(0).x(), 0.0);
  BOOST_CHECK_EQUAL(curve.controlPointAt(1).x(), 1.0);
  BOOST_CHECK_EQUAL(curve.controlPointAt(2).x(), 3.0);
}

/// const Point& controlPointAt(size_t index) const ;
/// Point& controlPointAt(size_t index) ;
/// void setControlPoint(size_t index, const Point &point) ;
BOOST_AUTO_TEST_CASE(testControlPointAccess)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);

  BSplineCurve curve(controlPoints, 1);

  // Test const access
  const BSplineCurve &constCurve = curve;
  BOOST_CHECK_EQUAL(constCurve.controlPointAt(0).x(), 0.0);
  BOOST_CHECK_EQUAL(constCurve.controlPointAt(1).y(), 1.0);

  // Test mutable access
  curve.controlPointAt(0) = Point(5.0, 5.0);
  BOOST_CHECK_EQUAL(curve.controlPointAt(0).x(), 5.0);

  // Test setControlPoint
  curve.setControlPoint(1, Point(10.0, 10.0));
  BOOST_CHECK_EQUAL(curve.controlPointAt(1).x(), 10.0);
  BOOST_CHECK_EQUAL(curve.controlPointAt(1).y(), 10.0);
}

/// std::vector<Point> controlPoints() const ;
BOOST_AUTO_TEST_CASE(testControlPoints)
{
  std::vector<Point> originalPoints;
  originalPoints.emplace_back(0.0, 0.0);
  originalPoints.emplace_back(1.0, 1.0);
  originalPoints.emplace_back(2.0, 0.0);

  BSplineCurve       curve(originalPoints, 2);
  std::vector<Point> retrievedPoints = curve.controlPoints();

  BOOST_CHECK_EQUAL(retrievedPoints.size(), 3U);
  BOOST_CHECK_EQUAL(retrievedPoints[0].x(), 0.0);
  BOOST_CHECK_EQUAL(retrievedPoints[1].x(), 1.0);
  BOOST_CHECK_EQUAL(retrievedPoints[2].x(), 2.0);
}

//-- Geometry interface tests

/// GeometryType geometryTypeId() const ;
BOOST_AUTO_TEST_CASE(testGeometryTypeId)
{
  BSplineCurve const curve;
  BOOST_CHECK_EQUAL(curve.geometryTypeId(), TYPE_BSPLINECURVE);
}

/// std::string geometryType() const ;
BOOST_AUTO_TEST_CASE(testGeometryType)
{
  BSplineCurve const curve;
  BOOST_CHECK_EQUAL(curve.geometryType(), "BSplineCurve");
}

/// BSplineCurve* clone() const ;
BOOST_AUTO_TEST_CASE(testClone)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);

  BSplineCurve                  curve(controlPoints, 1);
  std::unique_ptr<BSplineCurve> clone(curve.clone());

  BOOST_REQUIRE(clone != nullptr);
  BOOST_CHECK_EQUAL(clone->numControlPoints(), 2U);
  BOOST_CHECK_EQUAL(clone->degree(), 1U);
  BOOST_CHECK_EQUAL(clone->controlPointAt(0).x(), 0.0);
  BOOST_CHECK_EQUAL(clone->controlPointAt(1).y(), 1.0);
}

/// bool isEmpty() const ;
BOOST_AUTO_TEST_CASE(testIsEmpty)
{
  BSplineCurve emptyCurve;
  BOOST_CHECK(emptyCurve.isEmpty());

  BSplineCurve nonEmptyCurve;
  nonEmptyCurve.addControlPoint(Point(0.0, 0.0));
  BOOST_CHECK(!nonEmptyCurve.isEmpty());
}

/// bool is3D() const ;
BOOST_AUTO_TEST_CASE(testIs3D)
{
  BSplineCurve curve2D;
  curve2D.addControlPoint(Point(0.0, 0.0));
  BOOST_CHECK(!curve2D.is3D());

  BSplineCurve curve3D;
  curve3D.addControlPoint(Point(0.0, 0.0, 0.0));
  BOOST_CHECK(curve3D.is3D());

  BSplineCurve emptyCurve;
  BOOST_CHECK(!emptyCurve.is3D());
}

/// bool isMeasured() const ;
BOOST_AUTO_TEST_CASE(testIsMeasured)
{
  BSplineCurve curveRegular;
  curveRegular.addControlPoint(Point(0.0, 0.0));
  BOOST_CHECK(!curveRegular.isMeasured());

  BSplineCurve curveMeasured;
  curveMeasured.addControlPoint(Point(0.0, 0.0, Kernel::FT(NaN()), 1.0));
  BOOST_CHECK(curveMeasured.isMeasured());

  BSplineCurve emptyCurve;
  BOOST_CHECK(!emptyCurve.isMeasured());
}

//-- Curve interface tests

/// Point evaluate(double parameter) const ;
BOOST_AUTO_TEST_CASE(testEvaluate_linear)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(2.0, 2.0);

  BSplineCurve curve(controlPoints, 1);

  Point start = curve.evaluate(0.0);
  BOOST_CHECK_EQUAL(start.x(), 0.0);
  BOOST_CHECK_EQUAL(start.y(), 0.0);

  Point end = curve.evaluate(1.0);
  BOOST_CHECK_EQUAL(end.x(), 2.0);
  BOOST_CHECK_EQUAL(end.y(), 2.0);

  Point mid = curve.evaluate(0.5);
  BOOST_CHECK_EQUAL(mid.x(), 1.0);
  BOOST_CHECK_EQUAL(mid.y(), 1.0);
}

BOOST_AUTO_TEST_CASE(testEvaluate_quadratic)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 2.0);
  controlPoints.emplace_back(2.0, 0.0);

  BSplineCurve curve(controlPoints, 2);

  Point start = curve.evaluate(0.0);
  BOOST_CHECK_EQUAL(start.x(), 0.0);
  BOOST_CHECK_EQUAL(start.y(), 0.0);

  Point end = curve.evaluate(1.0);
  BOOST_CHECK_EQUAL(end.x(), 2.0);
  BOOST_CHECK_EQUAL(end.y(), 0.0);

  // B-spline with degree 2 and 3 control points is equivalent to a Bezier curve
  Point mid = curve.evaluate(0.5);
  BOOST_CHECK_EQUAL(mid.x(), 1.0);
  BOOST_CHECK_EQUAL(mid.y(), 1.0);
}

BOOST_AUTO_TEST_CASE(testEvaluate_empty)
{
  BSplineCurve emptyCurve;
  BOOST_CHECK_THROW(auto result = emptyCurve.evaluate(0.5), Exception);
}

BOOST_AUTO_TEST_CASE(testEvaluate_parameterBounds)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);
  controlPoints.emplace_back(2.0, 0.0);

  BSplineCurve curve(controlPoints, 1);
  auto         bounds = curve.parameterBounds();

  BOOST_CHECK_THROW(auto result1 = curve.evaluate(bounds.first - 0.1),
                    Exception);
  BOOST_CHECK_THROW(auto result2 = curve.evaluate(bounds.second + 0.1),
                    Exception);
  BOOST_CHECK_NO_THROW(auto result3 = curve.evaluate(bounds.first));
  BOOST_CHECK_NO_THROW(auto result4 = curve.evaluate(bounds.second));
}

/// Point derivative(double parameter, unsigned int order = 1) const ;
BOOST_AUTO_TEST_CASE(testDerivative_linear)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(2.0, 4.0);

  BSplineCurve curve(controlPoints, 1);

  Point derivative = curve.derivative(0.5, 1);
  BOOST_CHECK_EQUAL(derivative.x(), 2.0);
  BOOST_CHECK_EQUAL(derivative.y(), 4.0);
}

BOOST_AUTO_TEST_CASE(testDerivative_orderZero)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);

  BSplineCurve curve(controlPoints, 1);

  Point point         = curve.derivative(0.5, 0);
  Point evaluatePoint = curve.evaluate(0.5);
  BOOST_CHECK_EQUAL(point.x(), evaluatePoint.x());
  BOOST_CHECK_EQUAL(point.y(), evaluatePoint.y());
}

BOOST_AUTO_TEST_CASE(testDerivative_highOrder)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);

  BSplineCurve curve(controlPoints, 1);

  // Order higher than degree should return zero vector
  Point highDerivative = curve.derivative(0.5, 5);
  BOOST_CHECK_EQUAL(highDerivative.x(), 0.0);
  BOOST_CHECK_EQUAL(highDerivative.y(), 0.0);
}

/// std::unique_ptr<LineString> toLineString(unsigned int numSegments = 32)
/// const ;
BOOST_AUTO_TEST_CASE(testToLineString_empty)
{
  BSplineCurve                emptyCurve;
  std::unique_ptr<LineString> lineString = emptyCurve.toLineString(10);
  BOOST_CHECK(lineString->isEmpty());
}

BOOST_AUTO_TEST_CASE(testToLineString_linear)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(2.0, 2.0);

  BSplineCurve                curve(controlPoints, 1);
  std::unique_ptr<LineString> lineString = curve.toLineString(4);

  BOOST_CHECK_EQUAL(lineString->numPoints(), 5U); // numSegments + 1
  BOOST_CHECK_EQUAL(lineString->pointN(0).x(), 0.0);
  BOOST_CHECK_EQUAL(lineString->pointN(0).y(), 0.0);
  BOOST_CHECK_EQUAL(lineString->pointN(4).x(), 2.0);
  BOOST_CHECK_EQUAL(lineString->pointN(4).y(), 2.0);
}

BOOST_AUTO_TEST_CASE(testToLineString_invalidSegments)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);

  BSplineCurve curve(controlPoints, 1);
  BOOST_CHECK_THROW(auto result = curve.toLineString(1), Exception);
}

/// std::pair<double, double> parameterBounds() const ;
BOOST_AUTO_TEST_CASE(testParameterBounds_empty)
{
  BSplineCurve              curve;
  std::pair<double, double> bounds = curve.parameterBounds();
  BOOST_CHECK_EQUAL(bounds.first, 0.0);
  BOOST_CHECK_EQUAL(bounds.second, 1.0);
}

BOOST_AUTO_TEST_CASE(testParameterBounds_uniform)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);
  controlPoints.emplace_back(2.0, 0.0);
  controlPoints.emplace_back(3.0, 1.0);

  BSplineCurve              curve(controlPoints, 2);
  std::pair<double, double> bounds = curve.parameterBounds();

  // For degree 2 with 4 control points, bounds should be from knot[2] to
  // knot[4]
  BOOST_CHECK_EQUAL(bounds.first, 0.0);
  BOOST_CHECK_EQUAL(bounds.second, 1.0);
}

/// bool isValid() const ;
BOOST_AUTO_TEST_CASE(testIsValid_empty)
{
  BSplineCurve emptyCurve;
  BOOST_CHECK(emptyCurve.isValid());
}

BOOST_AUTO_TEST_CASE(testIsValid_consistent)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);
  controlPoints.emplace_back(2.0, 0.0);

  BSplineCurve curve(controlPoints, 2);
  BOOST_CHECK(curve.isValid());
}

//-- BSplineCurve specific tests

/// unsigned int degree() const ;
BOOST_AUTO_TEST_CASE(testDegree)
{
  BSplineCurve emptyCurve;
  BOOST_CHECK_EQUAL(emptyCurve.degree(), 0U);

  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);
  controlPoints.emplace_back(2.0, 0.0);

  BSplineCurve linear(controlPoints, 1);
  BOOST_CHECK_EQUAL(linear.degree(), 1U);

  BSplineCurve quadratic(controlPoints, 2);
  BOOST_CHECK_EQUAL(quadratic.degree(), 2U);
}

/// size_t numControlPoints() const ;
BOOST_AUTO_TEST_CASE(testNumControlPoints)
{
  BSplineCurve curve;
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 0U);

  curve.addControlPoint(Point(0.0, 0.0));
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 1U);

  curve.addControlPoint(Point(1.0, 1.0));
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 2U);
}

/// const std::vector<double>& knotVector() const ;
/// void setKnotVector(const std::vector<double> &knots) ;
BOOST_AUTO_TEST_CASE(testKnotVector)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);
  controlPoints.emplace_back(2.0, 0.0);
  controlPoints.emplace_back(3.0, 1.0);

  BSplineCurve curve(controlPoints, 2);

  // Check generated uniform knot vector
  auto knots = curve.knotVector();
  BOOST_CHECK_EQUAL(knots.size(), 7U); // n + p + 1 = 4 + 2 + 1

  // Check multiplicity at ends
  BOOST_CHECK_EQUAL(knots[0], 0.0);
  BOOST_CHECK_EQUAL(knots[1], 0.0);
  BOOST_CHECK_EQUAL(knots[2], 0.0);
  BOOST_CHECK_EQUAL(knots[4], 1.0);
  BOOST_CHECK_EQUAL(knots[5], 1.0);
  BOOST_CHECK_EQUAL(knots[6], 1.0);

  // Set custom knot vector
  std::vector<double> customKnots = {0, 0, 0, 0.3, 1, 1, 1};
  curve.setKnotVector(customKnots);
  BOOST_CHECK_EQUAL(curve.knotVector()[3], 0.3);

  // Invalid knot vector should throw
  std::vector<double> invalidKnots = {0, 0, 1, 1}; // Wrong size
  BOOST_CHECK_THROW(curve.setKnotVector(invalidKnots), Exception);
}

/// void generateUniformKnotVector() ;
BOOST_AUTO_TEST_CASE(testGenerateUniformKnotVector)
{
  std::vector<Point> controlPoints;
  for (int i = 0; i <= 5; ++i) {
    controlPoints.emplace_back(static_cast<double>(i), 0.0);
  }

  BSplineCurve curve(controlPoints, 3);
  curve.generateUniformKnotVector();

  auto knots = curve.knotVector();
  BOOST_CHECK_EQUAL(knots.size(), 10U); // 6 + 3 + 1

  // Check clamped ends
  for (size_t i = 0; i <= 3; ++i) {
    BOOST_CHECK_EQUAL(knots[i], 0.0);
    BOOST_CHECK_EQUAL(knots[9 - i], 1.0);
  }
}

//-- 3D and Measured coordinate tests

BOOST_AUTO_TEST_CASE(testEvaluate_3D)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0, 2.0);

  BSplineCurve curve(controlPoints, 1);
  BOOST_CHECK(curve.is3D());

  Point mid = curve.evaluate(0.5);
  BOOST_CHECK(mid.is3D());
  BOOST_CHECK_EQUAL(mid.x(), 0.5);
  BOOST_CHECK_EQUAL(mid.y(), 0.5);
  BOOST_CHECK_EQUAL(mid.z(), 1.0);
}

BOOST_AUTO_TEST_CASE(testEvaluate_measured)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0, 0.0, 0.0,
                             CoordinateType::COORDINATE_XYM);
  controlPoints.emplace_back(2.0, 2.0, 0.0, 4.0,
                             CoordinateType::COORDINATE_XYM);

  BSplineCurve curve(controlPoints, 1);
  BOOST_CHECK(curve.isMeasured());

  Point mid = curve.evaluate(0.5);
  BOOST_CHECK(mid.isMeasured());
  BOOST_CHECK_EQUAL(mid.x(), 1.0);
  BOOST_CHECK_EQUAL(mid.y(), 1.0);
  BOOST_CHECK_EQUAL(mid.m(), 2.0);
}

BOOST_AUTO_TEST_CASE(testEvaluate_3D_measured)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0, 0.0, 0.0);
  controlPoints.emplace_back(2.0, 2.0, 4.0, 8.0);

  BSplineCurve curve(controlPoints, 1);
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

//-- Test B-spline specific properties

BOOST_AUTO_TEST_CASE(testLocalControl)
{
  // Test that B-splines have local control property
  // Moving one control point should only affect nearby portions of the curve
  std::vector<Point> controlPoints;
  for (int i = 0; i <= 6; ++i) {
    controlPoints.emplace_back(static_cast<double>(i), 0.0);
  }

  BSplineCurve curve(controlPoints, 2);

  // Evaluate at the beginning
  Point startBefore = curve.evaluate(0.0);

  // Modify a control point near the end
  curve.setControlPoint(5, Point(5.0, 5.0));

  // Beginning should be unaffected due to local control
  Point startAfter = curve.evaluate(0.0);
  BOOST_CHECK_EQUAL(startBefore.x(), startAfter.x());
  BOOST_CHECK_EQUAL(startBefore.y(), startAfter.y());
}

BOOST_AUTO_TEST_CASE(testConvexHullProperty)
{
  // B-spline curve lies within convex hull of control points
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 2.0);
  controlPoints.emplace_back(3.0, 1.0);
  controlPoints.emplace_back(4.0, 0.0);

  BSplineCurve curve(controlPoints, 2);

  // Sample curve at many points
  for (int i = 0; i <= 100; ++i) {
    auto   bounds = curve.parameterBounds();
    double t      = bounds.first + (bounds.second - bounds.first) * i / 100.0;
    Point  p      = curve.evaluate(t);

    // Check that point is within bounding box of control points
    BOOST_CHECK(p.x() >= 0.0 && p.x() <= 4.0);
    BOOST_CHECK(p.y() >= 0.0 && p.y() <= 2.0);
  }
}

//-- Template tests (similar to LineString)

template <typename Derived>
void
testIsInstanceOf(const Geometry &geom)
{
  BOOST_CHECK(geom.is<Derived>());
}

BOOST_AUTO_TEST_CASE(isBSplineCurve)
{
  BSplineCurve const curve;
  testIsInstanceOf<BSplineCurve>(curve);
}

//-- Envelope tests

BOOST_AUTO_TEST_CASE(testEnvelope_empty)
{
  BSplineCurve emptyCurve;
  BOOST_CHECK(emptyCurve.envelope().isEmpty());
}

BOOST_AUTO_TEST_CASE(testEnvelope_2D)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(1.0, 5.0);
  controlPoints.emplace_back(3.0, 9.0);
  controlPoints.emplace_back(2.0, 7.0);

  BSplineCurve   curve(controlPoints, 2);
  Envelope const box = curve.envelope();

  BOOST_CHECK(!box.isEmpty());
  BOOST_CHECK(!box.is3D());

  // B-spline stays within convex hull
  BOOST_CHECK(box.xMin() >= 1.0);
  BOOST_CHECK(box.xMax() <= 3.0);
  BOOST_CHECK(box.yMin() >= 5.0);
  BOOST_CHECK(box.yMax() <= 9.0);
}

BOOST_AUTO_TEST_CASE(testEnvelope_3D)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(1.0, 5.0, 11.0);
  controlPoints.emplace_back(3.0, 9.0, 17.0);
  controlPoints.emplace_back(2.0, 7.0, 15.0);

  BSplineCurve   curve(controlPoints, 2);
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

//-- Basis functions test

BOOST_AUTO_TEST_CASE(testEvaluateWithBasisFunctions)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 2.0);
  controlPoints.emplace_back(3.0, 2.0);
  controlPoints.emplace_back(4.0, 0.0);

  BSplineCurve curve(controlPoints, 2);

  auto [point, basis] = curve.evaluateWithBasisFunctions(0.5);

  // Check that basis functions sum to 1 (partition of unity)
  double sum = 0.0;
  for (double b : basis) {
    sum += b;
  }
  BOOST_CHECK_CLOSE(sum, 1.0, 0.0001);

  // Check that all basis functions are non-negative
  for (double b : basis) {
    BOOST_CHECK(b >= 0.0);
  }

  // Check that we have the right number of basis functions
  BOOST_CHECK_EQUAL(basis.size(), 3U); // degree + 1
}

//-- String representation tests
BOOST_AUTO_TEST_CASE(testWkt_2D)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 2.0);
  controlPoints.emplace_back(3.0, 2.0);
  controlPoints.emplace_back(4.0, 0.0);

  BSplineCurve const curve(controlPoints, 2);
  std::string        expectedWkt = "BSPLINECURVE ((0 0,1 2,3 2,4 0),2)";
  BOOST_CHECK_EQUAL(curve.asText(0), expectedWkt);

  std::unique_ptr<Geometry> readCurve = io::readWkt(expectedWkt);
  BOOST_REQUIRE(readCurve->is<BSplineCurve>());
  BOOST_CHECK_EQUAL(readCurve->as<BSplineCurve>().numControlPoints(), 4U);
  BOOST_CHECK_EQUAL(readCurve->as<BSplineCurve>().degree(), 2U);
}

BOOST_AUTO_TEST_CASE(testWkt_3D)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0, 0.0);
  controlPoints.emplace_back(1.0, 2.0, 1.0);
  controlPoints.emplace_back(3.0, 2.0, 1.0);
  controlPoints.emplace_back(4.0, 0.0, 0.0);

  BSplineCurve const curve(controlPoints, 2);
  std::string expectedWkt = "BSPLINECURVE Z ((0 0 0,1 2 1,3 2 1,4 0 0),2)";
  BOOST_CHECK_EQUAL(curve.asText(0), expectedWkt);

  std::unique_ptr<Geometry> readCurve = io::readWkt(expectedWkt);
  BOOST_REQUIRE(readCurve->is<BSplineCurve>());
  BOOST_CHECK(readCurve->is3D());
  BOOST_CHECK_EQUAL(readCurve->as<BSplineCurve>().numControlPoints(), 4U);
  BOOST_CHECK_EQUAL(readCurve->as<BSplineCurve>().degree(), 2U);
}

BOOST_AUTO_TEST_CASE(testWkt_Empty)
{
  BSplineCurve const emptyCurve;
  std::string        expectedWkt = "BSPLINECURVE EMPTY";
  BOOST_CHECK_EQUAL(emptyCurve.asText(0), expectedWkt);

  std::unique_ptr<Geometry> readCurve = io::readWkt(expectedWkt);
  BOOST_REQUIRE(readCurve->is<BSplineCurve>());
  BOOST_CHECK(readCurve->isEmpty());
}

BOOST_AUTO_TEST_SUITE_END()
