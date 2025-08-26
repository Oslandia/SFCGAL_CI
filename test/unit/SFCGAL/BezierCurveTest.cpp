// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "SFCGAL/BezierCurve.h"
#include "SFCGAL/Envelope.h"
#include "SFCGAL/Exception.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/io/wkt.h"

using namespace SFCGAL;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_BezierCurveTest)

//-- Constructor tests

/// BezierCurve() ;
BOOST_AUTO_TEST_CASE(defaultConstructor)
{
  BezierCurve const curve;
  BOOST_CHECK(curve.isEmpty());
  BOOST_CHECK(!curve.is3D());
  BOOST_CHECK(!curve.isMeasured());
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 0U);
  BOOST_CHECK_EQUAL(curve.degree(), 0U);
}

/// BezierCurve(const std::vector<Point> &controlPoints) ;
BOOST_AUTO_TEST_CASE(constructorFromVector)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 2.0);
  controlPoints.emplace_back(3.0, 2.0);
  controlPoints.emplace_back(4.0, 0.0);

  BezierCurve curve(controlPoints);
  BOOST_REQUIRE_EQUAL(curve.numControlPoints(), 4U);
  BOOST_CHECK_EQUAL(curve.degree(), 3U);
  BOOST_CHECK(!curve.isEmpty());
  BOOST_CHECK_EQUAL(curve.controlPointAt(0).x(), 0.0);
  BOOST_CHECK_EQUAL(curve.controlPointAt(0).y(), 0.0);
  BOOST_CHECK_EQUAL(curve.controlPointAt(3).x(), 4.0);
  BOOST_CHECK_EQUAL(curve.controlPointAt(3).y(), 0.0);
}

/// BezierCurve(const BezierCurve& other) ;
/// BezierCurve& operator=(const BezierCurve& other) ;
BOOST_AUTO_TEST_CASE(copyConstructorAndAssignment)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);

  BezierCurve original(controlPoints);
  BezierCurve copy(original);
  BezierCurve assigned;
  assigned = original;

  BOOST_CHECK_EQUAL(copy.numControlPoints(), 2U);
  BOOST_CHECK_EQUAL(copy.degree(), 1U);
  BOOST_CHECK_EQUAL(assigned.numControlPoints(), 2U);
  BOOST_CHECK_EQUAL(assigned.degree(), 1U);

  // Verify independence
  copy.addControlPoint(Point(2.0, 2.0));
  BOOST_CHECK_EQUAL(original.numControlPoints(), 2U);
  BOOST_CHECK_EQUAL(copy.numControlPoints(), 3U);
}

//-- Basic operations tests

/// void clear() ;
BOOST_AUTO_TEST_CASE(testClear)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);

  BezierCurve curve(controlPoints);
  BOOST_CHECK(!curve.isEmpty());

  curve.clear();
  BOOST_CHECK(curve.isEmpty());
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 0U);
  BOOST_CHECK_EQUAL(curve.degree(), 0U);
}

/// void addControlPoint(const Point &point) ;
BOOST_AUTO_TEST_CASE(testAddControlPoint)
{
  BezierCurve curve;
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 0U);

  curve.addControlPoint(Point(0.0, 0.0));
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 1U);
  BOOST_CHECK_EQUAL(curve.degree(), 0U);

  curve.addControlPoint(Point(1.0, 1.0));
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 2U);
  BOOST_CHECK_EQUAL(curve.degree(), 1U);

  curve.addControlPoint(Point(2.0, 0.0));
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 3U);
  BOOST_CHECK_EQUAL(curve.degree(), 2U);
}

BOOST_AUTO_TEST_CASE(testAddControlPoint_dimensionalConsistency)
{
  BezierCurve curve;
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

  BezierCurve curve(controlPoints);
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 3U);

  curve.removeControlPoint(1);
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 2U);
  BOOST_CHECK_EQUAL(curve.controlPointAt(0).x(), 0.0);
  BOOST_CHECK_EQUAL(curve.controlPointAt(1).x(), 2.0);
}

/// const Point& controlPointAt(size_t index) const ;
/// Point& controlPointAt(size_t index) ;
/// void setControlPoint(size_t index, const Point &point) ;
BOOST_AUTO_TEST_CASE(testControlPointAccess)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);

  BezierCurve curve(controlPoints);

  // Test const access
  const BezierCurve &constCurve = curve;
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

  BezierCurve        curve(originalPoints);
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
  BezierCurve const curve;
  BOOST_CHECK_EQUAL(curve.geometryTypeId(), TYPE_BEZIERCURVE);
}

/// std::string geometryType() const ;
BOOST_AUTO_TEST_CASE(testGeometryType)
{
  BezierCurve const curve;
  BOOST_CHECK_EQUAL(curve.geometryType(), "BezierCurve");
}

/// BezierCurve* clone() const ;
BOOST_AUTO_TEST_CASE(testClone)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);

  BezierCurve                  curve(controlPoints);
  std::unique_ptr<BezierCurve> clone(curve.clone());

  BOOST_REQUIRE(clone != nullptr);
  BOOST_CHECK_EQUAL(clone->numControlPoints(), 2U);
  BOOST_CHECK_EQUAL(clone->controlPointAt(0).x(), 0.0);
  BOOST_CHECK_EQUAL(clone->controlPointAt(1).y(), 1.0);
}

/// bool isEmpty() const ;
BOOST_AUTO_TEST_CASE(testIsEmpty)
{
  BezierCurve emptyCurve;
  BOOST_CHECK(emptyCurve.isEmpty());

  BezierCurve nonEmptyCurve;
  nonEmptyCurve.addControlPoint(Point(0.0, 0.0));
  BOOST_CHECK(!nonEmptyCurve.isEmpty());
}

/// bool is3D() const ;
BOOST_AUTO_TEST_CASE(testIs3D)
{
  BezierCurve curve2D;
  curve2D.addControlPoint(Point(0.0, 0.0));
  BOOST_CHECK(!curve2D.is3D());

  BezierCurve curve3D;
  curve3D.addControlPoint(Point(0.0, 0.0, 0.0));
  BOOST_CHECK(curve3D.is3D());

  BezierCurve emptyCurve;
  BOOST_CHECK(!emptyCurve.is3D());
}

/// bool isMeasured() const ;
BOOST_AUTO_TEST_CASE(testIsMeasured)
{
  BezierCurve curveRegular;
  curveRegular.addControlPoint(Point(0.0, 0.0));
  BOOST_CHECK(!curveRegular.isMeasured());

  BezierCurve curveMeasured;
  curveMeasured.addControlPoint(Point(0.0, 0.0, Kernel::FT(NaN()), 1.0));
  BOOST_CHECK(curveMeasured.isMeasured());

  BezierCurve emptyCurve;
  BOOST_CHECK(!emptyCurve.isMeasured());
}

//-- Curve interface tests

/// Point evaluate(double parameter) const ;
BOOST_AUTO_TEST_CASE(testEvaluate_linear)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(2.0, 2.0);

  BezierCurve curve(controlPoints);

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

  BezierCurve curve(controlPoints);

  Point start = curve.evaluate(0.0);
  BOOST_CHECK_EQUAL(start.x(), 0.0);
  BOOST_CHECK_EQUAL(start.y(), 0.0);

  Point end = curve.evaluate(1.0);
  BOOST_CHECK_EQUAL(end.x(), 2.0);
  BOOST_CHECK_EQUAL(end.y(), 0.0);

  Point mid = curve.evaluate(0.5);
  BOOST_CHECK_EQUAL(mid.x(), 1.0);
  BOOST_CHECK_EQUAL(mid.y(), 1.0); // Should be at control point height
}

BOOST_AUTO_TEST_CASE(testEvaluate_empty)
{
  BezierCurve emptyCurve;
  BOOST_CHECK_THROW(auto result = emptyCurve.evaluate(0.5), Exception);
}

BOOST_AUTO_TEST_CASE(testEvaluate_parameterBounds)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);

  BezierCurve curve(controlPoints);

  BOOST_CHECK_THROW(auto result1 = curve.evaluate(-0.1), Exception);
  BOOST_CHECK_THROW(auto result2 = curve.evaluate(1.1), Exception);
  BOOST_CHECK_NO_THROW(auto result3 = curve.evaluate(0.0));
  BOOST_CHECK_NO_THROW(auto result4 = curve.evaluate(1.0));
}

/// Point derivative(double parameter, unsigned int order = 1) const ;
BOOST_AUTO_TEST_CASE(testDerivative_linear)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(2.0, 4.0);

  BezierCurve curve(controlPoints);

  Point derivative = curve.derivative(0.5, 1);
  BOOST_CHECK_EQUAL(derivative.x(), 2.0);
  BOOST_CHECK_EQUAL(derivative.y(), 4.0);
}

BOOST_AUTO_TEST_CASE(testDerivative_quadratic)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 2.0);
  controlPoints.emplace_back(2.0, 0.0);

  BezierCurve curve(controlPoints);

  Point firstDerivative = curve.derivative(0.0, 1);
  // First derivative at t=0 for quadratic Bezier
  BOOST_CHECK_EQUAL(firstDerivative.x(), 2.0); // 2 * (P1 - P0)
  BOOST_CHECK_EQUAL(firstDerivative.y(), 4.0);

  Point secondDerivative = curve.derivative(0.5, 2);
  // Second derivative should be constant for quadratic
  BOOST_CHECK_EQUAL(secondDerivative.x(), 0.0); // 2 * (P0 - 2*P1 + P2)
  BOOST_CHECK_EQUAL(secondDerivative.y(), -8.0);
}

BOOST_AUTO_TEST_CASE(testDerivative_orderZero)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);

  BezierCurve curve(controlPoints);

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

  BezierCurve curve(controlPoints);

  // Order higher than degree should return zero vector
  Point highDerivative = curve.derivative(0.5, 5);
  BOOST_CHECK_EQUAL(highDerivative.x(), 0.0);
  BOOST_CHECK_EQUAL(highDerivative.y(), 0.0);
}

/// std::unique_ptr<LineString> toLineString(unsigned int numSegments = 32)
/// const ;
BOOST_AUTO_TEST_CASE(testToLineString_empty)
{
  BezierCurve                 emptyCurve;
  std::unique_ptr<LineString> lineString = emptyCurve.toLineString(10);
  BOOST_CHECK(lineString->isEmpty());
}

BOOST_AUTO_TEST_CASE(testToLineString_linear)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(2.0, 2.0);

  BezierCurve                 curve(controlPoints);
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

  BezierCurve curve(controlPoints);
  BOOST_CHECK_THROW(auto result = curve.toLineString(1), Exception);
}

/// std::pair<double, double> parameterBounds() const ;
BOOST_AUTO_TEST_CASE(testParameterBounds)
{
  BezierCurve               curve;
  std::pair<double, double> bounds = curve.parameterBounds();
  BOOST_CHECK_EQUAL(bounds.first, 0.0);
  BOOST_CHECK_EQUAL(bounds.second, 1.0);
}

/// bool isValid() const ;
BOOST_AUTO_TEST_CASE(testIsValid_empty)
{
  BezierCurve emptyCurve;
  BOOST_CHECK(emptyCurve.isValid());
}

BOOST_AUTO_TEST_CASE(testIsValid_consistent)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);
  controlPoints.emplace_back(2.0, 0.0);

  BezierCurve curve(controlPoints);
  BOOST_CHECK(curve.isValid());
}

BOOST_AUTO_TEST_CASE(testIsValid_inconsistentDimensions)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);      // 2D
  controlPoints.emplace_back(1.0, 1.0, 1.0); // 3D

  BezierCurve curve(controlPoints);
  BOOST_CHECK(!curve.isValid());
}

//-- BezierCurve specific tests

/// unsigned int degree() const ;
BOOST_AUTO_TEST_CASE(testDegree)
{
  BezierCurve emptyCurve;
  BOOST_CHECK_EQUAL(emptyCurve.degree(), 0U);

  BezierCurve linear;
  linear.addControlPoint(Point(0.0, 0.0));
  linear.addControlPoint(Point(1.0, 1.0));
  BOOST_CHECK_EQUAL(linear.degree(), 1U);

  BezierCurve quadratic;
  quadratic.addControlPoint(Point(0.0, 0.0));
  quadratic.addControlPoint(Point(1.0, 2.0));
  quadratic.addControlPoint(Point(2.0, 0.0));
  BOOST_CHECK_EQUAL(quadratic.degree(), 2U);
}

/// size_t numControlPoints() const ;
BOOST_AUTO_TEST_CASE(testNumControlPoints)
{
  BezierCurve curve;
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 0U);

  curve.addControlPoint(Point(0.0, 0.0));
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 1U);

  curve.addControlPoint(Point(1.0, 1.0));
  BOOST_CHECK_EQUAL(curve.numControlPoints(), 2U);
}

/// std::vector<std::vector<Point>> evaluateWithSteps(double parameter) const ;
BOOST_AUTO_TEST_CASE(testEvaluateWithSteps_empty)
{
  BezierCurve                     emptyCurve;
  std::vector<std::vector<Point>> steps = emptyCurve.evaluateWithSteps(0.5);
  BOOST_CHECK(steps.empty());
}

BOOST_AUTO_TEST_CASE(testEvaluateWithSteps_linear)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(2.0, 2.0);

  BezierCurve                     curve(controlPoints);
  std::vector<std::vector<Point>> steps = curve.evaluateWithSteps(0.5);

  BOOST_CHECK_EQUAL(steps.size(), 2U);    // Initial points + final point
  BOOST_CHECK_EQUAL(steps[0].size(), 2U); // Two control points
  BOOST_CHECK_EQUAL(steps[1].size(), 1U); // Final evaluated point

  // Final point should be at (1,1)
  BOOST_CHECK_EQUAL(steps[1][0].x(), 1.0);
  BOOST_CHECK_EQUAL(steps[1][0].y(), 1.0);
}

BOOST_AUTO_TEST_CASE(testEvaluateWithSteps_quadratic)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 2.0);
  controlPoints.emplace_back(2.0, 0.0);

  BezierCurve                     curve(controlPoints);
  std::vector<std::vector<Point>> steps = curve.evaluateWithSteps(0.5);

  BOOST_CHECK_EQUAL(steps.size(),
                    3U); // 3 levels: original, intermediate, final
  BOOST_CHECK_EQUAL(steps[0].size(), 3U); // Three original control points
  BOOST_CHECK_EQUAL(steps[1].size(), 2U); // Two intermediate points
  BOOST_CHECK_EQUAL(steps[2].size(), 1U); // One final point
}

//-- 3D and Measured coordinate tests

BOOST_AUTO_TEST_CASE(testEvaluate_3D)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0, 2.0);

  BezierCurve curve(controlPoints);
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

  BezierCurve curve(controlPoints);
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

  BezierCurve curve(controlPoints);
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

//-- Error handling tests

BOOST_AUTO_TEST_CASE(testEvaluate_singlePoint)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(1.0, 2.0);

  BezierCurve curve(controlPoints);
  Point       result = curve.evaluate(0.0);
  BOOST_CHECK_EQUAL(result.x(), 1.0);
  BOOST_CHECK_EQUAL(result.y(), 2.0);

  Point resultMid = curve.evaluate(0.5);
  BOOST_CHECK_EQUAL(resultMid.x(), 1.0);
  BOOST_CHECK_EQUAL(resultMid.y(), 2.0);

  Point resultEnd = curve.evaluate(1.0);
  BOOST_CHECK_EQUAL(resultEnd.x(), 1.0);
  BOOST_CHECK_EQUAL(resultEnd.y(), 2.0);
}

//-- Template tests (similar to LineString)

template <typename Derived>
void
testIsInstanceOf(const Geometry &geom)
{
  BOOST_CHECK(geom.is<Derived>());
}

BOOST_AUTO_TEST_CASE(isBezierCurve)
{
  BezierCurve const curve;
  testIsInstanceOf<BezierCurve>(curve);
}

//-- Envelope tests

BOOST_AUTO_TEST_CASE(testEnvelope_empty)
{
  BezierCurve emptyCurve;
  BOOST_CHECK(emptyCurve.envelope().isEmpty());
}

BOOST_AUTO_TEST_CASE(testEnvelope_2D)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(1.0, 5.0);
  controlPoints.emplace_back(3.0, 9.0);
  controlPoints.emplace_back(2.0, 7.0);

  BezierCurve    curve(controlPoints);
  Envelope const box = curve.envelope();

  BOOST_CHECK(!box.isEmpty());
  BOOST_CHECK(!box.is3D());

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

  BezierCurve    curve(controlPoints);
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

//-- String representation tests
BOOST_AUTO_TEST_CASE(testWkt_2D)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 2.0);
  controlPoints.emplace_back(3.0, 2.0);
  controlPoints.emplace_back(4.0, 0.0);

  BezierCurve const curve(controlPoints);
  std::string       expectedWkt = "BEZIERCURVE (0 0,1 2,3 2,4 0)";
  BOOST_CHECK_EQUAL(curve.asText(0), expectedWkt);

  std::unique_ptr<Geometry> readCurve = io::readWkt(expectedWkt);
  BOOST_REQUIRE(readCurve->is<BezierCurve>());
  BOOST_CHECK_EQUAL(readCurve->as<BezierCurve>().numControlPoints(), 4U);
}

BOOST_AUTO_TEST_CASE(testWkt_3D)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0, 0.0);
  controlPoints.emplace_back(1.0, 2.0, 1.0);
  controlPoints.emplace_back(3.0, 2.0, 1.0);
  controlPoints.emplace_back(4.0, 0.0, 0.0);

  BezierCurve const curve(controlPoints);
  std::string       expectedWkt = "BEZIERCURVE Z (0 0 0,1 2 1,3 2 1,4 0 0)";
  BOOST_CHECK_EQUAL(curve.asText(0), expectedWkt);

  std::unique_ptr<Geometry> readCurve = io::readWkt(expectedWkt);
  BOOST_REQUIRE(readCurve->is<BezierCurve>());
  BOOST_CHECK(readCurve->is3D());
  BOOST_CHECK_EQUAL(readCurve->as<BezierCurve>().numControlPoints(), 4U);
}

BOOST_AUTO_TEST_CASE(testWkt_Empty)
{
  BezierCurve const emptyCurve;
  std::string       expectedWkt = "BEZIERCURVE EMPTY";
  BOOST_CHECK_EQUAL(emptyCurve.asText(0), expectedWkt);

  std::unique_ptr<Geometry> readCurve = io::readWkt(expectedWkt);
  BOOST_REQUIRE(readCurve->is<BezierCurve>());
  BOOST_CHECK(readCurve->isEmpty());
}

BOOST_AUTO_TEST_SUITE_END()
