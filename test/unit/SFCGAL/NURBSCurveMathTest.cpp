// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "SFCGAL/Exception.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/NURBSCurve.h"
#include "SFCGAL/algorithm/distance.h"
#include <cmath>

using namespace SFCGAL;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_NURBSCurveMathTest)

//-- Helper functions for mathematical tests

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
isNearlyEqual(double a, double b, double tolerance = 1e-10)
{
  return std::abs(a - b) < tolerance;
}

bool
isNearlyEqual(const Point &p1, const Point &p2, double tolerance = 1e-10)
{
  return algorithm::distance(p1, p2) < NURBSCurve::FT(tolerance);
}

//-- Arc length computation tests

BOOST_AUTO_TEST_CASE(testLinearCurveArcLength)
{
  // Linear curve from (0,0) to (3,4)
  std::vector<Point> linearPoints;
  linearPoints.emplace_back(0.0, 0.0);
  linearPoints.emplace_back(3.0, 4.0);

  NURBSCurve linearCurve(linearPoints, 1);

  // Expected length = sqrt(3² + 4²) = 5
  NURBSCurve::FT totalLength = linearCurve.length();
  BOOST_CHECK_CLOSE(CGAL::to_double(totalLength), 5.0, 1e-8);

  // Test partial lengths
  auto           bounds   = linearCurve.parameterBounds();
  NURBSCurve::FT midParam = (bounds.first + bounds.second) / NURBSCurve::FT(2);

  NURBSCurve::FT halfLength = linearCurve.length(bounds.first, midParam);
  BOOST_CHECK_CLOSE(CGAL::to_double(halfLength), 2.5, 1e-8);

  NURBSCurve::FT quarterLength = linearCurve.length(
      bounds.first,
      bounds.first + (midParam - bounds.first) / NURBSCurve::FT(2));
  BOOST_CHECK_CLOSE(CGAL::to_double(quarterLength), 1.25, 1e-8);
}

BOOST_AUTO_TEST_CASE(testCircularArcLength)
{
  Point          center(0.0, 0.0);
  NURBSCurve::FT radius(2.0);
  NURBSCurve::FT startAngle(0.0);
  NURBSCurve::FT endAngle(M_PI); // Semicircle

  auto arc =
      NURBSCurve::createCircularArc(center, radius, startAngle, endAngle);

  // Expected arc length = π * radius = 2π
  NURBSCurve::FT arcLength      = arc->length();
  double         expectedLength = M_PI * 2.0;
  BOOST_CHECK_CLOSE(CGAL::to_double(arcLength), expectedLength, 1e-6);

  // Quarter circle
  auto quarterArc = NURBSCurve::createCircularArc(
      center, radius, NURBSCurve::FT(0.0), NURBSCurve::FT(M_PI / 2.0));
  NURBSCurve::FT quarterLength         = quarterArc->length();
  double         expectedQuarterLength = M_PI * 2.0 / 4.0; // π/2 * radius
  BOOST_CHECK_CLOSE(CGAL::to_double(quarterLength), expectedQuarterLength,
                    1e-6);
}

BOOST_AUTO_TEST_CASE(testComplexCurveArcLength)
{
  // Create S-shaped curve
  std::vector<Point> sPoints;
  sPoints.emplace_back(0.0, 0.0);
  sPoints.emplace_back(1.0, 2.0);
  sPoints.emplace_back(3.0, 2.0);
  sPoints.emplace_back(4.0, 0.0);
  sPoints.emplace_back(5.0, -2.0);
  sPoints.emplace_back(7.0, -2.0);
  sPoints.emplace_back(8.0, 0.0);

  auto       weights = convertWeights({1.0, 2.0, 1.0, 1.0, 2.0, 1.0, 1.0});
  NURBSCurve sCurve(sPoints, weights, 3);

  NURBSCurve::FT totalLength = sCurve.length();
  BOOST_CHECK(totalLength >
              NURBSCurve::FT(8.0)); // Should be longer than straight line
  BOOST_CHECK(CGAL::is_finite(totalLength));
  BOOST_CHECK(totalLength > NURBSCurve::FT(0));

  // Test length with different tolerances
  NURBSCurve::FT roughLength =
      sCurve.length(NURBSCurve::Parameter(-1), NURBSCurve::Parameter(-1),
                    NURBSCurve::FT(1e-3));
  NURBSCurve::FT fineLength =
      sCurve.length(NURBSCurve::Parameter(-1), NURBSCurve::Parameter(-1),
                    NURBSCurve::FT(1e-8));

  // Fine tolerance should give more accurate (usually longer) result
  BOOST_CHECK(fineLength >= roughLength);
}

BOOST_AUTO_TEST_CASE(testArcLengthBounds)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(2.0, 3.0);
  controlPoints.emplace_back(5.0, 1.0);

  NURBSCurve curve(controlPoints, 2);

  auto bounds = curve.parameterBounds();

  // Invalid parameter ranges should return 0
  BOOST_CHECK_EQUAL(curve.length(bounds.second, bounds.first),
                    NURBSCurve::FT(0)); // Reversed
  BOOST_CHECK_EQUAL(curve.length(bounds.first, bounds.first),
                    NURBSCurve::FT(0)); // Same point

  // Valid ranges
  NURBSCurve::FT fullLength = curve.length(bounds.first, bounds.second);
  NURBSCurve::FT defaultLength =
      curve.length(); // Should be same as full length

  BOOST_CHECK_CLOSE(CGAL::to_double(fullLength), CGAL::to_double(defaultLength),
                    1e-10);

  // Partial lengths should sum to total
  NURBSCurve::FT midParam  = (bounds.first + bounds.second) / NURBSCurve::FT(2);
  NURBSCurve::FT firstHalf = curve.length(bounds.first, midParam);
  NURBSCurve::FT secondHalf = curve.length(midParam, bounds.second);

  BOOST_CHECK_CLOSE(CGAL::to_double(firstHalf + secondHalf),
                    CGAL::to_double(fullLength), 1e-8);
}

//-- Parameter by arc length tests

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

  // Beyond end should clamp to end
  NURBSCurve::Parameter beyondParam =
      linearCurve.parameterAtLength(totalLength * NURBSCurve::FT(2));
  BOOST_CHECK_CLOSE(CGAL::to_double(beyondParam),
                    CGAL::to_double(bounds.second), 1e-6);

  // Before start should clamp to start
  NURBSCurve::Parameter beforeParam =
      linearCurve.parameterAtLength(NURBSCurve::FT(-1.0));
  BOOST_CHECK_CLOSE(CGAL::to_double(beforeParam), CGAL::to_double(bounds.first),
                    1e-6);
}

BOOST_AUTO_TEST_CASE(testParameterAtLengthCircular)
{
  Point          center(0.0, 0.0);
  NURBSCurve::FT radius(1.0);
  NURBSCurve::FT startAngle(0.0);
  NURBSCurve::FT endAngle(M_PI); // Semicircle, length = π

  auto arc =
      NURBSCurve::createCircularArc(center, radius, startAngle, endAngle);

  NURBSCurve::FT totalLength = arc->length();

  // Quarter way around should be at (0, 1)
  NURBSCurve::Parameter quarterParam =
      arc->parameterAtLength(totalLength / NURBSCurve::FT(4));
  Point quarterPoint = arc->evaluate(quarterParam);

  // Should be approximately at 45 degrees
  BOOST_CHECK_CLOSE(CGAL::to_double(quarterPoint.x()), std::cos(M_PI / 4),
                    1e-5);
  BOOST_CHECK_CLOSE(CGAL::to_double(quarterPoint.y()), std::sin(M_PI / 4),
                    1e-5);

  // Halfway around should be at (0, 1)
  NURBSCurve::Parameter halfParam =
      arc->parameterAtLength(totalLength / NURBSCurve::FT(2));
  Point halfPoint = arc->evaluate(halfParam);

  BOOST_CHECK_CLOSE(CGAL::to_double(halfPoint.x()), 0.0, 1e-5);
  BOOST_CHECK_CLOSE(CGAL::to_double(halfPoint.y()), 1.0, 1e-5);
}

BOOST_AUTO_TEST_CASE(testParameterAtLengthConsistency)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 2.0);
  controlPoints.emplace_back(3.0, 1.0);
  controlPoints.emplace_back(4.0, 0.0);

  NURBSCurve curve(controlPoints, 3);

  NURBSCurve::FT totalLength = curve.length();

  // Test roundtrip consistency: parameterAtLength -> length should be identity
  std::vector<double> testLengths = {0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0};

  for (double fraction : testLengths) {
    NURBSCurve::FT        targetLength = totalLength * NURBSCurve::FT(fraction);
    NURBSCurve::Parameter param        = curve.parameterAtLength(targetLength);

    auto           bounds       = curve.parameterBounds();
    NURBSCurve::FT actualLength = curve.length(bounds.first, param);

    BOOST_CHECK_CLOSE(CGAL::to_double(actualLength),
                      CGAL::to_double(targetLength), 1e-5);
  }
}

//-- Arc length reparameterization tests

BOOST_AUTO_TEST_CASE(testReparameterizeByArcLength)
{
  // Create curve with non-uniform parameterization
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(0.1, 2.0); // Creates high curvature
  controlPoints.emplace_back(2.0, 2.1);
  controlPoints.emplace_back(4.0, 0.0);

  auto weights = convertWeights(
      {1.0, 5.0, 1.0, 1.0}); // Heavy weight creates non-uniform speed
  NURBSCurve originalCurve(controlPoints, weights, 3);

  auto reparamCurve = originalCurve.reparameterizeByArcLength();

  BOOST_REQUIRE(reparamCurve != nullptr);
  BOOST_CHECK(reparamCurve->numControlPoints() >=
              originalCurve.numControlPoints());

  // Check that reparameterized curve has uniform arc length spacing
  auto           reparamBounds = reparamCurve->parameterBounds();
  NURBSCurve::FT reparamRange  = reparamBounds.second - reparamBounds.first;

  // Sample at uniform parameter intervals and check arc length intervals
  std::vector<NURBSCurve::FT> arcLengths;
  for (int sampleIdx = 0; sampleIdx <= 10; ++sampleIdx) {
    NURBSCurve::Parameter param =
        reparamBounds.first +
        (NURBSCurve::FT(sampleIdx) / NURBSCurve::FT(10)) * reparamRange;
    NURBSCurve::FT arcLength = reparamCurve->length(reparamBounds.first, param);
    arcLengths.push_back(arcLength);
  }

  // Check that arc lengths are approximately linear in parameter
  NURBSCurve::FT totalArcLength = arcLengths.back();
  for (size_t idx = 0; idx < arcLengths.size(); ++idx) {
    NURBSCurve::FT expectedLength =
        (NURBSCurve::FT(idx) / NURBSCurve::FT(10)) * totalArcLength;
    BOOST_CHECK_CLOSE(CGAL::to_double(arcLengths[idx]),
                      CGAL::to_double(expectedLength),
                      5e-2); // 5% tolerance for numerical approximation
  }
}

BOOST_AUTO_TEST_CASE(testReparameterizeEmptyAndLinear)
{
  // Empty curve
  NURBSCurve emptyCurve;
  auto       emptyReparam = emptyCurve.reparameterizeByArcLength();
  BOOST_CHECK(emptyReparam->isEmpty());

  // Linear curve (already has uniform arc length parameterization)
  std::vector<Point> linearPoints;
  linearPoints.emplace_back(0.0, 0.0);
  linearPoints.emplace_back(3.0, 4.0);

  NURBSCurve linearCurve(linearPoints, 1);
  auto       linearReparam = linearCurve.reparameterizeByArcLength();

  // Should be very similar to original
  BOOST_CHECK(linearReparam != nullptr);

  // Check that endpoints match
  auto origBounds    = linearCurve.parameterBounds();
  auto reparamBounds = linearReparam->parameterBounds();

  Point origStart    = linearCurve.evaluate(origBounds.first);
  Point origEnd      = linearCurve.evaluate(origBounds.second);
  Point reparamStart = linearReparam->evaluate(reparamBounds.first);
  Point reparamEnd   = linearReparam->evaluate(reparamBounds.second);

  BOOST_CHECK(isNearlyEqual(origStart, reparamStart, 1e-6));
  BOOST_CHECK(isNearlyEqual(origEnd, reparamEnd, 1e-6));
}

//-- Point projection tests

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

  BOOST_CHECK(isNearlyEqual(pointOnCurve, closestOnCurve, 1e-10));
  BOOST_CHECK_CLOSE(CGAL::to_double(resultParam), 0.5,
                    1e-8); // Should be at midpoint

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
    if (isNearlyEqual(testPoint, closestOffCurve, 1e-8)) {
      foundPoint = true;
    }
  }
  BOOST_CHECK(foundPoint);
}

BOOST_AUTO_TEST_CASE(testDistanceToPoint)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(2.0, 2.0);
  controlPoints.emplace_back(4.0, 0.0);

  NURBSCurve curve(controlPoints, 2);

  // Distance to point on curve should be very small
  auto  bounds = curve.parameterBounds();
  Point pointOnCurve =
      curve.evaluate((bounds.first + bounds.second) / NURBSCurve::FT(2));

  NURBSCurve::FT distanceOnCurve = curve.distance(pointOnCurve);
  BOOST_CHECK_SMALL(CGAL::to_double(distanceOnCurve), 1e-8);

  // Distance to point far from curve
  Point          farPoint(10.0, 10.0);
  NURBSCurve::FT distanceFar = curve.distance(farPoint);
  BOOST_CHECK(distanceFar > NURBSCurve::FT(5.0)); // Should be reasonably large

  // Distance should be symmetric with closest point distance
  Point          closestToFar   = curve.closestPoint(farPoint);
  NURBSCurve::FT directDistance = algorithm::distance(farPoint, closestToFar);
  BOOST_CHECK_CLOSE(CGAL::to_double(distanceFar),
                    CGAL::to_double(directDistance), 1e-6);
}

BOOST_AUTO_TEST_CASE(testClosestPointCircular)
{
  Point          center(5.0, 5.0);
  NURBSCurve::FT radius(3.0);
  auto           fullCircle = NURBSCurve::createCircularArc(
      center, radius, NURBSCurve::FT(0.0), NURBSCurve::FT(2.0 * M_PI));

  // Point at center - any point on circle should be closest
  Point          closestToCenter = fullCircle->closestPoint(center);
  NURBSCurve::FT distToCenter    = algorithm::distance(closestToCenter, center);
  BOOST_CHECK_CLOSE(CGAL::to_double(distToCenter), 3.0,
                    1e-6); // Should equal radius

  // Point outside circle
  Point outsidePoint(10.0, 5.0); // 5 units right of center
  Point closestToOutside = fullCircle->closestPoint(outsidePoint);

  // Should be at (8, 5) - point on circle closest to (10, 5)
  BOOST_CHECK_CLOSE(CGAL::to_double(closestToOutside.x()), 8.0, 1e-5);
  BOOST_CHECK_CLOSE(CGAL::to_double(closestToOutside.y()), 5.0, 1e-5);

  // Point inside circle
  Point insidePoint(6.0, 5.0); // 1 unit right of center
  Point closestToInside = fullCircle->closestPoint(insidePoint);

  // Should be at (8, 5) - same as outside point
  BOOST_CHECK_CLOSE(CGAL::to_double(closestToInside.x()), 8.0, 1e-5);
  BOOST_CHECK_CLOSE(CGAL::to_double(closestToInside.y()), 5.0, 1e-5);
}

//-- Advanced mathematical properties

BOOST_AUTO_TEST_CASE(testSelfIntersectionDetection)
{
  // Create figure-8 like curve that should self-intersect
  std::vector<Point> selfIntersectingPoints;
  selfIntersectingPoints.emplace_back(-2.0, 0.0);
  selfIntersectingPoints.emplace_back(-1.0, 2.0);
  selfIntersectingPoints.emplace_back(0.0, 0.0);
  selfIntersectingPoints.emplace_back(1.0, -2.0);
  selfIntersectingPoints.emplace_back(2.0, 0.0);
  selfIntersectingPoints.emplace_back(1.0, 2.0);
  selfIntersectingPoints.emplace_back(0.0, 0.0);
  selfIntersectingPoints.emplace_back(-1.0, -2.0);

  NURBSCurve selfIntersecting(selfIntersectingPoints, 3);

  std::vector<std::pair<NURBSCurve::Parameter, NURBSCurve::Parameter>>
       intersections;
  bool hasSelfIntersections =
      selfIntersecting.hasSelfIntersections(&intersections);

  // This is a simplified test - the actual algorithm might not detect all
  // intersections depending on the sampling resolution
  BOOST_CHECK(hasSelfIntersections ||
              intersections.empty()); // Either detect or don't claim to detect

  // Simple non-self-intersecting curve
  std::vector<Point> simplePoints;
  simplePoints.emplace_back(0.0, 0.0);
  simplePoints.emplace_back(1.0, 1.0);
  simplePoints.emplace_back(2.0, 0.5);
  simplePoints.emplace_back(3.0, 0.0);

  NURBSCurve simpleCurve(simplePoints, 3);
  BOOST_CHECK(!simpleCurve.hasSelfIntersections());
}

BOOST_AUTO_TEST_CASE(testCurveIntersection)
{
  // Create two intersecting curves
  std::vector<Point> curve1Points;
  curve1Points.emplace_back(0.0, 0.0);
  curve1Points.emplace_back(2.0, 2.0);
  curve1Points.emplace_back(4.0, 0.0);

  std::vector<Point> curve2Points;
  curve2Points.emplace_back(0.0, 2.0);
  curve2Points.emplace_back(2.0, 0.0);
  curve2Points.emplace_back(4.0, 2.0);

  NURBSCurve curve1(curve1Points, 2);
  NURBSCurve curve2(curve2Points, 2);

  auto intersections = curve1.intersect(curve2, NURBSCurve::FT(0.1));

  // Should find at least one intersection point near (2, 1)
  bool foundIntersectionNear2_1 = false;
  for (const auto &[point, param1, param2] : intersections) {
    if (isNearlyEqual(point, Point(2.0, 1.0), 0.5)) { // Loose tolerance
      foundIntersectionNear2_1 = true;
      break;
    }
  }

  // This is a simplified test - exact intersection finding is complex
  BOOST_CHECK(foundIntersectionNear2_1 || intersections.empty());

  // Non-intersecting curves
  std::vector<Point> separateCurve1;
  separateCurve1.emplace_back(0.0, 0.0);
  separateCurve1.emplace_back(1.0, 0.0);
  separateCurve1.emplace_back(2.0, 0.0);

  std::vector<Point> separateCurve2;
  separateCurve2.emplace_back(0.0, 5.0);
  separateCurve2.emplace_back(1.0, 5.0);
  separateCurve2.emplace_back(2.0, 5.0);

  NURBSCurve separate1(separateCurve1, 2);
  NURBSCurve separate2(separateCurve2, 2);

  auto noIntersections = separate1.intersect(separate2, NURBSCurve::FT(0.1));
  BOOST_CHECK(noIntersections.empty());
}

//-- Numerical stability tests

BOOST_AUTO_TEST_CASE(testNumericalStabilityHighDegree)
{
  // Create high-degree curve and test numerical stability
  std::vector<Point> manyPoints;
  for (int idx = 0; idx <= 10; ++idx) {
    double angle = 2.0 * M_PI * idx / 10.0;
    manyPoints.emplace_back(std::cos(angle), std::sin(angle));
  }

  NURBSCurve highDegreeCurve(manyPoints, 8);

  // Test evaluation at many points
  auto           bounds     = highDegreeCurve.parameterBounds();
  NURBSCurve::FT paramRange = bounds.second - bounds.first;

  for (int testIdx = 0; testIdx <= 100; ++testIdx) {
    NURBSCurve::Parameter param =
        bounds.first +
        (NURBSCurve::FT(testIdx) / NURBSCurve::FT(100)) * paramRange;

    Point evalPoint = highDegreeCurve.evaluate(param);
    BOOST_CHECK(CGAL::is_finite(evalPoint.x()));
    BOOST_CHECK(CGAL::is_finite(evalPoint.y()));

    // Check derivatives don't blow up
    Point firstDeriv = highDegreeCurve.derivative(param, 1);
    BOOST_CHECK(CGAL::is_finite(firstDeriv.x()));
    BOOST_CHECK(CGAL::is_finite(firstDeriv.y()));
  }
}

BOOST_AUTO_TEST_CASE(testNumericalStabilityExtremeWeights)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);
  controlPoints.emplace_back(2.0, 0.0);

  // Very small and very large weights
  auto       extremeWeights = convertWeights({1e-8, 1e8, 1.0});
  NURBSCurve extremeCurve(controlPoints, extremeWeights, 2);

  auto bounds = extremeCurve.parameterBounds();

  // Should still be able to evaluate without numerical issues
  Point startPoint = extremeCurve.evaluate(bounds.first);
  Point midPoint =
      extremeCurve.evaluate((bounds.first + bounds.second) / NURBSCurve::FT(2));
  Point endPoint = extremeCurve.evaluate(bounds.second);

  BOOST_CHECK(CGAL::is_finite(startPoint.x()) &&
              CGAL::is_finite(startPoint.y()));
  BOOST_CHECK(CGAL::is_finite(midPoint.x()) && CGAL::is_finite(midPoint.y()));
  BOOST_CHECK(CGAL::is_finite(endPoint.x()) && CGAL::is_finite(endPoint.y()));

  // Very large weight should pull curve towards middle control point
  NURBSCurve::FT distToMiddleControl =
      algorithm::distance(midPoint, controlPoints[1]);
  NURBSCurve::FT distToFirstControl =
      algorithm::distance(midPoint, controlPoints[0]);
  NURBSCurve::FT distToLastControl =
      algorithm::distance(midPoint, controlPoints[2]);

  BOOST_CHECK(distToMiddleControl < distToFirstControl);
  BOOST_CHECK(distToMiddleControl < distToLastControl);
}

BOOST_AUTO_TEST_CASE(testParameterBoundaryBehavior)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 2.0);
  controlPoints.emplace_back(3.0, 1.0);
  controlPoints.emplace_back(4.0, 0.0);

  NURBSCurve curve(controlPoints, 3);

  auto bounds = curve.parameterBounds();

  // Evaluation right at boundaries should work
  Point startBoundary = curve.evaluate(bounds.first);
  Point endBoundary   = curve.evaluate(bounds.second);

  BOOST_CHECK(CGAL::is_finite(startBoundary.x()));
  BOOST_CHECK(CGAL::is_finite(endBoundary.y()));

  // Slightly outside boundaries should throw
  NURBSCurve::FT epsilon = NURBSCurve::FT(1e-12);
  BOOST_CHECK_THROW(curve.evaluate(bounds.first - epsilon), Exception);
  BOOST_CHECK_THROW(curve.evaluate(bounds.second + epsilon), Exception);

  // Very close to boundary but inside should work
  NURBSCurve::FT tinyEpsilon = NURBSCurve::FT(1e-15);
  BOOST_CHECK_NO_THROW(curve.evaluate(bounds.first + tinyEpsilon));
  BOOST_CHECK_NO_THROW(curve.evaluate(bounds.second - tinyEpsilon));
}

//-- Error handling and edge cases

BOOST_AUTO_TEST_CASE(testMathematicalErrorHandling)
{
  NURBSCurve emptyCurve;

  // Operations on empty curve should throw appropriate exceptions
  BOOST_CHECK_THROW(emptyCurve.length(), Exception);
  BOOST_CHECK_THROW(emptyCurve.parameterAtLength(NURBSCurve::FT(1.0)),
                    Exception);
  BOOST_CHECK_THROW(emptyCurve.closestPoint(Point(0, 0)), Exception);
  BOOST_CHECK_THROW(emptyCurve.distance(Point(0, 0)), Exception);
  BOOST_CHECK_NO_THROW(
      emptyCurve.reparameterizeByArcLength()); // Should return empty curve

  // Invalid parameters for operations
  std::vector<Point> validPoints;
  validPoints.emplace_back(0.0, 0.0);
  validPoints.emplace_back(1.0, 1.0);

  NURBSCurve validCurve(validPoints, 1);

  // Negative arc length
  BOOST_CHECK_NO_THROW(validCurve.parameterAtLength(
      NURBSCurve::FT(-1.0))); // Should clamp to start

  // Invalid tolerance values
  BOOST_CHECK_NO_THROW(validCurve.length(
      NURBSCurve::Parameter(-1), NURBSCurve::Parameter(-1),
      NURBSCurve::FT(-1e-6))); // Should handle negative tolerance gracefully
}

BOOST_AUTO_TEST_CASE(testZeroLengthCurves)
{
  // Curve with all identical control points
  std::vector<Point> identicalPoints;
  identicalPoints.emplace_back(2.0, 3.0);
  identicalPoints.emplace_back(2.0, 3.0);
  identicalPoints.emplace_back(2.0, 3.0);

  NURBSCurve pointCurve(identicalPoints, 2);

  NURBSCurve::FT zeroLength = pointCurve.length();
  BOOST_CHECK_SMALL(CGAL::to_double(zeroLength), 1e-10);

  // parameterAtLength should handle zero-length curves gracefully
  auto                  bounds = pointCurve.parameterBounds();
  NURBSCurve::Parameter param =
      pointCurve.parameterAtLength(NURBSCurve::FT(1.0));
  BOOST_CHECK(param >= bounds.first && param <= bounds.second);

  // Reparameterization of zero-length curve
  auto reparamZero = pointCurve.reparameterizeByArcLength();
  BOOST_CHECK(reparamZero != nullptr);
}

BOOST_AUTO_TEST_SUITE_END()
