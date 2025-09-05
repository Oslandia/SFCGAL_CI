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

BOOST_AUTO_TEST_SUITE(SFCGAL_NURBSCurveManipulationTest)

//-- Helper functions for manipulation tests

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
isNearlyEqual(const Point &p1, const Point &p2, double tolerance = 1e-10)
{
  return algorithm::distance(p1, p2) < NURBSCurve::FT(tolerance);
}

std::vector<Point>
createTestCurvePoints()
{
  std::vector<Point> points;
  points.emplace_back(0.0, 0.0);
  points.emplace_back(2.0, 4.0);
  points.emplace_back(6.0, 4.0);
  points.emplace_back(8.0, 0.0);
  return points;
}

std::vector<Point>
createCircularPoints(double radius = 1.0, size_t numPoints = 8)
{
  std::vector<Point> points;
  for (size_t pointIdx = 0; pointIdx < numPoints; ++pointIdx) {
    double angle = 2.0 * M_PI * pointIdx / numPoints;
    points.emplace_back(radius * std::cos(angle), radius * std::sin(angle));
  }
  return points;
}

//-- Curve splitting tests

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

BOOST_AUTO_TEST_CASE(testSplitQuadraticCurve)
{
  std::vector<Point> quadraticPoints;
  quadraticPoints.emplace_back(0.0, 0.0);
  quadraticPoints.emplace_back(2.0, 4.0);
  quadraticPoints.emplace_back(4.0, 0.0);

  auto       weights = convertWeights({1.0, 2.0, 1.0});
  NURBSCurve quadraticCurve(quadraticPoints, weights, 2);

  auto                  bounds = quadraticCurve.parameterBounds();
  NURBSCurve::Parameter splitParam =
      bounds.first +
      (bounds.second - bounds.first) * NURBSCurve::FT(0.3); // Split at 30%

  auto [firstPart, secondPart] = quadraticCurve.split(splitParam);

  BOOST_REQUIRE(firstPart != nullptr);
  BOOST_REQUIRE(secondPart != nullptr);

  // Check that split curves cover the original range
  Point originalSplitPoint = quadraticCurve.evaluate(splitParam);

  auto firstPartBounds  = firstPart->parameterBounds();
  auto secondPartBounds = secondPart->parameterBounds();

  Point firstPartEnd    = firstPart->evaluate(firstPartBounds.second);
  Point secondPartStart = secondPart->evaluate(secondPartBounds.first);

  std::cout << originalSplitPoint.asText(0) << "\n";
  std::cout << firstPartEnd.asText(0) << "\n";
  std::cout << secondPartStart.asText(0) << "\n";
  BOOST_CHECK(isNearlyEqual(originalSplitPoint, firstPartEnd, 1e-6));
  BOOST_CHECK(isNearlyEqual(originalSplitPoint, secondPartStart, 1e-6));
}

BOOST_AUTO_TEST_CASE(testSplitAtBoundaries)
{
  auto       controlPoints = createTestCurvePoints();
  NURBSCurve curve(controlPoints, 3);

  auto bounds = curve.parameterBounds();

  // Split at start should give empty first curve and full second curve
  auto [emptyFirst, fullSecond] = curve.split(bounds.first);

  BOOST_REQUIRE(emptyFirst != nullptr);
  BOOST_REQUIRE(fullSecond != nullptr);
  BOOST_CHECK(emptyFirst->isEmpty());
  BOOST_CHECK(!fullSecond->isEmpty());

  // Split at end should give full first curve and empty second curve
  auto [fullFirst, emptySecond] = curve.split(bounds.second);

  BOOST_REQUIRE(fullFirst != nullptr);
  BOOST_REQUIRE(emptySecond != nullptr);
  BOOST_CHECK(!fullFirst->isEmpty());
  BOOST_CHECK(emptySecond->isEmpty());

  // Split outside bounds should clamp to boundaries
  auto [clampedFirst, clampedSecond] =
      curve.split(bounds.second + NURBSCurve::FT(1.0));
  BOOST_CHECK(!clampedFirst->isEmpty());
  BOOST_CHECK(clampedSecond->isEmpty());
}

BOOST_AUTO_TEST_CASE(testSplitEmpty)
{
  NURBSCurve emptyCurve;

  auto [first, second] = emptyCurve.split(NURBSCurve::Parameter(0.5));

  BOOST_REQUIRE(first != nullptr);
  BOOST_REQUIRE(second != nullptr);
  BOOST_CHECK(first->isEmpty());
  BOOST_CHECK(second->isEmpty());
}

//-- Subcurve extraction tests

BOOST_AUTO_TEST_CASE(testSubcurveLinear)
{
  std::vector<Point> linearPoints;
  linearPoints.emplace_back(0.0, 0.0);
  linearPoints.emplace_back(6.0, 8.0);

  NURBSCurve linearCurve(linearPoints, 1);

  auto                  bounds = linearCurve.parameterBounds();
  NURBSCurve::Parameter startSub =
      bounds.first + (bounds.second - bounds.first) * NURBSCurve::FT(0.25);
  NURBSCurve::Parameter endSub =
      bounds.first + (bounds.second - bounds.first) * NURBSCurve::FT(0.75);

  auto subcurve = linearCurve.subcurve(startSub, endSub);

  BOOST_REQUIRE(subcurve != nullptr);
  BOOST_CHECK(!subcurve->isEmpty());

  // Check subcurve endpoints
  auto  subBounds = subcurve->parameterBounds();
  Point subStart  = subcurve->evaluate(subBounds.first);
  Point subEnd    = subcurve->evaluate(subBounds.second);

  Point originalStart = linearCurve.evaluate(startSub);
  Point originalEnd   = linearCurve.evaluate(endSub);

  BOOST_CHECK(isNearlyEqual(subStart, originalStart, 1e-8));
  BOOST_CHECK(isNearlyEqual(subEnd, originalEnd, 1e-8));

  // Subcurve of linear should be linear
  BOOST_CHECK(subcurve->isLinear());
}

BOOST_AUTO_TEST_CASE(testSubcurveComplex)
{
  auto       controlPoints = createTestCurvePoints();
  auto       weights       = convertWeights({1.0, 2.5, 1.5, 1.0});
  NURBSCurve complexCurve(controlPoints, weights, 3);

  auto                  bounds = complexCurve.parameterBounds();
  NURBSCurve::Parameter startSub =
      bounds.first + (bounds.second - bounds.first) * NURBSCurve::FT(0.2);
  NURBSCurve::Parameter endSub =
      bounds.first + (bounds.second - bounds.first) * NURBSCurve::FT(0.8);

  auto subcurve = complexCurve.subcurve(startSub, endSub);

  BOOST_REQUIRE(subcurve != nullptr);
  BOOST_CHECK(!subcurve->isEmpty());

  // Check endpoint correspondence
  auto  subBounds = subcurve->parameterBounds();
  Point subStart  = subcurve->evaluate(subBounds.first);
  Point subEnd    = subcurve->evaluate(subBounds.second);

  Point originalStart = complexCurve.evaluate(startSub);
  Point originalEnd   = complexCurve.evaluate(endSub);

  BOOST_CHECK(isNearlyEqual(subStart, originalStart, 1e-6));
  BOOST_CHECK(isNearlyEqual(subEnd, originalEnd, 1e-6));
}

BOOST_AUTO_TEST_CASE(testSubcurveInvalidRange)
{
  auto       controlPoints = createTestCurvePoints();
  NURBSCurve curve(controlPoints, 2);

  auto bounds = curve.parameterBounds();

  // Reversed range should return empty curve
  auto reversedSub = curve.subcurve(bounds.second, bounds.first);
  BOOST_REQUIRE(reversedSub != nullptr);
  BOOST_CHECK(reversedSub->isEmpty());

  // Same parameter should return empty curve
  auto sameSub = curve.subcurve(bounds.first, bounds.first);
  BOOST_REQUIRE(sameSub != nullptr);
  BOOST_CHECK(sameSub->isEmpty());

  // Full range should be similar to original
  auto fullSub = curve.subcurve(bounds.first, bounds.second);
  BOOST_REQUIRE(fullSub != nullptr);
  BOOST_CHECK_EQUAL(fullSub->degree(), curve.degree());

  // Parameters outside bounds should be clamped
  auto clampedSub = curve.subcurve(bounds.first - NURBSCurve::FT(1.0),
                                   bounds.second + NURBSCurve::FT(1.0));
  BOOST_REQUIRE(clampedSub != nullptr);
  BOOST_CHECK(!clampedSub->isEmpty());
}

//-- Curve reversal tests

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

BOOST_AUTO_TEST_CASE(testReverseWeightedCurve)
{
  auto       controlPoints = createTestCurvePoints();
  auto       weights       = convertWeights({1.0, 3.0, 0.5, 2.0});
  NURBSCurve originalCurve(controlPoints, weights, 3);

  auto reversedCurve = originalCurve.reverse();

  BOOST_REQUIRE(reversedCurve != nullptr);
  BOOST_CHECK(reversedCurve->isRational());
  BOOST_CHECK_EQUAL(reversedCurve->numControlPoints(),
                    originalCurve.numControlPoints());

  // Weights should be reversed
  for (size_t controlIdx = 0; controlIdx < originalCurve.numControlPoints();
       ++controlIdx) {
    size_t reversedIdx = originalCurve.numControlPoints() - 1 - controlIdx;
    BOOST_CHECK_EQUAL(originalCurve.weight(controlIdx),
                      reversedCurve->weight(reversedIdx));
  }

  // Control points should be reversed
  for (size_t controlIdx = 0; controlIdx < originalCurve.numControlPoints();
       ++controlIdx) {
    size_t reversedIdx   = originalCurve.numControlPoints() - 1 - controlIdx;
    Point  originalPoint = originalCurve.controlPointN(controlIdx);
    Point  reversedPoint = reversedCurve->controlPointN(reversedIdx);
    BOOST_CHECK(isNearlyEqual(originalPoint, reversedPoint, 1e-10));
  }
}

BOOST_AUTO_TEST_CASE(testReverseDoubleReverse)
{
  auto       controlPoints = createTestCurvePoints();
  auto       weights       = convertWeights({1.0, 2.0, 1.5, 1.0});
  NURBSCurve originalCurve(controlPoints, weights, 3);

  auto reversed       = originalCurve.reverse();
  auto doubleReversed = reversed->reverse();

  BOOST_REQUIRE(doubleReversed != nullptr);

  // Double reversal should give back original curve (approximately)
  BOOST_CHECK_EQUAL(originalCurve.numControlPoints(),
                    doubleReversed->numControlPoints());
  BOOST_CHECK_EQUAL(originalCurve.degree(), doubleReversed->degree());

  // Check control points are back to original order
  for (size_t controlIdx = 0; controlIdx < originalCurve.numControlPoints();
       ++controlIdx) {
    Point originalPoint       = originalCurve.controlPointN(controlIdx);
    Point doubleReversedPoint = doubleReversed->controlPointN(controlIdx);
    BOOST_CHECK(isNearlyEqual(originalPoint, doubleReversedPoint, 1e-10));

    BOOST_CHECK_CLOSE(CGAL::to_double(originalCurve.weight(controlIdx)),
                      CGAL::to_double(doubleReversed->weight(controlIdx)),
                      1e-10);
  }
}

BOOST_AUTO_TEST_CASE(testReverseEmptyCurve)
{
  NURBSCurve emptyCurve;

  auto reversedEmpty = emptyCurve.reverse();

  BOOST_REQUIRE(reversedEmpty != nullptr);
  BOOST_CHECK(reversedEmpty->isEmpty());
}

//-- Curve joining tests

BOOST_AUTO_TEST_CASE(testJoinLinearCurves)
{
  // First curve: (0,0) to (2,2)
  std::vector<Point> firstPoints;
  firstPoints.emplace_back(0.0, 0.0);
  firstPoints.emplace_back(2.0, 2.0);
  NURBSCurve firstCurve(firstPoints, 1);

  // Second curve: (2,2) to (4,0) - starts where first ends
  std::vector<Point> secondPoints;
  secondPoints.emplace_back(2.0, 2.0);
  secondPoints.emplace_back(4.0, 0.0);
  NURBSCurve secondCurve(secondPoints, 1);

  auto joinedCurve = firstCurve.join(secondCurve);

  BOOST_REQUIRE(joinedCurve != nullptr);
  BOOST_CHECK(!joinedCurve->isEmpty());

  // Joined curve should have 3 control points (overlapping point removed)
  BOOST_CHECK_EQUAL(joinedCurve->numControlPoints(), 3U);

  // Check endpoints
  auto  joinedBounds = joinedCurve->parameterBounds();
  Point joinedStart  = joinedCurve->evaluate(joinedBounds.first);
  Point joinedEnd    = joinedCurve->evaluate(joinedBounds.second);

  BOOST_CHECK(isNearlyEqual(joinedStart, Point(0.0, 0.0), 1e-10));
  BOOST_CHECK(isNearlyEqual(joinedEnd, Point(4.0, 0.0), 1e-10));
}

BOOST_AUTO_TEST_CASE(testJoinNonAdjacentCurves)
{
  // First curve
  std::vector<Point> firstPoints;
  firstPoints.emplace_back(0.0, 0.0);
  firstPoints.emplace_back(1.0, 1.0);
  NURBSCurve firstCurve(firstPoints, 1);

  // Second curve - not adjacent (gap)
  std::vector<Point> secondPoints;
  secondPoints.emplace_back(2.0, 2.0); // Gap from (1,1) to (2,2)
  secondPoints.emplace_back(3.0, 1.0);
  NURBSCurve secondCurve(secondPoints, 1);

  // Join with default tolerance should fail
  BOOST_CHECK_THROW((void)firstCurve.join(secondCurve), Exception);

  // Join with large tolerance should succeed
  auto joinedCurve = firstCurve.join(secondCurve, NURBSCurve::Continuity::C0,
                                     NURBSCurve::FT(2.0)); // Large tolerance

  BOOST_REQUIRE(joinedCurve != nullptr);
  BOOST_CHECK(!joinedCurve->isEmpty());
}

BOOST_AUTO_TEST_CASE(testJoinDifferentDegrees)
{
  // Linear curve
  std::vector<Point> linearPoints;
  linearPoints.emplace_back(0.0, 0.0);
  linearPoints.emplace_back(1.0, 1.0);
  NURBSCurve linearCurve(linearPoints, 1);

  // Quadratic curve starting at (1,1)
  std::vector<Point> quadraticPoints;
  quadraticPoints.emplace_back(1.0, 1.0);
  quadraticPoints.emplace_back(2.0, 3.0);
  quadraticPoints.emplace_back(3.0, 1.0);
  NURBSCurve quadraticCurve(quadraticPoints, 2);

  auto joinedCurve = linearCurve.join(quadraticCurve);

  BOOST_REQUIRE(joinedCurve != nullptr);
  BOOST_CHECK(!joinedCurve->isEmpty());

  // Should use maximum degree
  BOOST_CHECK_EQUAL(joinedCurve->degree(), 2U);
}

BOOST_AUTO_TEST_CASE(testJoinEmptyCurves)
{
  NURBSCurve emptyCurve1, emptyCurve2;

  auto joinedEmpty = emptyCurve1.join(emptyCurve2);
  BOOST_REQUIRE(joinedEmpty != nullptr);
  BOOST_CHECK(joinedEmpty->isEmpty());

  // Empty with non-empty
  auto       controlPoints = createTestCurvePoints();
  NURBSCurve nonEmptyCurve(controlPoints, 2);

  auto joinedWithEmpty = emptyCurve1.join(nonEmptyCurve);
  BOOST_REQUIRE(joinedWithEmpty != nullptr);
  BOOST_CHECK_EQUAL(joinedWithEmpty->numControlPoints(),
                    nonEmptyCurve.numControlPoints());

  auto joinedEmptyWith = nonEmptyCurve.join(emptyCurve1);
  BOOST_REQUIRE(joinedEmptyWith != nullptr);
  BOOST_CHECK_EQUAL(joinedEmptyWith->numControlPoints(),
                    nonEmptyCurve.numControlPoints());
}

BOOST_AUTO_TEST_CASE(testJoinWithDifferentTypes)
{
  auto       controlPoints = createTestCurvePoints();
  NURBSCurve nurbsCurve(controlPoints, 2);

  std::vector<Point> incompatiblePoints;
  incompatiblePoints.emplace_back(100.0, 100.0);
  incompatiblePoints.emplace_back(102.0, 102.0);
  NURBSCurve incompatibleCurve(incompatiblePoints, 1);
  BOOST_CHECK_THROW(static_cast<void>(nurbsCurve.join(incompatibleCurve)),
                    Exception);
  // WHen LineString will inherited Curve
  // Create a LineString
  // std::vector<Point> linePoints;
  // linePoints.emplace_back(8.0, 0.0); // Start where NURBS ends
  // linePoints.emplace_back(10.0, 2.0);
  // LineString lineString;
  // for (const auto &point : linePoints) {
  //   lineString.addPoint(point);
  // }
  //
  // // Joining with different curve type should throw
  // BOOST_CHECK_THROW(nurbsCurve.join(lineString), Exception);
}

//-- Curve offset tests

BOOST_AUTO_TEST_CASE(testOffsetLinearCurve2D)
{
  std::vector<Point> linearPoints;
  linearPoints.emplace_back(0.0, 0.0);
  linearPoints.emplace_back(4.0, 0.0); // Horizontal line

  NURBSCurve horizontalCurve(linearPoints, 1);

  // Offset upward
  NURBSCurve::FT offsetDistance(1.0);
  auto           offsetCurve = horizontalCurve.offset(offsetDistance);

  BOOST_REQUIRE(offsetCurve != nullptr);
  BOOST_CHECK(!offsetCurve->isEmpty());

  // Check that offset curve is parallel and displaced
  auto  offsetBounds = offsetCurve->parameterBounds();
  Point offsetStart  = offsetCurve->evaluate(offsetBounds.first);
  Point offsetEnd    = offsetCurve->evaluate(offsetBounds.second);

  // Should be offset upward by 1 unit
  BOOST_CHECK_CLOSE(CGAL::to_double(offsetStart.x()), 0.0, 1e-6);
  BOOST_CHECK_CLOSE(CGAL::to_double(offsetStart.y()), 1.0, 1e-6);
  BOOST_CHECK_CLOSE(CGAL::to_double(offsetEnd.x()), 4.0, 1e-6);
  BOOST_CHECK_CLOSE(CGAL::to_double(offsetEnd.y()), 1.0, 1e-6);

  // Negative offset should go in opposite direction
  auto negativeOffset = horizontalCurve.offset(NURBSCurve::FT(-1.0));
  BOOST_REQUIRE(negativeOffset != nullptr);

  auto  negBounds = negativeOffset->parameterBounds();
  Point negStart  = negativeOffset->evaluate(negBounds.first);

  BOOST_CHECK_CLOSE(CGAL::to_double(negStart.y()), -1.0, 1e-6);
}

BOOST_AUTO_TEST_CASE(testOffsetCurvedCurve2D)
{
  // Create semi-circular curve
  std::vector<Point> circularPoints;
  circularPoints.emplace_back(0.0, 0.0);
  circularPoints.emplace_back(1.0, 1.0);
  circularPoints.emplace_back(2.0, 0.0);

  auto weights =
      convertWeights({1.0, std::sqrt(2.0) / 2.0, 1.0}); // Circular arc weights
  NURBSCurve circularCurve(circularPoints, weights, 2);

  NURBSCurve::FT offsetDistance(0.5);
  auto           offsetCurve = circularCurve.offset(offsetDistance);

  BOOST_REQUIRE(offsetCurve != nullptr);
  BOOST_CHECK(!offsetCurve->isEmpty());

  // Offset curve should maintain the same general shape but be displaced
  auto origBounds   = circularCurve.parameterBounds();
  auto offsetBounds = offsetCurve->parameterBounds();

  Point origMid = circularCurve.evaluate(
      (origBounds.first + origBounds.second) / NURBSCurve::FT(2));
  Point offsetMid = offsetCurve->evaluate(
      (offsetBounds.first + offsetBounds.second) / NURBSCurve::FT(2));

  // Offset midpoint should be higher than original (for positive offset)
  BOOST_CHECK(offsetMid.y() > origMid.y());
}

BOOST_AUTO_TEST_CASE(testOffsetZeroDistance)
{
  auto       controlPoints = createTestCurvePoints();
  NURBSCurve originalCurve(controlPoints, 2);

  // Zero offset should return identical curve
  auto zeroOffset = originalCurve.offset(NURBSCurve::FT(0.0));

  BOOST_REQUIRE(zeroOffset != nullptr);
  BOOST_CHECK_EQUAL(zeroOffset->numControlPoints(),
                    originalCurve.numControlPoints());
  BOOST_CHECK_EQUAL(zeroOffset->degree(), originalCurve.degree());

  // Should be essentially the same curve
  auto origBounds = originalCurve.parameterBounds();
  auto zeroBounds = zeroOffset->parameterBounds();

  Point origStart = originalCurve.evaluate(origBounds.first);
  Point zeroStart = zeroOffset->evaluate(zeroBounds.first);

  BOOST_CHECK(isNearlyEqual(origStart, zeroStart, 1e-10));
}

BOOST_AUTO_TEST_CASE(testOffset3DCurve)
{
  std::vector<Point> points3D;
  points3D.emplace_back(0.0, 0.0, 0.0);
  points3D.emplace_back(1.0, 1.0, 1.0);
  points3D.emplace_back(2.0, 0.0, 0.0);

  NURBSCurve curve3D(points3D, 2);

  // 3D offset should throw exception (not implemented for 3D)
  BOOST_CHECK_THROW((void)curve3D.offset(NURBSCurve::FT(1.0)), Exception);
}

BOOST_AUTO_TEST_CASE(testOffsetEmptyCurve)
{
  NURBSCurve emptyCurve;

  // Offset of empty curve should return empty curve
  auto offsetEmpty = emptyCurve.offset(NURBSCurve::FT(1.0));
  BOOST_REQUIRE(offsetEmpty != nullptr);
  BOOST_CHECK(offsetEmpty->isEmpty());
}

//-- Advanced NURBS operations tests

BOOST_AUTO_TEST_CASE(testKnotInsertionBasic)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 2.0);
  controlPoints.emplace_back(2.0, 0.0);

  NURBSCurve originalCurve(controlPoints, 2);

  auto             bounds = originalCurve.parameterBounds();
  NURBSCurve::Knot insertKnot =
      bounds.first + (bounds.second - bounds.first) * NURBSCurve::FT(0.5);

  auto insertedCurve = originalCurve.insertKnot(insertKnot);

  BOOST_REQUIRE(insertedCurve != nullptr);
  // Note: Current implementation is stub, so this just tests the interface
  BOOST_CHECK_EQUAL(insertedCurve->degree(), originalCurve.degree());
}

BOOST_AUTO_TEST_CASE(testMultipleKnotInsertion)
{
  auto       controlPoints = createTestCurvePoints();
  NURBSCurve originalCurve(controlPoints, 3);

  auto             bounds = originalCurve.parameterBounds();
  NURBSCurve::Knot insertKnot =
      bounds.first + (bounds.second - bounds.first) * NURBSCurve::FT(0.4);

  // Insert knot multiple times
  auto onceInserted  = originalCurve.insertKnot(insertKnot, 1);
  auto twiceInserted = originalCurve.insertKnot(insertKnot, 2);

  BOOST_REQUIRE(onceInserted != nullptr);
  BOOST_REQUIRE(twiceInserted != nullptr);

  // Both should maintain curve degree (stub implementation)
  BOOST_CHECK_EQUAL(onceInserted->degree(), originalCurve.degree());
  BOOST_CHECK_EQUAL(twiceInserted->degree(), originalCurve.degree());
}

BOOST_AUTO_TEST_CASE(testKnotVectorRefinement)
{
  auto       controlPoints = createTestCurvePoints();
  NURBSCurve originalCurve(controlPoints, 2);

  auto                          bounds = originalCurve.parameterBounds();
  std::vector<NURBSCurve::Knot> newKnots;
  newKnots.push_back(bounds.first +
                     (bounds.second - bounds.first) * NURBSCurve::FT(0.3));
  newKnots.push_back(bounds.first +
                     (bounds.second - bounds.first) * NURBSCurve::FT(0.7));

  auto refinedCurve = originalCurve.refineKnotVector(newKnots);

  BOOST_REQUIRE(refinedCurve != nullptr);
  // Stub implementation returns copy
  BOOST_CHECK_EQUAL(refinedCurve->numControlPoints(),
                    originalCurve.numControlPoints());
}

BOOST_AUTO_TEST_CASE(testDegreeElevation)
{
  std::vector<Point> controlPoints;
  controlPoints.emplace_back(0.0, 0.0);
  controlPoints.emplace_back(1.0, 1.0);
  controlPoints.emplace_back(2.0, 0.0);

  NURBSCurve originalCurve(controlPoints, 1); // Degree 1 (linear)

  auto elevatedCurve = originalCurve.elevateDegree(2); // Elevate by 2

  BOOST_REQUIRE(elevatedCurve != nullptr);
  // Stub implementation returns copy
  BOOST_CHECK_EQUAL(elevatedCurve->degree(), originalCurve.degree());
}

BOOST_AUTO_TEST_CASE(testDegreeReduction)
{
  auto       controlPoints = createTestCurvePoints();
  NURBSCurve originalCurve(controlPoints, 3); // Degree 3

  auto reducedCurve = originalCurve.reduceDegree(NURBSCurve::FT(1e-6));

  BOOST_REQUIRE(reducedCurve != nullptr);
  // Stub implementation returns copy
  BOOST_CHECK_EQUAL(reducedCurve->degree(), originalCurve.degree());
}

BOOST_AUTO_TEST_CASE(testKnotRemoval)
{
  auto       controlPoints = createTestCurvePoints();
  NURBSCurve originalCurve(controlPoints, 3);

  auto simplifiedCurve = originalCurve.removeKnots(NURBSCurve::FT(1e-6));

  BOOST_REQUIRE(simplifiedCurve != nullptr);
  // Stub implementation returns copy
  BOOST_CHECK_EQUAL(simplifiedCurve->numControlPoints(),
                    originalCurve.numControlPoints());
}

//-- Error handling for manipulation operations

BOOST_AUTO_TEST_CASE(testManipulationErrorHandling)
{
  NURBSCurve emptyCurve;

  // Operations on empty curves should handle gracefully
  BOOST_CHECK_NO_THROW(auto result =
                           emptyCurve.split(NURBSCurve::Parameter(0.5)));
  BOOST_CHECK_NO_THROW(auto result = emptyCurve.reverse());
  BOOST_CHECK_NO_THROW(auto result =
                           emptyCurve.subcurve(NURBSCurve::Parameter(0.0),
                                               NURBSCurve::Parameter(1.0)));

  auto emptyOffset = emptyCurve.offset(NURBSCurve::FT(1.0));
  BOOST_CHECK(emptyOffset->isEmpty());

  // Advanced operations on empty curves
  BOOST_CHECK_NO_THROW(auto result =
                           emptyCurve.insertKnot(NURBSCurve::Parameter(0.5)));
  BOOST_CHECK_NO_THROW(auto result = emptyCurve.elevateDegree(1));
}

BOOST_AUTO_TEST_CASE(testInvalidParametersHandling)
{
  auto       controlPoints = createTestCurvePoints();
  NURBSCurve validCurve(controlPoints, 2);

  auto bounds = validCurve.parameterBounds();

  // Invalid split parameters should be handled gracefully (clamping)
  auto [beforeStart1, afterStart1] =
      validCurve.split(bounds.first - NURBSCurve::FT(1.0));
  auto [beforeEnd1, afterEnd1] =
      validCurve.split(bounds.second + NURBSCurve::FT(1.0));

  BOOST_CHECK(beforeStart1->isEmpty());
  BOOST_CHECK(!afterStart1->isEmpty());
  BOOST_CHECK(!beforeEnd1->isEmpty());
  BOOST_CHECK(afterEnd1->isEmpty());

  // Invalid subcurve ranges
  auto invalidSub =
      validCurve.subcurve(bounds.second, bounds.first); // Reversed
  BOOST_CHECK(invalidSub->isEmpty());
}

//-- Performance tests for manipulation operations

BOOST_AUTO_TEST_CASE(testManipulationPerformance)
{
  // Create complex curve with many control points
  auto       circularPoints = createCircularPoints(2.0, 50);
  auto       weights        = convertWeights(std::vector<double>(50, 1.0));
  NURBSCurve complexCurve(circularPoints, weights, 5);

  // All operations should complete without timing out
  BOOST_CHECK_NO_THROW(auto result =
                           complexCurve.split(NURBSCurve::Parameter(0.5)));
  BOOST_CHECK_NO_THROW(auto result = complexCurve.reverse());
  BOOST_CHECK_NO_THROW(auto result =
                           complexCurve.subcurve(NURBSCurve::Parameter(0.2),
                                                 NURBSCurve::Parameter(0.8)));
  BOOST_CHECK_NO_THROW(auto result = complexCurve.offset(NURBSCurve::FT(0.5)));

  // Advanced operations
  BOOST_CHECK_NO_THROW(auto result =
                           complexCurve.insertKnot(NURBSCurve::Parameter(0.6)));
  BOOST_CHECK_NO_THROW(auto result = complexCurve.elevateDegree(1));
}

BOOST_AUTO_TEST_SUITE_END()
