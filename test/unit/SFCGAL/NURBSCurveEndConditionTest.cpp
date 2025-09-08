// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "SFCGAL/Exception.h"
#include "SFCGAL/NURBSCurve.h"
#include "SFCGAL/Point.h"

using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_NURBSCurveEndConditionTest)

BOOST_AUTO_TEST_CASE(testClampedEndCondition)
{
  std::vector<Point> points = {Point(0.0, 0.0), Point(1.0, 2.0),
                               Point(2.0, 1.0), Point(3.0, 0.0)};

  auto curve = NURBSCurve::interpolateCurve(points, 3,
                                            NURBSCurve::KnotMethod::CENTRIPETAL,
                                            NURBSCurve::EndCondition::CLAMPED);

  BOOST_CHECK(!curve->isEmpty());
  BOOST_CHECK_EQUAL(curve->numPoints(), points.size());
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
  BOOST_CHECK_EQUAL(curve->numPoints(), points.size());
  BOOST_CHECK(curve->isValid());
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
  BOOST_CHECK(curve->isValid());
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

BOOST_AUTO_TEST_CASE(testTangentEndConditionFallback)
{
  std::vector<Point> points = {Point(0.0, 0.0), Point(1.0, 2.0),
                               Point(2.0, 1.0), Point(3.0, 0.0)};

  // TANGENT should fall back to CLAMPED for now
  auto curve = NURBSCurve::interpolateCurve(points, 3,
                                            NURBSCurve::KnotMethod::CENTRIPETAL,
                                            NURBSCurve::EndCondition::TANGENT);

  BOOST_CHECK(!curve->isEmpty());
  BOOST_CHECK_EQUAL(curve->numPoints(), points.size());
}

BOOST_AUTO_TEST_CASE(testFitMethodInterface)
{
  std::vector<Point> points = {Point(0.0, 0.0), Point(1.0, 1.0),
                               Point(2.0, 0.5), Point(3.0, 0.0)};

  // Test INTERPOLATE mode
  auto interpolatedCurve = NURBSCurve::fitCurve(
      points, 3, NURBSCurve::FitMethod::INTERPOLATE,
      NURBSCurve::KnotMethod::CENTRIPETAL, NURBSCurve::EndCondition::NATURAL);

  BOOST_CHECK(!interpolatedCurve->isEmpty());

  // Test APPROXIMATE mode
  auto approximatedCurve = NURBSCurve::fitCurve(
      points, 2, NURBSCurve::FitMethod::APPROXIMATE,
      NURBSCurve::KnotMethod::CHORD_LENGTH, NURBSCurve::EndCondition::CLAMPED,
      NURBSCurve::FT(0.1));

  BOOST_CHECK(!approximatedCurve->isEmpty());
}

BOOST_AUTO_TEST_CASE(testKnotMethodConsistency)
{
  std::vector<Point> points = {Point(0.0, 0.0), Point(1.0, 1.0),
                               Point(2.0, 0.0), Point(3.0, 1.0),
                               Point(4.0, 0.0)};

  // All knot methods should produce valid curves
  auto uniformCurve =
      NURBSCurve::interpolateCurve(points, 3, NURBSCurve::KnotMethod::UNIFORM);
  BOOST_CHECK(uniformCurve->isValid());

  auto chordCurve = NURBSCurve::interpolateCurve(
      points, 3, NURBSCurve::KnotMethod::CHORD_LENGTH);
  BOOST_CHECK(chordCurve->isValid());

  auto centripetalCurve = NURBSCurve::interpolateCurve(
      points, 3, NURBSCurve::KnotMethod::CENTRIPETAL);
  BOOST_CHECK(centripetalCurve->isValid());
}

BOOST_AUTO_TEST_CASE(testNaturalVsClampedDifference)
{
  std::vector<Point> points = {Point(0.0, 0.0), Point(1.0, 2.0),
                               Point(2.0, -1.0), Point(3.0, 0.0),
                               Point(4.0, 1.0)};

  auto clampedCurve = NURBSCurve::interpolateCurve(
      points, 3, NURBSCurve::KnotMethod::CENTRIPETAL,
      NURBSCurve::EndCondition::CLAMPED);

  auto naturalCurve = NURBSCurve::interpolateCurve(
      points, 3, NURBSCurve::KnotMethod::CENTRIPETAL,
      NURBSCurve::EndCondition::NATURAL);

  // Both should be valid but potentially different
  BOOST_CHECK(clampedCurve->isValid());
  BOOST_CHECK(naturalCurve->isValid());
  BOOST_CHECK_EQUAL(clampedCurve->numPoints(), naturalCurve->numPoints());
}

BOOST_AUTO_TEST_CASE(testBackwardCompatibility)
{
  std::vector<Point> points = {Point(0.0, 0.0), Point(1.0, 1.0),
                               Point(2.0, 0.0)};

  // Old-style calls should still work
  auto oldStyleCurve = NURBSCurve::interpolateCurve(points);
  BOOST_CHECK(!oldStyleCurve->isEmpty());

  auto oldStyleApprox =
      NURBSCurve::approximateCurve(points, 2, NURBSCurve::FT(0.1));
  BOOST_CHECK(!oldStyleApprox->isEmpty());
}

BOOST_AUTO_TEST_SUITE_END()