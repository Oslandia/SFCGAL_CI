// Copyright (c) 2025-2025, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later
#include <SFCGAL/Exception.h>
#include <SFCGAL/Segment.h>

#include <boost/test/unit_test.hpp>

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_SegmentTest)

// Test constructors
BOOST_AUTO_TEST_CASE(testDefaultConstructor)
{
  Segment segment;
  BOOST_CHECK(segment.isEmpty());
}

BOOST_AUTO_TEST_CASE(testConstructorWithPoints)
{
  // 2D points
  {
    Point   p1(0.0, 0.0);
    Point   p2(3.0, 4.0);
    Segment segment(p1, p2);

    BOOST_CHECK(!segment.isEmpty());
    BOOST_CHECK(!segment.is3D());
    BOOST_CHECK(!segment.isMeasured());
    BOOST_CHECK_EQUAL(segment.source().x(), 0.0);
    BOOST_CHECK_EQUAL(segment.source().y(), 0.0);
    BOOST_CHECK_EQUAL(segment.target().x(), 3.0);
    BOOST_CHECK_EQUAL(segment.target().y(), 4.0);
    BOOST_CHECK_EQUAL(CGAL::to_double(segment.length()), 5.0);
  }

  // 3D points
  {
    Point   p1(0.0, 0.0, 0.0);
    Point   p2(3.0, 0.0, 4.0);
    Segment segment(p1, p2);

    BOOST_CHECK(!segment.isEmpty());
    BOOST_CHECK(segment.is3D());
    BOOST_CHECK(!segment.isMeasured());
    BOOST_CHECK_EQUAL(segment.source().z(), 0.0);
    BOOST_CHECK_EQUAL(segment.target().z(), 4.0);
    BOOST_CHECK_EQUAL(CGAL::to_double(segment.length()), 5.0);
  }

  // Measured points (2D with M)
  {
    Point p1(1.0, 2.0);
    p1.setM(10.0);
    Point p2(5.0, 7.0);
    p2.setM(20.0);

    Segment segment(p1, p2);
    BOOST_CHECK(!segment.isEmpty());
    BOOST_CHECK(!segment.is3D());
    BOOST_CHECK(segment.isMeasured());
    BOOST_CHECK_EQUAL(segment.source().m(), 10.0);
    BOOST_CHECK_EQUAL(segment.target().m(), 20.0);
  }

  // Measured points (3D with M)
  {
    Point p1(1.0, 2.0, 3.0);
    p1.setM(10.0);
    Point p2(5.0, 7.0, 9.0);
    p2.setM(20.0);

    Segment segment(p1, p2);
    BOOST_CHECK(!segment.isEmpty());
    BOOST_CHECK(segment.is3D());
    BOOST_CHECK(segment.isMeasured());
    BOOST_CHECK_EQUAL(segment.source().m(), 10.0);
    BOOST_CHECK_EQUAL(segment.target().m(), 20.0);
  }
}

BOOST_AUTO_TEST_CASE(testConstructorWithCGALPoints)
{
  // 2D CGAL points
  {
    Kernel::Point_2 p1(0.0, 0.0);
    Kernel::Point_2 p2(3.0, 4.0);
    Segment         segment(p1, p2);

    BOOST_CHECK(!segment.isEmpty());
    BOOST_CHECK(!segment.is3D());
    BOOST_CHECK(!segment.isMeasured());
    BOOST_CHECK_EQUAL(segment.source().x(), 0.0);
    BOOST_CHECK_EQUAL(segment.source().y(), 0.0);
    BOOST_CHECK_EQUAL(segment.target().x(), 3.0);
    BOOST_CHECK_EQUAL(segment.target().y(), 4.0);
    BOOST_CHECK_EQUAL(CGAL::to_double(segment.length()), 5.0);
  }

  // 3D CGAL points
  {
    Kernel::Point_3 p1(0.0, 0.0, 0.0);
    Kernel::Point_3 p2(3.0, 0.0, 4.0);
    Segment         segment(p1, p2);

    BOOST_CHECK(!segment.isEmpty());
    BOOST_CHECK(segment.is3D());
    BOOST_CHECK(!segment.isMeasured());
    BOOST_CHECK_EQUAL(segment.source().z(), 0.0);
    BOOST_CHECK_EQUAL(segment.target().z(), 4.0);
    BOOST_CHECK_EQUAL(CGAL::to_double(segment.length()), 5.0);
  }
}

BOOST_AUTO_TEST_CASE(testConstructorWithCGALSegments)
{
  // 2D CGAL segment
  {
    Kernel::Segment_2 cgalSegment(Kernel::Point_2(0.0, 0.0),
                                  Kernel::Point_2(3.0, 4.0));
    Segment           segment(cgalSegment);

    BOOST_CHECK(!segment.isEmpty());
    BOOST_CHECK(!segment.is3D());
    BOOST_CHECK(!segment.isMeasured());
    BOOST_CHECK_EQUAL(segment.source().x(), 0.0);
    BOOST_CHECK_EQUAL(segment.source().y(), 0.0);
    BOOST_CHECK_EQUAL(segment.target().x(), 3.0);
    BOOST_CHECK_EQUAL(segment.target().y(), 4.0);
    BOOST_CHECK_EQUAL(CGAL::to_double(segment.length()), 5.0);
  }

  // 3D CGAL segment
  {
    Kernel::Segment_3 cgalSegment(Kernel::Point_3(0.0, 0.0, 0.0),
                                  Kernel::Point_3(3.0, 0.0, 4.0));
    Segment           segment(cgalSegment);

    BOOST_CHECK(!segment.isEmpty());
    BOOST_CHECK(segment.is3D());
    BOOST_CHECK(!segment.isMeasured());
    BOOST_CHECK_EQUAL(segment.source().z(), 0.0);
    BOOST_CHECK_EQUAL(segment.target().z(), 4.0);
    BOOST_CHECK_EQUAL(CGAL::to_double(segment.length()), 5.0);
  }
}

BOOST_AUTO_TEST_CASE(testConstructorWithInvalidPoints)
{
  // Mixed 2D/3D points
  {
    Point p1(0.0, 0.0);
    Point p2(3.0, 4.0, 5.0);
    BOOST_CHECK_THROW(Segment segment(p1, p2), Exception);
  }

  // Mixed measured/unmeasured points
  {
    Point p1(0.0, 0.0);
    Point p2(3.0, 4.0);
    p2.setM(10.0);
    BOOST_CHECK_THROW(Segment segment(p1, p2), Exception);
  }

  // One empty point
  {
    Point   p1;
    Point   p2(3.0, 4.0);
    Segment segment(p1, p2);
    BOOST_CHECK(segment.isEmpty());
  }

  // Two empty points
  {
    Point   p1;
    Point   p2;
    Segment segment(p1, p2);
    BOOST_CHECK(segment.isEmpty());
  }
}

// Test emptiness and dimension checking
BOOST_AUTO_TEST_CASE(testEmptiness)
{
  Segment empty;
  BOOST_CHECK(empty.isEmpty());

  Point   p1(0.0, 0.0);
  Point   p2(3.0, 4.0);
  Segment nonEmpty(p1, p2);
  BOOST_CHECK(!nonEmpty.isEmpty());

  // Make segment empty
  nonEmpty.clear();
  BOOST_CHECK(nonEmpty.isEmpty());
}

BOOST_AUTO_TEST_CASE(testDimensions)
{
  Point   p1(0.0, 0.0);
  Point   p2(3.0, 4.0);
  Segment segment2D(p1, p2);
  BOOST_CHECK(!segment2D.is3D());
  BOOST_CHECK(!segment2D.isMeasured());

  Point   p3(0.0, 0.0, 0.0);
  Point   p4(3.0, 4.0, 5.0);
  Segment segment3D(p3, p4);
  BOOST_CHECK(segment3D.is3D());
  BOOST_CHECK(!segment3D.isMeasured());

  p1.setM(10.0);
  p2.setM(20.0);
  Segment segment2DM(p1, p2);
  BOOST_CHECK(!segment2DM.is3D());
  BOOST_CHECK(segment2DM.isMeasured());

  p3.setM(10.0);
  p4.setM(20.0);
  Segment segment3DM(p3, p4);
  BOOST_CHECK(segment3DM.is3D());
  BOOST_CHECK(segment3DM.isMeasured());
}

BOOST_AUTO_TEST_CASE(testHasSameDimension)
{
  // Normally, always true!
  // Empty segment
  {
    Segment empty;
    BOOST_CHECK(empty.hasSameDimension());
  }

  // Valid segments
  {
    Point   p1(0.0, 0.0);
    Point   p2(3.0, 4.0);
    Segment segment(p1, p2);
    BOOST_CHECK(segment.hasSameDimension());
  }
}

// Test geometric properties
BOOST_AUTO_TEST_CASE(testLength)
{
  // Empty segment
  {
    Segment empty;
    BOOST_CHECK_EQUAL(empty.length(), 0.0);
  }

  // Horizontal segment
  {
    Segment horizontal(Point(0.0, 0.0), Point(5.0, 0.0));
    BOOST_CHECK_EQUAL(CGAL::to_double(horizontal.length()), 5.0);
  }

  // Vertical segment
  {
    Segment vertical(Point(0.0, 0.0), Point(0.0, 7.0));
    BOOST_CHECK_EQUAL(CGAL::to_double(vertical.length()), 7.0);
  }

  // 45-degree segment
  {
    Segment diagonal(Point(0.0, 0.0), Point(1.0, 1.0));
    BOOST_CHECK_CLOSE(CGAL::to_double(diagonal.length()), std::sqrt(2.0),
                      0.0001);
  }

  // 3D segment
  {
    Segment seg3D(Point(0.0, 0.0, 0.0), Point(1.0, 1.0, 1.0));
    BOOST_CHECK_CLOSE(CGAL::to_double(seg3D.length()), std::sqrt(3.0), 0.0001);
  }
}

BOOST_AUTO_TEST_CASE(testIsDegenerate)
{
  // Degenerate segment (zero length)
  {
    Segment degenerate(Point(1.0, 1.0), Point(1.0, 1.0));
    BOOST_CHECK(degenerate.isDegenerate());
  }

  // Non-degenerate segment
  {
    Segment nonDegenerate(Point(1.0, 1.0), Point(2.0, 2.0));
    BOOST_CHECK(!nonDegenerate.isDegenerate());
  }

  // Empty segment
  {
    Segment empty;
    BOOST_CHECK(empty.isDegenerate());
  }
}

// Test distance calculation
BOOST_AUTO_TEST_CASE(testDistanceToPoint)
{
  Segment segment(Point(0.0, 0.0), Point(3.0, 0.0));

  // Point on the segment
  BOOST_CHECK_SMALL(segment.distanceToPoint(Point(1.0, 0.0)), EPSILON);

  // Point perpendicular to the segment
  BOOST_CHECK_CLOSE(segment.distanceToPoint(Point(1.0, 1.0)), 1.0, EPSILON);

  // Point before the segment
  BOOST_CHECK_CLOSE(segment.distanceToPoint(Point(-1.0, 0.0)), 1.0, EPSILON);

  // Point after the segment
  BOOST_CHECK_CLOSE(segment.distanceToPoint(Point(4.0, 0.0)), 1.0, EPSILON);

  // Distance with double coordinates
  BOOST_CHECK_CLOSE(segment.distanceToPoint(1.5, 2.0), 2.0, EPSILON);

  // Distance with CGAL points
  BOOST_CHECK_CLOSE(segment.distanceToPoint(Kernel::Point_2(1.5, 2.0)), 2.0,
                    EPSILON);

  // 3D segment and point
  Segment segment3D(Point(0.0, 0.0, 0.0), Point(3.0, 0.0, 0.0));
  BOOST_CHECK_CLOSE(segment3D.distanceToPoint(Point(1.0, 1.0, 1.0)),
                    std::sqrt(2.0), EPSILON);

  // Empty segment
  Segment empty;
  BOOST_CHECK_EQUAL(empty.distanceToPoint(Point(0.0, 0.0)), 0.0);
}

// Test interpolation
BOOST_AUTO_TEST_CASE(testInterpolationParameter)
{
  Segment segment(Point(0.0, 0.0), Point(10.0, 0.0));

  // Parameter at source
  BOOST_CHECK_SMALL(
      CGAL::to_double(segment.exactInterpolationParameter(Point(0.0, 0.0))),
      EPSILON);

  // Parameter at target
  BOOST_CHECK_CLOSE(
      CGAL::to_double(segment.exactInterpolationParameter(Point(10.0, 0.0))),
      1.0, EPSILON);

  // Parameter in the middle
  BOOST_CHECK_CLOSE(
      CGAL::to_double(segment.exactInterpolationParameter(Point(5.0, 0.0))),
      0.5, EPSILON);

  // Parameter for point not on the segment (should project)
  BOOST_CHECK_CLOSE(
      CGAL::to_double(segment.exactInterpolationParameter(Point(5.0, 2.0))),
      0.5, EPSILON);

  // Parameter for point before the segment (should clamp to 0)
  BOOST_CHECK_CLOSE(
      CGAL::to_double(segment.exactInterpolationParameter(Point(-5.0, 0.0))),
      0.0, EPSILON);

  // Parameter for point after the segment (should clamp to 1)
  BOOST_CHECK_CLOSE(
      CGAL::to_double(segment.exactInterpolationParameter(Point(15.0, 0.0))),
      1.0, EPSILON);

  // Double version
  BOOST_CHECK_CLOSE(segment.interpolationParameter(5.0, 2.0), 0.5, EPSILON);

  // 3D version
  Segment segment3D(Point(0.0, 0.0, 0.0), Point(10.0, 0.0, 0.0));
  BOOST_CHECK_CLOSE(CGAL::to_double(segment3D.exactInterpolationParameter(
                        Point(7.0, 2.0, 3.0))),
                    0.7, EPSILON);

  // Empty segment
  Segment empty;
  BOOST_CHECK_EQUAL(empty.exactInterpolationParameter(Point(0.0, 0.0)), 0);
}

BOOST_AUTO_TEST_CASE(testInterpolate)
{
  // 2D segment
  {
    Segment segment(Point(0.0, 0.0), Point(10.0, 0.0));

    // Interpolate at parameter 0
    Point p0 = segment.interpolate(0.0);
    BOOST_CHECK_EQUAL(p0.x(), 0.0);
    BOOST_CHECK_EQUAL(p0.y(), 0.0);

    // Interpolate at parameter 1
    Point p1 = segment.interpolate(1.0);
    BOOST_CHECK_EQUAL(p1.x(), 10.0);
    BOOST_CHECK_EQUAL(p1.y(), 0.0);

    // Interpolate at parameter 0.25
    Point p25 = segment.interpolate(0.25);
    BOOST_CHECK_EQUAL(p25.x(), 2.5);
    BOOST_CHECK_EQUAL(p25.y(), 0.0);

    // Test clamping - parameter below 0
    Point pNeg = segment.interpolate(-0.5);
    BOOST_CHECK_EQUAL(pNeg.x(), 0.0);
    BOOST_CHECK_EQUAL(pNeg.y(), 0.0);

    // Test clamping - parameter above 1
    Point pBig = segment.interpolate(1.5);
    BOOST_CHECK_EQUAL(pBig.x(), 10.0);
    BOOST_CHECK_EQUAL(pBig.y(), 0.0);
  }

  // 3D segment
  {
    Segment segment(Point(0.0, 0.0, 0.0), Point(10.0, 0.0, 10.0));

    Point p05 = segment.interpolate(0.5);
    BOOST_CHECK_EQUAL(p05.x(), 5.0);
    BOOST_CHECK_EQUAL(p05.y(), 0.0);
    BOOST_CHECK_EQUAL(p05.z(), 5.0);
  }

  // Segment with M values
  {
    Point p1(0.0, 0.0);
    p1.setM(10.0);
    Point p2(10.0, 0.0);
    p2.setM(20.0);
    Segment segment(p1, p2);

    Point p05 = segment.interpolate(0.5);
    BOOST_CHECK_EQUAL(p05.x(), 5.0);
    BOOST_CHECK_EQUAL(p05.y(), 0.0);
    BOOST_CHECK_EQUAL(p05.m(), 15.0);
  }

  // 3D segment with M values
  {
    Point p1(0.0, 0.0, 0.0);
    p1.setM(10.0);
    Point p2(10.0, 0.0, 10.0);
    p2.setM(20.0);
    Segment segment(p1, p2);

    Point p05 = segment.interpolate(0.5);
    BOOST_CHECK_EQUAL(p05.x(), 5.0);
    BOOST_CHECK_EQUAL(p05.y(), 0.0);
    BOOST_CHECK_EQUAL(p05.z(), 5.0);
    BOOST_CHECK_EQUAL(p05.m(), 15.0);
  }

  // Empty segment
  {
    Segment empty;
    BOOST_CHECK(empty.interpolate(0.5).isEmpty());
  }
}

BOOST_AUTO_TEST_CASE(testMidpoint)
{
  // 2D segment
  {
    Segment segment(Point(0.0, 0.0), Point(10.0, 0.0));
    Point   midpoint = segment.midpoint();

    BOOST_CHECK_EQUAL(midpoint.x(), 5.0);
    BOOST_CHECK_EQUAL(midpoint.y(), 0.0);
  }

  // 3D segment
  {
    Segment segment(Point(0.0, 0.0, 0.0), Point(10.0, 0.0, 10.0));
    Point   midpoint = segment.midpoint();

    BOOST_CHECK_EQUAL(midpoint.x(), 5.0);
    BOOST_CHECK_EQUAL(midpoint.y(), 0.0);
    BOOST_CHECK_EQUAL(midpoint.z(), 5.0);
  }

  // Segment with M values
  {
    Point p1(0.0, 0.0);
    p1.setM(10.0);
    Point p2(10.0, 0.0);
    p2.setM(20.0);
    Segment segment(p1, p2);

    Point midpoint = segment.midpoint();
    BOOST_CHECK_EQUAL(midpoint.x(), 5.0);
    BOOST_CHECK_EQUAL(midpoint.y(), 0.0);
    BOOST_CHECK_EQUAL(midpoint.m(), 15.0);
  }

  // Empty segment
  {
    Segment empty;
    BOOST_CHECK(empty.midpoint().isEmpty());
  }
}

// Test point on segment
BOOST_AUTO_TEST_CASE(testHasOn)
{
  Segment segment(Point(0.0, 0.0), Point(10.0, 0.0));

  // Test points on the segment
  BOOST_CHECK(segment.hasOn(Point(0.0, 0.0)));
  BOOST_CHECK(segment.hasOn(Point(5.0, 0.0)));
  BOOST_CHECK(segment.hasOn(Point(10.0, 0.0)));

  // Test points not on the segment
  BOOST_CHECK(!segment.hasOn(Point(5.0, 1.0)));
  BOOST_CHECK(!segment.hasOn(Point(-1.0, 0.0)));
  BOOST_CHECK(!segment.hasOn(Point(11.0, 0.0)));

  // Test with tolerance
  BOOST_CHECK(segment.hasOn(Point(5.0, 0.001), 0.01));
  BOOST_CHECK(!segment.hasOn(Point(5.0, 0.1), 0.01));

  // 3D segment
  Segment segment3D(Point(0.0, 0.0, 0.0), Point(10.0, 0.0, 0.0));

  // Test points on the segment
  BOOST_CHECK(segment3D.hasOn(Point(0.0, 0.0, 0.0)));
  BOOST_CHECK(segment3D.hasOn(Point(5.0, 0.0, 0.0)));

  // Test points not on the segment
  BOOST_CHECK(!segment3D.hasOn(Point(5.0, 0.0, 1.0)));

  // Test with tolerance
  BOOST_CHECK(segment3D.hasOn(Point(5.0, 0.0, 0.001), 0.01));

  // Empty segment
  Segment empty;
  BOOST_CHECK(!empty.hasOn(Point(0.0, 0.0)));
}

// Test setter methods
BOOST_AUTO_TEST_CASE(testSetters)
{
  // setSource/setTarget on non-empty segment
  {
    Point   p1(0.0, 0.0);
    Point   p2(10.0, 0.0);
    Segment segment(p1, p2);

    Point newSource(1.0, 1.0);
    segment.setSource(newSource);

    BOOST_CHECK_EQUAL(segment.source().x(), 1.0);
    BOOST_CHECK_EQUAL(segment.source().y(), 1.0);

    Point newTarget(9.0, 1.0);
    segment.setTarget(newTarget);

    BOOST_CHECK_EQUAL(segment.target().x(), 9.0);
    BOOST_CHECK_EQUAL(segment.target().y(), 1.0);
  }

  // setSource/setTarget with dimension mismatch
  {
    Point   p1(0.0, 0.0);
    Point   p2(10.0, 0.0);
    Segment segment(p1, p2);

    Point source3D(1.0, 1.0, 1.0);
    BOOST_CHECK_THROW(segment.setSource(source3D), Exception);

    Point p1_m(0.0, 0.0);
    p1_m.setM(1.0);
    BOOST_CHECK_THROW(segment.setSource(p1_m), Exception);
  }

  // setSource/setTarget on empty segment
  {
    Segment empty;
    Point   p1(0.0, 0.0);

    BOOST_CHECK_THROW(empty.setSource(p1), Exception);
    BOOST_CHECK_THROW(empty.setTarget(p1), Exception);
  }

  // setPoints on empty segment
  {
    Segment empty;
    Point   p1(0.0, 0.0);
    Point   p2(10.0, 0.0);

    empty.setPoints(p1, p2);
    BOOST_CHECK(!empty.isEmpty());
    BOOST_CHECK_EQUAL(empty.source().x(), 0.0);
    BOOST_CHECK_EQUAL(empty.target().x(), 10.0);
  }

  // setPoints with dimension mismatch
  {
    Segment empty;
    Point   p1(0.0, 0.0);
    Point   p2(10.0, 0.0, 1.0);

    BOOST_CHECK_THROW(empty.setPoints(p1, p2), Exception);
  }

  // setPoints with CGAL points
  {
    Segment         segment;
    Kernel::Point_2 p1(0.0, 0.0);
    Kernel::Point_2 p2(10.0, 0.0);

    segment.setPoints(p1, p2);
    BOOST_CHECK(!segment.isEmpty());
    BOOST_CHECK_EQUAL(segment.source().x(), 0.0);
    BOOST_CHECK_EQUAL(segment.target().x(), 10.0);
  }

  // clear
  {
    Point   p1(0.0, 0.0);
    Point   p2(10.0, 0.0);
    Segment segment(p1, p2);

    BOOST_CHECK(!segment.isEmpty());
    segment.clear();
    BOOST_CHECK(segment.isEmpty());
  }
}

// Test reverse
BOOST_AUTO_TEST_CASE(testReverse)
{
  // Normal segment
  {
    Segment segment(Point(1.0, 2.0), Point(3.0, 4.0));
    segment.reverse();

    BOOST_CHECK_EQUAL(segment.source().x(), 3.0);
    BOOST_CHECK_EQUAL(segment.source().y(), 4.0);
    BOOST_CHECK_EQUAL(segment.target().x(), 1.0);
    BOOST_CHECK_EQUAL(segment.target().y(), 2.0);
  }
}

// Test conversion methods
BOOST_AUTO_TEST_CASE(testConversion)
{
  // 2D segment to CGAL::Segment_2
  {
    Segment           segment(Point(1.0, 2.0), Point(3.0, 4.0));
    Kernel::Segment_2 seg2 = segment.toSegment_2();

    BOOST_CHECK_EQUAL(seg2.source().x(), 1.0);
    BOOST_CHECK_EQUAL(seg2.source().y(), 2.0);
    BOOST_CHECK_EQUAL(seg2.target().x(), 3.0);
    BOOST_CHECK_EQUAL(seg2.target().y(), 4.0);
  }

  // 3D segment to CGAL::Segment_3
  {
    Segment           segment3D(Point(1.0, 2.0, 3.0), Point(4.0, 5.0, 6.0));
    Kernel::Segment_3 seg3 = segment3D.toSegment_3();

    BOOST_CHECK_EQUAL(seg3.source().x(), 1.0);
    BOOST_CHECK_EQUAL(seg3.source().y(), 2.0);
    BOOST_CHECK_EQUAL(seg3.source().z(), 3.0);
    BOOST_CHECK_EQUAL(seg3.target().x(), 4.0);
    BOOST_CHECK_EQUAL(seg3.target().y(), 5.0);
    BOOST_CHECK_EQUAL(seg3.target().z(), 6.0);
  }
}

BOOST_AUTO_TEST_SUITE_END()
