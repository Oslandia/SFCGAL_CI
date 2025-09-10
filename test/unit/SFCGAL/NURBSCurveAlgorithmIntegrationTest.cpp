#include <boost/test/unit_test.hpp>

#include "SFCGAL/NURBSCurve.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/algorithm/area.h"
#include "SFCGAL/algorithm/length.h"
#include "SFCGAL/algorithm/centroid.h"
#include "SFCGAL/algorithm/distance.h"
#include "SFCGAL/algorithm/distance3d.h"
#include "SFCGAL/algorithm/volume.h"
#include "SFCGAL/algorithm/convexHull.h"
#include "SFCGAL/algorithm/intersects.h"
#include "SFCGAL/algorithm/intersection.h"
#include "SFCGAL/algorithm/BoundaryVisitor.h"

using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_NURBSCurveAlgorithmIntegrationTest)

BOOST_AUTO_TEST_CASE(testBasicAlgorithms)
{
  // Create a simple NURBS curve
  std::vector<Point> points = {
    Point(0.0, 0.0),
    Point(1.0, 1.0),
    Point(2.0, 0.0),
    Point(3.0, 1.0)
  };
  
  auto curve = NURBSCurve::interpolateCurve(points, 3);
  
  // Test area (should return 0 for curves)
  double area = algorithm::area(*curve);
  BOOST_CHECK_EQUAL(area, 0.0);
  
  // Test volume (should return 0 for curves)
  auto vol = algorithm::volume(*curve);
  BOOST_CHECK_EQUAL(vol, 0.0);
  
  // Test length (should return positive value)
  double length = algorithm::length(*curve);
  BOOST_CHECK(length > 0.0);
  
  // Test centroid (should return a point)
  auto centroidPoint = algorithm::centroid(*curve);
  BOOST_CHECK(!centroidPoint->isEmpty());
  
  // Test convex hull (should work via GetPointsVisitor)
  auto hull = algorithm::convexHull(*curve);
  BOOST_CHECK(!hull->isEmpty());
}

BOOST_AUTO_TEST_CASE(testDistanceAlgorithms)
{
  std::vector<Point> points1 = {Point(0.0, 0.0), Point(1.0, 0.0), Point(2.0, 0.0)};
  std::vector<Point> points2 = {Point(0.0, 1.0), Point(1.0, 1.0), Point(2.0, 1.0)};
  
  auto curve1 = NURBSCurve::interpolateCurve(points1, 2);
  auto curve2 = NURBSCurve::interpolateCurve(points2, 2);
  
  // Test 2D distance
  double dist2d = algorithm::distance(*curve1, *curve2);
  BOOST_CHECK(dist2d >= 0.0);
  BOOST_CHECK(dist2d <= 1.1); // Should be approximately 1.0
  
  // Test 3D distance
  double dist3d = algorithm::distance3D(*curve1, *curve2);
  BOOST_CHECK(dist3d >= 0.0);
  BOOST_CHECK(dist3d <= 1.1); // Should be approximately 1.0
}

BOOST_AUTO_TEST_CASE(testIntersectionAlgorithms)
{
  // Create two intersecting curves
  std::vector<Point> points1 = {Point(0.0, 0.0), Point(2.0, 2.0)};
  std::vector<Point> points2 = {Point(0.0, 2.0), Point(2.0, 0.0)};
  
  auto curve1 = NURBSCurve::interpolateCurve(points1, 1);
  auto curve2 = NURBSCurve::interpolateCurve(points2, 1);
  
  // Test intersects
  bool intersect = algorithm::intersects(*curve1, *curve2);
  BOOST_CHECK(intersect == true);
  
  // Test intersection (should return non-empty geometry)
  auto result = algorithm::intersection(*curve1, *curve2);
  BOOST_CHECK(!result->isEmpty());
}

BOOST_AUTO_TEST_CASE(testBoundaryVisitor)
{
  std::vector<Point> points = {Point(0.0, 0.0), Point(1.0, 1.0), Point(2.0, 0.0)};
  auto curve = NURBSCurve::interpolateCurve(points, 2);
  
  algorithm::BoundaryVisitor visitor;
  curve->accept(visitor);
  
  auto boundary = visitor.releaseBoundary();
  BOOST_CHECK(boundary != nullptr);
  // For an open curve, boundary should be a MultiPoint with start and end points
  BOOST_CHECK(!boundary->isEmpty());
}

BOOST_AUTO_TEST_SUITE_END()
