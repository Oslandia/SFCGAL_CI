/**
 *   SFCGAL
 *
 *   Copyright (C) 2012-2013 Oslandia <infos@oslandia.com>
 *   Copyright (C) 2012-2013 IGN (http://www.ign.fr)
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Library General Public
 *   License as published by the Free Software Foundation; either
 *   version 2 of the License, or (at your option) any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Library General Public License for more details.

 *   You should have received a copy of the GNU Library General Public
 *   License along with this library; if not, see
 <http://www.gnu.org/licenses/>.
 */
#include <boost/test/unit_test.hpp>

#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/algorithm/covers.h>
#include <SFCGAL/algorithm/visibility.h>

#include <SFCGAL/io/wkt.h>
using namespace SFCGAL;

#include "../../../test_config.h"
// always after CGAL
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_VisibilityTest)

// algorithm::alphaShapes

BOOST_AUTO_TEST_CASE(testVisibility_PointInPolygon)
{
  std::vector<Point> points;
  points.push_back(Point(0.0, 4.0));
  points.push_back(Point(0.0, 0.0));
  points.push_back(Point(3.0, 2.0));
  points.push_back(Point(4.0, 0.0));
  points.push_back(Point(4.0, 4.0));
  points.push_back(Point(1.0, 2.0));
  points.push_back(Point(0.0, 4.0));

  LineString lineString(points);
  Polygon    poly(lineString);

  Point queryPoint(0.5, 2.0);

  std::unique_ptr<Polygon> result(algorithm::visibility(poly, queryPoint));
  std::string              expectedWkt =
      "POLYGON((3.0 2.0,1.0 2.0,0.0 4.0,0.0 0.0,3.0 2.0))";
  BOOST_CHECK_EQUAL(result->asText(1), expectedWkt);
}

BOOST_AUTO_TEST_CASE(testVisibility_PointInPolygonHole)
{
  std::vector<Point> points;
  points.push_back(Point(0.0, 4.0));
  points.push_back(Point(0.0, 0.0));
  points.push_back(Point(3.0, 2.0));
  points.push_back(Point(4.0, 0.0));
  points.push_back(Point(4.0, 4.0));
  points.push_back(Point(1.0, 2.0));
  points.push_back(Point(0.0, 4.0));

  std::vector<Point> points_hole;
  points_hole.push_back(Point(0.2, 1.75));
  points_hole.push_back(Point(0.9, 1.8));
  points_hole.push_back(Point(0.7, 1.2));
  points_hole.push_back(Point(0.2, 1.75));

  LineString lineString(points);
  LineString hole(points_hole);

  Polygon poly(lineString);
  poly.addInteriorRing(hole);

  Point queryPoint(0.5, 2.0);

  std::unique_ptr<Polygon> result(algorithm::visibility(poly, queryPoint));
  std::string expectedWkt = "POLYGON((0.0 1.6,0.2 1.8,0.9 1.8,1.9 1.3,3.0 "
                            "2.0,1.0 2.0,0.0 4.0,0.0 1.6))";
  BOOST_CHECK_EQUAL(result->asText(1), expectedWkt);
}

BOOST_AUTO_TEST_CASE(testVisibility_SegmentInPolygon)
{
  std::vector<Point> points;
  points.push_back(Point(0.0, 4.0));
  points.push_back(Point(0.0, 0.0));
  points.push_back(Point(3.0, 2.0));
  points.push_back(Point(4.0, 0.0));
  points.push_back(Point(4.0, 4.0));
  points.push_back(Point(1.0, 2.0));
  points.push_back(Point(0.0, 4.0));

  LineString lineString(points);
  Polygon    poly(lineString);

  Point startPoint(1.0, 2.0);
  Point endPoint(4.0, 4.0);

  std::unique_ptr<Polygon> result(
      algorithm::visibility(poly, startPoint, endPoint));
  std::string expectedWkt =
      "POLYGON((4.0 0.0,4.0 4.0,1.0 2.0,0.0 1.3,0.0 0.0,3.0 2.0,4.0 0.0))";
  BOOST_CHECK_EQUAL(result->asText(1), expectedWkt);
}

BOOST_AUTO_TEST_CASE(testVisibility_SegmentInPolygonHole)
{
  std::vector<Point> points;
  points.push_back(Point(1.0, 2.0));
  points.push_back(Point(12.0, 3.0));
  points.push_back(Point(19.0, -2.0));
  points.push_back(Point(12.0, 6.0));
  points.push_back(Point(14.0, 14.0));
  points.push_back(Point(9.0, 5.0));
  points.push_back(Point(1.0, 2.0));

  std::vector<Point> points_hole1;
  points_hole1.push_back(Point(8.0, 3.0));
  points_hole1.push_back(Point(8.0, 4.0));
  points_hole1.push_back(Point(10.0, 3.0));
  points_hole1.push_back(Point(8.0, 3.0));
  std::vector<Point> points_hole2;
  points_hole2.push_back(Point(10.0, 6.0));
  points_hole2.push_back(Point(11.0, 7.0));
  points_hole2.push_back(Point(11.0, 6.0));
  points_hole2.push_back(Point(10.0, 6.0));

  LineString lineString(points);
  LineString hole1(points_hole1);
  LineString hole2(points_hole2);

  Polygon poly(lineString);
  poly.addInteriorRing(hole1);
  poly.addInteriorRing(hole2);

  Point startPoint(19.0, -2.0);
  Point endPoint(12.0, 6.0);

  std::unique_ptr<Polygon> result(
      algorithm::visibility(poly, startPoint, endPoint));
  std::string expectedWkt =
      "POLYGON((19.0 -2.0,12.0 6.0,14.0 14.0,10.4 7.6,11.0 7.0,11.0 6.0,10.0 "
      "6.0,9.6 6.0,9.0 5.0,1.0 2.0,4.7 2.3,8.0 4.0,10.0 3.0,9.9 2.8,12.0 "
      "3.0,19.0 -2.0))";
  BOOST_CHECK_EQUAL(result->asText(1), expectedWkt);
}
BOOST_AUTO_TEST_SUITE_END()
