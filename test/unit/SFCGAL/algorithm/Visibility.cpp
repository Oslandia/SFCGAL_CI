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
  points.emplace_back(0.0, 4.0);
  points.emplace_back(0.0, 0.0);
  points.emplace_back(3.0, 2.0);
  points.emplace_back(4.0, 0.0);
  points.emplace_back(4.0, 4.0);
  points.emplace_back(1.0, 2.0);
  points.emplace_back(0.0, 4.0);

  LineString const lineString(points);
  Polygon const    poly(lineString);

  Point const queryPoint(0.5, 2.0);

  std::unique_ptr<Polygon> result(algorithm::visibility(poly, queryPoint));
  std::string const        expectedWkt =
      "POLYGON((3.0 2.0,1.0 2.0,0.0 4.0,0.0 0.0,3.0 2.0))";
  BOOST_CHECK_EQUAL(result->asText(1), expectedWkt);
}

BOOST_AUTO_TEST_CASE(testVisibility_PointOnPolygon)
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

  Point queryPoint(0.0, 2.0);

  std::unique_ptr<Polygon> result(algorithm::visibility(poly, queryPoint));
  std::string              expectedWkt =
      "POLYGON((1.0 2.0,0.0 4.0,0.0 0.0,3.0 2.0,1.0 2.0))";
  BOOST_CHECK_EQUAL(result->asText(1), expectedWkt);
}

BOOST_AUTO_TEST_CASE(testVisibility_PointVertexOnPolygon)
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

  Point queryPoint(3.0, 2.0);

  std::unique_ptr<Polygon> result(algorithm::visibility(poly, queryPoint));
  std::string              expectedWkt =
      "POLYGON((0.0 0.0,3.0 2.0,4.0 0.0,4.0 4.0,1.0 2.0,0.0 2.0,0.0 0.0))";
  BOOST_CHECK_EQUAL(result->asText(1), expectedWkt);
}

BOOST_AUTO_TEST_CASE(testVisibility_PointInPolygonHole)
{
  std::vector<Point> points;
  points.emplace_back(0.0, 4.0);
  points.emplace_back(0.0, 0.0);
  points.emplace_back(3.0, 2.0);
  points.emplace_back(4.0, 0.0);
  points.emplace_back(4.0, 4.0);
  points.emplace_back(1.0, 2.0);
  points.emplace_back(0.0, 4.0);

  std::vector<Point> points_hole;
  points_hole.emplace_back(0.2, 1.75);
  points_hole.emplace_back(0.9, 1.8);
  points_hole.emplace_back(0.7, 1.2);
  points_hole.emplace_back(0.2, 1.75);

  LineString const lineString(points);
  LineString const hole(points_hole);

  Polygon poly(lineString);
  poly.addInteriorRing(hole);

  Point const queryPoint(0.5, 2.0);

  std::unique_ptr<Polygon> result(algorithm::visibility(poly, queryPoint));
  std::string const        expectedWkt =
      "POLYGON((0.0 1.6,0.2 1.8,0.9 1.8,1.9 1.3,3.0 "
      "2.0,1.0 2.0,0.0 4.0,0.0 1.6))";
  BOOST_CHECK_EQUAL(result->asText(1), expectedWkt);
}

BOOST_AUTO_TEST_CASE(testVisibility_PointOnPolygonHole)
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

  Point queryPoint(0.0, 2.0);

  std::unique_ptr<Polygon> result(algorithm::visibility(poly, queryPoint));
  std::string expectedWkt = "POLYGON((1.0 2.0,0.0 4.0,0.0 0.0,1.0 0.7,0.2 "
                            "1.8,0.9 1.8,2.2 1.5,3.0 2.0,1.0 2.0))";
  BOOST_CHECK_EQUAL(result->asText(1), expectedWkt);
}

BOOST_AUTO_TEST_CASE(testVisibility_PointVertexOnPolygonHole)
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

  Point queryPoint(3.0, 2.0);

  std::unique_ptr<Polygon> result(algorithm::visibility(poly, queryPoint));
  std::string              expectedWkt =
      "POLYGON((0.0 0.0,3.0 2.0,4.0 0.0,4.0 4.0,1.0 2.0,0.0 2.0,0.0 1.7,0.2 "
      "1.8,0.9 1.8,0.7 1.2,0.0 1.0,0.0 0.0))";
  BOOST_CHECK_EQUAL(result->asText(1), expectedWkt);
}

BOOST_AUTO_TEST_CASE(testVisibility_PointOnHolePolygonHole)
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

  Point queryPoint(0.550, 1.775);

  std::unique_ptr<Polygon> result(algorithm::visibility(poly, queryPoint));
  std::string expectedWkt = "POLYGON((0.7 1.2,0.9 1.8,0.2 1.8,0.7 1.2))";
  BOOST_CHECK_EQUAL(result->asText(1), expectedWkt);
}

BOOST_AUTO_TEST_CASE(testVisibility_PointVertexOnHolePolygonHole)
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

  Point queryPoint(0.9, 1.8);

  std::unique_ptr<Polygon> result(algorithm::visibility(poly, queryPoint));
  std::string expectedWkt = "POLYGON((0.7 1.2,0.9 1.8,0.2 1.8,0.7 1.2))";
  BOOST_CHECK_EQUAL(result->asText(1), expectedWkt);
}

BOOST_AUTO_TEST_CASE(testVisibility_SegmentInPolygon)
{
  std::vector<Point> points;
  points.emplace_back(0.0, 4.0);
  points.emplace_back(0.0, 0.0);
  points.emplace_back(3.0, 2.0);
  points.emplace_back(4.0, 0.0);
  points.emplace_back(4.0, 4.0);
  points.emplace_back(1.0, 2.0);
  points.emplace_back(0.0, 4.0);

  LineString const lineString(points);
  Polygon const    poly(lineString);

  Point const startPoint(1.0, 2.0);
  Point const endPoint(4.0, 4.0);

  std::unique_ptr<Polygon> result(
      algorithm::visibility(poly, startPoint, endPoint));
  std::string const expectedWkt =
      "POLYGON((4.0 0.0,4.0 4.0,1.0 2.0,0.0 1.3,0.0 0.0,3.0 2.0,4.0 0.0))";
  BOOST_CHECK_EQUAL(result->asText(1), expectedWkt);
}

BOOST_AUTO_TEST_CASE(testVisibility_SegmentInPolygonHole)
{
  std::vector<Point> points;
  points.emplace_back(1.0, 2.0);
  points.emplace_back(12.0, 3.0);
  points.emplace_back(19.0, -2.0);
  points.emplace_back(12.0, 6.0);
  points.emplace_back(14.0, 14.0);
  points.emplace_back(9.0, 5.0);
  points.emplace_back(1.0, 2.0);

  std::vector<Point> points_hole1;
  points_hole1.emplace_back(8.0, 3.0);
  points_hole1.emplace_back(8.0, 4.0);
  points_hole1.emplace_back(10.0, 3.0);
  points_hole1.emplace_back(8.0, 3.0);
  std::vector<Point> points_hole2;
  points_hole2.emplace_back(10.0, 6.0);
  points_hole2.emplace_back(11.0, 7.0);
  points_hole2.emplace_back(11.0, 6.0);
  points_hole2.emplace_back(10.0, 6.0);

  LineString const lineString(points);
  LineString const hole1(points_hole1);
  LineString const hole2(points_hole2);

  Polygon poly(lineString);
  poly.addInteriorRing(hole1);
  poly.addInteriorRing(hole2);

  Point const startPoint(19.0, -2.0);
  Point const endPoint(12.0, 6.0);

  std::unique_ptr<Polygon> result(
      algorithm::visibility(poly, startPoint, endPoint));
  std::string const expectedWkt =
      "POLYGON((19.0 -2.0,12.0 6.0,14.0 14.0,10.4 7.6,11.0 7.0,11.0 6.0,10.0 "
      "6.0,9.6 6.0,9.0 5.0,1.0 2.0,4.7 2.3,8.0 4.0,10.0 3.0,9.9 2.8,12.0 "
      "3.0,19.0 -2.0))";
  BOOST_CHECK_EQUAL(result->asText(1), expectedWkt);
}

BOOST_AUTO_TEST_CASE(testVisibility_PointOutPolygon)
{
  std::vector<Point> points;

  points.emplace_back(24.2222222, 40);
  points.emplace_back(24.183792760806462, 39.609819355967744);
  points.emplace_back(24.069981265022573, 39.23463313526982);
  points.emplace_back(23.88516142460509, 38.8888595339608);
  points.emplace_back(23.636435762373097, 38.58578643762691);
  points.emplace_back(23.333362666039207, 38.33706077539491);
  points.emplace_back(22.98758906473018, 38.15224093497743);
  points.emplace_back(22.612402844032257, 38.038429439193536);
  points.emplace_back(22.2222222, 38);
  points.emplace_back(21.832041555967745, 38.038429439193536);
  points.emplace_back(21.45685533526982, 38.15224093497743);
  points.emplace_back(21.1110817339608, 38.33706077539491);
  points.emplace_back(20.808008637626905, 38.58578643762691);
  points.emplace_back(20.55928297539491, 38.8888595339608);
  points.emplace_back(20.37446313497743, 39.23463313526982);
  points.emplace_back(20.26065163919354, 39.609819355967744);
  points.emplace_back(20.2222222, 40);
  points.emplace_back(20.26065163919354, 40.390180644032256);
  points.emplace_back(20.37446313497743, 40.76536686473018);
  points.emplace_back(20.55928297539491, 41.1111404660392);
  points.emplace_back(20.808008637626905, 41.41421356237309);
  points.emplace_back(21.111081733960795, 41.66293922460509);
  points.emplace_back(21.45685533526982, 41.84775906502257);
  points.emplace_back(21.832041555967745, 41.96157056080646);
  points.emplace_back(22.2222222, 42);
  points.emplace_back(22.612402844032257, 41.961570560806464);
  points.emplace_back(22.98758906473018, 41.84775906502257);
  points.emplace_back(23.333362666039204, 41.66293922460509);
  points.emplace_back(23.636435762373097, 41.41421356237309);
  points.emplace_back(23.88516142460509, 41.1111404660392);
  points.emplace_back(24.069981265022573, 40.76536686473018);
  points.emplace_back(24.183792760806462, 40.390180644032256);
  points.emplace_back(24.2222222, 40);

  LineString lineString(points);
  lineString.reverse();
  Polygon const poly(lineString);

  Point const queryPoint(-10, 60);

  try {
    std::unique_ptr<Polygon> const result(
        algorithm::visibility(poly, queryPoint));
  } catch (std::exception &e) {
    BOOST_CHECK_EQUAL(e.what(), "Can not find corresponding face.");
  }
}

BOOST_AUTO_TEST_SUITE_END()
