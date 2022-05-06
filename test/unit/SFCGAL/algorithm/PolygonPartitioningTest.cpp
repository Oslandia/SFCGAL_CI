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

#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/LineString.h>
#include <SFCGAL/MultiLineString.h>
#include <SFCGAL/MultiPoint.h>
#include <SFCGAL/MultiPolygon.h>
#include <SFCGAL/MultiSolid.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Solid.h>
#include <SFCGAL/Triangle.h>
#include <SFCGAL/TriangulatedSurface.h>
#include <SFCGAL/algorithm/polygonPartitioning.h>
#include <SFCGAL/io/wkt.h>

using namespace SFCGAL;

#include "../../../test_config.h"
// always after CGAL
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_PolygonPartitioningTest)

// algorithm::greene_approx_convex_partition

BOOST_AUTO_TEST_CASE(testPolygonPartitioning2D_Empty)
{
  GeometryCollection collect;
  collect.addGeometry(Polygon());
  collect.addGeometry(Polygon());
  std::unique_ptr<Geometry> polygonPartitioning(
      algorithm::greene_approx_convex_partition(collect));
  BOOST_CHECK(polygonPartitioning->isEmpty());
}

// Fail: CGAL error: precondition violation!
// Expression : is_simple_2(first, last, traits)
// BOOST_AUTO_TEST_CASE(testPolygonPartitioning2D_Example) {
//   std::vector<Point> points;
//   points.push_back(Point(0.0, 0.0));
//   points.push_back(Point(0.5, 0.5));
//   points.push_back(Point(1.0, 0.0));
//   points.push_back(Point(0.0, 1.0));
//
//   LineString lineString(points);
//   std::unique_ptr<Geometry>
//   polygonPartitioning(algorithm::greene_approx_convex_partition(lineString));
//   //BOOST_CHECK(polygonPartitioning->is<Polygon>());
//   std::string expectedWkt =
//       "POLYGON((0.0 0.0,0.0 1.0,0.5 0.5,1.0 0.0,0.0 0.0))";
//   //BOOST_CHECK_EQUAL(polygonPartitioning->asText(1), expectedWkt);
//   std::cout << polygonPartitioning->asText(1) << "\n";
// }

BOOST_AUTO_TEST_CASE(testPolygonPartitioning2D_Polygon)
{
  std::vector<Point> points;
  points.push_back(Point(391, 374));
  points.push_back(Point(240, 431));
  points.push_back(Point(252, 340));
  points.push_back(Point(374, 320));
  points.push_back(Point(289, 214));
  points.push_back(Point(134, 390));
  points.push_back(Point(68, 186));
  points.push_back(Point(154, 259));
  points.push_back(Point(161, 107));
  points.push_back(Point(435, 108));
  points.push_back(Point(208, 148));
  points.push_back(Point(295, 160));
  points.push_back(Point(421, 212));
  points.push_back(Point(441, 303));

  LineString                lineString(points);
  std::unique_ptr<Geometry> polygonPartitioning(
      algorithm::greene_approx_convex_partition(lineString));
  std::string expectedWkt =
      "MULTIPOLYGON(((134.0 390.0,68.0 186.0,154.0 259.0,134.0 390.0)),"
      "((161.0 107.0,435.0 108.0,208.0 148.0,161.0 107.0)),"
      "((208.0 148.0,295.0 160.0,421.0 212.0,289.0 214.0,208.0 148.0)),"
      "((154.0 259.0,161.0 107.0,208.0 148.0,154.0 259.0)),"
      "((289.0 214.0,134.0 390.0,154.0 259.0,208.0 148.0,289.0 214.0)),"
      "((374.0 320.0,289.0 214.0,421.0 212.0,374.0 320.0)),"
      "((374.0 320.0,421.0 212.0,441.0 303.0,391.0 374.0,374.0 320.0)),"
      "((391.0 374.0,240.0 431.0,252.0 340.0,374.0 320.0,391.0 374.0)))";
  BOOST_CHECK_EQUAL(polygonPartitioning->asText(1), expectedWkt);
}

BOOST_AUTO_TEST_CASE(testPolygonPartitioning2D_Greenland_polygon)
{

  std::string inputData(SFCGAL_TEST_DIRECTORY);
  inputData += "/data/greenland.txt";
  std::ifstream ifs(inputData.c_str());
  BOOST_REQUIRE(ifs.good());

  std::string expectedData(SFCGAL_TEST_DIRECTORY);
  expectedData += "/data/greenland_expected.txt";
  std::ifstream efs(expectedData.c_str());
  BOOST_REQUIRE(efs.good());

  std::string inputWkt;
  std::string expectedWkt;

  std::getline(ifs, inputWkt);

  std::unique_ptr<Geometry> g(io::readWkt(inputWkt));
  BOOST_CHECK(g->is<Polygon>());

  // expectedWkt
  std::getline(efs, expectedWkt);

  std::unique_ptr<Geometry> polygonPartitioning(
      algorithm::greene_approx_convex_partition(
          g->as<const SFCGAL::Geometry>()));

  BOOST_CHECK_EQUAL(polygonPartitioning->asText(1), expectedWkt);
}

BOOST_AUTO_TEST_SUITE_END()
