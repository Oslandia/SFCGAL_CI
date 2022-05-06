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

// algorithm::optimal_convex_partition

BOOST_AUTO_TEST_CASE(testPolygonPartitioning2D_Empty) {
  GeometryCollection collect;
  collect.addGeometry(Polygon());
  collect.addGeometry(Polygon());
  std::unique_ptr<Geometry> polygonPartitioning(algorithm::optimal_convex_partition(collect));
  BOOST_CHECK(polygonPartitioning->isEmpty());
}

// BOOST_AUTO_TEST_CASE(testPolygonPartitioning2D_Example) {
//   std::vector<Point> points;
//   points.push_back(Point(0.0, 0.0));
//   points.push_back(Point(0.5, 0.5));
//   points.push_back(Point(1.0, 0.0));
//   points.push_back(Point(0.0, 1.0));
//
//   LineString lineString(points);
//   std::unique_ptr<Geometry> polygonPartitioning(algorithm::optimal_convex_partition(lineString));
//   //BOOST_CHECK(polygonPartitioning->is<Polygon>());
//   std::string expectedWkt =
//       "POLYGON((0.0 0.0,0.0 1.0,0.5 0.5,1.0 0.0,0.0 0.0))";
//   //BOOST_CHECK_EQUAL(polygonPartitioning->asText(1), expectedWkt);
//   std::cout << polygonPartitioning->asText(1) << "\n";
// }

// BOOST_AUTO_TEST_CASE(testPolygonPartitioning2D_Polygon) {
//   std::vector<Point> points;
//      points.push_back(Point(391, 374));
//    points.push_back(Point(240, 431));
//    points.push_back(Point(252, 340));
//    points.push_back(Point(374, 320));
//    points.push_back(Point(289, 214));
//    points.push_back(Point(134, 390));
//    points.push_back(Point( 68, 186));
//    points.push_back(Point(154, 259));
//    points.push_back(Point(161, 107));
//    points.push_back(Point(435, 108));
//    points.push_back(Point(208, 148));
//    points.push_back(Point(295, 160));
//    points.push_back(Point(421, 212));
//    points.push_back(Point(441, 303));
//
//   LineString lineString(points);
//   std::unique_ptr<Geometry> polygonPartitioning(algorithm::optimal_convex_partition(lineString));
//   std::string expectedWkt =
//       "POLYGON((0.0 0.0,0.0 1.0,1.0 1.0,1.0 0.0,0.0 0.0))";
//   std::cout << polygonPartitioning->asText(1) << "\n";
//  //BOOST_CHECK_EQUAL(polygonPartitioning->asText(1), expectedWkt);
// }

BOOST_AUTO_TEST_CASE(testPolygonPartitioning2D_Greenland) {
  std::vector<Point> points;

  points.push_back(Point(-47, 83));
  points.push_back(Point(-45, 82));
  points.push_back(Point(-63, 82));
  points.push_back(Point(-68, 80));
  points.push_back(Point(-66, 79));
  points.push_back(Point(-73, 78));
  points.push_back(Point(-67, 77));
  points.push_back(Point(-71, 77));
  points.push_back(Point(-69, 76));
  points.push_back(Point(-59, 76));
  points.push_back(Point(-55, 73));
  points.push_back(Point(-56, 72));
  points.push_back(Point(-51, 71));
  points.push_back(Point(-54, 71));
  points.push_back(Point(-55, 70));
  points.push_back(Point(-51, 70));
  points.push_back(Point(-54, 67));
  points.push_back(Point(-52, 64));
  points.push_back(Point(-48, 61));
  points.push_back(Point(-43, 60));
  points.push_back(Point(-40, 65));
  points.push_back(Point(-22, 70));
  points.push_back(Point(-26, 70));
  points.push_back(Point(-26, 71));
  points.push_back(Point(-22, 71));
  points.push_back(Point(-25, 72));
  points.push_back(Point(-19, 74));
  points.push_back(Point(-22, 77));
  points.push_back(Point(-18, 77));
  points.push_back(Point(-20, 79));
  points.push_back(Point(-18, 80));
  points.push_back(Point(-20, 80));
  points.push_back(Point(-12, 81));
  points.push_back(Point(-32, 82));
  points.push_back(Point(-21, 83));
  points.push_back(Point(-27, 84));


  LineString lineString(points);
  std::unique_ptr<Geometry> polygonPartitioning(algorithm::optimal_convex_partition(lineString));
  //BOOST_CHECK(polygonPartitioning->is<Polygon>());
  std::string expectedWkt =
      "POLYGON((0.0 0.0,0.0 1.0,0.5 0.5,1.0 0.0,0.0 0.0))";
  //BOOST_CHECK_EQUAL(polygonPartitioning->asText(1), expectedWkt);
  std::cout << polygonPartitioning->asText(1) << "\n";
}

// BOOST_AUTO_TEST_CASE(testPolygonPartitioning2D_Greenland) {
//
//   std::cout << "ici\n";
//   std::string inputData(SFCGAL_TEST_DIRECTORY);
//   inputData += "/data/greenland.txt";
//   std::ifstream ifs(inputData.c_str());
//   BOOST_REQUIRE(ifs.good());
//
//   std::cout << "là\n";
//   std::string expectedData(SFCGAL_TEST_DIRECTORY);
//   expectedData += "/data/greenland_expected.txt";
//   std::ifstream efs(expectedData.c_str());
//   BOOST_REQUIRE(efs.good());
//
//   std::string inputWkt;
//   std::string expectedWkt;
//   
//   std::cout << "prout\n";
//   std::getline(ifs, inputWkt);
//   std::cout << "pas prout\n";
//   
//     std::unique_ptr<Geometry> g(io::readWkt(inputWkt));
//   BOOST_CHECK(g->is<Polygon>());
//
//     std::cout << "debug\n";
//     // expectedWkt
// //    std::getline(efs, expectedWkt);
// //
//     std::unique_ptr<Geometry> polygonPartitioning(
//         algorithm::optimal_convex_partition(g->as<const SFCGAL::Geometry>()));
//
//     std::cout << "algo done\n";
//   std::cout << polygonPartitioning->asText(1) << "\n";
//     // BOOST_CHECK_EQUAL(polygonPartitioning->asText(1), expectedWkt);
//
// }

BOOST_AUTO_TEST_SUITE_END()
