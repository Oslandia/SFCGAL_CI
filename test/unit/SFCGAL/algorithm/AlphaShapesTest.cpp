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

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/algorithm/alphaShapes.h"
#include "SFCGAL/io/wkt.h"

using namespace SFCGAL;

#include "../../../test_config.h"
// always after CGAL
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_AlphaShapesTest)

// algorithm::alphaShapes

BOOST_AUTO_TEST_CASE(testAlphaShapes2D_ComputeAlpha)
{
  std::vector<Point> points;
  points.emplace_back(0.0, 0.0);
  points.emplace_back(1.0, 0.0);
  points.emplace_back(1.0, 1.0);
  points.emplace_back(0.0, 1.0);

  LineString const                lineString(points);
  size_t const                    nb_comp = 3;
  std::unique_ptr<Geometry> const alphaShapes(
      algorithm::optimal_alpha_shapes(lineString, false, nb_comp));
}

BOOST_AUTO_TEST_CASE(testAlphaShapes2D_Empty)
{
  GeometryCollection collect;
  collect.addGeometry(Polygon());
  collect.addGeometry(Polygon());
  std::unique_ptr<Geometry> alphaShapes(algorithm::alphaShapes(collect));
  BOOST_CHECK(alphaShapes->isEmpty());
}

BOOST_AUTO_TEST_CASE(testAlphaShapes2D_ColinearProduceEmpty)
{
  LineString lineString;
  lineString.addPoint(Point(0.0, 0.0));
  lineString.addPoint(Point(1.0, 1.0));
  lineString.addPoint(Point(2.0, 2.0));

  std::unique_ptr<Geometry> alphaShapes(algorithm::alphaShapes(lineString));
  BOOST_CHECK(alphaShapes->isEmpty());
}

BOOST_AUTO_TEST_CASE(testAlphaShapes2D_Triangle)
{
  std::vector<Point> points;
  points.emplace_back(0.0, 0.0);
  points.emplace_back(0.5, 0.5);
  points.emplace_back(1.0, 0.0);
  points.emplace_back(0.0, 1.0);

  LineString const          lineString(points);
  std::unique_ptr<Geometry> alphaShapes(algorithm::alphaShapes(lineString));
  BOOST_CHECK(alphaShapes->is<Polygon>());
  std::string const expectedWkt =
      "POLYGON ((0.0 0.0,0.0 1.0,0.5 0.5,1.0 0.0,0.0 0.0))";
  BOOST_CHECK_EQUAL(alphaShapes->asText(1), expectedWkt);
}

BOOST_AUTO_TEST_CASE(testAlphaShapes2D_Polygon)
{
  std::vector<Point> points;
  points.emplace_back(0.0, 0.0);
  points.emplace_back(1.0, 0.0);
  points.emplace_back(1.0, 1.0);
  points.emplace_back(0.0, 1.0);

  LineString const          lineString(points);
  std::unique_ptr<Geometry> alphaShapes(algorithm::alphaShapes(lineString));
  BOOST_CHECK(alphaShapes->is<Polygon>());
  std::string const expectedWkt =
      "POLYGON ((0.0 0.0,0.0 1.0,1.0 1.0,1.0 0.0,0.0 0.0))";
  BOOST_CHECK_EQUAL(alphaShapes->asText(1), expectedWkt);
}

BOOST_AUTO_TEST_CASE(testAlphaShapes2D_MultiPoint)
{
  std::string inputData(SFCGAL_TEST_DIRECTORY);
  inputData += "/data/AlphaShapesWkt.txt";
  std::ifstream ifs(inputData.c_str());
  BOOST_REQUIRE(ifs.good());

  std::string expectedData(SFCGAL_TEST_DIRECTORY);
  expectedData += "/data/AlphaShapesWkt_expected.txt";
  std::ifstream efs(expectedData.c_str());
  BOOST_REQUIRE(efs.good());

  std::string expectedDataOptimal(SFCGAL_TEST_DIRECTORY);
  expectedDataOptimal += "/data/AlphaShapesWkt_expected_optimal.txt";
  std::ifstream efsOptimal(expectedDataOptimal.c_str());
  BOOST_REQUIRE(efsOptimal.good());

  std::string expectedDataOptimalHoles(SFCGAL_TEST_DIRECTORY);
  expectedDataOptimalHoles += "/data/AlphaShapesWkt_expected_optimal_holes.txt";
  std::ifstream efsOptimalHoles(expectedDataOptimalHoles.c_str());
  BOOST_REQUIRE(efsOptimalHoles.good());

  std::string inputWkt;
  std::string expectedWkt;
  std::string expectedWkt_optimal;
  std::string expectedWkt_optimal_holes;

  while (std::getline(ifs, inputWkt)) {
    std::unique_ptr<Geometry> g(io::readWkt(inputWkt));

    // expectedWkt
    std::getline(efs, expectedWkt);
    std::unique_ptr<Geometry> alphaShapes(
        algorithm::alphaShapes(g->as<const SFCGAL::Geometry>(), 1000));
    BOOST_CHECK_EQUAL(alphaShapes->asText(1), expectedWkt);

    // expectedWktOptimal
    std::getline(efsOptimal, expectedWkt_optimal);
    std::unique_ptr<Geometry> alphaShapesOptim(
        algorithm::optimal_alpha_shapes(g->as<const SFCGAL::Geometry>()));
    BOOST_CHECK_EQUAL(alphaShapesOptim->asText(1), expectedWkt_optimal);

    // expectedWktOptimalHoles
    std::getline(efsOptimalHoles, expectedWkt_optimal_holes);
    std::unique_ptr<Geometry> alphaShapesOptimHoles(
        algorithm::optimal_alpha_shapes(g->as<const SFCGAL::Geometry>(), true));
    BOOST_CHECK_EQUAL(alphaShapesOptimHoles->asText(1),
                      expectedWkt_optimal_holes);
  }
}

// https://gitlab.com/sfcgal/SFCGAL/-/issues/254
BOOST_AUTO_TEST_CASE(testAlphaShapes2D_InvalidPolygon_Issue254)
{
  std::vector<Point> points;
  points.emplace_back(1.0, 2.0);
  points.emplace_back(1.0, 2.0);
  points.emplace_back(1.0, 2.0);
  points.emplace_back(1.0, 2.0);

  LineString const lineString(points);
  BOOST_CHECK_THROW(algorithm::alphaShapes(lineString, 20.1, false),
                    std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
