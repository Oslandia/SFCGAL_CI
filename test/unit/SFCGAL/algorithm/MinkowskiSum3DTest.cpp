// Copyright (c) 2012-2024, SFCGAL Contributors and Oslandia
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Sphere.h"
#include "SFCGAL/algorithm/area.h"
#include "SFCGAL/algorithm/minkowskiSum3D.h"
#include "SFCGAL/algorithm/volume.h"
#include "SFCGAL/detail/tools/Registry.h"
#include "SFCGAL/io/wkt.h"

using namespace SFCGAL;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_MinkowskiSum3DTest)

BOOST_AUTO_TEST_CASE(testMinkowskiSum3D_Sphere_Polyline)
{
  // Create a sphere
  Kernel::FT      radius(10.0);
  int             num_vertical   = 4;
  int             num_horizontal = 8;
  Kernel::Point_3 center(0, 0, 0);
  Sphere          sphere(radius, center, num_vertical, num_horizontal);
  std::unique_ptr<Geometry> sphereGeom(
      new PolyhedralSurface(sphere.generatePolyhedron()));

  // Create a polyline
  std::vector<Point> polyline_points = {Point(-100, 0, 0), Point(40, -70, 0),
                                        Point(40, 50, 40), Point(-90, -60, 60),
                                        Point(0, 0, -100), Point(30, 0, 150)};
  LineString         polyline(polyline_points);

  // Perform Minkowski sum
  std::unique_ptr<Geometry> result =
      algorithm::minkowskiSum3D(*sphereGeom, polyline);

  // Check that the result is not empty
  BOOST_CHECK(!result->isEmpty());

  // Check that the result is a PolyhedralSurface
  BOOST_CHECK_EQUAL(result->geometryTypeId(), TYPE_POLYHEDRALSURFACE);

  // Print result for visual inspection
  std::cout << "Minkowski sum result: " << result->asText(1) << std::endl;
}

BOOST_AUTO_TEST_CASE(testMinkowskiSum3D_Cube_Point)
{
  // Create a cube
  std::string cubeWkt =
      "SOLID((((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)), ((0 0 0, 0 0 1, 0 1 1, 0 "
      "1 0, 0 0 0)), ((0 0 0, 1 0 0, 1 0 1, 0 0 1, 0 0 0)), ((1 1 1, 0 1 1, 0 "
      "0 1, 1 0 1, 1 1 1)), ((1 1 1, 1 0 1, 1 0 0, 1 1 0, 1 1 1)), ((1 1 1, 1 "
      "1 0, 0 1 0, 0 1 1, 1 1 1))))";
  std::unique_ptr<Geometry> cube(io::readWkt(cubeWkt));

  // Create a point
  Point point(5, 5, 5);

  // Perform Minkowski sum
  std::unique_ptr<Geometry> result = algorithm::minkowskiSum3D(*cube, point);

  // Check that the result is not empty
  BOOST_CHECK(!result->isEmpty());

  // Check that the result is a PolyhedralSurface
  BOOST_CHECK_EQUAL(result->geometryTypeId(), TYPE_POLYHEDRALSURFACE);

  // Print result for visual inspection
  std::cout << "Minkowski sum result (Cube + Point): " << result->asText(1)
            << std::endl;
}

BOOST_AUTO_TEST_CASE(testMinkowskiSum3D_EmptyGeometries)
{
  // Create empty geometries
  GeometryCollection emptyGeom1;
  GeometryCollection emptyGeom2;

  // Perform Minkowski sum
  std::unique_ptr<Geometry> result =
      algorithm::minkowskiSum3D(emptyGeom1, emptyGeom2);

  // Check that the result is an empty GeometryCollection
  BOOST_CHECK(result->isEmpty());
  BOOST_CHECK_EQUAL(result->geometryTypeId(), TYPE_GEOMETRYCOLLECTION);

  // Print result for visual inspection
  std::cout << "Minkowski sum result (Empty geometries): " << result->asText(1)
            << std::endl;
}

BOOST_AUTO_TEST_CASE(testMinkowskiSum3D_Square_Polyline_2D)
{
  // Create a square
  std::string               squareWkt = "POLYGON((0 0, 0 1, 1 1, 1 0, 0 0))";
  std::unique_ptr<Geometry> square(io::readWkt(squareWkt));

  // Create a polyline
  std::string               polylineWkt = "LINESTRING(0 0, 1 1, 2 0)";
  std::unique_ptr<Geometry> polyline(io::readWkt(polylineWkt));

  // Perform Minkowski sum
  std::unique_ptr<Geometry> result =
      algorithm::minkowskiSum3D(*square, *polyline);

  // Check that the result is not empty
  BOOST_CHECK(!result->isEmpty());

  // Check that the result is a PolyhedralSurface (3D representation of the 2D
  // result)
  BOOST_CHECK_EQUAL(result->geometryTypeId(), TYPE_POLYHEDRALSURFACE);

  // Print result for visual inspection
  std::cout << "Minkowski sum result (Square + Polyline 2D): "
            << result->asText(1) << std::endl;
}

BOOST_AUTO_TEST_CASE(testMinkowskiSum3D_WKT_Cube_Point)
{
  std::string cubeWkt =
      "SOLID((((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)), ((0 0 0, 0 0 1, 0 1 1, 0 "
      "1 0, 0 0 0)), ((0 0 0, 1 0 0, 1 0 1, 0 0 1, 0 0 0)), ((1 1 1, 0 1 1, 0 "
      "0 1, 1 0 1, 1 1 1)), ((1 1 1, 1 0 1, 1 0 0, 1 1 0, 1 1 1)), ((1 1 1, 1 "
      "1 0, 0 1 0, 0 1 1, 1 1 1))))";
  std::string pointWkt = "POINT(2 2 2)";

  std::unique_ptr<Geometry> cube(io::readWkt(cubeWkt));
  std::unique_ptr<Geometry> point(io::readWkt(pointWkt));

  std::unique_ptr<Geometry> result = algorithm::minkowskiSum3D(*cube, *point);

  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK_EQUAL(result->geometryTypeId(), TYPE_POLYHEDRALSURFACE);

  std::cout << "Minkowski sum result (WKT Cube + Point): " << result->asText(1)
            << std::endl;
}

BOOST_AUTO_TEST_CASE(testMinkowskiSum3D_WKT_Tetrahedron_Segment)
{
  std::string tetrahedronWkt =
      "SOLID((((0 0 0, 1 0 0, 0 1 0, 0 0 0)), ((0 0 0, 0 1 0, 0 0 1, 0 0 0)), "
      "((0 0 0, 0 0 1, 1 0 0, 0 0 0)), ((1 0 0, 0 0 1, 0 1 0, 1 0 0))))";
  std::string segmentWkt = "LINESTRING(0 0 0, 1 1 1)";

  std::unique_ptr<Geometry> tetrahedron(io::readWkt(tetrahedronWkt));
  std::unique_ptr<Geometry> segment(io::readWkt(segmentWkt));

  BOOST_CHECK(tetrahedron != nullptr);
  BOOST_CHECK(segment != nullptr);

  BOOST_CHECK(!tetrahedron->isEmpty());
  BOOST_CHECK(!segment->isEmpty());

  BOOST_CHECK_EQUAL(tetrahedron->geometryTypeId(), TYPE_SOLID);
  BOOST_CHECK_EQUAL(segment->geometryTypeId(), TYPE_LINESTRING);

  std::unique_ptr<Geometry> result;
  BOOST_CHECK_NO_THROW(result =
                           algorithm::minkowskiSum3D(*tetrahedron, *segment));

  BOOST_CHECK(result != nullptr);
  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK_EQUAL(result->geometryTypeId(), TYPE_POLYHEDRALSURFACE);

  std::cout << "Minkowski sum result (WKT Tetrahedron + Segment): "
            << result->asText(1) << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()
