// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Kernel.h"
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
#include "SFCGAL/algorithm/normal.h"
#include "SFCGAL/io/wkt.h"

using namespace SFCGAL;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_NormalTest)

BOOST_AUTO_TEST_CASE(testNormal1)
{
  Point_3 const a(0.0, 0.0, 0.0);
  Point_3 const b(1.0, 0.0, 0.0);
  Point_3 const c(1.0, 1.0, 0.0);

  Vector_3 const normal = algorithm::normal3D(a, b, c);
  BOOST_CHECK_EQUAL(normal.x(), 0.0);
  BOOST_CHECK_EQUAL(normal.y(), 0.0);
  BOOST_CHECK_EQUAL(normal.z(), 1.0);
}

BOOST_AUTO_TEST_CASE(testNormal2)
{
  // a square ccw
  std::unique_ptr<Geometry> gA(io::readWkt("POLYGON ((0 0,1 0,1 1,0 1,0 0))"));
  // a square cw oriented
  std::unique_ptr<Geometry> gB(io::readWkt("POLYGON ((0 0,0 1,1 1,1 0,0 0))"));

  // a pseudo-square ccw oriented, with a concave part
  std::unique_ptr<Geometry> gC(
      io::readWkt("POLYGON ((0 0,0.5 0.5,1 0,1 1,0 1,0 0))"));

  {
    CGAL::Vector_3<Kernel> const normal =
        algorithm::normal3D<Kernel>(gA->as<Polygon>());
    BOOST_CHECK_EQUAL(normal.x(), 0.0);
    BOOST_CHECK_EQUAL(normal.y(), 0.0);
    BOOST_CHECK_EQUAL(normal.z(), 2.0);
  }

  {
    CGAL::Vector_3<Kernel> const normal =
        algorithm::normal3D<Kernel>(gB->as<Polygon>());
    BOOST_CHECK_EQUAL(normal.x(), 0.0);
    BOOST_CHECK_EQUAL(normal.y(), 0.0);
    BOOST_CHECK_EQUAL(normal.z(), -2.0);
  }

  {
    CGAL::Vector_3<Kernel> const normal =
        algorithm::normal3D<Kernel>(gC->as<Polygon>());
    BOOST_CHECK_EQUAL(normal.x(), 0.0);
    BOOST_CHECK_EQUAL(normal.y(), 0.0);
    // ok, the normal is pointing up (z > 0)
    BOOST_CHECK_EQUAL(normal.z(), 1.5);
  }
}

BOOST_AUTO_TEST_CASE(testNormal3)
{
  std::unique_ptr<Geometry> gA(
      io::readWkt("POLYGON ((0 1 0,0 1 1,1 1 1,1 1 0,0 1 0))"));
  // exact
  {
    CGAL::Vector_3<Kernel> const normal =
        algorithm::normal3D<Kernel>(gA->as<Polygon>(), true);
    // std::cout << CGAL::exact(normal) << std::endl;
    CGAL::Plane_3<Kernel> const plane(
        gA->as<Polygon>().exteriorRing().startPoint().toPoint_3(), normal);
    // std::cout << CGAL::exact(plane) << std::endl;
    BOOST_CHECK(!plane.is_degenerate());
  }
  // round
  {
    CGAL::Vector_3<Kernel> const normal =
        algorithm::normal3D<Kernel>(gA->as<Polygon>(), false);
    // std::cout << CGAL::exact(normal) << std::endl;
    CGAL::Plane_3<Kernel> const plane(
        gA->as<Polygon>().exteriorRing().startPoint().toPoint_3(), normal);
    // std::cout << CGAL::exact(plane) << std::endl;
    BOOST_CHECK(!plane.is_degenerate());
  }
}

BOOST_AUTO_TEST_SUITE_END()
