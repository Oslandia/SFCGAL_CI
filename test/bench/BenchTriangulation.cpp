// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <fstream>

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
#include "SFCGAL/detail/generator/disc.h"
#include "SFCGAL/detail/generator/hoch.h"
#include "SFCGAL/io/wkt.h"
#include "SFCGAL/triangulate/triangulate2DZ.h"
#include "SFCGAL/triangulate/triangulatePolygon.h"

#include "../test_config.h"
#include "Bench.h"

#include <boost/format.hpp>
#include <boost/test/unit_test.hpp>

#include <CGAL/point_generators_2.h>

using namespace boost::unit_test;
using namespace SFCGAL;
using namespace SFCGAL::triangulate;

#define N_POINTS 100000

typedef CGAL::Creator_uniform_2<double, Kernel::Point_2> Creator;

BOOST_AUTO_TEST_SUITE(SFCGAL_BenchTriangulation)

BOOST_AUTO_TEST_CASE(testMultiPointTriangulation)
{
  CGAL::Random_points_in_disc_2<Kernel::Point_2, Creator> g(150.0);

  MultiPoint multiPoint;

  for (int i = 0; i < N_POINTS; i++) {
    Kernel::Point_2 p = *(++g);
    multiPoint.addGeometry(new Point(p));
  }

  bench().start(boost::format("triangulate2DZ %s points") % N_POINTS);
  ConstraintDelaunayTriangulation cdt = triangulate2DZ(multiPoint);
  bench().stop();
}

BOOST_AUTO_TEST_CASE(testPolygonTriangulationHoch)
{
  const int                N = 7;
  std::unique_ptr<Polygon> fractal(generator::hoch(N));
  BOOST_CHECK_EQUAL(fractal->exteriorRing().numPoints(), 49153U);

  bench().start(boost::format("triangulate hoch(%s)") % N);
  TriangulatedSurface triangulatedSurface;
  SFCGAL::triangulate::triangulatePolygon3D(*fractal, triangulatedSurface);
  bench().stop();
}

BOOST_AUTO_TEST_CASE(testPolygonTriangulationHoch_roundingSixDecimal)
{
  const int                N = 7;
  std::unique_ptr<Polygon> fractal(generator::hoch(N));
  BOOST_CHECK_EQUAL(fractal->exteriorRing().numPoints(), 49153U);

  fractal->round(100000);

  bench().start(boost::format("triangulate hoch(%s) (round 6 six decimals)") %
                N);
  TriangulatedSurface triangulatedSurface;
  SFCGAL::triangulate::triangulatePolygon3D(*fractal, triangulatedSurface);
  bench().stop();
}

BOOST_AUTO_TEST_CASE(testPolygonTriangulationDisc)
{
  const int N = 20000;

  std::unique_ptr<Polygon> disc(generator::disc(Point(0.0, 0.0), 1.0, 8U));
  //	std::cout << fractal->asText(5) << std::endl ;

  bench().start(boost::format("triangulate disc x %s") % N);

  for (int i = 0; i < N; i++) {
    TriangulatedSurface triangulatedSurface;
    SFCGAL::triangulate::triangulatePolygon3D(*disc, triangulatedSurface);
  }

  bench().stop();
}

BOOST_AUTO_TEST_CASE(testMultiPointTriangulation2D)
{
  CGAL::Random_points_in_disc_2<Kernel::Point_2, Creator> g(150.0);

  MultiPoint multiPoint;

  for (int i = 0; i < N_POINTS; i++) {
    Kernel::Point_2 p = *(++g);
    multiPoint.addGeometry(new Point(p));
  }

  bench().start(boost::format("triangulate2D %s points") % N_POINTS);
  SFCGAL::triangulate::ConstraintDelaunayTriangulation cdt;
  SFCGAL::triangulate::triangulate2DZ(multiPoint, cdt);
  bench().stop();
}

BOOST_AUTO_TEST_CASE(testPolygonTriangulationHoch2D)
{
  const int                N = 7;
  std::unique_ptr<Polygon> fractal(generator::hoch(N));
  //	std::cout << fractal->asText(5) << std::endl ;

  bench().start(boost::format("triangulate2D hoch(%s)") % N);
  SFCGAL::triangulate::ConstraintDelaunayTriangulation cdt;
  SFCGAL::triangulate::triangulate2DZ(*fractal, cdt);
  bench().stop();
}

BOOST_AUTO_TEST_CASE(testPolygonTriangulationDisc2D)
{
  const int N = 20000;

  std::unique_ptr<Polygon> disc(generator::disc(Point(0.0, 0.0), 1.0, 8U));
  //	std::cout << fractal->asText(5) << std::endl ;

  bench().start(boost::format("triangulate2D disc x %s") % N);

  for (int i = 0; i < N; i++) {
    SFCGAL::triangulate::ConstraintDelaunayTriangulation cdt;
    SFCGAL::triangulate::triangulate2DZ(*disc, cdt);
  }

  bench().stop();
}

BOOST_AUTO_TEST_CASE(testPolygonTriangulationDisc2D_roundingSixDecimal)
{
  const int N = 20000;

  std::unique_ptr<Polygon> disc(generator::disc(Point(0.0, 0.0), 1.0, 8U));
  disc->round(100000);

  bench().start(
      boost::format("triangulate2D disc x %s (round 6 six decimals)") % N);

  for (int i = 0; i < N; i++) {
    SFCGAL::triangulate::ConstraintDelaunayTriangulation cdt;
    SFCGAL::triangulate::triangulate2DZ(*disc, cdt);
  }

  bench().stop();
}

BOOST_AUTO_TEST_SUITE_END()
