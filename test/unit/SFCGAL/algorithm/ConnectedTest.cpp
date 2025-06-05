// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

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
#include "SFCGAL/algorithm/connection.h"
#include "SFCGAL/io/wkt.h"

using namespace boost::unit_test;
using namespace SFCGAL;
using namespace SFCGAL::algorithm;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_Connected)

BOOST_AUTO_TEST_CASE(allFine)
{
  std::unique_ptr<Geometry> geom(
      io::readWkt("POLYHEDRALSURFACE (((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)),\
                                   ((0 0 0, 0 0 1, 0 1 1, 0 1 0, 0 0 0)),\
                                   ((0 0 0, 1 0 0, 1 0 1, 0 0 1, 0 0 0)),\
                                   ((1 1 1, 0 1 1, 0 0 1, 1 0 1, 1 1 1)),\
                                   ((1 1 1, 1 0 1, 1 0 0, 1 1 0, 1 1 1)),\
                                   ((1 1 1, 1 1 0, 0 1 0, 0 1 1, 1 1 1)))"));

  SurfaceGraph const graph(geom->as<PolyhedralSurface>());
  BOOST_CHECK_MESSAGE(isConnected(graph), "not connected");
  BOOST_CHECK_MESSAGE(isClosed(graph), "not closed");
}

BOOST_AUTO_TEST_CASE(notConnected)
{
  std::unique_ptr<Geometry> geom(io::readWkt(
      "POLYHEDRALSURFACE (((0 0 -1, 0 1 -1, 1 1 -1, 1 0 -1, 0 0 -1)),\
                                   ((0 0 0, 0 0 1, 0 1 1, 0 1 0, 0 0 0)),\
                                   ((0 0 0, 1 0 0, 1 0 1, 0 0 1, 0 0 0)),\
                                   ((1 1 1, 0 1 1, 0 0 1, 1 0 1, 1 1 1)),\
                                   ((1 1 1, 1 0 1, 1 0 0, 1 1 0, 1 1 1)),\
                                   ((1 1 1, 1 1 0, 0 1 0, 0 1 1, 1 1 1)))"));

  SurfaceGraph const graph(geom->as<PolyhedralSurface>());
  BOOST_CHECK_MESSAGE(!isConnected(graph), "connected");
  BOOST_CHECK_MESSAGE(!isClosed(graph), "closed");
}

BOOST_AUTO_TEST_CASE(notClosed)
{
  std::unique_ptr<Geometry> geom(
      io::readWkt("POLYHEDRALSURFACE (((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)),\
                                   ((0 0 0, 0 0 1, 0 1 1, 0 1 0, 0 0 0)),\
                                   ((0 0 0, 1 0 0, 1 0 1, 0 0 1, 0 0 0)),\
                                   ((1 1 1, 0 1 1, 0 0 1, 1 0 1, 1 1 1)),\
                                   ((1 1 1, 1 0 1, 1 0 0, 1 1 0, 1 1 1)))"));

  SurfaceGraph const graph(geom->as<PolyhedralSurface>());
  BOOST_CHECK_MESSAGE(isConnected(graph), "not connected");
  BOOST_CHECK_MESSAGE(!isClosed(graph), "closed");
}

BOOST_AUTO_TEST_SUITE_END()
