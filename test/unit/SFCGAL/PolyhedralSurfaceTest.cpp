// Copyright (c) 2025-2025, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/Geometry.h"
#include "SFCGAL/Kernel.h"

#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/io/wkt.h"

#include <boost/test/unit_test.hpp>
#include <memory>
using namespace boost::unit_test;

using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_PolyhedralSurfaceTest)

BOOST_AUTO_TEST_CASE(geometryNTest)
{
  std::unique_ptr<Geometry> emptyGeom(io::readWkt("POLYHEDRALSURFACE EMPTY"));
  BOOST_CHECK(emptyGeom->is<PolyhedralSurface>());
  BOOST_CHECK(emptyGeom->isEmpty());
  BOOST_CHECK_EQUAL(emptyGeom->numGeometries(), 0);
  BOOST_CHECK_THROW(emptyGeom->geometryN(0), Exception);

  std::string const polyhedralStr =
    "POLYHEDRALSURFACE Z ("
    "((0 0 0, 10 0 0, 10 10 0, 0 10 0, 0 0 0)),"
    "((0 0 0, 10 0 0, 5 0 5, 0 0 0)),"
    "((0 0 0, 0 10 0, 5 5 5, 0 0 0))"
    ")";

  std::unique_ptr<Geometry> geom(io::readWkt(polyhedralStr));
  BOOST_CHECK(!geom->isEmpty());
  BOOST_CHECK_EQUAL(geom->numGeometries(), 3);
  BOOST_CHECK_EQUAL(geom->geometryN(0).asText(0), "POLYGON Z ((0 0 0,10 0 0,10 10 0,0 10 0,0 0 0))");
  BOOST_CHECK_EQUAL(geom->geometryN(1).asText(0), "POLYGON Z ((0 0 0,10 0 0,5 0 5,0 0 0))");
  BOOST_CHECK_EQUAL(geom->geometryN(2).asText(0), "POLYGON Z ((0 0 0,0 10 0,5 5 5,0 0 0))");
  BOOST_CHECK_THROW(geom->geometryN(3), Exception);
}

BOOST_AUTO_TEST_SUITE_END()
