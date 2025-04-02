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

BOOST_AUTO_TEST_CASE(setGeometryNTest)
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

  // set new Polygon at index 1 from a Geometry object
  std::string const newPolygonStr = "POLYGON Z ((0 0 0, 0 5 0, 5 5 5, 5 0 5, 0 0 0))";
  std::unique_ptr<Geometry> newGeom(io::readWkt(newPolygonStr));
  geom->setGeometryN(newGeom->clone(), 1);

  BOOST_CHECK_EQUAL(geom->numGeometries(), 3);
  BOOST_CHECK_EQUAL(geom->geometryN(0).asText(0), "POLYGON Z ((0 0 0,10 0 0,10 10 0,0 10 0,0 0 0))");
  BOOST_CHECK_EQUAL(geom->geometryN(1).asText(), newGeom->asText());
  BOOST_CHECK_EQUAL(geom->geometryN(2).asText(0), "POLYGON Z ((0 0 0,0 10 0,5 5 5,0 0 0))");
  BOOST_CHECK_THROW(geom->geometryN(3), Exception);

  // set New Polygon at index 2 from a Polygon
  std::string const newPolygonStr2 = "POLYGON Z ((0 0 0, 0 15 0, 15 15 15, 15 0 15, 0 0 0))";
  std::unique_ptr<Geometry> newGeom2(io::readWkt(newPolygonStr2));
  Polygon *newPolygon2 = dynamic_cast<Polygon *>(newGeom2.get());
  BOOST_CHECK(newPolygon2);
  geom->setGeometryN(newPolygon2->clone(), 2);

  BOOST_CHECK_EQUAL(geom->numGeometries(), 3);
  BOOST_CHECK_EQUAL(geom->geometryN(0).asText(0), "POLYGON Z ((0 0 0,10 0 0,10 10 0,0 10 0,0 0 0))");
  BOOST_CHECK_EQUAL(geom->geometryN(1).asText(), newGeom->asText());
  BOOST_CHECK_EQUAL(geom->geometryN(2).asText(), newGeom2->asText());
  BOOST_CHECK_THROW(geom->geometryN(3), Exception);
}

BOOST_AUTO_TEST_CASE(dropZTest)
{
  std::unique_ptr<Geometry> emptyGeom(io::readWkt("POLYHEDRALSURFACE EMPTY"));
  BOOST_CHECK(!emptyGeom->is3D());
  BOOST_CHECK(!emptyGeom->dropZ());

  std::string const polyhedralStr =
    "POLYHEDRALSURFACE Z ("
    "((0 0 0, 10 0 0, 10 10 0, 0 10 0, 0 0 0)),"
    "((0 0 0, 10 0 0, 5 0 5, 0 0 0)),"
    "((0 0 0, 0 10 0, 5 5 5, 0 0 0))"
    ")";

  std::unique_ptr<Geometry> geom(io::readWkt(polyhedralStr));
  BOOST_CHECK(geom->is3D());
  BOOST_CHECK(geom->dropZ());
  BOOST_CHECK_EQUAL(geom->asText(0),
                    "POLYHEDRALSURFACE (((0 0,10 0,10 10,0 10,0 0)),"
                    "((0 0,10 0,5 0,0 0)),((0 0,0 10,5 5,0 0)))");
  BOOST_CHECK(!geom->is3D());
  BOOST_CHECK(!geom->dropZ());
}

BOOST_AUTO_TEST_SUITE_END()
