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

  // set new Polygon at index 1 from a Geometry object
  std::string const newPolygonStr = "POLYGON Z ((0 0 0, 0 5 0, 5 5 5, 5 0 5, 0 0 0))";
  std::unique_ptr<Geometry> newGeom(io::readWkt(newPolygonStr));
  geom->setGeometryN(newGeom->clone(), 1);

  BOOST_CHECK_EQUAL(geom->numGeometries(), 3);
  BOOST_CHECK_EQUAL(geom->geometryN(0).asText(0), "POLYGON Z ((0 0 0,10 0 0,10 10 0,0 10 0,0 0 0))");
  BOOST_CHECK_EQUAL(geom->geometryN(1).asText(), newGeom->asText());
  BOOST_CHECK_EQUAL(geom->geometryN(2).asText(0), "POLYGON Z ((0 0 0,0 10 0,5 5 5,0 0 0))");

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
}

BOOST_AUTO_TEST_SUITE_END()
