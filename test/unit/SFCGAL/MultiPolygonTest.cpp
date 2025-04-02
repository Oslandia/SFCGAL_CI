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

#include <exception>

#include "SFCGAL/Envelope.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/Polygon.h"

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_MultiPolygonTest)

BOOST_AUTO_TEST_CASE(defaultConstructor)
{
  MultiPolygon const g;
  BOOST_CHECK(g.isEmpty());
  BOOST_CHECK(!g.is3D());
  BOOST_CHECK_EQUAL(g.numGeometries(), 0U);
}

BOOST_AUTO_TEST_CASE(testGeometryTypeId)
{
  MultiPolygon const g;
  BOOST_CHECK_EQUAL(g.geometryTypeId(), TYPE_MULTIPOLYGON);
}

//-- addAllowedGeometry
BOOST_AUTO_TEST_CASE(addPolygon)
{
  MultiPolygon g;
  g.addGeometry(new Polygon());
  BOOST_CHECK_EQUAL(g.numGeometries(), 1U);
  g.addGeometry(new Polygon());
  BOOST_CHECK_EQUAL(g.numGeometries(), 2U);
}
//-- addForbidenGeometry
BOOST_AUTO_TEST_CASE(addLineStringThrow)
{
  MultiPolygon g;
  BOOST_CHECK_THROW(g.addGeometry(LineString()), std::exception);
}

//-- asText

BOOST_AUTO_TEST_CASE(asTextEmpty)
{
  MultiPolygon const g;
  BOOST_CHECK_EQUAL(g.asText(1), "MULTIPOLYGON EMPTY");
}
BOOST_AUTO_TEST_CASE(asText2d)
{
  MultiPolygon g;
  g.addGeometry(Envelope(0.0, 1.0, 0.0, 1.0).toPolygon().release());
  g.addGeometry(Envelope(2.0, 3.0, 4.0, 5.0).toPolygon().release());
  BOOST_CHECK_EQUAL(g.asText(3),
                    "MULTIPOLYGON (((0.000 0.000,1.000 0.000,1.000 1.000,0.000 "
                    "1.000,0.000 0.000)),((2.000 4.000,3.000 4.000,3.000 "
                    "5.000,2.000 5.000,2.000 4.000)))");
}

BOOST_AUTO_TEST_CASE(dropZ)
{
  MultiPolygon emptyMultiPolygon;
  BOOST_CHECK(!emptyMultiPolygon.is3D());
  BOOST_CHECK(!emptyMultiPolygon.dropZ());

  MultiPolygon multiPolygon;
  Polygon poly1;
  poly1.exteriorRing().addPoint(Point(0.0, 0.0, 2.0));
  poly1.exteriorRing().addPoint(Point(1.0, 0.0, 2.0));
  poly1.exteriorRing().addPoint(Point(1.0, 1.0, 2.0));
  poly1.exteriorRing().addPoint(Point(0.0, 1.0, 2.0));
  poly1.exteriorRing().addPoint(Point(0.0, 0.0, 2.0));

  Polygon poly2;
  poly2.exteriorRing().addPoint(Point(3.0, 5.0, 3.0));
  poly2.exteriorRing().addPoint(Point(4.0, 5.0, 3.0));
  poly2.exteriorRing().addPoint(Point(4.0, 6.0, 3.0));
  poly2.exteriorRing().addPoint(Point(3.0, 6.0, 3.0));
  poly2.exteriorRing().addPoint(Point(3.0, 5.0, 3.0));

  multiPolygon.addGeometry(poly1);
  multiPolygon.addGeometry(poly2);
  BOOST_CHECK(multiPolygon.is3D());
  BOOST_CHECK(multiPolygon.dropZ());

  BOOST_CHECK_EQUAL(multiPolygon.asText(1),
                    "MULTIPOLYGON (((0.0 0.0,1.0 0.0,1.0 1.0,0.0 1.0,0.0 0.0)),"
                    "((3.0 5.0,4.0 5.0,4.0 6.0,3.0 6.0,3.0 5.0)))");

  BOOST_CHECK(!multiPolygon.is3D());
  BOOST_CHECK(!multiPolygon.dropZ());
}

//-- is< T >

BOOST_AUTO_TEST_CASE(isGeometryCollection)
{
  MultiPolygon const g;
  BOOST_CHECK(g.is<GeometryCollection>());
}

BOOST_AUTO_TEST_CASE(isMultiPolygon)
{
  MultiPolygon const g;
  BOOST_CHECK(g.is<MultiPolygon>());
}

BOOST_AUTO_TEST_SUITE_END()
