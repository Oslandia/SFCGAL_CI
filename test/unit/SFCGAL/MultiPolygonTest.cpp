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
#include "SFCGAL/io/wkt.h"

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

BOOST_AUTO_TEST_CASE(dropZM)
{
  MultiPolygon emptyMultiPolygon;
  BOOST_CHECK(!emptyMultiPolygon.is3D());
  BOOST_CHECK(!emptyMultiPolygon.isMeasured());
  BOOST_CHECK(!emptyMultiPolygon.dropM());
  BOOST_CHECK(!emptyMultiPolygon.dropZ());

  MultiPolygon multiPolygon2D;
  Polygon poly12D;
  poly12D.exteriorRing().addPoint(Point(0.0, 0.0));
  poly12D.exteriorRing().addPoint(Point(1.0, 0.0));
  poly12D.exteriorRing().addPoint(Point(1.0, 1.0));
  poly12D.exteriorRing().addPoint(Point(0.0, 1.0));
  poly12D.exteriorRing().addPoint(Point(0.0, 0.0));

  Polygon poly22D;
  poly22D.exteriorRing().addPoint(Point(3.0, 5.0));
  poly22D.exteriorRing().addPoint(Point(4.0, 5.0));
  poly22D.exteriorRing().addPoint(Point(4.0, 6.0));
  poly22D.exteriorRing().addPoint(Point(3.0, 6.0));
  poly22D.exteriorRing().addPoint(Point(3.0, 5.0));

  multiPolygon2D.addGeometry(poly12D);
  multiPolygon2D.addGeometry(poly22D);
  BOOST_CHECK(!multiPolygon2D.is3D());
  BOOST_CHECK(!multiPolygon2D.isMeasured());
  BOOST_CHECK(!multiPolygon2D.dropM());
  BOOST_CHECK(!multiPolygon2D.dropZ());

  MultiPolygon multiPolygon3D;
  Polygon poly13D;
  poly13D.exteriorRing().addPoint(Point(0.0, 0.0, 2.0));
  poly13D.exteriorRing().addPoint(Point(1.0, 0.0, 2.0));
  poly13D.exteriorRing().addPoint(Point(1.0, 1.0, 2.0));
  poly13D.exteriorRing().addPoint(Point(0.0, 1.0, 2.0));
  poly13D.exteriorRing().addPoint(Point(0.0, 0.0, 2.0));

  Polygon poly23D;
  poly23D.exteriorRing().addPoint(Point(3.0, 5.0, 3.0));
  poly23D.exteriorRing().addPoint(Point(4.0, 5.0, 3.0));
  poly23D.exteriorRing().addPoint(Point(4.0, 6.0, 3.0));
  poly23D.exteriorRing().addPoint(Point(3.0, 6.0, 3.0));
  poly23D.exteriorRing().addPoint(Point(3.0, 5.0, 3.0));

  multiPolygon3D.addGeometry(poly13D);
  multiPolygon3D.addGeometry(poly23D);
  BOOST_CHECK(multiPolygon3D.is3D());
  BOOST_CHECK(!multiPolygon3D.isMeasured());
  BOOST_CHECK(!multiPolygon3D.dropM());
  BOOST_CHECK(multiPolygon3D.dropZ());
  BOOST_CHECK_EQUAL(multiPolygon3D.asText(1),
                    "MULTIPOLYGON (((0.0 0.0,1.0 0.0,1.0 1.0,0.0 1.0,0.0 0.0)),"
                    "((3.0 5.0,4.0 5.0,4.0 6.0,3.0 6.0,3.0 5.0)))");
  BOOST_CHECK(!multiPolygon3D.isMeasured());
  BOOST_CHECK(!multiPolygon3D.dropM());
  BOOST_CHECK(!multiPolygon3D.is3D());
  BOOST_CHECK(!multiPolygon3D.dropZ());

  MultiPolygon multiPolygonM;
  multiPolygonM.addGeometry(io::readWkt("POLYGON M ((0 0 1, 0 3 2, 3 3 3, 3 0 4, 0 0 1))").release());
  multiPolygonM.addGeometry(io::readWkt("POLYGON M ((4 5 1, 1 4 2, 4 0 4, 4 5 1))").release());
  BOOST_CHECK(!multiPolygonM.is3D());
  BOOST_CHECK(multiPolygonM.isMeasured());
  BOOST_CHECK(!multiPolygonM.dropZ());
  BOOST_CHECK(multiPolygonM.dropM());
  BOOST_CHECK_EQUAL(multiPolygonM.asText(0), "MULTIPOLYGON (((0 0,0 3,3 3,3 0,0 0)),((4 5,1 4,4 0,4 5)))");
  BOOST_CHECK(!multiPolygonM.is3D());
  BOOST_CHECK(!multiPolygonM.isMeasured());
  BOOST_CHECK(!multiPolygonM.dropM());
  BOOST_CHECK(!multiPolygonM.dropZ());

  MultiPolygon multiPolygonZM;
  multiPolygonZM.addGeometry(io::readWkt("POLYGON ZM ((0 0 1 4, 0 3 2 5, 3 3 3 6, 3 0 4 7, 0 0 1 4))").release());
  multiPolygonZM.addGeometry(io::readWkt("POLYGON ZM ((2 4 1 5, 7 7 7 6, 4 2 5 7, 2 4 1 5))").release());
  BOOST_CHECK(multiPolygonZM.is3D());
  BOOST_CHECK(multiPolygonZM.isMeasured());
  BOOST_CHECK(multiPolygonZM.dropM());
  BOOST_CHECK(multiPolygonZM.is3D());
  BOOST_CHECK(!multiPolygonZM.isMeasured());
  BOOST_CHECK_EQUAL(multiPolygonZM.asText(0), "MULTIPOLYGON Z (((0 0 1,0 3 2,3 3 3,3 0 4,0 0 1)),((2 4 1,7 7 7,4 2 5,2 4 1)))");
  BOOST_CHECK(!multiPolygonZM.dropM());
  BOOST_CHECK(multiPolygonZM.dropZ());
  BOOST_CHECK_EQUAL(multiPolygonZM.asText(0), "MULTIPOLYGON (((0 0,0 3,3 3,3 0,0 0)),((2 4,7 7,4 2,2 4)))");
  BOOST_CHECK(!multiPolygonZM.is3D());
  BOOST_CHECK(!multiPolygonZM.isMeasured());
  BOOST_CHECK(!multiPolygonZM.dropM());
  BOOST_CHECK(!multiPolygonZM.dropZ());
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
