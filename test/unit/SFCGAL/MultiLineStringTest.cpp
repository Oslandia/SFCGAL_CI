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

#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/io/wkt.h"

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_MultiLineStringTest)

BOOST_AUTO_TEST_CASE(defaultConstructor)
{
  MultiLineString const g;
  BOOST_CHECK(g.isEmpty());
  BOOST_CHECK(!g.is3D());
  BOOST_CHECK_EQUAL(g.numGeometries(), 0U);
}

//-- addAllowedGeometry
BOOST_AUTO_TEST_CASE(addLineString)
{
  MultiLineString g;
  g.addGeometry(new LineString());
  BOOST_CHECK_EQUAL(g.numGeometries(), 1U);

  g.addGeometry(LineString(Point(0.0, 0.0), Point(1.0, 1.0)));
  BOOST_CHECK_EQUAL(g.numGeometries(), 2U);
}
//-- addForbidenGeometry
BOOST_AUTO_TEST_CASE(addLineStringThrow)
{
  MultiLineString g;
  BOOST_CHECK_THROW(g.addGeometry(Point()), std::exception);
}

//-- asText

BOOST_AUTO_TEST_CASE(asTextEmpty)
{
  MultiLineString const g;
  BOOST_CHECK_EQUAL(g.asText(1), "MULTILINESTRING EMPTY");
}

BOOST_AUTO_TEST_CASE(asText2d)
{
  MultiLineString g;
  g.addGeometry(LineString(Point(0.0, 0.0), Point(1.0, 1.0)));
  g.addGeometry(LineString(Point(1.0, 1.0), Point(2.0, 2.0)));
  BOOST_CHECK_EQUAL(g.asText(1),
                    "MULTILINESTRING ((0.0 0.0,1.0 1.0),(1.0 1.0,2.0 2.0))");
}

//-- is< T >

BOOST_AUTO_TEST_CASE(isGeometryCollection)
{
  MultiLineString const g;
  BOOST_CHECK(g.is<GeometryCollection>());
}

BOOST_AUTO_TEST_CASE(isMultiLineString)
{
  MultiLineString const g;
  BOOST_CHECK(g.is<MultiLineString>());
}

BOOST_AUTO_TEST_CASE(dropZM)
{
  MultiLineString multiLineStringEmpty;
  BOOST_CHECK(multiLineStringEmpty.isEmpty());
  BOOST_CHECK(!multiLineStringEmpty.is3D());
  BOOST_CHECK(!multiLineStringEmpty.isMeasured());
  BOOST_CHECK(!multiLineStringEmpty.dropZ());
  BOOST_CHECK(!multiLineStringEmpty.dropM());

  MultiLineString multiLineString2D;
  multiLineString2D.addGeometry(new LineString(Point(2.0, 3.0), Point(4.0, 5.0)));
  multiLineString2D.addGeometry(new LineString(Point(6.0, 7.0), Point(9.0, 10.0)));
  BOOST_CHECK(!multiLineString2D.is3D());
  BOOST_CHECK(!multiLineString2D.isMeasured());
  BOOST_CHECK(!multiLineString2D.dropM());
  BOOST_CHECK(!multiLineString2D.dropZ());

  MultiLineString multiLineString3D;
  multiLineString3D.addGeometry(new LineString(Point(2.0, 3.0, 5.0), Point(4.0, 5.0, 5.0)));
  multiLineString3D.addGeometry(new LineString(Point(6.0, 7.0, 5.0), Point(9.0, 10.0, 5.0)));
  BOOST_CHECK(multiLineString3D.is3D());
  BOOST_CHECK(!multiLineString3D.isMeasured());
  BOOST_CHECK(!multiLineString3D.dropM());
  BOOST_CHECK(multiLineString3D.dropZ());
  BOOST_CHECK_EQUAL(multiLineString3D.asText(0), "MULTILINESTRING ((2 3,4 5),(6 7,9 10))");
  BOOST_CHECK(!multiLineString3D.is3D());
  BOOST_CHECK(!multiLineString3D.isMeasured());
  BOOST_CHECK(!multiLineString3D.dropM());
  BOOST_CHECK(!multiLineString3D.dropZ());

  MultiLineString multiLineStringM;
  multiLineStringM.addGeometry(io::readWkt("LINESTRING M (0 0 4, 1 1 5, 2 2 6)").release());
  multiLineStringM.addGeometry(io::readWkt("LINESTRING M (3 2 4, 4 2 5)").release());
  BOOST_CHECK(!multiLineStringM.is3D());
  BOOST_CHECK(multiLineStringM.isMeasured());
  BOOST_CHECK(!multiLineStringM.dropZ());
  BOOST_CHECK(multiLineStringM.dropM());
  BOOST_CHECK_EQUAL(multiLineStringM.asText(0), "MULTILINESTRING ((0 0,1 1,2 2),(3 2,4 2))");
  BOOST_CHECK(!multiLineStringM.is3D());
  BOOST_CHECK(!multiLineStringM.isMeasured());
  BOOST_CHECK(!multiLineStringM.dropM());
  BOOST_CHECK(!multiLineStringM.dropZ());

  MultiLineString multiLineStringZM;
  multiLineStringZM.addGeometry(new LineString(Point(2.0, 3.0, 5.0, 2.0), Point(4.0, 5.0, 5.0, 2.0)));
  multiLineStringZM.addGeometry(new LineString(Point(6.0, 7.0, 5.0, 1.0), Point(9.0, 10.0, 5.0, 1.0)));
  BOOST_CHECK(multiLineStringZM.is3D());
  BOOST_CHECK(multiLineStringZM.isMeasured());
  BOOST_CHECK(multiLineStringZM.dropM());
  BOOST_CHECK(multiLineStringZM.is3D());
  BOOST_CHECK(!multiLineStringZM.isMeasured());
  BOOST_CHECK_EQUAL(multiLineStringZM.asText(0), "MULTILINESTRING Z ((2 3 5,4 5 5),(6 7 5,9 10 5))");
  BOOST_CHECK(!multiLineStringZM.dropM());
  BOOST_CHECK(multiLineStringZM.dropZ());
  BOOST_CHECK_EQUAL(multiLineStringZM.asText(0), "MULTILINESTRING ((2 3,4 5),(6 7,9 10))");
  BOOST_CHECK(!multiLineStringZM.is3D());
  BOOST_CHECK(!multiLineStringZM.isMeasured());
  BOOST_CHECK(!multiLineStringZM.dropM());
  BOOST_CHECK(!multiLineStringZM.dropZ());
}

BOOST_AUTO_TEST_SUITE_END()
