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
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/io/wkt.h"

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_MultiPointTest)

BOOST_AUTO_TEST_CASE(defaultConstructor)
{
  MultiPoint const g;
  BOOST_CHECK(g.isEmpty());
  BOOST_CHECK(!g.is3D());
  BOOST_CHECK_EQUAL(g.numGeometries(), 0U);
}

BOOST_AUTO_TEST_CASE(testGeometryTypeId)
{
  MultiPoint const g;
  BOOST_CHECK_EQUAL(g.geometryTypeId(), TYPE_MULTIPOINT);
}

//-- addAllowedGeometry
BOOST_AUTO_TEST_CASE(addPoint)
{
  MultiPoint g;
  g.addGeometry(new Point(2.0, 3.0));
  BOOST_CHECK_EQUAL(g.numGeometries(), 1U);
  g.addGeometry(new Point(4.0, 5.0));
  BOOST_CHECK_EQUAL(g.numGeometries(), 2U);
}
//-- addForbidenGeometry
BOOST_AUTO_TEST_CASE(addLineStringThrow)
{
  MultiPoint g;
  BOOST_CHECK_THROW(g.addGeometry(LineString()), std::exception);
}

//-- asText

BOOST_AUTO_TEST_CASE(asTextEmpty)
{
  MultiPoint const g;
  BOOST_CHECK_EQUAL(g.asText(1), "MULTIPOINT EMPTY");
}

BOOST_AUTO_TEST_CASE(asText2d)
{
  MultiPoint g;
  g.addGeometry(Point(2.0, 3.0));
  g.addGeometry(Point(3.0, 4.0));
  BOOST_CHECK_EQUAL(g.asText(3), "MULTIPOINT ((2.000 3.000),(3.000 4.000))");
}

BOOST_AUTO_TEST_CASE(dropZM)
{
  MultiPoint multiPointEmpty;
  BOOST_CHECK(multiPointEmpty.isEmpty());
  BOOST_CHECK(!multiPointEmpty.is3D());
  BOOST_CHECK(!multiPointEmpty.isMeasured());
  BOOST_CHECK(!multiPointEmpty.dropZ());
  BOOST_CHECK(!multiPointEmpty.dropM());

  MultiPoint multiPoint2D;
  multiPoint2D.addGeometry(new Point(2.0, 3.0));
  multiPoint2D.addGeometry(new Point(4.0, 5.0));
  BOOST_CHECK(!multiPoint2D.is3D());
  BOOST_CHECK(!multiPoint2D.isMeasured());
  BOOST_CHECK(!multiPoint2D.dropM());
  BOOST_CHECK(!multiPoint2D.dropZ());

  MultiPoint multiPoint3D;
  multiPoint3D.addGeometry(new Point(2.0, 3.0, 5.0));
  multiPoint3D.addGeometry(new Point(4.0, 5.0, 7.0));
  BOOST_CHECK(multiPoint3D.is3D());
  BOOST_CHECK(!multiPoint3D.isMeasured());
  BOOST_CHECK(!multiPoint3D.dropM());
  BOOST_CHECK(multiPoint3D.dropZ());
  BOOST_CHECK_EQUAL(multiPoint3D.asText(1), "MULTIPOINT ((2.0 3.0),(4.0 5.0))");
  BOOST_CHECK(!multiPoint3D.is3D());
  BOOST_CHECK(!multiPoint3D.isMeasured());
  BOOST_CHECK(!multiPoint3D.dropM());
  BOOST_CHECK(!multiPoint3D.dropZ());

  MultiPoint multiPointM;
  multiPointM.addGeometry(io::readWkt("POINT M (2 3 4)").release());
  multiPointM.addGeometry(io::readWkt("POINT M (4 5 7)").release());
  BOOST_CHECK(!multiPointM.is3D());
  BOOST_CHECK(multiPointM.isMeasured());
  BOOST_CHECK(!multiPointM.dropZ());
  BOOST_CHECK(multiPointM.dropM());
  BOOST_CHECK_EQUAL(multiPointM.asText(0), "MULTIPOINT ((2 3),(4 5))");
  BOOST_CHECK(!multiPointM.is3D());
  BOOST_CHECK(!multiPointM.isMeasured());
  BOOST_CHECK(!multiPointM.dropM());
  BOOST_CHECK(!multiPointM.dropZ());

  MultiPoint multiPointZM;
  multiPointZM.addGeometry(new Point(2.0, 3.0, 5.0, 6.0));
  multiPointZM.addGeometry(new Point(4.0, 5.0, 7.0, 8.0));
  BOOST_CHECK(multiPointZM.is3D());
  BOOST_CHECK(multiPointZM.isMeasured());
  BOOST_CHECK(multiPointZM.dropM());
  BOOST_CHECK(multiPointZM.is3D());
  BOOST_CHECK(!multiPointZM.isMeasured());
  BOOST_CHECK_EQUAL(multiPointZM.asText(0), "MULTIPOINT Z ((2 3 5),(4 5 7))");
  BOOST_CHECK(!multiPointZM.dropM());
  BOOST_CHECK(multiPointZM.dropZ());
  BOOST_CHECK_EQUAL(multiPointZM.asText(0), "MULTIPOINT ((2 3),(4 5))");
  BOOST_CHECK(!multiPointZM.is3D());
  BOOST_CHECK(!multiPointZM.isMeasured());
  BOOST_CHECK(!multiPointZM.dropM());
  BOOST_CHECK(!multiPointZM.dropZ());
}

BOOST_AUTO_TEST_CASE(swapXY)
{
  MultiPoint multiPointEmpty;
  BOOST_CHECK(multiPointEmpty.isEmpty());
  multiPointEmpty.swapXY();
  BOOST_CHECK(multiPointEmpty.isEmpty());

  MultiPoint multiPoint2D;
  multiPoint2D.addGeometry(new Point(2.0, 3.0));
  multiPoint2D.addGeometry(new Point(4.0, 5.0));
  multiPoint2D.swapXY();
  BOOST_CHECK_EQUAL(multiPoint2D.asText(0), "MULTIPOINT ((3 2),(5 4))");

  MultiPoint multiPoint3D;
  multiPoint3D.addGeometry(new Point(9.0, 3.0, 5.0));
  multiPoint3D.addGeometry(new Point(12.0, 5.0, 7.0));
  multiPoint3D.swapXY();
  BOOST_CHECK_EQUAL(multiPoint3D.asText(0), "MULTIPOINT Z ((3 9 5),(5 12 7))");

  MultiPoint multiPointM;
  multiPointM.addGeometry(io::readWkt("POINT M (20 7 4)").release());
  multiPointM.addGeometry(io::readWkt("POINT M (14 9 7)").release());
  multiPointM.swapXY();
  BOOST_CHECK_EQUAL(multiPointM.asText(0), "MULTIPOINT M ((7 20 4),(9 14 7))");

  MultiPoint multiPointZM;
  multiPointZM.addGeometry(new Point(-2.0, -3.0, 5.0, 6.0));
  multiPointZM.addGeometry(new Point(42.0, -5.0, 7.0, 8.0));
  multiPointZM.swapXY();
  BOOST_CHECK_EQUAL(multiPointZM.asText(0),
                    "MULTIPOINT ZM ((-3 -2 5 6),(-5 42 7 8))");
}
//-- is< T >

BOOST_AUTO_TEST_CASE(isGeometryCollection)
{
  MultiPoint const g;
  BOOST_CHECK(g.is<GeometryCollection>());
}

BOOST_AUTO_TEST_CASE(isMultiPoint)
{
  MultiPoint const g;
  BOOST_CHECK(g.is<MultiPoint>());
}

BOOST_AUTO_TEST_CASE(getCoordinateType)
{
  BOOST_CHECK_EQUAL(
      io::readWkt("MULTIPOINT((0 0), (1 1))")->getCoordinateType(),
      CoordinateType::COORDINATE_XY);
  BOOST_CHECK_EQUAL(
      io::readWkt("MULTIPOINT Z((0 0 1), (1 1 1))")->getCoordinateType(),
      CoordinateType::COORDINATE_XYZ);
  BOOST_CHECK_EQUAL(
      io::readWkt("MULTIPOINT M((0 0 2), (1 1 2))")->getCoordinateType(),
      CoordinateType::COORDINATE_XYM);
  BOOST_CHECK_EQUAL(
      io::readWkt("MULTIPOINT ZM((0 0 1 2), (1 1 1 2))")->getCoordinateType(),
      CoordinateType::COORDINATE_XYZM);
}
BOOST_AUTO_TEST_SUITE_END()
