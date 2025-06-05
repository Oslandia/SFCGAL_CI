// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include <exception>

#include "SFCGAL/Envelope.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/io/wkt.h"

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_MultiSolidTest)

BOOST_AUTO_TEST_CASE(defaultConstructor)
{
  MultiSolid const g;
  BOOST_CHECK(g.isEmpty());
  BOOST_CHECK(!g.is3D());
  BOOST_CHECK_EQUAL(g.numGeometries(), 0U);
}

BOOST_AUTO_TEST_CASE(testGeometryTypeId)
{
  MultiSolid const g;
  BOOST_CHECK_EQUAL(g.geometryTypeId(), TYPE_MULTISOLID);
}

//-- addAllowedGeometry
BOOST_AUTO_TEST_CASE(addSolid)
{
  MultiSolid g;
  g.addGeometry(new Solid());
  BOOST_CHECK_EQUAL(g.numGeometries(), 1U);
  g.addGeometry(new Solid());
  BOOST_CHECK_EQUAL(g.numGeometries(), 2U);
}
//-- addForbidenGeometry
BOOST_AUTO_TEST_CASE(addLineStringThrow)
{
  MultiSolid g;
  BOOST_CHECK_THROW(g.addGeometry(LineString()), std::exception);
}

//-- asText

BOOST_AUTO_TEST_CASE(asTextEmpty)
{
  MultiSolid const g;
  BOOST_CHECK_EQUAL(g.asText(1), "MULTISOLID EMPTY");
}
BOOST_AUTO_TEST_CASE(asText2d)
{
  MultiSolid g;
  g.addGeometry(Envelope(0.0, 1.0, 0.0, 1.0, 0.0, 1.0).toSolid().release());
  g.addGeometry(Envelope(2.0, 3.0, 4.0, 5.0, 6.0, 7.0).toSolid().release());
  BOOST_CHECK_EQUAL(
      g.asText(0),
      "MULTISOLID Z (((((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)),((0 0 1,1 0 1,1 1 1,0 "
      "1 1,0 0 1)),((0 0 0,1 0 0,1 0 1,0 0 1,0 0 0)),((1 1 0,0 1 0,0 1 1,1 1 "
      "1,1 1 0)),((1 0 0,1 1 0,1 1 1,1 0 1,1 0 0)),((0 0 0,0 0 1,0 1 1,0 1 0,0 "
      "0 0)))),((((2 4 6,2 5 6,3 5 6,3 4 6,2 4 6)),((2 4 7,3 4 7,3 5 7,2 5 7,2 "
      "4 7)),((2 4 6,3 4 6,3 4 7,2 4 7,2 4 6)),((3 5 6,2 5 6,2 5 7,3 5 7,3 5 "
      "6)),((3 4 6,3 5 6,3 5 7,3 4 7,3 4 6)),((2 4 6,2 4 7,2 5 7,2 5 6,2 4 "
      "6)))))");
}

//-- is< T >

BOOST_AUTO_TEST_CASE(isGeometryCollection)
{
  MultiSolid const g;
  BOOST_CHECK(g.is<GeometryCollection>());
}

BOOST_AUTO_TEST_CASE(isMultiSolid)
{
  MultiSolid const g;
  BOOST_CHECK(g.is<MultiSolid>());
}

BOOST_AUTO_TEST_CASE(dropZ)
{
  MultiSolid geom;
  BOOST_CHECK(geom.isEmpty());
  BOOST_CHECK(!geom.dropZ());

  geom.addGeometry(Envelope(0.0, 1.0, 0.0, 1.0, 0.0, 1.0).toSolid().release());
  geom.addGeometry(Envelope(2.0, 3.0, 4.0, 5.0, 6.0, 7.0).toSolid().release());
  BOOST_CHECK(geom.dropZ());
  BOOST_CHECK_EQUAL(geom.asText(0),
                    "MULTISOLID (((((0 0,0 1,1 1,1 0,0 0)),((0 0,1 0,1 1,0 "
                    "1,0 0)),((0 0,1 0,1 0,0 0,0 0)),((1 1,0 1,0 1,1 1,"
                    "1 1)),((1 0,1 1,1 1,1 0,1 0)),((0 0,0 0,0 1,0 1,0 "
                    "0)))),((((2 4,2 5,3 5,3 4,2 4)),((2 4,3 4,3 5,2 5,2 "
                    "4)),((2 4,3 4,3 4,2 4,2 4)),((3 5,2 5,2 5,3 5,3 5))"
                    ",((3 4,3 5,3 5,3 4,3 4)),((2 4,2 4,2 5,2 5,2 4)))))");

  BOOST_CHECK(!geom.dropZ());
}

BOOST_AUTO_TEST_CASE(swapXY)
{
  MultiSolid geom;

  geom.addGeometry(Envelope(0.0, 1.0, 0.0, 1.0, 0.0, 1.0).toSolid().release());
  geom.addGeometry(Envelope(2.0, 3.0, 4.0, 5.0, 6.0, 7.0).toSolid().release());
  geom.swapXY();
  BOOST_CHECK_EQUAL(geom.asText(0), "MULTISOLID Z ("
                                    "((((0 0 0,1 0 0,1 1 0,0 1 0,0 0 0)),"
                                    "((0 0 1,0 1 1,1 1 1,1 0 1,0 0 1)),"
                                    "((0 0 0,0 1 0,0 1 1,0 0 1,0 0 0)),"
                                    "((1 1 0,1 0 0,1 0 1,1 1 1,1 1 0)),"
                                    "((0 1 0,1 1 0,1 1 1,0 1 1,0 1 0)),"
                                    "((0 0 0,0 0 1,1 0 1,1 0 0,0 0 0)))),"
                                    "((((4 2 6,5 2 6,5 3 6,4 3 6,4 2 6)),"
                                    "((4 2 7,4 3 7,5 3 7,5 2 7,4 2 7)),"
                                    "((4 2 6,4 3 6,4 3 7,4 2 7,4 2 6)),"
                                    "((5 3 6,5 2 6,5 2 7,5 3 7,5 3 6)),"
                                    "((4 3 6,5 3 6,5 3 7,4 3 7,4 3 6)),"
                                    "((4 2 6,4 2 7,5 2 7,5 2 6,4 2 6)))))");
}

BOOST_AUTO_TEST_CASE(getCoordinateType)
{
  std::unique_ptr<Geometry> geom{
      io::readWkt("MULTISOLID ZM ((("
                  "((0 0 0 1, 1 0 0 1, 1 1 0 1, 0 1 0 1, 0 0 0 1)),"
                  "((0 0 1 1, 1 0 1 1, 1 1 1 1, 0 1 1 1, 0 0 1 1)),"
                  "((0 0 0 1, 0 1 0 1, 0 1 1 1, 0 0 1 1, 0 0 0 1)),"
                  "((1 0 0 1, 1 1 0 1, 1 1 1 1, 1 0 1 1, 1 0 0 1)),"
                  "((0 0 0 1, 1 0 0 1, 1 0 1 1, 0 0 1 1, 0 0 0 1)),"
                  "((0 1 0 1, 1 1 0 1, 1 1 1 1, 0 1 1 1, 0 1 0 1)) "
                  ")))")
          .release()};
  BOOST_CHECK_EQUAL(geom->getCoordinateType(), CoordinateType::COORDINATE_XYZM);
  std::unique_ptr<Geometry> cube{geom->clone()};
  cube->dropM();
  BOOST_CHECK_EQUAL(cube->getCoordinateType(), CoordinateType::COORDINATE_XYZ);
  geom->dropZ();
  BOOST_CHECK_EQUAL(geom->getCoordinateType(), CoordinateType::COORDINATE_XYM);
  geom->dropM();
  BOOST_CHECK_EQUAL(geom->getCoordinateType(), CoordinateType::COORDINATE_XY);
}
BOOST_AUTO_TEST_SUITE_END()
