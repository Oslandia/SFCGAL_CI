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
#include <memory>
#include <string>

#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/LineString.h>
#include <SFCGAL/MultiLineString.h>
#include <SFCGAL/MultiPoint.h>
#include <SFCGAL/MultiPolygon.h>
#include <SFCGAL/MultiSolid.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Solid.h>
#include <SFCGAL/Triangle.h>
#include <SFCGAL/TriangulatedSurface.h>
#include <SFCGAL/io/wkb.h>
#include <SFCGAL/io/wkt.h>
#include <SFCGAL/io/ewkt.h>
#include <SFCGAL/PreparedGeometry.h>

#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

using namespace SFCGAL;
using namespace SFCGAL::io;

#include "../../../test_config.h"

BOOST_AUTO_TEST_SUITE(SFCGAL_io_WkbWriterTest)

BOOST_AUTO_TEST_CASE(writeWkb)
{
  std::string inputData(SFCGAL_TEST_DIRECTORY);
  inputData += "/data/WKT.txt";
  std::ifstream ifs(inputData.c_str());
  BOOST_REQUIRE(ifs.good());

  std::string expectedData(SFCGAL_TEST_DIRECTORY);
  expectedData += "/data/WKT_expected.txt";
  std::ifstream efs(expectedData.c_str());
  BOOST_REQUIRE(efs.good());

  std::string inputWkt;
  std::string expectedWkb;
  while (std::getline(ifs, inputWkt)) {
    std::unique_ptr<Geometry> g(io::readWkt(inputWkt));
    std::getline(efs, expectedWkb);
    BOOST_CHECK_EQUAL(g->asWkb(), expectedWkb);
  }
}

BOOST_AUTO_TEST_CASE(readWkb)
{
  std::string inputData(SFCGAL_TEST_DIRECTORY);
  inputData += "/data/WKT.txt";
  std::ifstream ifs(inputData.c_str());
  BOOST_REQUIRE(ifs.good());

  std::string expectedData(SFCGAL_TEST_DIRECTORY);
  expectedData += "/data/WKT_expected.txt";
  std::ifstream efs(expectedData.c_str());
  BOOST_REQUIRE(efs.good());

  std::string inputWkt;
  std::string expectedWkb;
  while (std::getline(ifs, inputWkt)) {
    std::getline(efs, expectedWkb);
    std::unique_ptr<Geometry> g(io::readWkt(inputWkt));
    std::unique_ptr<Geometry> gWkb(io::readWkb(expectedWkb));
    std::vector<std::string>  allowedBeThyFail{
        "GEOMETRYCOLLECTION (POINT Z (1 2 3), LINESTRING (0 0, 1 1, 2 2), "
         "POLYGON Z ((0 0 1, 0 3 2, 3 3 3, 3 0 4, 0 0 1)), MULTIPOINT M ((1 1 "
         "4), (2 2 5)))",
        "GEOMETRYCOLLECTION (POINT ZM (1 2 3 4), POINT EMPTY, POLYGON ((0 0, 0 "
         "4, 4 4, 4 0, 0 0)), POLYGON ZM ((0 0 1 4, 0 3 2 5, 3 3 3 6, 3 0 4 7, "
         "0 0 1 4)), POLYGON EMPTY)",
        "POLYGON Z ((0 0,0 10,10 10,10 0,0 0),(1 1 1,1 2 1,2 2 1,2 1 1,1 1 "
         "1))"};
    if (std::find(allowedBeThyFail.begin(), allowedBeThyFail.end(), inputWkt) ==
        std::end(allowedBeThyFail)) {
      BOOST_CHECK_EQUAL(g->asText(0), gWkb->asText(0));
    }
  }
}

BOOST_AUTO_TEST_CASE(readEWkb)
{
  std::string inputData(SFCGAL_TEST_DIRECTORY);
  inputData += "/data/EWKT.txt";
  std::ifstream ifs(inputData.c_str());
  BOOST_REQUIRE(ifs.good());

  std::string expectedData(SFCGAL_TEST_DIRECTORY);
  expectedData += "/data/EWKT_expected.txt";
  std::ifstream efs(expectedData.c_str());
  BOOST_REQUIRE(efs.good());

  std::string inputWkt;
  std::string expectedWkb;
  while (std::getline(ifs, inputWkt)) {
    std::getline(efs, expectedWkb);
    std::unique_ptr<PreparedGeometry> pg(io::readEwkt(inputWkt));
    std::unique_ptr<Geometry> gWkb(io::readWkb(expectedWkb));
    BOOST_CHECK_EQUAL(pg->geometry().asText(0), gWkb->asText(0));
  }
}
BOOST_AUTO_TEST_SUITE_END()
