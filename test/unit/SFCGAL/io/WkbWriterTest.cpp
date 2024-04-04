// Copyright (c) 2023-2023, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <memory>
#include <string>

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/PreparedGeometry.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/io/ewkt.h"
#include "SFCGAL/io/wkb.h"
#include "SFCGAL/io/wkt.h"

#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

using namespace SFCGAL;
using namespace SFCGAL::io;

#include "../../../test_config.h"

BOOST_AUTO_TEST_SUITE(SFCGAL_io_WkbWriterTest)

const std::vector<std::string> allowedBeThyFail{
    "GEOMETRYCOLLECTION (POINT Z (1 2 3), LINESTRING (0 0, 1 1, 2 2), "
    "POLYGON Z ((0 0 1, 0 3 2, 3 3 3, 3 0 4, 0 0 1)), MULTIPOINT M ((1 1 "
    "4), (2 2 5)))",
    "GEOMETRYCOLLECTION (POINT ZM (1 2 3 4), POINT EMPTY, POLYGON ((0 0, 0 "
    "4, 4 4, 4 0, 0 0)), POLYGON ZM ((0 0 1 4, 0 3 2 5, 3 3 3 6, 3 0 4 7, "
    "0 0 1 4)), POLYGON EMPTY)",
    "POLYGON Z ((0 0,0 10,10 10,10 0,0 0),(1 1 1,1 2 1,2 2 1,2 1 1,1 1 "
    "1))"};
BOOST_AUTO_TEST_CASE(writeWkb)
{
  std::string inputData(SFCGAL_TEST_DIRECTORY);
  inputData += "/data/WKT.txt";
  std::ifstream ifs(inputData.c_str());
  BOOST_REQUIRE(ifs.good());

  std::string expectedData(SFCGAL_TEST_DIRECTORY);
  expectedData += "/data/WKB_expected.txt";
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
  expectedData += "/data/WKB_expected.txt";
  std::ifstream efs(expectedData.c_str());
  BOOST_REQUIRE(efs.good());

  std::string inputWkt;
  std::string expectedWkb;
  while (std::getline(ifs, inputWkt)) {
    std::getline(efs, expectedWkb);
    std::unique_ptr<Geometry> g(io::readWkt(inputWkt));
    std::unique_ptr<Geometry> gWkb(io::readWkb(expectedWkb));
    if (std::find(allowedBeThyFail.begin(), allowedBeThyFail.end(), inputWkt) ==
        std::end(allowedBeThyFail)) {
      BOOST_CHECK_EQUAL(g->asText(0), gWkb->asText(0));
    }
  }
}

BOOST_AUTO_TEST_CASE(PostgisEWkb)
{
  std::string inputData(SFCGAL_TEST_DIRECTORY);
  inputData += "/data/EWKB_postgis.txt";
  std::ifstream ifs(inputData.c_str());
  BOOST_REQUIRE(ifs.good());

  std::string expectedData(SFCGAL_TEST_DIRECTORY);
  expectedData += "/data/WKT.txt";
  std::ifstream efs(expectedData.c_str());
  BOOST_REQUIRE(efs.good());

  std::string inputWkb;
  std::string expectedWkt;
  while (std::getline(ifs, inputWkb)) {
    std::getline(efs, expectedWkt);
    std::unique_ptr<Geometry>         gWkt(io::readWkt(expectedWkt));
    std::string const                 ewkt = "SRID=3946;" + expectedWkt;
    std::unique_ptr<PreparedGeometry> gEwkt(io::readEwkt(ewkt));
    if (!(expectedWkt.find("EMPTY") != std::string::npos) &&
        !inputWkb.empty()) {
      BOOST_CHECK_EQUAL(gEwkt->asEWKB(), inputWkb);
    }

    std::vector allowedBeThyFailFull = allowedBeThyFail;
    allowedBeThyFailFull.emplace_back(
        "GEOMETRYCOLLECTION (LINESTRING (0 0, 1 1), POLYGON EMPTY)");

    allowedBeThyFailFull.emplace_back(
        "GEOMETRYCOLLECTION (POINT (1 2), POINT EMPTY)");
    allowedBeThyFailFull.emplace_back("MULTIPOINT((1 1), EMPTY)");
    allowedBeThyFailFull.emplace_back("MULTIPOINT(EMPTY,  EMPTY)");
    allowedBeThyFailFull.emplace_back("MULTIPOINT(EMPTY, (1 1))");
    if (std::find(allowedBeThyFailFull.begin(), allowedBeThyFailFull.end(),
                  expectedWkt) == std::end(allowedBeThyFailFull)) {
      if (!inputWkb.empty()) {
        std::unique_ptr<PreparedGeometry> gEwkbFile(io::readEwkb(inputWkb));
        BOOST_CHECK_EQUAL(gEwkbFile->geometry().asText(0), gWkt->asText(0));
        BOOST_CHECK_EQUAL(3946, gEwkbFile->SRID());
      }
    }
  }
}
BOOST_AUTO_TEST_SUITE_END()
