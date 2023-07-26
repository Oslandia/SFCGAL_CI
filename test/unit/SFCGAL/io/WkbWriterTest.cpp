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
#include <SFCGAL/io/wkt.h>

#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

using namespace SFCGAL;
using namespace SFCGAL::io;

#include "../../../test_config.h"

BOOST_AUTO_TEST_SUITE(SFCGAL_io_WkbWriterTest)

//-- WKB POINT

BOOST_AUTO_TEST_CASE(wktFiles)
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

BOOST_AUTO_TEST_SUITE_END()
