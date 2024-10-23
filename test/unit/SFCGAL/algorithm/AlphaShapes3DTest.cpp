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
#include <boost/format/parsing.hpp>
#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/algorithm/alphaShapes3D.h"
#include "SFCGAL/io/wkt.h"

using namespace SFCGAL;

#include "../../../test_config.h"
// always after CGAL
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_AlphaShapesTest)

// algorithm::alphaShapes3D


BOOST_AUTO_TEST_CASE(testAlphaShapes3D_Empty)
{
  GeometryCollection emptyCollection;
  emptyCollection.addGeometry(Polygon());
  emptyCollection.addGeometry(Polygon());
  std::unique_ptr<Geometry> emptyAlphaShape3D (algorithm::alphaShapes3D(emptyCollection));
  BOOST_CHECK(emptyAlphaShape3D->isEmpty());
}

BOOST_AUTO_TEST_CASE(testAlphaShapes3D_MultiPoint)
{
  std::string inputData(SFCGAL_TEST_DIRECTORY);
  inputData += "/data/bunny1000Wkt.txt";
  std::ifstream bunnyFS(inputData.c_str());
  BOOST_REQUIRE(bunnyFS.good());

  std::ostringstream inputWkt;
  inputWkt << bunnyFS.rdbuf();

  std::unique_ptr<Geometry> geomInput(io::readWkt(inputWkt.str()));

  std::unique_ptr<Geometry> alphaShapesGeneral(algorithm::alphaShapes3D(geomInput->as<const SFCGAL::Geometry>()));

  // check
  try {
    const SFCGAL::Validity validity = SFCGAL::algorithm::isValid(*alphaShapesGeneral);
    bool is_valid = validity;
    if (!is_valid) {
      std::cout << validity.reason() << "\n";
      exit(1);
    }
  } catch (SFCGAL::Exception &e) {
    std::cout << e.what() << "\n";
    exit(1);
  }

  exit(1);
}


BOOST_AUTO_TEST_SUITE_END()
