// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/format/parsing.hpp>
#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/algorithm/alphaWrapping3D.h"
#include "SFCGAL/algorithm/covers.h"
#include "SFCGAL/io/wkt.h"

using namespace SFCGAL;

#include "../../../test_config.h"
// always after CGAL
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_AlphaWrapping3DTest)

// algorithm::alphaWrapping3D

BOOST_AUTO_TEST_CASE(testAlphaWrapping3D_Empty)
{
  GeometryCollection emptyCollection;
  emptyCollection.addGeometry(Polygon());
  emptyCollection.addGeometry(Polygon());
  std::unique_ptr<Geometry> emptyAlphaWrapping3D(
      algorithm::alphaWrapping3D(emptyCollection, 300, 5000));
  BOOST_CHECK(emptyAlphaWrapping3D->isEmpty());
}

BOOST_AUTO_TEST_CASE(testAlphaWrapping3D_MultiPoint)
{
  std::string inputData(SFCGAL_TEST_DIRECTORY);
  inputData += "/data/bunny1000Wkt.txt";
  std::ifstream bunnyFSInput(inputData.c_str());
  BOOST_REQUIRE(bunnyFSInput.good());
  std::ostringstream inputWkt;
  inputWkt << bunnyFSInput.rdbuf();

  std::unique_ptr<Geometry> inputGeom(io::readWkt(inputWkt.str()));
  BOOST_REQUIRE(inputGeom->is3D());

  std::unique_ptr<Geometry> alphaWrappingResult(
      algorithm::alphaWrapping3D(inputGeom->as<const SFCGAL::Geometry>(), 20));
#if CGAL_VERSION_MAJOR < 6
  // 2304 on Linux
  // 2306 on Mac/FreeBSD maybe other
  BOOST_CHECK_GE(alphaWrappingResult->as<PolyhedralSurface>().numPatches(),
                 2304);
#else
  BOOST_CHECK_EQUAL(alphaWrappingResult->as<PolyhedralSurface>().numPatches(),
                    2386);
#endif

  BOOST_REQUIRE(alphaWrappingResult->is3D());
}

BOOST_AUTO_TEST_SUITE_END()
