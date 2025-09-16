// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <fstream>

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/algorithm/intersection.h"
#include "SFCGAL/detail/GetPointsVisitor.h"
#include "SFCGAL/detail/transform/AffineTransform3.h"
#include "SFCGAL/io/GeometryStreams.h"
#include "SFCGAL/io/wkt.h"

#include "../../../test_config.h"

#include <boost/ptr_container/ptr_map.hpp>
#include <boost/test/unit_test.hpp>

using namespace SFCGAL;
using namespace boost::unit_test;
// namespace SFCGAL{
// inline Geometry* new_clone( const Geometry& g )
//{
//     return g.clone();
// }
// }

void
insertOrReplace(boost::ptr_map<std::string, Geometry> &map, std::string key,
                Geometry *value)
{
  boost::ptr_map<std::string, Geometry>::iterator const found = map.find(key);

  if (found != map.end()) {
    map.erase(found);
  }

  map.insert(key, value);
}
BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_IntersectionTest)

BOOST_AUTO_TEST_CASE(testFileIntersectionTest)
{
  int const argc = framework::master_test_suite().argc;
  char    **argv = framework::master_test_suite().argv;

  // look for options
  int  test_one_line     = -1;
  bool print_line_number = false;

  for (int i = 0; i < argc; ++i) {
    std::string const argi(argv[i]);

    if (argi == "--line") {
      // only test one line
      if (argc >= i + 1) {
        sscanf(argv[i + 1], "%d", &test_one_line);
        ++i;
        continue;
      }
    }

    if (argi == "--line-number") {
      print_line_number = true;
    }
  }

  boost::ptr_map<std::string, Geometry> storedGeom;

  // logger().setLogLevel( Logger::Debug );

  std::string filename(SFCGAL_TEST_DIRECTORY);
  filename += "/data/IntersectionTest.txt";

  std::ifstream ifs(filename.c_str());
  BOOST_REQUIRE(ifs.good());

  int         numLine = 0;
  std::string line;

  while (std::getline(ifs, line)) {
    numLine++;

    if (line[0] == '#' || line.empty()) {
      continue;
    }

    if (print_line_number) {
      std::cout << numLine << '\n';
    }

    if (test_one_line == -1) {
      BOOST_TEST_MESSAGE(boost::format("line#%s:%s") % numLine % line);
    }

    std::istringstream iss(line);

    std::string dimension;
    std::string wktGA;
    std::string wktGB;
    std::string wktOut;

    std::getline(iss, dimension, '|');

    if (dimension == "S") {
      // store the geometry
      std::string geomName;
      std::getline(iss, geomName, '|');
      std::getline(iss, wktGA, '|');
      insertOrReplace(storedGeom, geomName, io::readWkt(wktGA).release());
      continue;
    }

    // we need to read the lines until here to store reference
    if (test_one_line != -1) {
      if (numLine != test_one_line) {
        continue;
      }
      BOOST_TEST_MESSAGE(boost::format("line#%s:%s") % numLine % line);
    }

    std::unique_ptr<Geometry> gA;
    std::unique_ptr<Geometry> gB;
    std::unique_ptr<Geometry> gOut;

    std::getline(iss, wktGA, '|');

    if (wktGA[0] == '@') {
      // stored geometry reference
      std::string const name = wktGA.substr(1);
      const boost::ptr_map<std::string, Geometry>::const_iterator found =
          storedGeom.find(name);

      if (found == storedGeom.end()) {
        BOOST_REQUIRE_MESSAGE(
            false, numLine << ": can't find the geometry named " << name);
      }

      gA = found->second->clone();
    } else {
      gA = io::readWkt(wktGA);
    }

    insertOrReplace(storedGeom, "A", gA->clone().release());

    std::getline(iss, wktGB, '|');

    if (wktGB[0] == '@') {
      // stored geometry reference
      std::string const name = wktGB.substr(1);
      const boost::ptr_map<std::string, Geometry>::const_iterator found =
          storedGeom.find(name);

      if (found == storedGeom.end()) {
        BOOST_REQUIRE_MESSAGE(
            false, numLine << ": can't find the geometry named " << name);
      }

      gB = found->second->clone();
    } else {
      gB = io::readWkt(wktGB);
    }

    insertOrReplace(storedGeom, "B", gB->clone().release());

    // exception management
    bool expectException      = false;
    bool expectNotImplemented = false;

    std::getline(iss, wktOut, '|');

    if (wktOut[0] == '@') {
      // stored geometry reference
      std::string const name = wktOut.substr(1);
      const boost::ptr_map<std::string, Geometry>::const_iterator found =
          storedGeom.find(name);

      if (found == storedGeom.end()) {
        BOOST_CHECK_MESSAGE(false, numLine << ": can't find the geometry named "
                                           << name);
      }

      gOut = found->second->clone();
    }
    // expect an exception
    else if (wktOut[0] == '!') {
      if (wktOut == "!NotImplemented") {
        expectNotImplemented = true;
      }

      expectException = true;
    } else {
      gOut = io::readWkt(wktOut);
    }

    try {
      if (dimension == "2") {
        std::unique_ptr<Geometry> result = algorithm::intersection(*gA, *gB);
        BOOST_CHECK_MESSAGE(*result == *gOut,
                            numLine << ": intersection(" << gA->asText() << ", "
                                    << gB->asText() << ") is "
                                    << result->asText() << " and should be "
                                    << gOut->asText());
      } else if (dimension == "3") {
        std::unique_ptr<Geometry> result = algorithm::intersection3D(*gA, *gB);
        BOOST_CHECK_MESSAGE(
            result.get()->almostEqual(
                *gOut.get(), 0.0,
                algorithm::EqualityStrictness::coverSubGeomNonOrdered()),
            numLine << ": intersection3D(" << gA->asText() << ", "
                    << gB->asText() << ") is " << result->asText()
                    << " and should be " << gOut->asText());
      } else {
        BOOST_CHECK(false);
      }
    } catch (SFCGAL::NotImplementedException &e) {
      BOOST_CHECK_MESSAGE(expectNotImplemented, numLine << ": " << e.what());
    } catch (std::exception &e) {
      BOOST_CHECK_MESSAGE(expectException, numLine << ": " << e.what());
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
