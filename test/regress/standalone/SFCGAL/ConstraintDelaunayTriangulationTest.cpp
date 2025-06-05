// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <fstream>

#include "SFCGAL/GeometryCollection.h"
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
#include "SFCGAL/algorithm/area.h"
#include "SFCGAL/detail/triangulate/ConstraintDelaunayTriangulation.h"
#include <cmath>

#include "../../../test_config.h"

#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

using namespace SFCGAL;
using namespace SFCGAL::triangulate;

BOOST_AUTO_TEST_SUITE(SFCGAL_ConstraintDelaunayTriangulationTest)

BOOST_AUTO_TEST_CASE(testTriangulateRGC)
{
  ConstraintDelaunayTriangulation triangulation;

  /*
   * read points from file
   */
  std::string filename(SFCGAL_TEST_DIRECTORY);
  filename += "/data/rgc-france-ign.xyz";
  std::ifstream ifs(filename.c_str());
  BOOST_REQUIRE(ifs.good());

  double x = NAN;
  double y = NAN;
  double z = NAN;

  while (ifs >> x >> y >> z) {
    triangulation.addVertex(Coordinate(x, y, z));
  }

  ifs.close();

  // std::string wkt = triangulation.getTriangulatedSurface()->asText(5.0) ;
  // std::cerr << "INSERT INTO draw (geometry) VALUES ( '" << "MULTIPOLYGON" <<
  // wkt.substr(3) << "'::geometry );" << std::endl;

  BOOST_CHECK_EQUAL(triangulation.numVertices(), 36566U);
  BOOST_CHECK_EQUAL(triangulation.numTriangles(), 73114U);

  std::unique_ptr<TriangulatedSurface> triangulatedSurface =
      triangulation.getTriangulatedSurface();
  BOOST_CHECK_EQUAL(triangulatedSurface->numPatches(), 73114U);
  BOOST_CHECK_CLOSE(algorithm::area(*triangulatedSurface), 818056610000.0, 0.1);
}

BOOST_AUTO_TEST_SUITE_END()
