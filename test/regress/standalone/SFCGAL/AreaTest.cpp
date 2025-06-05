// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

#include "../../../test_config.h"

#include <boost/format.hpp>
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
#include "SFCGAL/io/wkt.h"

#include "SFCGAL/Transform.h"
#include "SFCGAL/algorithm/area.h"
#include <cmath>

using namespace boost::unit_test;
using namespace SFCGAL;

class RotateCoordinate : public Transform {
public:
  void
  transform(Point &p) override
  {
    BOOST_ASSERT(!p.isEmpty());
    p = Point(p.is3D() ? p.z() : 0.0, p.x(), p.y());
  }
};

BOOST_AUTO_TEST_SUITE(SFCGAL_AreaTest)

/**
 * Triangulate polygon and make some checks
 * @todo Check inPolygon.area3D() == outPolygon.area3D();
 */
BOOST_AUTO_TEST_CASE(testComputeArea)
{
  std::string filename(SFCGAL_TEST_DIRECTORY);
  filename += "/data/AreaTest.txt";

  std::ifstream ifs(filename.c_str());
  BOOST_REQUIRE(ifs.good());

  std::string line;

  while (std::getline(ifs, line)) {
    if (line[0] == '#' || line.empty()) {
      continue;
    }

    std::istringstream iss(line);

    std::string id;
    iss >> id;

    double expectedArea = NAN;
    iss >> expectedArea;

    std::string inputWkt;
    std::getline(iss, inputWkt);

    std::unique_ptr<Geometry> g(io::readWkt(inputWkt));
    double                    area = algorithm::area3D(*g);
    BOOST_TEST_MESSAGE(boost::format("area( '%1%' ) = %2%") % inputWkt % area);

    RotateCoordinate rotateCoordinate;
    g->accept(rotateCoordinate);
    double const areaRotate = algorithm::area3D(*g);

    // check area == areaRotate
    BOOST_CHECK_CLOSE(area, areaRotate, 0.5);
    // check area == expectedArea
    BOOST_CHECK_CLOSE(area, expectedArea, 0.5);
  }
}

BOOST_AUTO_TEST_SUITE_END()
