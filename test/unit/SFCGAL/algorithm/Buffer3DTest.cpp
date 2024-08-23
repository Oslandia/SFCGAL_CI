#include <boost/test/unit_test.hpp>

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
#include "SFCGAL/algorithm/buffer3D.h"
#include "SFCGAL/io/wkt.h"

using namespace SFCGAL;

#include "../../../test_config.h"
// always after CGAL
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_Buffer3DTest)

// algorithm::buffer3D

BOOST_AUTO_TEST_CASE(testBuffer3D_Point)
{
  double radius   = 10.0;
  int    segments = 16;

  Point               point(0, 0, 0);
  algorithm::Buffer3D buffer3d(point, radius, segments);

  std::vector<algorithm::Buffer3D::BufferType> bufferTypes = {
      algorithm::Buffer3D::ROUND, algorithm::Buffer3D::CYLSPHERE,
      algorithm::Buffer3D::FLAT};

  for (auto bufferType : bufferTypes) {
    PolyhedralSurface buffer = buffer3d.compute(bufferType);
    BOOST_CHECK(buffer.is3D());
    BOOST_CHECK(buffer.numGeometries() > 0);

    std::string typeName;
    switch (bufferType) {
    case algorithm::Buffer3D::ROUND:
      typeName = "ROUND";
      break;
    case algorithm::Buffer3D::CYLSPHERE:
      typeName = "CYLSPHERE";
      break;
    case algorithm::Buffer3D::FLAT:
      typeName = "FLAT";
      break;
    }

    std::cout << "Point " << typeName << " Buffer 3D WKT: " << buffer.asText(1)
              << std::endl;
    SFCGAL::io::OBJ::save(buffer, "/tmp/point_" + typeName + "_buffer_3d.obj");
  }
}

BOOST_AUTO_TEST_CASE(testBuffer3D_LineString)
{
  double radius   = 10.0;
  int    segments = 16;

  std::vector<Point> points = {Point(-100, 0, 0), Point(40, -70, 0),
                               Point(40, 50, 40), Point(-90, -60, 60),
                               Point(0, 0, -100), Point(30, 0, 150)};
  LineString         lineString(points);

  algorithm::Buffer3D buffer3d(lineString, radius, segments);

  std::vector<algorithm::Buffer3D::BufferType> bufferTypes = {
      algorithm::Buffer3D::ROUND, algorithm::Buffer3D::CYLSPHERE,
      algorithm::Buffer3D::FLAT};

  for (auto bufferType : bufferTypes) {
    PolyhedralSurface buffer = buffer3d.compute(bufferType);
    BOOST_CHECK(buffer.is3D());
    BOOST_CHECK(buffer.numGeometries() > 0);

    std::string typeName;
    switch (bufferType) {
    case algorithm::Buffer3D::ROUND:
      typeName = "ROUND";
      break;
    case algorithm::Buffer3D::CYLSPHERE:
      typeName = "CYLSPHERE";
      break;
    case algorithm::Buffer3D::FLAT:
      typeName = "FLAT";
      break;
    }

    std::cout << "LineString " << typeName
              << " Buffer 3D WKT: " << buffer.asText(1) << std::endl;
    SFCGAL::io::OBJ::save(buffer,
                          "/tmp/linestring_" + typeName + "_buffer_3d.obj");
  }
}

BOOST_AUTO_TEST_CASE(testBuffer3D_InvalidGeometry)
{
  double radius   = 10.0;
  int    segments = 16;

  Triangle triangle(Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0));
  BOOST_CHECK_THROW(algorithm::Buffer3D(triangle, radius, segments),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
