// scale_tests.cpp
#include <SFCGAL/LineString.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Triangle.h>
#include <SFCGAL/algorithm/scale.h>
#include <SFCGAL/io/wkt.h>
#include <boost/test/unit_test.hpp>

using namespace SFCGAL;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_ScaleTest)

// Helper function to check if two geometries are equal, ignoring formatting
// differences
auto
geometriesEqual(const std::string &g1, const std::string &g2) -> bool
{
  // Remove all spaces and convert to lowercase for comparison
  auto normalize = [](std::string s) {
    s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
  };
  return normalize(g1) == normalize(g2);
}

BOOST_AUTO_TEST_CASE(testScaleUniform2D)
{
  std::unique_ptr<Geometry> g(io::readWkt("POINT (1 2)"));
  algorithm::scale(*g, 2.0);
  BOOST_CHECK(geometriesEqual(g->asText(1), "POINT (2.0 4.0)"));
}

BOOST_AUTO_TEST_CASE(testScaleUniform3D)
{
  std::unique_ptr<Geometry> g(io::readWkt("POINT Z (1 2 3)"));
  algorithm::scale(*g, 2.0);
  BOOST_CHECK(geometriesEqual(g->asText(1), "POINT Z (2.0 4.0 6.0)"));
}

BOOST_AUTO_TEST_CASE(testScaleUniformZM)
{
  std::unique_ptr<Geometry> g(io::readWkt("POINT ZM (1 2 3 4)"));
  algorithm::scale(*g, 2.0);
  BOOST_CHECK(geometriesEqual(g->asText(1), "POINT ZM (2.0 4.0 6.0 4.0)"));
}

BOOST_AUTO_TEST_CASE(testScaleNonUniform2DDefaultZ)
{
  std::unique_ptr<Geometry> g(io::readWkt("POINT (1 2)"));
  algorithm::scale(*g, 2.0, 3.0);
  BOOST_CHECK(geometriesEqual(g->asText(1), "POINT (2.0 6.0)"));
}

BOOST_AUTO_TEST_CASE(testScaleNonUniform2D)
{
  std::unique_ptr<Geometry> g(io::readWkt("POINT (1 2)"));
  algorithm::scale(*g, 2.0, 3.0, 1.0);
  BOOST_CHECK(geometriesEqual(g->asText(1), "POINT (2.0 6.0)"));
}

BOOST_AUTO_TEST_CASE(testScaleNonUniform3D)
{
  std::unique_ptr<Geometry> g(io::readWkt("POINT Z (1 2 3)"));
  algorithm::scale(*g, 2.0, 3.0, 4.0);
  BOOST_CHECK(geometriesEqual(g->asText(1), "POINT Z (2.0 6.0 12.0)"));
}

BOOST_AUTO_TEST_CASE(testScaleNonUniformZM)
{
  std::unique_ptr<Geometry> g(io::readWkt("POINT ZM (1 2 3 4)"));
  algorithm::scale(*g, 2.0, 3.0, 4.0);
  BOOST_CHECK(geometriesEqual(g->asText(1), "POINT ZM (2.0 6.0 12.0 4.0)"));
}

BOOST_AUTO_TEST_CASE(testScaleAroundCenter2D)
{
  std::unique_ptr<Geometry> g(io::readWkt("POINT (3 4)"));
  algorithm::scale(*g, 2.0, 2.0, 1.0, 1, 1, 0);
  BOOST_CHECK(geometriesEqual(g->asText(1), "POINT (5.0 7.0)"));
}

BOOST_AUTO_TEST_CASE(testScaleAroundCenter3D)
{
  std::unique_ptr<Geometry> g(io::readWkt("POINT Z (3 4 5)"));
  algorithm::scale(*g, 2.0, 2.0, 2.0, 1, 1, 1);
  BOOST_CHECK(geometriesEqual(g->asText(1), "POINT Z (5.0 7.0 9.0)"));
}

BOOST_AUTO_TEST_CASE(testScaleAroundCenterZM)
{
  std::unique_ptr<Geometry> g(io::readWkt("POINT ZM (3 4 5 6)"));
  algorithm::scale(*g, 2.0, 2.0, 2.0, 1, 1, 1);
  BOOST_CHECK(geometriesEqual(g->asText(1), "POINT ZM (5.0 7.0 9.0 6.0)"));
}

BOOST_AUTO_TEST_CASE(testScaleLineString2D)
{
  std::unique_ptr<Geometry> g(io::readWkt("LINESTRING (0 0, 1 1, 2 0)"));
  algorithm::scale(*g, 2.0);
  BOOST_CHECK(
      geometriesEqual(g->asText(1), "LINESTRING (0.0 0.0,2.0 2.0,4.0 0.0)"));
}

BOOST_AUTO_TEST_CASE(testScaleLineString3D)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("LINESTRING Z (0 0 0, 1 1 1, 2 0 2)"));
  algorithm::scale(*g, 2.0);
  BOOST_CHECK(geometriesEqual(
      g->asText(1), "LINESTRING Z (0.0 0.0 0.0,2.0 2.0 2.0,4.0 0.0 4.0)"));
}

BOOST_AUTO_TEST_CASE(testScalePolygon2D)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))"));
  algorithm::scale(*g, 2.0, 3.0, 1.0);
  BOOST_CHECK(geometriesEqual(
      g->asText(1), "POLYGON ((0.0 0.0,2.0 0.0,2.0 3.0,0.0 3.0,0.0 0.0))"));
}

BOOST_AUTO_TEST_CASE(testScalePolygon2DDefaultZ)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))"));
  algorithm::scale(*g, 2.0, 3.0);
  BOOST_CHECK(geometriesEqual(
      g->asText(1), "POLYGON ((0.0 0.0,2.0 0.0,2.0 3.0,0.0 3.0,0.0 0.0))"));
}

BOOST_AUTO_TEST_CASE(testScalePolygon3D)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("POLYGON Z ((0 0 0, 1 0 0, 1 1 1, 0 1 1, 0 0 0))"));
  algorithm::scale(*g, 2.0, 3.0, 4.0);
  BOOST_CHECK(geometriesEqual(g->asText(1),
                              "POLYGON Z ((0.0 0.0 0.0,2.0 0.0 0.0,2.0 "
                              "3.0 4.0,0.0 3.0 4.0,0.0 0.0 0.0))"));
}

BOOST_AUTO_TEST_CASE(testScaleTriangle2D)
{
  std::unique_ptr<Geometry> g(io::readWkt("TRIANGLE ((0 0, 1 0, 0 1, 0 0))"));
  algorithm::scale(*g, 2.0);
  BOOST_CHECK(geometriesEqual(g->asText(1),
                              "TRIANGLE ((0.0 0.0,2.0 0.0,0.0 2.0,0.0 0.0))"));
}

BOOST_AUTO_TEST_CASE(testScaleTriangle3D)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("TRIANGLE Z ((0 0 0, 1 0 0, 0 1 0, 0 0 0))"));
  algorithm::scale(*g, 2.0);
  BOOST_CHECK(geometriesEqual(
      g->asText(1),
      "TRIANGLE Z ((0.0 0.0 0.0,2.0 0.0 0.0,0.0 2.0 0.0,0.0 0.0 0.0))"));
}

BOOST_AUTO_TEST_CASE(testScalePolyhedralSurface)
{
  std::string wkt =
      "POLYHEDRALSURFACE Z (((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)), "
      "((0 0 0, 0 1 0, 0 1 1, 0 0 1, 0 0 0)), "
      "((0 0 0, 1 0 0, 1 0 1, 0 0 1, 0 0 0)), "
      "((1 1 1, 1 0 1, 0 0 1, 0 1 1, 1 1 1)), "
      "((1 1 1, 1 0 1, 1 0 0, 1 1 0, 1 1 1)), "
      "((1 1 1, 0 1 1, 0 1 0, 1 1 0, 1 1 1)))";
  std::unique_ptr<Geometry> g(io::readWkt(wkt));
  algorithm::scale(*g, 2.0, 2.0, 2.0);

  auto *surface = dynamic_cast<PolyhedralSurface *>(g.get());
  BOOST_CHECK(surface != nullptr);
  BOOST_CHECK(
      geometriesEqual(surface->polygonN(0).exteriorRing().pointN(2).asText(1),
                      "POINT Z (2.0 2.0 0.0)"));
  BOOST_CHECK(
      geometriesEqual(surface->polygonN(3).exteriorRing().pointN(0).asText(1),
                      "POINT Z (2.0 2.0 2.0)"));
}

BOOST_AUTO_TEST_CASE(testScaleZero)
{
  std::unique_ptr<Geometry> g(io::readWkt("POINT (1 2 3)"));
  algorithm::scale(*g, 0.0);
  BOOST_CHECK(geometriesEqual(g->asText(1), "POINT Z (0.0 0.0 0.0)"));
}

BOOST_AUTO_TEST_CASE(testScaleNegative)
{
  std::unique_ptr<Geometry> g(io::readWkt("POINT (1 2 3)"));
  algorithm::scale(*g, -1.0);
  BOOST_CHECK(geometriesEqual(g->asText(1), "POINT Z (-1.0 -2.0 -3.0)"));
}

BOOST_AUTO_TEST_CASE(testPreserveDimensionality2D)
{
  std::unique_ptr<Geometry> g(io::readWkt("POINT (1 2)"));
  algorithm::scale(*g, 2.0, 2.0, 2.0);
  BOOST_CHECK(geometriesEqual(g->asText(1), "POINT (2.0 4.0)"));
}

BOOST_AUTO_TEST_CASE(testPreserveDimensionality3D)
{
  std::unique_ptr<Geometry> g(io::readWkt("POINT Z (1 2 3)"));
  algorithm::scale(*g, 2.0, 2.0, 2.0);
  BOOST_CHECK(geometriesEqual(g->asText(1), "POINT Z (2.0 4.0 6.0)"));
}

BOOST_AUTO_TEST_CASE(testPreserveDimensionality3DM)
{
  std::unique_ptr<Geometry> g(io::readWkt("POINT M (1 2 3)"));
  algorithm::scale(*g, 2.0, 2.0, 2.0);
  BOOST_CHECK(geometriesEqual(g->asText(1), "POINT M (2.0 4.0 3.0)"));
}

BOOST_AUTO_TEST_CASE(testPreserveDimensionalityZM)
{
  std::unique_ptr<Geometry> g(io::readWkt("POINT ZM (1 2 3 4)"));
  algorithm::scale(*g, 2.0, 2.0, 2.0);
  BOOST_CHECK(geometriesEqual(g->asText(1), "POINT ZM (2.0 4.0 6.0 4.0)"));
}

BOOST_AUTO_TEST_SUITE_END()
