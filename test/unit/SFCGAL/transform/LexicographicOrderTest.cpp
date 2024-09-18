#include "SFCGAL/algorithm/lexicographicOrder.h"
#include "SFCGAL/io/wkt.h"
#include <boost/test/unit_test.hpp>

using namespace SFCGAL;
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_LexicographicOrderTest)

#define CHECK_WKT(g, expected_wkt) BOOST_CHECK_EQUAL(g->asText(0), expected_wkt)

BOOST_AUTO_TEST_CASE(testPoint)
{
  std::unique_ptr<Geometry> g(io::readWkt("POINT(1 2)"));
  transform::lexicographicOrder(*g, false);
  CHECK_WKT(g, "POINT(1 2)");
}

BOOST_AUTO_TEST_CASE(testPoint_reorder)
{
  std::unique_ptr<Geometry> g(io::readWkt("POINT(1 2)"));
  transform::lexicographicOrder(*g, true);
  CHECK_WKT(g, "POINT(1 2)");
}

BOOST_AUTO_TEST_CASE(testMultiPoint)
{
  std::unique_ptr<Geometry> g(io::readWkt("MULTIPOINT((1 1),(0 0))"));
  transform::lexicographicOrder(*g, false);
  CHECK_WKT(g, "MULTIPOINT((1 1),(0 0))");
}

BOOST_AUTO_TEST_CASE(testMultiPoint_reorder)
{
  std::unique_ptr<Geometry> g(io::readWkt("MULTIPOINT((1 1),(0 0))"));
  transform::lexicographicOrder(*g, true);
  CHECK_WKT(g, "MULTIPOINT((0 0),(1 1))");
}

BOOST_AUTO_TEST_CASE(testLineString)
{
  std::unique_ptr<Geometry> g(io::readWkt("LINESTRING(0 0,-1 -1,-2 -2)"));
  transform::lexicographicOrder(*g, false);
  CHECK_WKT(g, "LINESTRING(0 0,-1 -1,-2 -2)");
}

BOOST_AUTO_TEST_CASE(testLineString_reorder)
{
  std::unique_ptr<Geometry> g(io::readWkt("LINESTRING(0 0,-1 -1,-2 -2)"));
  transform::lexicographicOrder(*g, true);
  CHECK_WKT(g, "LINESTRING(-2 -2,-1 -1,0 0)");
}

BOOST_AUTO_TEST_CASE(testMultiLineString)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("MULTILINESTRING((0 0,1 1),(2 2,1 0))"));
  transform::lexicographicOrder(*g, false);
  CHECK_WKT(g, "MULTILINESTRING((0 0,1 1),(2 2,1 0))");
}

BOOST_AUTO_TEST_CASE(testMultiLineString_reorder)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("MULTILINESTRING((0 0,1 1),(2 2,1 0))"));
  transform::lexicographicOrder(*g, true);
  CHECK_WKT(g, "MULTILINESTRING((0 0,1 1),(1 0,2 2))");
}

BOOST_AUTO_TEST_CASE(testPolygon)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("POLYGON((0 0,0 5,5 5,5 0,0 0),(1 1,2 1,2 2,1 2,1 1),(3 3,4 "
                  "3,4 4,3 4,3 3))"));
  transform::lexicographicOrder(*g, false);
  CHECK_WKT(g, "POLYGON((0 0,0 5,5 5,5 0,0 0),(1 1,2 1,2 2,1 2,1 1),(3 3,4 3,4 "
               "4,3 4,3 3))");

  std::unique_ptr<Geometry> g2(
      io::readWkt("POLYGON((0 0,0 5,5 5,5 0,0 0),(2 2,1 2,1 1,2 1,2 2),(3 3,4 "
                  "3,4 4,3 4,3 3))"));
  transform::lexicographicOrder(*g2, false);
  CHECK_WKT(g2, "POLYGON((0 0,0 5,5 5,5 0,0 0),(2 2,1 2,1 1,2 1,2 2),(3 3,4 "
                "3,4 4,3 4,3 3))");
}

BOOST_AUTO_TEST_CASE(testPolygon_reorder)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("POLYGON((0 0,0 5,5 5,5 0,0 0),(1 1,2 1,2 2,1 2,1 1),(3 3,4 "
                  "3,4 4,3 4,3 3))"));
  transform::lexicographicOrder(*g, true);
  CHECK_WKT(g, "POLYGON((0 0,0 5,5 5,5 0,0 0),(1 1,2 1,2 2,1 2,1 1),(3 3,4 3,4 "
               "4,3 4,3 3))");

  std::unique_ptr<Geometry> g2(
      io::readWkt("POLYGON((0 0,0 5,5 5,5 0,0 0),(2 2,1 2,1 1,2 1,2 2),(3 3,4 "
                  "3,4 4,3 4,3 3))"));
  transform::lexicographicOrder(*g2, true);
  CHECK_WKT(g2, "POLYGON((0 0,0 5,5 5,5 0,0 0),(1 1,2 1,2 2,1 2,1 1),(3 3,4 "
                "3,4 4,3 4,3 3))");
}

BOOST_AUTO_TEST_CASE(testMultiPolygon)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("MULTIPOLYGON(((0 0,0 1,1 1,1 0,0 0)),((2 2,2 3,3 3,3 2,2 "
                  "2)),((3 3,4 3,4 4,3 4,3 3)))"));
  transform::lexicographicOrder(*g, false);
  CHECK_WKT(g, "MULTIPOLYGON(((0 0,0 1,1 1,1 0,0 0)),((2 2,2 3,3 3,3 2,2 "
               "2)),((3 3,4 3,4 4,3 4,3 3)))");

  std::unique_ptr<Geometry> g2(
      io::readWkt("MULTIPOLYGON(((2 2,2 3,3 3,3 2,2 2)),((0 0,0 1,1 1,1 0,0 "
                  "0)),((3 3,4 3,4 4,3 4,3 3)))"));
  transform::lexicographicOrder(*g2, false);
  CHECK_WKT(g2, "MULTIPOLYGON(((0 0,0 1,1 1,1 0,0 0)),((2 2,2 3,3 3,3 2,2 "
                "2)),((3 3,4 3,4 4,3 4,3 3)))");

  std::unique_ptr<Geometry> g3(
      io::readWkt("MULTIPOLYGON Z(((2 2 2,2 3 3,3 3 3,3 2 3,2 2 2)),((0 0 3,0 "
                  "1 2,1 1 2,1 0 2,0 0 3)))"));
  transform::lexicographicOrder(*g3, false);
  CHECK_WKT(g3, "MULTIPOLYGON Z(((0 0 3,0 1 2,1 1 2,1 0 2,0 0 3)),((2 2 2,2 3 "
                "3,3 3 3,3 2 3,2 2 2)))");
}

BOOST_AUTO_TEST_CASE(testMultiPolygon_reorder)
{
  std::unique_ptr<Geometry> g(io::readWkt(
      "MULTIPOLYGON(((0 0,0 1,1 1,1 0,0 0)),((2 2,2 3,3 3,3 2,2 2)))"));
  transform::lexicographicOrder(*g, true);
  CHECK_WKT(g, "MULTIPOLYGON(((0 0,0 1,1 1,1 0,0 0)),((2 2,2 3,3 3,3 2,2 2)))");

  std::unique_ptr<Geometry> g2(io::readWkt(
      "MULTIPOLYGON(((1 1,2 1,2 2,1 2,1 1)),((0 0,1 0,1 1,0 1,0 0)))"));
  transform::lexicographicOrder(*g2, true);
  CHECK_WKT(g2,
            "MULTIPOLYGON(((0 0,1 0,1 1,0 1,0 0)),((1 1,2 1,2 2,1 2,1 1)))");
}

BOOST_AUTO_TEST_CASE(testTriangle)
{
  std::unique_ptr<Geometry> g(io::readWkt("TRIANGLE((0 0,1 0,0 1,0 0))"));
  transform::lexicographicOrder(*g, true);
  CHECK_WKT(g, "TRIANGLE((0 0,1 0,0 1,0 0))");
}

BOOST_AUTO_TEST_CASE(testTriangle_no_reorder)
{
  std::unique_ptr<Geometry> g(io::readWkt("TRIANGLE((0 0,1 0,0 1,0 0))"));
  transform::lexicographicOrder(*g, false);
  CHECK_WKT(g, "TRIANGLE((0 0,1 0,0 1,0 0))");
}

BOOST_AUTO_TEST_CASE(testPolyhedralSurface)
{
  std::unique_ptr<Geometry> g(io::readWkt(
      "POLYHEDRALSURFACE Z(((5.0 5.0 5.0,5.0 5.0 6.0,5.0 6.0 5.0,5.0 5.0 "
      "5.0)),((5.0 6.0 5.0,5.0 5.0 6.0,5.0 6.0 6.0,5.0 6.0 5.0)),((5.0 5.0 "
      "6.0,5.0 5.0 5.0,6.0 5.0 5.0,5.0 5.0 6.0)),((6.0 5.0 6.0,5.0 5.0 6.0,6.0 "
      "5.0 5.0,6.0 5.0 6.0)),((6.0 6.0 5.0,6.0 5.0 5.0,5.0 6.0 5.0,6.0 6.0 "
      "5.0)),((5.0 6.0 5.0,6.0 5.0 5.0,5.0 5.0 5.0,5.0 6.0 5.0)),((6.0 6.0 "
      "6.0,6.0 6.0 5.0,5.0 6.0 6.0,6.0 6.0 6.0)),((5.0 6.0 6.0,6.0 6.0 5.0,5.0 "
      "6.0 5.0,5.0 6.0 6.0)),((5.0 6.0 6.0,5.0 5.0 6.0,6.0 5.0 6.0,5.0 6.0 "
      "6.0)),((6.0 6.0 6.0,5.0 6.0 6.0,6.0 5.0 6.0,6.0 6.0 6.0)),((6.0 6.0 "
      "5.0,6.0 6.0 6.0,6.0 5.0 6.0,6.0 6.0 5.0)),((6.0 5.0 5.0,6.0 6.0 5.0,6.0 "
      "5.0 6.0,6.0 5.0 5.0)))"));
  transform::lexicographicOrder(*g, true);
  CHECK_WKT(g, "POLYHEDRALSURFACE Z(((5 5 5,5 5 6,5 6 5,5 5 5)),((5 5 5,6 5 "
               "5,5 5 6,5 5 5)),((5 5 5,5 6 5,6 5 5,5 5 5)),((5 5 6,5 6 6,5 6 "
               "5,5 5 6)),((5 5 6,6 5 5,6 5 6,5 5 6)),((5 5 6,6 5 6,5 6 6,5 5 "
               "6)),((5 6 5,6 6 5,6 5 5,5 6 5)),((5 6 5,5 6 6,6 6 5,5 6 "
               "5)),((5 6 6,6 6 6,6 6 5,5 6 6)),((5 6 6,6 5 6,6 6 6,5 6 "
               "6)),((6 5 5,6 6 5,6 5 6,6 5 5)),((6 5 6,6 6 5,6 6 6,6 5 6)))");
}

BOOST_AUTO_TEST_CASE(testPolyhedralSurface_no_reorder)
{
  std::unique_ptr<Geometry> g(io::readWkt(
      "POLYHEDRALSURFACE Z(((5.0 5.0 5.0,5.0 5.0 6.0,5.0 6.0 5.0,5.0 5.0 "
      "5.0)),((5.0 6.0 5.0,5.0 5.0 6.0,5.0 6.0 6.0,5.0 6.0 5.0)),((5.0 5.0 "
      "6.0,5.0 5.0 5.0,6.0 5.0 5.0,5.0 5.0 6.0)),((6.0 5.0 6.0,5.0 5.0 6.0,6.0 "
      "5.0 5.0,6.0 5.0 6.0)),((6.0 6.0 5.0,6.0 5.0 5.0,5.0 6.0 5.0,6.0 6.0 "
      "5.0)),((5.0 6.0 5.0,6.0 5.0 5.0,5.0 5.0 5.0,5.0 6.0 5.0)),((6.0 6.0 "
      "6.0,6.0 6.0 5.0,5.0 6.0 6.0,6.0 6.0 6.0)),((5.0 6.0 6.0,6.0 6.0 5.0,5.0 "
      "6.0 5.0,5.0 6.0 6.0)),((5.0 6.0 6.0,5.0 5.0 6.0,6.0 5.0 6.0,5.0 6.0 "
      "6.0)),((6.0 6.0 6.0,5.0 6.0 6.0,6.0 5.0 6.0,6.0 6.0 6.0)),((6.0 6.0 "
      "5.0,6.0 6.0 6.0,6.0 5.0 6.0,6.0 6.0 5.0)),((6.0 5.0 5.0,6.0 6.0 5.0,6.0 "
      "5.0 6.0,6.0 5.0 5.0)))"));
  transform::lexicographicOrder(*g, false);
  CHECK_WKT(g, "POLYHEDRALSURFACE Z(((5 5 5,5 5 6,5 6 5,5 5 5)),((5 5 6,5 5 "
               "5,6 5 5,5 5 6)),((5 6 5,6 5 5,5 5 5,5 6 5)),((5 6 5,5 5 6,5 6 "
               "6,5 6 5)),((6 5 6,5 5 6,6 5 5,6 5 6)),((5 6 6,5 5 6,6 5 6,5 6 "
               "6)),((6 6 5,6 5 5,5 6 5,6 6 5)),((5 6 6,6 6 5,5 6 5,5 6 "
               "6)),((6 6 6,6 6 5,5 6 6,6 6 6)),((6 6 6,5 6 6,6 5 6,6 6 "
               "6)),((6 5 5,6 6 5,6 5 6,6 5 5)),((6 6 5,6 6 6,6 5 6,6 6 5)))");
}

BOOST_AUTO_TEST_CASE(testSolid)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("SOLID((((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0),(0 0 1,1 0 1,1 1 "
                  "1,0 1 1,0 0 1)),((0 0 0,0 0 1,0 1 1,0 1 0,0 0 0)),((0 0 0,1 "
                  "0 0,1 0 1,0 0 1,0 0 0)),((1 0 0,1 1 0,1 1 1,1 0 1,1 0 "
                  "0)),((0 1 0,0 1 1,1 1 1,1 1 0,0 1 0))))"));
  transform::lexicographicOrder(*g, true);
  CHECK_WKT(g, "SOLID Z((((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0),(0 0 1,1 0 1,1 1 1,0 "
               "1 1,0 0 1)),((0 0 0,0 0 1,0 1 1,0 1 0,0 0 0)),((0 0 0,1 0 0,1 "
               "0 1,0 0 1,0 0 0)),((0 1 0,0 1 1,1 1 1,1 1 0,0 1 0)),((1 0 0,1 "
               "1 0,1 1 1,1 0 1,1 0 0))))");
}

BOOST_AUTO_TEST_CASE(testSolid_no_reorder)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("SOLID((((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0),(0 0 1,1 0 1,1 1 "
                  "1,0 1 1,0 0 1)),((0 0 0,0 0 1,0 1 1,0 1 0,0 0 0)),((0 0 0,1 "
                  "0 0,1 0 1,0 0 1,0 0 0)),((1 0 0,1 1 0,1 1 1,1 0 1,1 0 "
                  "0)),((0 1 0,0 1 1,1 1 1,1 1 0,0 1 0))))"));
  transform::lexicographicOrder(*g, false);
  CHECK_WKT(g, "SOLID Z((((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0),(0 0 1,1 0 1,1 1 1,0 "
               "1 1,0 0 1)),((0 0 0,0 0 1,0 1 1,0 1 0,0 0 0)),((0 0 0,1 0 0,1 "
               "0 1,0 0 1,0 0 0)),((0 1 0,0 1 1,1 1 1,1 1 0,0 1 0)),((1 0 0,1 "
               "1 0,1 1 1,1 0 1,1 0 0))))");
}

BOOST_AUTO_TEST_CASE(testGeometryCollection)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("GEOMETRYCOLLECTION(POINT(1 1),LINESTRING(0 0,1 "
                  "1),POLYGON((0 0,1 0,1 1,0 1,0 0)))"));
  transform::lexicographicOrder(*g, true);
  CHECK_WKT(g, "GEOMETRYCOLLECTION(POINT(1 1),LINESTRING(0 0,1 1),POLYGON((0 "
               "0,1 0,1 1,0 1,0 0)))");
}

BOOST_AUTO_TEST_CASE(testGeometryCollection_no_reorder)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("GEOMETRYCOLLECTION(POINT(1 1),LINESTRING(0 0,1 "
                  "1),POLYGON((0 0,1 0,1 1,0 1,0 0)))"));
  transform::lexicographicOrder(*g, false);
  CHECK_WKT(g, "GEOMETRYCOLLECTION(POINT(1 1),LINESTRING(0 0,1 1),POLYGON((0 "
               "0,1 0,1 1,0 1,0 0)))");
}

BOOST_AUTO_TEST_CASE(testSpecificPolyhedralSurfaceEquivalence)
{
  std::string wkt1 =
      "POLYHEDRALSURFACE Z(((5.0 5.0 5.0,5.0 5.0 6.0,5.0 6.0 5.0,5.0 5.0 "
      "5.0)),((5.0 6.0 5.0,5.0 5.0 6.0,5.0 6.0 6.0,5.0 6.0 5.0)),((5.0 5.0 "
      "6.0,5.0 5.0 5.0,6.0 5.0 5.0,5.0 5.0 6.0)),((6.0 5.0 6.0,5.0 5.0 6.0,6.0 "
      "5.0 5.0,6.0 5.0 6.0)),((6.0 6.0 5.0,6.0 5.0 5.0,5.0 6.0 5.0,6.0 6.0 "
      "5.0)),((5.0 6.0 5.0,6.0 5.0 5.0,5.0 5.0 5.0,5.0 6.0 5.0)),((6.0 6.0 "
      "6.0,6.0 6.0 5.0,5.0 6.0 6.0,6.0 6.0 6.0)),((5.0 6.0 6.0,6.0 6.0 5.0,5.0 "
      "6.0 5.0,5.0 6.0 6.0)),((5.0 6.0 6.0,5.0 5.0 6.0,6.0 5.0 6.0,5.0 6.0 "
      "6.0)),((6.0 6.0 6.0,5.0 6.0 6.0,6.0 5.0 6.0,6.0 6.0 6.0)),((6.0 6.0 "
      "5.0,6.0 6.0 6.0,6.0 5.0 6.0,6.0 6.0 5.0)),((6.0 5.0 5.0,6.0 6.0 5.0,6.0 "
      "5.0 6.0,6.0 5.0 5.0)))";
  std::string wkt2 =
      "POLYHEDRALSURFACE Z(((5.0 6.0 6.0,5.0 6.0 5.0,5.0 5.0 6.0,5.0 6.0 "
      "6.0)),((5.0 5.0 6.0,5.0 6.0 5.0,5.0 5.0 5.0,5.0 5.0 6.0)),((6.0 5.0 "
      "5.0,6.0 5.0 6.0,5.0 5.0 6.0,6.0 5.0 5.0)),((5.0 5.0 5.0,6.0 5.0 5.0,5.0 "
      "5.0 6.0,5.0 5.0 5.0)),((5.0 5.0 5.0,5.0 6.0 5.0,6.0 5.0 5.0,5.0 5.0 "
      "5.0)),((6.0 5.0 5.0,5.0 6.0 5.0,6.0 6.0 5.0,6.0 5.0 5.0)),((6.0 5.0 "
      "6.0,6.0 6.0 6.0,5.0 6.0 6.0,6.0 5.0 6.0)),((5.0 5.0 6.0,6.0 5.0 6.0,5.0 "
      "6.0 6.0,5.0 5.0 6.0)),((5.0 6.0 5.0,5.0 6.0 6.0,6.0 6.0 5.0,5.0 6.0 "
      "5.0)),((6.0 6.0 5.0,5.0 6.0 6.0,6.0 6.0 6.0,6.0 6.0 5.0)),((6.0 5.0 "
      "6.0,6.0 5.0 5.0,6.0 6.0 5.0,6.0 5.0 6.0)),((6.0 6.0 6.0,6.0 5.0 6.0,6.0 "
      "6.0 5.0,6.0 6.0 6.0)))";
  std::unique_ptr<Geometry> g1(io::readWkt(wkt1));
  std::unique_ptr<Geometry> g2(io::readWkt(wkt2));

  transform::lexicographicOrder(*g1, false);
  transform::lexicographicOrder(*g2, false);

  BOOST_CHECK_EQUAL(g1->asText(1), g2->asText(1));
}

BOOST_AUTO_TEST_CASE(testSpecificPolyhedralSurfaceEquivalence_order)
{
  std::string wkt1 =
      "POLYHEDRALSURFACE Z(((5.0 5.0 5.0,5.0 5.0 6.0,5.0 6.0 5.0,5.0 5.0 "
      "5.0)),((5.0 6.0 5.0,5.0 5.0 6.0,5.0 6.0 6.0,5.0 6.0 5.0)),((5.0 5.0 "
      "6.0,5.0 5.0 5.0,6.0 5.0 5.0,5.0 5.0 6.0)),((6.0 5.0 6.0,5.0 5.0 6.0,6.0 "
      "5.0 5.0,6.0 5.0 6.0)),((6.0 6.0 5.0,6.0 5.0 5.0,5.0 6.0 5.0,6.0 6.0 "
      "5.0)),((5.0 6.0 5.0,6.0 5.0 5.0,5.0 5.0 5.0,5.0 6.0 5.0)),((6.0 6.0 "
      "6.0,6.0 6.0 5.0,5.0 6.0 6.0,6.0 6.0 6.0)),((5.0 6.0 6.0,6.0 6.0 5.0,5.0 "
      "6.0 5.0,5.0 6.0 6.0)),((5.0 6.0 6.0,5.0 5.0 6.0,6.0 5.0 6.0,5.0 6.0 "
      "6.0)),((6.0 6.0 6.0,5.0 6.0 6.0,6.0 5.0 6.0,6.0 6.0 6.0)),((6.0 6.0 "
      "5.0,6.0 6.0 6.0,6.0 5.0 6.0,6.0 6.0 5.0)),((6.0 5.0 5.0,6.0 6.0 5.0,6.0 "
      "5.0 6.0,6.0 5.0 5.0)))";
  std::string wkt2 =
      "POLYHEDRALSURFACE Z(((5.0 6.0 6.0,5.0 6.0 5.0,5.0 5.0 6.0,5.0 6.0 "
      "6.0)),((5.0 5.0 6.0,5.0 6.0 5.0,5.0 5.0 5.0,5.0 5.0 6.0)),((6.0 5.0 "
      "5.0,6.0 5.0 6.0,5.0 5.0 6.0,6.0 5.0 5.0)),((5.0 5.0 5.0,6.0 5.0 5.0,5.0 "
      "5.0 6.0,5.0 5.0 5.0)),((5.0 5.0 5.0,5.0 6.0 5.0,6.0 5.0 5.0,5.0 5.0 "
      "5.0)),((6.0 5.0 5.0,5.0 6.0 5.0,6.0 6.0 5.0,6.0 5.0 5.0)),((6.0 5.0 "
      "6.0,6.0 6.0 6.0,5.0 6.0 6.0,6.0 5.0 6.0)),((5.0 5.0 6.0,6.0 5.0 6.0,5.0 "
      "6.0 6.0,5.0 5.0 6.0)),((5.0 6.0 5.0,5.0 6.0 6.0,6.0 6.0 5.0,5.0 6.0 "
      "5.0)),((6.0 6.0 5.0,5.0 6.0 6.0,6.0 6.0 6.0,6.0 6.0 5.0)),((6.0 5.0 "
      "6.0,6.0 5.0 5.0,6.0 6.0 5.0,6.0 5.0 6.0)),((6.0 6.0 6.0,6.0 5.0 6.0,6.0 "
      "6.0 5.0,6.0 6.0 6.0)))";
  std::unique_ptr<Geometry> g1(io::readWkt(wkt1));
  std::unique_ptr<Geometry> g2(io::readWkt(wkt2));

  transform::lexicographicOrder(*g1, true);
  transform::lexicographicOrder(*g2, true);

  BOOST_CHECK_EQUAL(g1->asText(1), g2->asText(1));
}

BOOST_AUTO_TEST_SUITE_END()
