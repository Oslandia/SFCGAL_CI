#include <boost/test/unit_test.hpp>

#include "SFCGAL/GeometryCollection.h"
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
#include "SFCGAL/algorithm/lexicographicOrder.h"
#include "SFCGAL/io/wkt.h"

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_LexicographicOrderTest)

BOOST_AUTO_TEST_CASE(testMultiPoint)
{
  std::unique_ptr<Geometry> g(io::readWkt("MULTIPOINT((1 1),(0 0))"));
  auto                      ordered = transform::lexicographicOrder(*g);
  BOOST_CHECK_EQUAL(ordered->asText(0), "MULTIPOINT((0 0),(1 1))");
}

BOOST_AUTO_TEST_CASE(testMultiLineString)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("MULTILINESTRING((0 0,1 1),(2 2,1 0))"));
  auto ordered = transform::lexicographicOrder(*g);
  BOOST_CHECK_EQUAL(ordered->asText(0), "MULTILINESTRING((0 0,1 1),(2 2,1 0))");

  std::unique_ptr<Geometry> g2(
      io::readWkt("MULTILINESTRING((2 2,1 0),(0 0,1 1))"));
  auto ordered2 = transform::lexicographicOrder(*g2);
  BOOST_CHECK_EQUAL(ordered2->asText(0),
                    "MULTILINESTRING((0 0,1 1),(2 2,1 0))");

  std::unique_ptr<Geometry> g3(
      io::readWkt("MULTILINESTRING Z((2 2 -1,1 0 -1),(0 0 0,1 1 0))"));
  auto ordered3 = transform::lexicographicOrder(*g3);
  BOOST_CHECK_EQUAL(ordered3->asText(0),
                    "MULTILINESTRING Z((0 0 0,1 1 0),(2 2 -1,1 0 -1))");
}

BOOST_AUTO_TEST_CASE(testPoint)
{
  std::unique_ptr<Geometry> g(io::readWkt("POINT(1 2)"));
  std::unique_ptr<Geometry> ordered = transform::lexicographicOrder(*g);
  BOOST_CHECK_EQUAL(ordered->asText(0), "POINT(1 2)");
}

BOOST_AUTO_TEST_CASE(testLineString)
{
  std::unique_ptr<Geometry> g(io::readWkt("LINESTRING(0 0,-1 -1,-2 -2)"));
  std::unique_ptr<Geometry> ordered = transform::lexicographicOrder(*g);
  BOOST_CHECK_EQUAL(ordered->asText(0), "LINESTRING(0 0,-1 -1,-2 -2)");
}

BOOST_AUTO_TEST_CASE(testPolygon)
{
  std::string wkt = "POLYGON((0 0,0 5,5 5,5 0,0 0),(2 2,1 2,1 1,2 1,2 2),(3 "
                    "3,4 3,4 4,3 4,3 3))";
  std::unique_ptr<Geometry> g(io::readWkt(wkt));

  std::unique_ptr<Geometry> ordered = transform::lexicographicOrder(*g);

  BOOST_CHECK_EQUAL(ordered->asText(0),
                    "POLYGON((0 0,0 5,5 5,5 0,0 0),(2 2,1 2,1 1,2 1,2 2),(3 "
                    "3,4 3,4 4,3 4,3 3))");

  std::string wkt2 = "POLYGON((3 3,4 3,4 4,3 4,3 3),(0 0,0 5,5 5,5 0,0 0),(2 "
                     "2,1 2,1 1,2 1,2 2))";
  std::unique_ptr<Geometry> g2(io::readWkt(wkt2));

  std::unique_ptr<Geometry> ordered2 = transform::lexicographicOrder(*g2);

  BOOST_CHECK_EQUAL(ordered2->asText(0),
                    "POLYGON((0 0,0 5,5 5,5 0,0 0),(2 2,1 2,1 1,2 1,2 2),(3 "
                    "3,4 3,4 4,3 4,3 3))");
}

BOOST_AUTO_TEST_CASE(testMultiPolygon)
{
  std::string wkt = "MULTIPOLYGON(((2 2,2 3,3 3,3 2,2 2)),((0 0,0 1,1 1,1 0,0 "
                    "0)),((3 3,4 3,4 4,3 4,3 3)))";
  std::unique_ptr<Geometry> g(io::readWkt(wkt));

  std::unique_ptr<Geometry> ordered = transform::lexicographicOrder(*g);

  BOOST_CHECK_EQUAL(ordered->asText(0),
                    "MULTIPOLYGON(((0 0,0 1,1 1,1 0,0 0)),((2 2,2 3,3 3,3 2,2 "
                    "2)),((3 3,4 3,4 4,3 4,3 3)))");
  std::unique_ptr<Geometry> g2(
      io::readWkt("MULTIPOLYGON(((2 2,3 2,3 3,2 3,2 2)),((0 0,1 0,1 1,0 1,0 0)))"));
  auto ordered2 = transform::lexicographicOrder(*g2);
  BOOST_CHECK_EQUAL(ordered2->asText(0),
                    "MULTIPOLYGON(((0 0,1 0,1 1,0 1,0 0)),((2 2,3 2,3 3,2 3,2 2)))");
}

BOOST_AUTO_TEST_CASE(testTriangle)
{
  std::unique_ptr<Geometry> g(io::readWkt("TRIANGLE((2 2,0 0,1 1,2 2))"));
  auto ordered = transform::lexicographicOrder(*g);
  BOOST_CHECK_EQUAL(ordered->asText(0), "TRIANGLE((0 0,1 1,2 2,0 0))");
}


BOOST_AUTO_TEST_CASE(testPolyhedralSurface)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("POLYHEDRALSURFACE Z(((2 0 0,3 0 0,3 1 0,2 0 0)),((0 0 0,1 0 0,1 1 0,0 0 0)))"));
  auto ordered = transform::lexicographicOrder(*g);
  BOOST_CHECK_EQUAL(ordered->asText(0),
                    "POLYHEDRALSURFACE Z(((0 0 0,1 0 0,1 1 0,0 0 0)),((2 0 0,3 0 0,3 1 0,2 0 0)))");
}

BOOST_AUTO_TEST_CASE(testGeometryCollection)
{
  std::unique_ptr<Geometry> g(io::readWkt("GEOMETRYCOLLECTION(LINESTRING(0 0,1 1),POINT(0 0))"));
  auto ordered = transform::lexicographicOrder(*g);
  BOOST_CHECK_EQUAL(ordered->asText(0), "GEOMETRYCOLLECTION(POINT(0 0),LINESTRING(0 0,1 1))");
}

BOOST_AUTO_TEST_CASE(testMultiSolid)
{
  std::string wkt = "MULTISOLID("
                    "((((2 2 2,3 2 2,3 3 2,2 2 2,2 2 2))))," 
                    "((((0 0 0,1 0 0,1 1 0,0 0 0,0 0 0))))"
                    ")";
  std::unique_ptr<Geometry> g(io::readWkt(wkt));
  auto ordered = transform::lexicographicOrder(*g);
  BOOST_CHECK_EQUAL(ordered->asText(0),
                    "MULTISOLID Z("
                    "((((0 0 0,1 0 0,1 1 0,0 0 0,0 0 0))))," 
                    "((((2 2 2,3 2 2,3 3 2,2 2 2,2 2 2))))"
                    ")");
}

BOOST_AUTO_TEST_CASE(testSpecificPolyhedralSurfaceEquivalence_Minko1)
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

  g1 = transform::lexicographicOrder(*g1);
  g2 = transform::lexicographicOrder(*g2);

  BOOST_CHECK_EQUAL(g1->asText(1), g2->asText(1));
}

BOOST_AUTO_TEST_CASE(testSpecificPolyhedralSurfaceEquivalence_Minko2)
{
  std::string wkt1 = "POLYHEDRALSURFACE Z(((2.0 2.0 2.0,2.0 2.0 3.0,2.0 3.0 2.0,2.0 2.0 2.0)),((2.0 3.0 2.0,2.0 2.0 3.0,2.0 3.0 3.0,2.0 3.0 2.0)),((2.0 2.0 3.0,2.0 2.0 2.0,3.0 2.0 2.0,2.0 2.0 3.0)),((3.0 2.0 3.0,2.0 2.0 3.0,3.0 2.0 2.0,3.0 2.0 3.0)),((3.0 3.0 2.0,3.0 2.0 2.0,2.0 3.0 2.0,3.0 3.0 2.0)),((2.0 3.0 2.0,3.0 2.0 2.0,2.0 2.0 2.0,2.0 3.0 2.0)),((3.0 3.0 3.0,3.0 3.0 2.0,2.0 3.0 3.0,3.0 3.0 3.0)),((2.0 3.0 3.0,3.0 3.0 2.0,2.0 3.0 2.0,2.0 3.0 3.0)),((2.0 3.0 3.0,2.0 2.0 3.0,3.0 2.0 3.0,2.0 3.0 3.0)),((3.0 3.0 3.0,2.0 3.0 3.0,3.0 2.0 3.0,3.0 3.0 3.0)),((3.0 3.0 2.0,3.0 3.0 3.0,3.0 2.0 3.0,3.0 3.0 2.0)),((3.0 2.0 2.0,3.0 3.0 2.0,3.0 2.0 3.0,3.0 2.0 2.0)))"; 
  std::string wkt2 = "POLYHEDRALSURFACE Z(((2.0 3.0 3.0,2.0 3.0 2.0,2.0 2.0 3.0,2.0 3.0 3.0)),((2.0 2.0 3.0,2.0 3.0 2.0,2.0 2.0 2.0,2.0 2.0 3.0)),((3.0 2.0 2.0,3.0 2.0 3.0,2.0 2.0 3.0,3.0 2.0 2.0)),((2.0 2.0 2.0,3.0 2.0 2.0,2.0 2.0 3.0,2.0 2.0 2.0)),((2.0 2.0 2.0,2.0 3.0 2.0,3.0 2.0 2.0,2.0 2.0 2.0)),((3.0 2.0 2.0,2.0 3.0 2.0,3.0 3.0 2.0,3.0 2.0 2.0)),((3.0 2.0 3.0,3.0 3.0 3.0,2.0 3.0 3.0,3.0 2.0 3.0)),((2.0 2.0 3.0,3.0 2.0 3.0,2.0 3.0 3.0,2.0 2.0 3.0)),((2.0 3.0 2.0,2.0 3.0 3.0,3.0 3.0 2.0,2.0 3.0 2.0)),((3.0 3.0 2.0,2.0 3.0 3.0,3.0 3.0 3.0,3.0 3.0 2.0)),((3.0 2.0 3.0,3.0 2.0 2.0,3.0 3.0 2.0,3.0 2.0 3.0)),((3.0 3.0 3.0,3.0 2.0 3.0,3.0 3.0 2.0,3.0 3.0 3.0)))";
  std::unique_ptr<Geometry> g1(io::readWkt(wkt1));
  std::unique_ptr<Geometry> g2(io::readWkt(wkt2));


  g1 = transform::lexicographicOrder(*g1);
  g2 = transform::lexicographicOrder(*g2);

  BOOST_CHECK_EQUAL(g1->asText(1), g2->asText(1));
}


BOOST_AUTO_TEST_CASE(testSpecificPolyhedralSurfaceEquivalence_Minko3)
{
  std::string wkt1 = "POLYHEDRALSURFACE Z(((0.0 1.0 0.0,0.0 0.0 0.0,0.0 0.0 1.0,0.0 1.0 0.0)),((0.0 0.0 0.0,0.0 0.0 0.0,0.0 0.0 1.0,0.0 0.0 0.0)),((0.0 0.0 1.0,0.0 0.0 0.0,0.0 0.0 1.0,0.0 0.0 1.0)),((0.0 0.0 0.0,0.0 0.0 0.0,0.0 0.0 0.0,0.0 0.0 0.0)),((0.0 0.0 0.0,0.0 0.0 0.0,0.0 0.0 0.0,0.0 0.0 0.0)),((0.0 0.0 0.0,0.0 1.0 0.0,0.0 0.0 0.0,0.0 0.0 0.0)),((0.0 0.0 0.0,0.0 1.0 0.0,0.0 1.0 0.0,0.0 0.0 0.0)),((0.0 1.0 0.0,1.0 2.0 1.0,0.0 1.0 0.0,0.0 1.0 0.0)),((0.0 1.0 0.0,1.0 2.0 1.0,1.0 2.0 1.0,0.0 1.0 0.0)),((0.0 0.0 1.0,1.0 1.0 2.0,0.0 1.0 0.0,0.0 0.0 1.0)),((0.0 1.0 0.0,1.0 1.0 2.0,1.0 2.0 1.0,0.0 1.0 0.0)),((0.0 0.0 1.0,0.0 0.0 1.0,1.0 1.0 2.0,0.0 0.0 1.0)),((1.0 1.0 2.0,0.0 0.0 1.0,1.0 1.0 2.0,1.0 1.0 2.0)),((1.0 1.0 2.0,0.0 0.0 1.0,0.0 0.0 1.0,1.0 1.0 2.0)),((1.0 1.0 2.0,1.0 1.0 2.0,0.0 0.0 1.0,1.0 1.0 2.0)),((0.0 0.0 1.0,0.0 0.0 0.0,0.0 0.0 0.0,0.0 0.0 1.0)),((0.0 0.0 1.0,0.0 0.0 1.0,0.0 0.0 0.0,0.0 0.0 1.0)),((0.0 0.0 1.0,0.0 0.0 0.0,1.0 0.0 0.0,0.0 0.0 1.0)),((0.0 0.0 0.0,1.0 0.0 0.0,0.0 0.0 0.0,0.0 0.0 0.0)),((0.0 0.0 0.0,1.0 0.0 0.0,1.0 0.0 0.0,0.0 0.0 0.0)),((1.0 0.0 0.0,0.0 0.0 0.0,0.0 1.0 0.0,1.0 0.0 0.0)),((0.0 1.0 0.0,1.0 2.0 1.0,1.0 0.0 0.0,0.0 1.0 0.0)),((1.0 0.0 0.0,1.0 2.0 1.0,2.0 1.0 1.0,1.0 0.0 0.0)),((1.0 1.0 2.0,1.0 1.0 2.0,1.0 2.0 1.0,1.0 1.0 2.0)),((1.0 2.0 1.0,1.0 1.0 2.0,1.0 1.0 2.0,1.0 2.0 1.0)),((1.0 2.0 1.0,1.0 1.0 2.0,2.0 1.0 1.0,1.0 2.0 1.0)),((1.0 2.0 1.0,2.0 1.0 1.0,2.0 1.0 1.0,1.0 2.0 1.0)),((1.0 2.0 1.0,2.0 1.0 1.0,1.0 2.0 1.0,1.0 2.0 1.0)),((1.0 1.0 2.0,0.0 0.0 1.0,1.0 0.0 0.0,1.0 1.0 2.0)),((2.0 1.0 1.0,1.0 1.0 2.0,1.0 0.0 0.0,2.0 1.0 1.0)),((1.0 0.0 0.0,1.0 0.0 0.0,2.0 1.0 1.0,1.0 0.0 0.0)),((2.0 1.0 1.0,1.0 0.0 0.0,2.0 1.0 1.0,2.0 1.0 1.0)))";
  std::string wkt2 = "POLYHEDRALSURFACE Z(((0.0 0.0 1.0,0.0 1.0 0.0,0.0 0.0 0.0,0.0 0.0 1.0)),((0.0 0.0 1.0,0.0 0.0 1.0,0.0 0.0 0.0,0.0 0.0 1.0)),((0.0 0.0 0.0,0.0 0.0 1.0,0.0 0.0 0.0,0.0 0.0 0.0)),((0.0 0.0 0.0,0.0 0.0 0.0,0.0 0.0 0.0,0.0 0.0 0.0)),((0.0 0.0 0.0,0.0 0.0 0.0,0.0 0.0 0.0,0.0 0.0 0.0)),((0.0 1.0 0.0,0.0 0.0 0.0,0.0 1.0 0.0,0.0 1.0 0.0)),((0.0 1.0 0.0,0.0 0.0 0.0,0.0 0.0 0.0,0.0 1.0 0.0)),((1.0 1.0 2.0,1.0 1.0 2.0,0.0 0.0 1.0,1.0 1.0 2.0)),((0.0 0.0 1.0,1.0 1.0 2.0,0.0 0.0 1.0,0.0 0.0 1.0)),((1.0 2.0 1.0,0.0 1.0 0.0,1.0 1.0 2.0,1.0 2.0 1.0)),((1.0 1.0 2.0,0.0 1.0 0.0,0.0 0.0 1.0,1.0 1.0 2.0)),((1.0 2.0 1.0,0.0 1.0 0.0,1.0 2.0 1.0,1.0 2.0 1.0)),((1.0 2.0 1.0,0.0 1.0 0.0,0.0 1.0 0.0,1.0 2.0 1.0)),((0.0 0.0 0.0,0.0 0.0 1.0,0.0 0.0 1.0,0.0 0.0 0.0)),((0.0 0.0 0.0,0.0 0.0 0.0,0.0 0.0 1.0,0.0 0.0 0.0)),((0.0 0.0 1.0,1.0 1.0 2.0,1.0 1.0 2.0,0.0 0.0 1.0)),((0.0 0.0 1.0,0.0 0.0 1.0,1.0 1.0 2.0,0.0 0.0 1.0)),((1.0 0.0 0.0,0.0 0.0 0.0,1.0 0.0 0.0,1.0 0.0 0.0)),((1.0 0.0 0.0,0.0 0.0 0.0,0.0 0.0 0.0,1.0 0.0 0.0)),((0.0 1.0 0.0,1.0 0.0 0.0,0.0 0.0 0.0,0.0 1.0 0.0)),((1.0 0.0 0.0,0.0 0.0 1.0,0.0 0.0 0.0,1.0 0.0 0.0)),((2.0 1.0 1.0,1.0 0.0 0.0,1.0 2.0 1.0,2.0 1.0 1.0)),((1.0 2.0 1.0,1.0 0.0 0.0,0.0 1.0 0.0,1.0 2.0 1.0)),((1.0 1.0 2.0,2.0 1.0 1.0,1.0 1.0 2.0,1.0 1.0 2.0)),((1.0 1.0 2.0,2.0 1.0 1.0,2.0 1.0 1.0,1.0 1.0 2.0)),((1.0 1.0 2.0,2.0 1.0 1.0,1.0 2.0 1.0,1.0 1.0 2.0)),((1.0 1.0 2.0,1.0 2.0 1.0,1.0 2.0 1.0,1.0 1.0 2.0)),((1.0 1.0 2.0,1.0 2.0 1.0,1.0 1.0 2.0,1.0 1.0 2.0)),((1.0 0.0 0.0,2.0 1.0 1.0,1.0 1.0 2.0,1.0 0.0 0.0)),((0.0 0.0 1.0,1.0 0.0 0.0,1.0 1.0 2.0,0.0 0.0 1.0)),((2.0 1.0 1.0,2.0 1.0 1.0,1.0 0.0 0.0,2.0 1.0 1.0)),((1.0 0.0 0.0,2.0 1.0 1.0,1.0 0.0 0.0,1.0 0.0 0.0)))";
  std::unique_ptr<Geometry> g1(io::readWkt(wkt1));
  std::unique_ptr<Geometry> g2(io::readWkt(wkt2));

  g1 = transform::lexicographicOrder(*g1);
  g2 = transform::lexicographicOrder(*g2);


  BOOST_CHECK_EQUAL(g1->asText(1), g2->asText(1));
}
BOOST_AUTO_TEST_SUITE_END()
