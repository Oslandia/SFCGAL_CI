#include <boost/test/unit_test.hpp>

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/algorithm/compare.h"

using namespace SFCGAL;
using namespace SFCGAL::algorithm::compare;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_CompareTest)

BOOST_AUTO_TEST_CASE(testComparePoint)
{
  Point p1(1, 2, 3);
  Point p2(1, 2, 3);
  Point p3(1, 2, 3.000001);
  Point p4(1, 2);

  BOOST_CHECK(strictCompare(p1, p2));
  BOOST_CHECK(sortedCompare(p1, p2));
  BOOST_CHECK(fuzzyCompare(p1, p2, 1e-6));
  BOOST_CHECK(fuzzyCompare(p1, p2, 1e-6, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(!strictCompare(p1, p3));
  BOOST_CHECK(!sortedCompare(p1, p3));
  BOOST_CHECK(fuzzyCompare(p1, p3, 1e-5));
  BOOST_CHECK(fuzzyCompare(p1, p3, 1e-5, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(!strictCompare(p1, p4));
  BOOST_CHECK(!sortedCompare(p1, p4));
  BOOST_CHECK(!fuzzyCompare(p1, p4, 1e-6));
  BOOST_CHECK(!fuzzyCompare(p1, p4, 1e-6, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(sortedCompare(p1, p3, true, 1e-5, FuzzyCompareMethod::Equal));
  BOOST_CHECK(
      sortedCompare(p1, p3, true, 1e-5, FuzzyCompareMethod::DistanceEqual));
}

BOOST_AUTO_TEST_CASE(testCompareLineString)
{
  LineString ls1;
  ls1.addPoint(Point(0, 0));
  ls1.addPoint(Point(1, 1));

  LineString ls2;
  ls2.addPoint(Point(0, 0));
  ls2.addPoint(Point(1, 1));

  LineString ls3;
  ls3.addPoint(Point(1, 1));
  ls3.addPoint(Point(0, 0));

  LineString ls4;
  ls4.addPoint(Point(0, 0));
  ls4.addPoint(Point(1.000001, 1.000001));

  BOOST_CHECK(strictCompare(ls1, ls2));
  BOOST_CHECK(sortedCompare(ls1, ls2));
  BOOST_CHECK(fuzzyCompare(ls1, ls2, 1e-6));
  BOOST_CHECK(fuzzyCompare(ls1, ls2, 1e-6, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(!strictCompare(ls1, ls3));
  BOOST_CHECK(!sortedCompare(ls1, ls3));
  BOOST_CHECK(!fuzzyCompare(ls1, ls3, 1e-6));
  BOOST_CHECK(!fuzzyCompare(ls1, ls3, 1e-6, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(!strictCompare(ls1, ls4));
  BOOST_CHECK(!sortedCompare(ls1, ls4));
  BOOST_CHECK(fuzzyCompare(ls1, ls4, 1e-5));
  BOOST_CHECK(fuzzyCompare(ls1, ls4, 1e-5, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(sortedCompare(ls1, ls4, true, 1e-5, FuzzyCompareMethod::Equal));
  BOOST_CHECK(
      sortedCompare(ls1, ls4, true, 1e-5, FuzzyCompareMethod::DistanceEqual));
}

BOOST_AUTO_TEST_CASE(testComparePolygon)
{
  Polygon poly1;
  poly1.exteriorRing().addPoint(Point(0, 0));
  poly1.exteriorRing().addPoint(Point(1, 0));
  poly1.exteriorRing().addPoint(Point(1, 1));
  poly1.exteriorRing().addPoint(Point(0, 1));
  poly1.exteriorRing().addPoint(Point(0, 0));

  Polygon poly2(poly1);

  Polygon    poly3(poly1);
  LineString hole;
  hole.addPoint(Point(0.25, 0.25));
  hole.addPoint(Point(0.75, 0.25));
  hole.addPoint(Point(0.75, 0.75));
  hole.addPoint(Point(0.25, 0.75));
  hole.addPoint(Point(0.25, 0.25));
  poly3.addRing(hole);

  Polygon poly4;
  poly4.exteriorRing().addPoint(Point(0, 0));
  poly4.exteriorRing().addPoint(Point(1.000001, 0));
  poly4.exteriorRing().addPoint(Point(1.000001, 1.000001));
  poly4.exteriorRing().addPoint(Point(0, 1.000001));
  poly4.exteriorRing().addPoint(Point(0, 0));

  BOOST_CHECK(strictCompare(poly1, poly2));
  BOOST_CHECK(sortedCompare(poly1, poly2));
  BOOST_CHECK(fuzzyCompare(poly1, poly2, 1e-6));
  BOOST_CHECK(
      fuzzyCompare(poly1, poly2, 1e-6, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(!strictCompare(poly1, poly3));
  BOOST_CHECK(!sortedCompare(poly1, poly3));
  BOOST_CHECK(!fuzzyCompare(poly1, poly3, 1e-6));
  BOOST_CHECK(
      !fuzzyCompare(poly1, poly3, 1e-6, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(!strictCompare(poly1, poly4));
  BOOST_CHECK(!sortedCompare(poly1, poly4));
  BOOST_CHECK(fuzzyCompare(poly1, poly4, 1e-5));
  BOOST_CHECK(
      fuzzyCompare(poly1, poly4, 1e-5, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(
      sortedCompare(poly1, poly4, true, 1e-5, FuzzyCompareMethod::Equal));
  BOOST_CHECK(sortedCompare(poly1, poly4, true, 1e-5,
                            FuzzyCompareMethod::DistanceEqual));
}

BOOST_AUTO_TEST_CASE(testCompareTriangle)
{
  Triangle tri1(Point(0, 0), Point(1, 0), Point(0, 1));
  Triangle tri2(Point(0, 0), Point(1, 0), Point(0, 1));
  Triangle tri3(Point(0, 0), Point(1.000001, 0), Point(0, 1.000001));

  BOOST_CHECK(strictCompare(tri1, tri2));
  BOOST_CHECK(sortedCompare(tri1, tri2));
  BOOST_CHECK(fuzzyCompare(tri1, tri2, 1e-6));
  BOOST_CHECK(
      fuzzyCompare(tri1, tri2, 1e-6, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(!strictCompare(tri1, tri3));
  BOOST_CHECK(!sortedCompare(tri1, tri3));
  BOOST_CHECK(fuzzyCompare(tri1, tri3, 1e-5));
  BOOST_CHECK(
      fuzzyCompare(tri1, tri3, 1e-5, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(sortedCompare(tri1, tri3, true, 1e-5, FuzzyCompareMethod::Equal));
  BOOST_CHECK(
      sortedCompare(tri1, tri3, true, 1e-5, FuzzyCompareMethod::DistanceEqual));
}

BOOST_AUTO_TEST_CASE(testCompareMultiPoint)
{
  MultiPoint mp1;
  mp1.addGeometry(new Point(0, 0));
  mp1.addGeometry(new Point(1, 1));

  MultiPoint mp2;
  mp2.addGeometry(new Point(0, 0));
  mp2.addGeometry(new Point(1, 1));

  MultiPoint mp3;
  mp3.addGeometry(new Point(1, 1));
  mp3.addGeometry(new Point(0, 0));

  MultiPoint mp4;
  mp4.addGeometry(new Point(0, 0));
  mp4.addGeometry(new Point(1.000001, 1.000001));

  BOOST_CHECK(strictCompare(mp1, mp2));
  BOOST_CHECK(sortedCompare(mp1, mp2));
  BOOST_CHECK(fuzzyCompare(mp1, mp2, 1e-6));
  BOOST_CHECK(fuzzyCompare(mp1, mp2, 1e-6, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(!strictCompare(mp1, mp3));
  BOOST_CHECK(sortedCompare(mp1, mp3));
  BOOST_CHECK(!fuzzyCompare(mp1, mp3, 1e-6));
  BOOST_CHECK(!fuzzyCompare(mp1, mp3, 1e-6, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(!strictCompare(mp1, mp4));
  BOOST_CHECK(!sortedCompare(mp1, mp4));
  BOOST_CHECK(fuzzyCompare(mp1, mp4, 1e-5));
  BOOST_CHECK(fuzzyCompare(mp1, mp4, 1e-5, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(sortedCompare(mp1, mp4, true, 1e-5, FuzzyCompareMethod::Equal));
  BOOST_CHECK(
      sortedCompare(mp1, mp4, true, 1e-5, FuzzyCompareMethod::DistanceEqual));
}

BOOST_AUTO_TEST_CASE(testCompareMultiLineString)
{
  MultiLineString mls1;
  LineString     *ls1 = new LineString();
  ls1->addPoint(Point(0, 0));
  ls1->addPoint(Point(1, 1));
  mls1.addGeometry(ls1);

  MultiLineString mls2(mls1);

  MultiLineString mls3;
  LineString     *ls3 = new LineString();
  ls3->addPoint(Point(1, 1));
  ls3->addPoint(Point(0, 0));
  mls3.addGeometry(ls3);

  MultiLineString mls4;
  LineString     *ls4 = new LineString();
  ls4->addPoint(Point(0, 0));
  ls4->addPoint(Point(1.000001, 1.000001));
  mls4.addGeometry(ls4);

  BOOST_CHECK(strictCompare(mls1, mls2));
  BOOST_CHECK(sortedCompare(mls1, mls2));
  BOOST_CHECK(fuzzyCompare(mls1, mls2, 1e-6));
  BOOST_CHECK(
      fuzzyCompare(mls1, mls2, 1e-6, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(!strictCompare(mls1, mls3));
  BOOST_CHECK(sortedCompare(mls1, mls3));
  BOOST_CHECK(!fuzzyCompare(mls1, mls3, 1e-6));
  BOOST_CHECK(
      !fuzzyCompare(mls1, mls3, 1e-6, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(!strictCompare(mls1, mls4));
  BOOST_CHECK(!sortedCompare(mls1, mls4));
  BOOST_CHECK(fuzzyCompare(mls1, mls4, 1e-5));
  BOOST_CHECK(
      fuzzyCompare(mls1, mls4, 1e-5, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(sortedCompare(mls1, mls4, true, 1e-5, FuzzyCompareMethod::Equal));
  BOOST_CHECK(
      sortedCompare(mls1, mls4, true, 1e-5, FuzzyCompareMethod::DistanceEqual));
}

BOOST_AUTO_TEST_CASE(testCompareMultiPolygon)
{
  MultiPolygon mp1;
  Polygon     *poly1 = new Polygon();
  poly1->exteriorRing().addPoint(Point(0, 0));
  poly1->exteriorRing().addPoint(Point(1, 0));
  poly1->exteriorRing().addPoint(Point(1, 1));
  poly1->exteriorRing().addPoint(Point(0, 1));
  poly1->exteriorRing().addPoint(Point(0, 0));
  mp1.addGeometry(poly1);

  MultiPolygon mp2(mp1);

  MultiPolygon mp3;
  Polygon     *poly3 = new Polygon();
  poly3->exteriorRing().addPoint(Point(1, 1));
  poly3->exteriorRing().addPoint(Point(0, 1));
  poly3->exteriorRing().addPoint(Point(0, 0));
  poly3->exteriorRing().addPoint(Point(1, 0));
  poly3->exteriorRing().addPoint(Point(1, 1));
  mp3.addGeometry(poly3);

  MultiPolygon mp4;
  Polygon     *poly4 = new Polygon();
  poly4->exteriorRing().addPoint(Point(0, 0));
  poly4->exteriorRing().addPoint(Point(1.000001, 0));
  poly4->exteriorRing().addPoint(Point(1.000001, 1.000001));
  poly4->exteriorRing().addPoint(Point(0, 1.000001));
  poly4->exteriorRing().addPoint(Point(0, 0));
  mp4.addGeometry(poly4);

  BOOST_CHECK(strictCompare(mp1, mp2));
  BOOST_CHECK(sortedCompare(mp1, mp2));
  BOOST_CHECK(fuzzyCompare(mp1, mp2, 1e-6));
  BOOST_CHECK(fuzzyCompare(mp1, mp2, 1e-6, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(!strictCompare(mp1, mp3));
  BOOST_CHECK(sortedCompare(mp1, mp3));
  BOOST_CHECK(!fuzzyCompare(mp1, mp3, 1e-6));
  BOOST_CHECK(!fuzzyCompare(mp1, mp3, 1e-6, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(!strictCompare(mp1, mp4));
  BOOST_CHECK(!sortedCompare(mp1, mp4));
  BOOST_CHECK(fuzzyCompare(mp1, mp4, 1e-5));
  BOOST_CHECK(fuzzyCompare(mp1, mp4, 1e-5, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(sortedCompare(mp1, mp4, true, 1e-5, FuzzyCompareMethod::Equal));
  BOOST_CHECK(
      sortedCompare(mp1, mp4, true, 1e-5, FuzzyCompareMethod::DistanceEqual));
}

BOOST_AUTO_TEST_CASE(testComparePolyhedralSurface)
{
  PolyhedralSurface ps1;
  Polygon          *poly1 = new Polygon();
  poly1->exteriorRing().addPoint(Point(0, 0, 0));
  poly1->exteriorRing().addPoint(Point(1, 0, 0));
  poly1->exteriorRing().addPoint(Point(1, 1, 0));
  poly1->exteriorRing().addPoint(Point(0, 1, 0));
  poly1->exteriorRing().addPoint(Point(0, 0, 0));
  ps1.addPolygon(*poly1);

  PolyhedralSurface ps2(ps1);

  PolyhedralSurface ps3;
  Polygon          *poly3 = new Polygon();
  poly3->exteriorRing().addPoint(Point(1, 1, 0));
  poly3->exteriorRing().addPoint(Point(0, 1, 0));
  poly3->exteriorRing().addPoint(Point(0, 0, 0));
  poly3->exteriorRing().addPoint(Point(1, 0, 0));
  poly3->exteriorRing().addPoint(Point(1, 1, 0));
  ps3.addPolygon(*poly3);

  PolyhedralSurface ps4;
  Polygon          *poly4 = new Polygon();
  poly4->exteriorRing().addPoint(Point(0, 0, 0));
  poly4->exteriorRing().addPoint(Point(1.000001, 0, 0));
  poly4->exteriorRing().addPoint(Point(1.000001, 1.000001, 0));
  poly4->exteriorRing().addPoint(Point(0, 1.000001, 0));
  poly4->exteriorRing().addPoint(Point(0, 0, 0));
  ps4.addPolygon(*poly4);

  BOOST_CHECK(strictCompare(ps1, ps2));
  BOOST_CHECK(sortedCompare(ps1, ps2));
  BOOST_CHECK(fuzzyCompare(ps1, ps2, 1e-6));
  BOOST_CHECK(fuzzyCompare(ps1, ps2, 1e-6, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(!strictCompare(ps1, ps3));
  BOOST_CHECK(sortedCompare(ps1, ps3));
  BOOST_CHECK(!fuzzyCompare(ps1, ps3, 1e-6));
  BOOST_CHECK(!fuzzyCompare(ps1, ps3, 1e-6, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(!strictCompare(ps1, ps4));
  BOOST_CHECK(!sortedCompare(ps1, ps4));
  BOOST_CHECK(fuzzyCompare(ps1, ps4, 1e-5));
  BOOST_CHECK(fuzzyCompare(ps1, ps4, 1e-5, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(sortedCompare(ps1, ps4, true, 1e-5, FuzzyCompareMethod::Equal));
  BOOST_CHECK(
      sortedCompare(ps1, ps4, true, 1e-5, FuzzyCompareMethod::DistanceEqual));

  delete poly1;
  delete poly3;
  delete poly4;
}

BOOST_AUTO_TEST_CASE(testCompareTriangulatedSurface)
{
  TriangulatedSurface ts1;
  ts1.addTriangle(Triangle(Point(0, 0), Point(1, 0), Point(0, 1)));

  TriangulatedSurface ts2(ts1);

  TriangulatedSurface ts3;
  ts3.addTriangle(Triangle(Point(1, 0), Point(0, 1), Point(0, 0)));

  TriangulatedSurface ts4;
  ts4.addTriangle(
      Triangle(Point(0, 0), Point(1.000001, 0), Point(0, 1.000001)));

  BOOST_CHECK(strictCompare(ts1, ts2));
  BOOST_CHECK(sortedCompare(ts1, ts2));
  BOOST_CHECK(fuzzyCompare(ts1, ts2, 1e-6));
  BOOST_CHECK(fuzzyCompare(ts1, ts2, 1e-6, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(!strictCompare(ts1, ts3));
  BOOST_CHECK(sortedCompare(ts1, ts3));
  BOOST_CHECK(!fuzzyCompare(ts1, ts3, 1e-6));
  BOOST_CHECK(!fuzzyCompare(ts1, ts3, 1e-6, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(!strictCompare(ts1, ts4));
  BOOST_CHECK(!sortedCompare(ts1, ts4));
  BOOST_CHECK(fuzzyCompare(ts1, ts4, 1e-5));
  BOOST_CHECK(fuzzyCompare(ts1, ts4, 1e-5, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(sortedCompare(ts1, ts4, true, 1e-5, FuzzyCompareMethod::Equal));
  BOOST_CHECK(
      sortedCompare(ts1, ts4, true, 1e-5, FuzzyCompareMethod::DistanceEqual));
}

BOOST_AUTO_TEST_CASE(testCompareSolid)
{
  PolyhedralSurface shell1;
  Polygon          *poly1 = new Polygon();
  poly1->exteriorRing().addPoint(Point(0, 0, 0));
  poly1->exteriorRing().addPoint(Point(1, 0, 0));
  poly1->exteriorRing().addPoint(Point(1, 1, 0));
  poly1->exteriorRing().addPoint(Point(0, 1, 0));
  poly1->exteriorRing().addPoint(Point(0, 0, 0));
  shell1.addPolygon(*poly1);

  Solid s1(shell1);
  Solid s2(s1);

  PolyhedralSurface shell3;
  Polygon          *poly3 = new Polygon();
  poly3->exteriorRing().addPoint(Point(1, 1, 0));
  poly3->exteriorRing().addPoint(Point(0, 1, 0));
  poly3->exteriorRing().addPoint(Point(0, 0, 0));
  poly3->exteriorRing().addPoint(Point(1, 0, 0));
  poly3->exteriorRing().addPoint(Point(1, 1, 0));
  shell3.addPolygon(*poly3);
  Solid s3(shell3);

  PolyhedralSurface shell4;
  Polygon          *poly4 = new Polygon();
  poly4->exteriorRing().addPoint(Point(0, 0, 0));
  poly4->exteriorRing().addPoint(Point(1.000001, 0, 0));
  poly4->exteriorRing().addPoint(Point(1.000001, 1.000001, 0));
  poly4->exteriorRing().addPoint(Point(0, 1.000001, 0));
  poly4->exteriorRing().addPoint(Point(0, 0, 0));
  shell4.addPolygon(*poly4);
  Solid s4(shell4);

  BOOST_CHECK(strictCompare(s1, s2));
  BOOST_CHECK(sortedCompare(s1, s2));
  BOOST_CHECK(fuzzyCompare(s1, s2, 1e-6));
  BOOST_CHECK(fuzzyCompare(s1, s2, 1e-6, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(!strictCompare(s1, s3));
  BOOST_CHECK(sortedCompare(s1, s3));
  BOOST_CHECK(!fuzzyCompare(s1, s3, 1e-6));
  BOOST_CHECK(!fuzzyCompare(s1, s3, 1e-6, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(!strictCompare(s1, s4));
  BOOST_CHECK(!sortedCompare(s1, s4));
  BOOST_CHECK(fuzzyCompare(s1, s4, 1e-5));
  BOOST_CHECK(fuzzyCompare(s1, s4, 1e-5, FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(sortedCompare(s1, s4, true, 1e-5, FuzzyCompareMethod::Equal));
  BOOST_CHECK(
      sortedCompare(s1, s4, true, 1e-5, FuzzyCompareMethod::DistanceEqual));

  delete poly1;
  delete poly3;
  delete poly4;
}

BOOST_AUTO_TEST_CASE(testCompareEmptyGeometries)
{
  Point               emptyPoint;
  LineString          emptyLineString;
  Polygon             emptyPolygon;
  Triangle            emptyTriangle;
  MultiPoint          emptyMultiPoint;
  MultiLineString     emptyMultiLineString;
  MultiPolygon        emptyMultiPolygon;
  GeometryCollection  emptyGeometryCollection;
  PolyhedralSurface   emptyPolyhedralSurface;
  TriangulatedSurface emptyTriangulatedSurface;
  Solid               emptySolid;

  BOOST_CHECK(strictCompare(emptyPoint, emptyPoint));
  BOOST_CHECK(sortedCompare(emptyLineString, emptyLineString));
  BOOST_CHECK(fuzzyCompare(emptyPolygon, emptyPolygon, 1e-6));
  BOOST_CHECK(fuzzyCompare(emptyTriangle, emptyTriangle, 1e-6,
                           FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(!strictCompare(emptyPoint, emptyLineString));
  BOOST_CHECK(!sortedCompare(emptyLineString, emptyPolygon));
  BOOST_CHECK(!fuzzyCompare(emptyPolygon, emptyMultiPoint, 1e-6));
  BOOST_CHECK(!fuzzyCompare(emptyMultiLineString, emptyMultiPolygon, 1e-6,
                            FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(strictCompare(emptyGeometryCollection, emptyGeometryCollection));
  BOOST_CHECK(sortedCompare(emptyPolyhedralSurface, emptyPolyhedralSurface));
  BOOST_CHECK(
      fuzzyCompare(emptyTriangulatedSurface, emptyTriangulatedSurface, 1e-6));
  BOOST_CHECK(fuzzyCompare(emptySolid, emptySolid, 1e-6,
                           FuzzyCompareMethod::DistanceEqual));
}

BOOST_AUTO_TEST_CASE(testCompareDifferentGeometryTypes)
{
  Point      point(0, 0);
  LineString lineString;
  lineString.addPoint(Point(0, 0));
  lineString.addPoint(Point(1, 1));
  Polygon polygon;
  polygon.exteriorRing().addPoint(Point(0, 0));
  polygon.exteriorRing().addPoint(Point(1, 0));
  polygon.exteriorRing().addPoint(Point(1, 1));
  polygon.exteriorRing().addPoint(Point(0, 1));
  polygon.exteriorRing().addPoint(Point(0, 0));

  BOOST_CHECK(!strictCompare(point, lineString));
  BOOST_CHECK(!sortedCompare(point, lineString));
  BOOST_CHECK(!fuzzyCompare(point, lineString, 1e-6));
  BOOST_CHECK(!fuzzyCompare(point, lineString, 1e-6,
                            FuzzyCompareMethod::DistanceEqual));

  BOOST_CHECK(!strictCompare(lineString, polygon));
  BOOST_CHECK(!sortedCompare(lineString, polygon));
  BOOST_CHECK(!fuzzyCompare(lineString, polygon, 1e-6));
  BOOST_CHECK(!fuzzyCompare(lineString, polygon, 1e-6,
                            FuzzyCompareMethod::DistanceEqual));
}

BOOST_AUTO_TEST_SUITE_END()
