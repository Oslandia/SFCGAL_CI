// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <boost/test/unit_test.hpp>

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
#include "SFCGAL/algorithm/covers.h"
#include "SFCGAL/algorithm/length.h"
#include "SFCGAL/algorithm/straightSkeleton.h"
#include "SFCGAL/io/wkt.h"

#include "SFCGAL/detail/tools/Registry.h"

using namespace boost::unit_test;
using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_StraightSkeletonTest)

BOOST_AUTO_TEST_CASE(testTriangle)
{
  std::unique_ptr<Geometry> g(io::readWkt("TRIANGLE ((1 1,2 1,2 2,1 1))"));

  std::string const expectedWKT("MULTILINESTRING ((1.0 1.0,1.7 1.3),(2.0 "
                                "1.0,1.7 1.3),(2.0 2.0,1.7 1.3))");
  {
    std::unique_ptr<MultiLineString> result(algorithm::straightSkeleton(*g));
    BOOST_CHECK_EQUAL(result->numGeometries(), 3U);
    BOOST_CHECK_EQUAL(result->asText(1), expectedWKT);
  }
  // check orientation insensitive
  {
    g->as<Triangle>().reverse();
    std::unique_ptr<MultiLineString> result(algorithm::straightSkeleton(*g));
    BOOST_CHECK_EQUAL(result->numGeometries(), 3U);
    BOOST_CHECK_EQUAL(result->asText(1), expectedWKT);
  }
}

BOOST_AUTO_TEST_CASE(testPolygon)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("POLYGON ((1 1,11 1,11 11,1 11,1 1))"));

  std::string const expectedWKT(
      "MULTILINESTRING ((1 1,6 6),(11 1,6 6),(11 11,6 6),(1 11,6 6))");
  {
    std::unique_ptr<MultiLineString> result(algorithm::straightSkeleton(*g));
    BOOST_CHECK_EQUAL(result->numGeometries(), 4U);
    BOOST_CHECK_EQUAL(result->asText(0), expectedWKT);
  }
  // check orientation insensitive
  {
    g->as<Polygon>().exteriorRing().reverse();
    std::unique_ptr<MultiLineString> result(algorithm::straightSkeleton(*g));
    BOOST_CHECK_EQUAL(result->numGeometries(), 4U);
    BOOST_CHECK_EQUAL(result->asText(0), expectedWKT);
  }
}

BOOST_AUTO_TEST_CASE(testPolygonWithHole)
{
  std::unique_ptr<Geometry> const g(
      io::readWkt("POLYGON ( (-1.0 -1.0,1.0 -1.0,1.0 1.0,-1.0 1.0,-1.0 -1.0)"
                  ", (-0.5 -0.5,-0.5 0.5,0.5 0.5,-0.5 -0.5)"
                  ")"));

  std::unique_ptr<MultiLineString> result(algorithm::straightSkeleton(*g));
  BOOST_CHECK_EQUAL(result->numGeometries(), 13);

  std::unique_ptr<Geometry> expected(
      io::readWkt("MULTILINESTRING ( (-1/1 -1/1,-3/4 -3/4)"
                  ", (1/1 -1/1,1865452045155277/4503599627370496 "
                  "-3730904090310553/9007199254740992)"
                  ", (1/1 1/1,3/4 3/4)"
                  ", (-1/1 1/1,-3/4 3/4)"
                  ", (-1/2 -1/2,-5822673418478105/9007199254740992 "
                  "-7688125463633383/9007199254740992)"
                  ", (-1/2 1/2,-3/4 3/4)"
                  ", (1/2 1/2,3844062731816691/4503599627370496 "
                  "727834177309763/1125899906842624)"
                  ", (-5822673418478105/9007199254740992 "
                  "-7688125463633383/9007199254740992,1865452045155277/"
                  "4503599627370496 -3730904090310553/9007199254740992)"
                  ", (-5822673418478105/9007199254740992 "
                  "-7688125463633383/9007199254740992,-3/4 -3/4)"
                  ", (3844062731816691/4503599627370496 "
                  "727834177309763/1125899906842624,3/4 3/4)"
                  ", (3844062731816691/4503599627370496 "
                  "727834177309763/1125899906842624,1865452045155277/"
                  "4503599627370496 -3730904090310553/9007199254740992)"
                  ", (-3/4 3/4,3/4 3/4)"
                  ", (-3/4 3/4,-3/4 -3/4)"
                  ")"));
  // results for gcc and clang differ due to rounding
  // To avoid rounding errors in the results (-3730904090310553/9007199254740992
  // vs -466363011288819/1125899906842624), a text comparison is used. This is
  // not optimal.
  std::unique_ptr<Geometry> const r(io::readWkt(result->asText(10)));
  std::unique_ptr<Geometry> const e(io::readWkt(expected->asText(10)));
  BOOST_CHECK(algorithm::covers(*r, *e));
}

BOOST_AUTO_TEST_CASE(testPolygonWithHoleTouchingShell)
{
  std::unique_ptr<Geometry> const g(
      io::readWkt("POLYGON ((-1.0 -1.0,1.0 -1.0,1.0 1.0,-1.0 1.0,-1.0 "
                  "-1.0),(-0.5 -0.5,-0.5 0.5,0.5 0.5,1.0 -0.5,-0.5 -0.5))"));
  BOOST_CHECK_THROW(algorithm::straightSkeleton(*g), NotImplementedException);
}

BOOST_AUTO_TEST_CASE(testPolygonWithTouchingHoles)
{
  std::unique_ptr<Geometry> const g(
      io::readWkt("POLYGON ((-1.0 -1.0,1.0 -1.0,1.0 1.0,-1.0 1.0,-1.0 "
                  "-1.0),(-0.5 -0.5,-0.5 0.5,-0.1 0.5,0.1 -0.5,-0.5 -0.5),(0.1 "
                  "-0.5,0.1 0.5,0.5 0.5,0.5 -0.5,0.1 -0.5))"));

  BOOST_CHECK_THROW(algorithm::straightSkeleton(*g), NotImplementedException);
}

BOOST_AUTO_TEST_CASE(testMultiPolygon)
{
  std::unique_ptr<Geometry> const  g(io::readWkt(
      "MULTIPOLYGON (((3.000000 0.000000,2.875000 0.484123,2.750000 "
       "0.661438,2.625000 0.780625,2.500000 0.866025,2.375000 0.927025,2.250000 "
       "0.968246,2.125000 0.992157,2.000000 1.000000,1.875000 1.484123,1.750000 "
       "1.661438,1.625000 1.780625,1.500000 1.866025,1.375000 1.927025,1.250000 "
       "1.968246,1.125000 1.992157,1.000000 2.000000,0.750000 2.661438,0.500000 "
       "2.866025,0.250000 2.968246,0.000000 3.000000,-0.250000 "
       "2.968246,-0.500000 2.866025,-0.750000 2.661438,-1.000000 "
       "2.000000,-1.125000 1.992157,-1.250000 1.968246,-1.375000 "
       "1.927025,-1.500000 1.866025,-1.625000 1.780625,-1.750000 "
       "1.661438,-1.875000 1.484123,-2.000000 1.000000,-2.125000 "
       "0.992157,-2.250000 0.968246,-2.375000 0.927025,-2.500000 "
       "0.866025,-2.625000 0.780625,-2.750000 0.661438,-2.875000 "
       "0.484123,-3.000000 0.000000,-2.875000 -0.484123,-2.750000 "
       "-0.661438,-2.625000 -0.780625,-2.500000 -0.866025,-2.375000 "
       "-0.927025,-2.250000 -0.968246,-2.125000 -0.992157,-2.000000 "
       "-1.000000,-1.875000 -1.484123,-1.750000 -1.661438,-1.625000 "
       "-1.780625,-1.500000 -1.866025,-1.375000 -1.927025,-1.250000 "
       "-1.968246,-1.125000 -1.992157,-1.000000 -2.000000,-0.750000 "
       "-2.661438,-0.500000 -2.866025,-0.250000 -2.968246,0.000000 "
       "-3.000000,0.250000 -2.968246,0.500000 -2.866025,0.750000 "
       "-2.661438,1.000000 -2.000000,1.125000 -1.992157,1.250000 "
       "-1.968246,1.375000 -1.927025,1.500000 -1.866025,1.625000 "
       "-1.780625,1.750000 -1.661438,1.875000 -1.484123,2.000000 "
       "-1.000000,2.125000 -0.992157,2.250000 -0.968246,2.375000 "
       "-0.927025,2.500000 -0.866025,2.625000 -0.780625,2.750000 "
       "-0.661438,2.875000 -0.484123,3.000000 0.000000),(0.000000 "
       "1.000000,0.125000 0.515877,0.250000 0.338562,0.375000 0.219375,0.500000 "
       "0.133975,0.625000 0.072975,0.750000 0.031754,0.875000 0.007843,1.000000 "
       "0.000000,0.875000 -0.007843,0.750000 -0.031754,0.625000 "
       "-0.072975,0.500000 -0.133975,0.375000 -0.219375,0.250000 "
       "-0.338562,0.125000 -0.515877,0.000000 -1.000000,-0.125000 "
       "-0.515877,-0.250000 -0.338562,-0.375000 -0.219375,-0.500000 "
       "-0.133975,-0.625000 -0.072975,-0.750000 -0.031754,-0.875000 "
       "-0.007843,-1.000000 0.000000,-0.875000 0.007843,-0.750000 "
       "0.031754,-0.625000 0.072975,-0.500000 0.133975,-0.375000 "
       "0.219375,-0.250000 0.338562,-0.125000 0.515877,0.000000 1.000000)))"));
  std::unique_ptr<MultiLineString> result(algorithm::straightSkeleton(*g));
  BOOST_CHECK_EQUAL(result->numGeometries(), 220U);
}

BOOST_AUTO_TEST_CASE(testInvalidTypes)
{
  std::vector<std::string> wkt;
  wkt.emplace_back("POINT (1 2)");
  wkt.emplace_back("LINESTRING (0 0,1 1)");

  for (auto &it : wkt) {
    std::unique_ptr<Geometry> const  g(io::readWkt(it));
    std::unique_ptr<MultiLineString> result(algorithm::straightSkeleton(*g));
    BOOST_CHECK_EQUAL(result->numGeometries(), 0U);
  }
}

// See https://github.com/Oslandia/SFCGAL/issues/75
BOOST_AUTO_TEST_CASE(testPostgisIssue3107)
{
  std::unique_ptr<Geometry> const g(io::readWkt(
      "POLYGON ((1347259.25 7184745.94,1347273.17 7184758.16,1347280.39 "
      "7184749.95,1347278.04 7184747.88,1347281.66 7184743.76,1347284.01 "
      "7184745.83,1347293.5 7184735.05,1347279.61 7184722.85,1347269.29 "
      "7184734.6,1347270.35 7184735.51,1347267.31 7184738.96,1347266.22 "
      "7184738.01,1347259.25 7184745.94),(1347267.31 7184738.96,1347269.31 "
      "7184736.7,1347272.57 7184739.55,1347270.56 7184741.83,1347267.31 "
      "7184738.96))"));
  BOOST_CHECK_THROW(algorithm::straightSkeleton(*g), NotImplementedException);
}

// See https://github.com/Oslandia/SFCGAL/issues/91
BOOST_AUTO_TEST_CASE(testMultiPolygonWithTouchingHoles)
{
  std::unique_ptr<Geometry> const g(io::readWkt(
      "MULTIPOLYGON (((1347259.25 7184745.94,1347273.17 7184758.16,1347280.39 "
      "7184749.95,1347278.04 7184747.88,1347281.66 7184743.76,1347284.01 "
      "7184745.83,1347293.5 7184735.05,1347279.61 7184722.85,1347269.29 "
      "7184734.6,1347270.35 7184735.51,1347267.31 7184738.96,1347266.22 "
      "7184738.01,1347259.25 7184745.94),(1347267.31 7184738.96,1347269.31 "
      "7184736.7,1347272.57 7184739.55,1347270.56 7184741.83,1347267.31 "
      "7184738.96)))"));
  BOOST_CHECK_THROW(algorithm::straightSkeleton(*g), NotImplementedException);
}

BOOST_AUTO_TEST_CASE(testDistanceInM)
{
  std::unique_ptr<Geometry> const g(
      io::readWkt("POLYGON ((0 0,1 0,1 1,0 1,0 0))"));
  std::unique_ptr<Geometry> out(algorithm::straightSkeleton(
      *g, /* autoOrientation */ true, /* innerOnly */ false,
      /* outputDistanceInM */ true));
  std::string const         expectedWKT(
      "MULTILINESTRING M ((0.0 0.0 0.0,0.5 0.5 0.5),(1.0 0.0 0.0,0.5 0.5 "
              "0.5),(1.0 1.0 0.0,0.5 0.5 0.5),(0.0 1.0 0.0,0.5 0.5 0.5))");
  BOOST_CHECK_EQUAL(out->asText(1), expectedWKT);
}

BOOST_AUTO_TEST_CASE(testMultiEmptyEmpty)
{
  std::unique_ptr<Geometry> const g(io::readWkt("MULTIPOLYGON (EMPTY,EMPTY)"));
  std::unique_ptr<Geometry>       out(algorithm::straightSkeleton(*g));
  std::string const               expectedWKT("MULTILINESTRING EMPTY");
  BOOST_CHECK_EQUAL(out->asText(1), expectedWKT);
}

// See https://gitlab.com/Oslandia/SFCGAL/-/issues/194
BOOST_AUTO_TEST_CASE(testDegenerateMultiLineString)
{
  std::unique_ptr<Geometry> const g(io::readWkt(
      "Polygon ((1294585.78643762995488942 200985.78643762698629871, 1294000 "
      "202400, 1294000 212400, 1294585.78643762995488942 "
      "213814.21356237301370129, 1296000 214400, 1297000 214400, "
      "1298414.21356237004511058 213814.21356237301370129, 1299000 212400, "
      "1299000 202400, 1298414.21356237004511058 200985.78643762698629871, "
      "1297000 200400, 1296000 200400, 1294585.78643762995488942 "
      "200985.78643762698629871),(1297000 202400, 1297000 212400, 1296000 "
      "212400, 1296000 202400, 1297000 202400))"));
  const double                    tolerance = EPSILON;
  std::unique_ptr<Geometry>       out(
      algorithm::straightSkeleton(*g, tolerance != 0.0));
  for (size_t i = 0; i < out->numGeometries(); i++) {
    BOOST_CHECK(algorithm::length(out->geometryN(i)) > tolerance);
  }
}

BOOST_AUTO_TEST_CASE(testEmptyExtrudeStraightSkeleton)
{
  std::unique_ptr<Geometry> const    geom(io::readWkt("POLYGON EMPTY"));
  std::unique_ptr<PolyhedralSurface> out(
      algorithm::extrudeStraightSkeleton(*geom, 2.0));
  std::string const expectedWKT("POLYHEDRALSURFACE EMPTY");
  BOOST_CHECK_EQUAL(out->asText(2), expectedWKT);
}

BOOST_AUTO_TEST_CASE(testEmptyExtrudeStraightSkeletonRoofTop)
{
  std::unique_ptr<Geometry> const    geom(io::readWkt("POLYGON EMPTY"));
  std::unique_ptr<PolyhedralSurface> out(
      algorithm::extrudeStraightSkeleton(*geom, 2.0, 10.0));
  std::string const expectedWKT("POLYHEDRALSURFACE EMPTY");
  BOOST_CHECK_EQUAL(out->asText(2), expectedWKT);
}

BOOST_AUTO_TEST_CASE(testExtrudeStraightSkeleton)
{

  std::unique_ptr<Geometry> const geom(
      io::readWkt("POLYGON (( 0 0, 5 0, 5 5, 4 5, 4 4, 0 4, 0 0 ))"));
  std::unique_ptr<PolyhedralSurface> out(
      algorithm::extrudeStraightSkeleton(*geom, 2.0));
  std::string const expectedWKT(
      "POLYHEDRALSURFACE Z (((4.00 5.00 0.00,5.00 5.00 0.00,4.00 4.00 "
      "0.00,4.00 5.00 0.00)),((0.00 4.00 0.00,4.00 4.00 0.00,0.00 0.00 "
      "0.00,0.00 4.00 0.00)),((4.00 4.00 0.00,5.00 0.00 0.00,0.00 0.00 "
      "0.00,4.00 4.00 0.00)),((5.00 5.00 0.00,5.00 0.00 0.00,4.00 4.00 "
      "0.00,5.00 5.00 0.00)),((0.00 4.00 0.00,0.00 0.00 0.00,2.00 2.00 "
      "2.00,0.00 4.00 0.00)),((0.00 0.00 0.00,5.00 0.00 0.00,3.00 2.00 "
      "2.00,0.00 0.00 0.00)),((2.00 2.00 2.00,0.00 0.00 0.00,3.00 2.00 "
      "2.00,2.00 2.00 2.00)),((4.50 3.50 0.50,5.00 5.00 0.00,4.50 4.50 "
      "0.50,4.50 3.50 0.50)),((3.00 2.00 2.00,5.00 0.00 0.00,4.50 3.50 "
      "0.50,3.00 2.00 2.00)),((4.50 3.50 0.50,5.00 0.00 0.00,5.00 5.00 "
      "0.00,4.50 3.50 0.50)),((5.00 5.00 0.00,4.00 5.00 0.00,4.50 4.50 "
      "0.50,5.00 5.00 0.00)),((4.50 4.50 0.50,4.00 4.00 0.00,4.50 3.50 "
      "0.50,4.50 4.50 0.50)),((4.50 4.50 0.50,4.00 5.00 0.00,4.00 4.00 "
      "0.00,4.50 4.50 0.50)),((4.00 4.00 0.00,0.00 4.00 0.00,2.00 2.00 "
      "2.00,4.00 4.00 0.00)),((4.50 3.50 0.50,4.00 4.00 0.00,3.00 2.00 "
      "2.00,4.50 3.50 0.50)),((3.00 2.00 2.00,4.00 4.00 0.00,2.00 2.00 "
      "2.00,3.00 2.00 2.00)))");
  BOOST_CHECK_EQUAL(out->asText(2), expectedWKT);
}

BOOST_AUTO_TEST_CASE(testExtrudeStraightSkeletonPolygonWithHole)
{

  std::unique_ptr<Geometry> const g(
      io::readWkt("POLYGON (( 0 0, 5 0, 5 5, 4 5, 4 4, 0 4, 0 0 ), (1 1, 1 2, "
                  "2 2, 2 1, 1 1))"));
  std::unique_ptr<PolyhedralSurface> out(
      algorithm::extrudeStraightSkeleton(*g, 2.0));
  std::string const expectedWKT(
      "POLYHEDRALSURFACE Z (((4.00 5.00 0.00,5.00 5.00 0.00,4.00 4.00 "
      "0.00,4.00 "
      "5.00 0.00)),((2.00 1.00 0.00,5.00 0.00 0.00,0.00 0.00 0.00,2.00 1.00 "
      "0.00)),((5.00 5.00 0.00,5.00 0.00 0.00,4.00 4.00 0.00,5.00 5.00 "
      "0.00)),((2.00 1.00 0.00,0.00 0.00 0.00,1.00 1.00 0.00,2.00 1.00 "
      "0.00)),((1.00 2.00 0.00,1.00 1.00 0.00,0.00 0.00 0.00,1.00 2.00 "
      "0.00)),((0.00 4.00 0.00,2.00 2.00 0.00,1.00 2.00 0.00,0.00 4.00 "
      "0.00)),((0.00 4.00 0.00,1.00 2.00 0.00,0.00 0.00 0.00,0.00 4.00 "
      "0.00)),((4.00 4.00 0.00,5.00 0.00 0.00,2.00 2.00 0.00,4.00 4.00 "
      "0.00)),((4.00 4.00 0.00,2.00 2.00 0.00,0.00 4.00 0.00,4.00 4.00 "
      "0.00)),((2.00 2.00 0.00,5.00 0.00 0.00,2.00 1.00 0.00,2.00 2.00 "
      "0.00)),((0.50 2.50 0.50,0.00 0.00 0.00,0.50 0.50 0.50,0.50 2.50 "
      "0.50)),((1.00 3.00 1.00,0.00 4.00 0.00,0.50 2.50 0.50,1.00 3.00 "
      "1.00)),((0.50 2.50 0.50,0.00 4.00 0.00,0.00 0.00 0.00,0.50 2.50 "
      "0.50)),((2.50 0.50 0.50,5.00 0.00 0.00,3.50 1.50 1.50,2.50 0.50 "
      "0.50)),((0.00 0.00 0.00,5.00 0.00 0.00,2.50 0.50 0.50,0.00 0.00 "
      "0.00)),((0.50 0.50 0.50,0.00 0.00 0.00,2.50 0.50 0.50,0.50 0.50 "
      "0.50)),((4.50 3.50 0.50,5.00 5.00 0.00,4.50 4.50 0.50,4.50 3.50 "
      "0.50)),((3.50 2.50 1.50,3.50 1.50 1.50,4.50 3.50 0.50,3.50 2.50 "
      "1.50)),((4.50 3.50 0.50,5.00 0.00 0.00,5.00 5.00 0.00,4.50 3.50 "
      "0.50)),((3.50 1.50 1.50,5.00 0.00 0.00,4.50 3.50 0.50,3.50 1.50 "
      "1.50)),((5.00 5.00 0.00,4.00 5.00 0.00,4.50 4.50 0.50,5.00 5.00 "
      "0.00)),((4.50 4.50 0.50,4.00 4.00 0.00,4.50 3.50 0.50,4.50 4.50 "
      "0.50)),((4.50 4.50 0.50,4.00 5.00 0.00,4.00 4.00 0.00,4.50 4.50 "
      "0.50)),((3.00 3.00 1.00,0.00 4.00 0.00,1.00 3.00 1.00,3.00 3.00 "
      "1.00)),((3.50 2.50 1.50,4.50 3.50 0.50,3.00 3.00 1.00,3.50 2.50 "
      "1.50)),((3.00 3.00 1.00,4.00 4.00 0.00,0.00 4.00 0.00,3.00 3.00 "
      "1.00)),((4.50 3.50 0.50,4.00 4.00 0.00,3.00 3.00 1.00,4.50 3.50 "
      "0.50)),((2.00 1.00 0.00,1.00 1.00 0.00,0.50 0.50 0.50,2.00 1.00 "
      "0.00)),((2.50 0.50 0.50,2.00 1.00 0.00,0.50 0.50 0.50,2.50 0.50 "
      "0.50)),((1.00 1.00 0.00,1.00 2.00 0.00,0.50 2.50 0.50,1.00 1.00 "
      "0.00)),((0.50 0.50 0.50,1.00 1.00 0.00,0.50 2.50 0.50,0.50 0.50 "
      "0.50)),((1.00 3.00 1.00,2.00 2.00 0.00,3.00 3.00 1.00,1.00 3.00 "
      "1.00)),((0.50 2.50 0.50,1.00 2.00 0.00,1.00 3.00 1.00,0.50 2.50 "
      "0.50)),((1.00 3.00 1.00,1.00 2.00 0.00,2.00 2.00 0.00,1.00 3.00 "
      "1.00)),((2.00 2.00 0.00,2.00 1.00 0.00,2.50 0.50 0.50,2.00 2.00 "
      "0.00)),((3.50 2.50 1.50,3.00 3.00 1.00,3.50 1.50 1.50,3.50 2.50 "
      "1.50)),((3.50 1.50 1.50,2.00 2.00 0.00,2.50 0.50 0.50,3.50 1.50 "
      "1.50)),((3.00 3.00 1.00,2.00 2.00 0.00,3.50 1.50 1.50,3.00 3.00 "
      "1.00)))");
  BOOST_CHECK_EQUAL(out->asText(2), expectedWKT);
}

BOOST_AUTO_TEST_CASE(testExtrudeStraightSkeletonGenerateBuilding)
{

  std::unique_ptr<Geometry> const geom(
      io::readWkt("POLYGON (( 0 0, 5 0, 5 5, 4 5, 4 4, 0 4, 0 0 ), (1 1, 1 2, "
                  "2 2, 2 1, 1 1))"));
  std::unique_ptr<Geometry> out(
      algorithm::extrudeStraightSkeleton(*geom, 9.0, 2.0));
  std::string const expectedWKT(
      "POLYHEDRALSURFACE Z (((0.00 0.00 0.00,0.00 4.00 0.00,4.00 4.00 "
      "0.00,4.00 5.00 0.00,5.00 5.00 0.00,5.00 0.00 0.00,0.00 0.00 0.00),(1.00 "
      "1.00 0.00,2.00 1.00 0.00,2.00 2.00 0.00,1.00 2.00 0.00,1.00 1.00 "
      "0.00)),((0.00 0.00 0.00,0.00 0.00 9.00,0.00 4.00 9.00,0.00 4.00 "
      "0.00,0.00 0.00 0.00)),((0.00 4.00 0.00,0.00 4.00 9.00,4.00 4.00 "
      "9.00,4.00 4.00 0.00,0.00 4.00 0.00)),((4.00 4.00 0.00,4.00 4.00 "
      "9.00,4.00 5.00 9.00,4.00 5.00 0.00,4.00 4.00 0.00)),((4.00 5.00 "
      "0.00,4.00 5.00 9.00,5.00 5.00 9.00,5.00 5.00 0.00,4.00 5.00 "
      "0.00)),((5.00 5.00 0.00,5.00 5.00 9.00,5.00 0.00 9.00,5.00 0.00 "
      "0.00,5.00 5.00 0.00)),((5.00 0.00 0.00,5.00 0.00 9.00,0.00 0.00 "
      "9.00,0.00 0.00 0.00,5.00 0.00 0.00)),((1.00 1.00 0.00,1.00 1.00 "
      "9.00,2.00 1.00 9.00,2.00 1.00 0.00,1.00 1.00 0.00)),((2.00 1.00 "
      "0.00,2.00 1.00 9.00,2.00 2.00 9.00,2.00 2.00 0.00,2.00 1.00 "
      "0.00)),((2.00 2.00 0.00,2.00 2.00 9.00,1.00 2.00 9.00,1.00 2.00 "
      "0.00,2.00 2.00 0.00)),((1.00 2.00 0.00,1.00 2.00 9.00,1.00 1.00 "
      "9.00,1.00 1.00 0.00,1.00 2.00 0.00)),((0.50 2.50 9.50,0.00 0.00 "
      "9.00,0.50 0.50 9.50,0.50 2.50 9.50)),((1.00 3.00 10.00,0.00 4.00 "
      "9.00,0.50 2.50 9.50,1.00 3.00 10.00)),((0.50 2.50 9.50,0.00 4.00 "
      "9.00,0.00 0.00 9.00,0.50 2.50 9.50)),((2.50 0.50 9.50,5.00 0.00 "
      "9.00,3.50 1.50 10.50,2.50 0.50 9.50)),((0.00 0.00 9.00,5.00 0.00 "
      "9.00,2.50 0.50 9.50,0.00 0.00 9.00)),((0.50 0.50 9.50,0.00 0.00 "
      "9.00,2.50 0.50 9.50,0.50 0.50 9.50)),((4.50 3.50 9.50,5.00 5.00 "
      "9.00,4.50 4.50 9.50,4.50 3.50 9.50)),((3.50 2.50 10.50,3.50 1.50 "
      "10.50,4.50 3.50 9.50,3.50 2.50 10.50)),((4.50 3.50 9.50,5.00 0.00 "
      "9.00,5.00 5.00 9.00,4.50 3.50 9.50)),((3.50 1.50 10.50,5.00 0.00 "
      "9.00,4.50 3.50 9.50,3.50 1.50 10.50)),((5.00 5.00 9.00,4.00 5.00 "
      "9.00,4.50 4.50 9.50,5.00 5.00 9.00)),((4.50 4.50 9.50,4.00 4.00 "
      "9.00,4.50 3.50 9.50,4.50 4.50 9.50)),((4.50 4.50 9.50,4.00 5.00 "
      "9.00,4.00 4.00 9.00,4.50 4.50 9.50)),((3.00 3.00 10.00,0.00 4.00 "
      "9.00,1.00 3.00 10.00,3.00 3.00 10.00)),((3.50 2.50 10.50,4.50 3.50 "
      "9.50,3.00 3.00 10.00,3.50 2.50 10.50)),((3.00 3.00 10.00,4.00 4.00 "
      "9.00,0.00 4.00 9.00,3.00 3.00 10.00)),((4.50 3.50 9.50,4.00 4.00 "
      "9.00,3.00 3.00 10.00,4.50 3.50 9.50)),((2.00 1.00 9.00,1.00 1.00 "
      "9.00,0.50 0.50 9.50,2.00 1.00 9.00)),((2.50 0.50 9.50,2.00 1.00 "
      "9.00,0.50 0.50 9.50,2.50 0.50 9.50)),((1.00 1.00 9.00,1.00 2.00 "
      "9.00,0.50 2.50 9.50,1.00 1.00 9.00)),((0.50 0.50 9.50,1.00 1.00 "
      "9.00,0.50 2.50 9.50,0.50 0.50 9.50)),((1.00 3.00 10.00,2.00 2.00 "
      "9.00,3.00 3.00 10.00,1.00 3.00 10.00)),((0.50 2.50 9.50,1.00 2.00 "
      "9.00,1.00 3.00 10.00,0.50 2.50 9.50)),((1.00 3.00 10.00,1.00 2.00 "
      "9.00,2.00 2.00 9.00,1.00 3.00 10.00)),((2.00 2.00 9.00,2.00 1.00 "
      "9.00,2.50 0.50 9.50,2.00 2.00 9.00)),((3.50 2.50 10.50,3.00 3.00 "
      "10.00,3.50 1.50 10.50,3.50 2.50 10.50)),((3.50 1.50 10.50,2.00 2.00 "
      "9.00,2.50 0.50 9.50,3.50 1.50 10.50)),((3.00 3.00 10.00,2.00 2.00 "
      "9.00,3.50 1.50 10.50,3.00 3.00 10.00)))");
  BOOST_CHECK_EQUAL(out->asText(2), expectedWKT);
}

BOOST_AUTO_TEST_CASE(testStraightSkeletonPartitionLShapedPolygon)
{

  std::unique_ptr<Geometry> const g(
      io::readWkt("POLYGON ((0 0, 0 2, 1 2, 1 1, 2 1, 2 0, 0 0))"));
  std::unique_ptr<Geometry> out(algorithm::straightSkeletonPartition(*g));
  std::string const         expectedWKT(
      "POLYHEDRALSURFACE (((0.00 0.00,0.50 0.50,0.50 1.50,0.00 2.00)),((2.00 "
              "0.00,1.50 0.50,0.50 0.50,0.00 0.00)),((2.00 1.00,1.50 0.50,2.00 "
              "0.00)),((1.00 1.00,0.50 0.50,1.50 0.50,2.00 1.00)),((1.00 2.00,0.50 "
              "1.50,0.50 0.50,1.00 1.00)),((0.00 2.00,0.50 1.50,1.00 2.00)))");
  BOOST_CHECK_EQUAL(out->asText(2), expectedWKT);
}

BOOST_AUTO_TEST_CASE(testStraightSkeletonPartitionSimpleRectangle)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("POLYGON ((0 0, 0 2, 3 2, 3 0, 0 0))"));
  std::unique_ptr<Geometry> out(algorithm::straightSkeletonPartition(*g));
  BOOST_CHECK(out->is<PolyhedralSurface>());
  BOOST_CHECK_EQUAL(out->as<PolyhedralSurface>().numGeometries(), 1);
  BOOST_CHECK_EQUAL(out->as<PolyhedralSurface>().numPolygons(), 4);
  std::string const expectedWKT(
      "POLYHEDRALSURFACE (((0.00 0.00,1.00 1.00,0.00 2.00)),((3.00 0.00,2.00 "
      "1.00,1.00 1.00,0.00 0.00)),((3.00 2.00,2.00 1.00,3.00 0.00)),((0.00 "
      "2.00,1.00 1.00,2.00 1.00,3.00 2.00)))");
  BOOST_CHECK_EQUAL(out->asText(2), expectedWKT);
}

BOOST_AUTO_TEST_CASE(testStraightSkeletonPartitionComplexPolygon)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("POLYGON ((0 0, 0 3, 1 3, 1 2, 2 2, 2 3, 3 3, 3 0, 0 0))"));
  std::unique_ptr<Geometry> out(algorithm::straightSkeletonPartition(*g));
  BOOST_CHECK(out->is<PolyhedralSurface>());
  BOOST_CHECK_EQUAL(out->as<PolyhedralSurface>().numGeometries(), 1);
  BOOST_CHECK_EQUAL(out->as<PolyhedralSurface>().numPolygons(), 8);
  std::string const expectedWKT(
      "POLYHEDRALSURFACE (((0.00 0.00,1.00 1.00,0.50 1.50,0.50 2.50,0.00 "
      "3.00)),((3.00 0.00,2.00 1.00,1.00 1.00,0.00 0.00)),((3.00 3.00,2.50 "
      "2.50,2.50 1.50,2.00 1.00,3.00 0.00)),((2.00 3.00,2.50 2.50,3.00 "
      "3.00)),((2.00 2.00,2.50 1.50,2.50 2.50,2.00 3.00)),((1.00 2.00,0.50 "
      "1.50,1.00 1.00,2.00 1.00,2.50 1.50,2.00 2.00)),((1.00 3.00,0.50 "
      "2.50,0.50 1.50,1.00 2.00)),((0.00 3.00,0.50 2.50,1.00 3.00)))");
  BOOST_CHECK_EQUAL(out->asText(2), expectedWKT);
}

BOOST_AUTO_TEST_CASE(testStraightSkeletonPartitionPolygonWithHole)
{
  std::unique_ptr<Geometry> g(io::readWkt(
      "POLYGON ((0 0, 0 4, 4 4, 4 0, 0 0), (1 1, 3 1, 3 3, 1 3, 1 1))"));
  std::unique_ptr<Geometry> out(algorithm::straightSkeletonPartition(*g));
  BOOST_CHECK(out->is<PolyhedralSurface>());
  std::string const expectedWKT(
      "POLYHEDRALSURFACE (((0.00 0.00,0.50 0.50,0.50 3.50,0.00 4.00)),((4.00 "
      "0.00,3.50 0.50,0.50 0.50,0.00 0.00)),((4.00 4.00,3.50 3.50,3.50 "
      "0.50,4.00 0.00)),((0.00 4.00,0.50 3.50,3.50 3.50,4.00 4.00)),((1.00 "
      "1.00,0.50 0.50,3.50 0.50,3.00 1.00)),((1.00 3.00,0.50 3.50,0.50 "
      "0.50,1.00 1.00)),((3.00 3.00,3.50 3.50,0.50 3.50,1.00 3.00)),((3.00 "
      "1.00,3.50 0.50,3.50 3.50,3.00 3.00)))");
  BOOST_CHECK_EQUAL(out->asText(2), expectedWKT);
}

BOOST_AUTO_TEST_CASE(testStraightSkeletonPartitionMultiPolygon)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("MULTIPOLYGON (((0 0, 0 2, 2 2, 2 0, 0 0)), ((3 3, 3 5, 5 5, "
                  "5 3, 3 3)))"));
  std::unique_ptr<Geometry> out(algorithm::straightSkeletonPartition(*g));
  BOOST_CHECK(out->is<PolyhedralSurface>());
  std::string const expectedWKT(
      "POLYHEDRALSURFACE (((0.00 0.00,1.00 1.00,0.00 2.00)),((2.00 0.00,1.00 "
      "1.00,0.00 0.00)),((2.00 2.00,1.00 1.00,2.00 0.00)),((0.00 2.00,1.00 "
      "1.00,2.00 2.00)),((3.00 3.00,4.00 4.00,3.00 5.00)),((5.00 3.00,4.00 "
      "4.00,3.00 3.00)),((5.00 5.00,4.00 4.00,5.00 3.00)),((3.00 5.00,4.00 "
      "4.00,5.00 5.00)))");
  BOOST_CHECK_EQUAL(out->asText(2), expectedWKT);
}

BOOST_AUTO_TEST_CASE(testStraightSkeletonPartitionEmptyPolygon)
{
  std::unique_ptr<Geometry> g(io::readWkt("POLYGON EMPTY"));
  std::unique_ptr<Geometry> out(algorithm::straightSkeletonPartition(*g));
  BOOST_CHECK(out->is<PolyhedralSurface>());
  BOOST_CHECK(out->isEmpty());
  std::string const expectedWKT("POLYHEDRALSURFACE EMPTY");
  BOOST_CHECK_EQUAL(out->asText(2), expectedWKT);
}

BOOST_AUTO_TEST_CASE(testStraightSkeletonPartitionNonPolygonGeometry)
{
  std::unique_ptr<Geometry> g(io::readWkt("LINESTRING (0 0, 1 1, 2 2)"));
  BOOST_CHECK_THROW(algorithm::straightSkeletonPartition(*g), std::exception);
}

BOOST_AUTO_TEST_CASE(testProjectMedialAxisToEdgesRectangle)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("POLYGON ((0 0, 10 0, 10 6, 0 6, 0 0))"));

  std::unique_ptr<MultiLineString> result(
      algorithm::projectMedialAxisToEdges(*g));

  // Should have one continuous line from left edge to right edge
  BOOST_CHECK_EQUAL(result->numGeometries(), 1U);

  const auto &line = result->geometryN(0).as<LineString>();
  BOOST_CHECK_EQUAL(line.numPoints(), 4U);

  // Check start point (left edge projection)
  BOOST_CHECK_CLOSE(CGAL::to_double(line.startPoint().x()), 0.0, 1e-6);
  BOOST_CHECK_CLOSE(CGAL::to_double(line.startPoint().y()), 3.0, 1e-6);

  // Check end point (right edge projection)
  BOOST_CHECK_CLOSE(CGAL::to_double(line.endPoint().x()), 10.0, 1e-6);
  BOOST_CHECK_CLOSE(CGAL::to_double(line.endPoint().y()), 3.0, 1e-6);

  // Check that medial axis points are included
  BOOST_CHECK_CLOSE(CGAL::to_double(line.pointN(1).x()), 3.0, 1e-6);
  BOOST_CHECK_CLOSE(CGAL::to_double(line.pointN(1).y()), 3.0, 1e-6);
  BOOST_CHECK_CLOSE(CGAL::to_double(line.pointN(2).x()), 7.0, 1e-6);
  BOOST_CHECK_CLOSE(CGAL::to_double(line.pointN(2).y()), 3.0, 1e-6);
}

BOOST_AUTO_TEST_CASE(testProjectMedialAxisToEdgesLShape)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("POLYGON ((0 0, 6 0, 6 4, 3 4, 3 6, 0 6, 0 0))"));

  std::unique_ptr<MultiLineString> result(
      algorithm::projectMedialAxisToEdges(*g));

  // Should have 3 segments (vertical, horizontal connectors, and horizontal)
  BOOST_CHECK_EQUAL(result->numGeometries(), 3U);

  // Verify each segment has extensions
  for (size_t i = 0; i < result->numGeometries(); ++i) {
    const auto &line = result->geometryN(i).as<LineString>();
    BOOST_CHECK_GE(line.numPoints(), 2U);
  }

  // Expected structure should contain:
  // - Vertical segment from top edge to junction
  // - Junction connection
  // - Horizontal segment from right edge to junction
  bool foundVerticalExtension   = false;
  bool foundHorizontalExtension = false;

  for (size_t i = 0; i < result->numGeometries(); ++i) {
    const auto &line = result->geometryN(i).as<LineString>();

    // Check for vertical extension (should reach y=6)
    if (line.numPoints() >= 3) {
      double maxY = std::max(CGAL::to_double(line.startPoint().y()),
                             CGAL::to_double(line.endPoint().y()));
      if (maxY > 5.9) {
        foundVerticalExtension = true;
      }
    }

    // Check for horizontal extension (should reach x=6)
    if (line.numPoints() >= 3) {
      double maxX = std::max(CGAL::to_double(line.startPoint().x()),
                             CGAL::to_double(line.endPoint().x()));
      if (maxX > 5.9) {
        foundHorizontalExtension = true;
      }
    }
  }

  BOOST_CHECK(foundVerticalExtension);
  BOOST_CHECK(foundHorizontalExtension);
}

BOOST_AUTO_TEST_CASE(testProjectMedialAxisToEdgesTShape)
{
  std::unique_ptr<Geometry> g(
      io::readWkt("POLYGON ((0 0, 10 0, 10 3, 6 3, 6 6, 3 6, 3 3, 0 3, 0 0))"));

  std::unique_ptr<MultiLineString> result(
      algorithm::projectMedialAxisToEdges(*g));

  // Should have 3 segments for T-shape
  BOOST_CHECK_EQUAL(result->numGeometries(), 3U);

  // Check that we have extensions to all 3 branches
  bool foundLeftExtension  = false;
  bool foundRightExtension = false;
  bool foundTopExtension   = false;

  for (size_t i = 0; i < result->numGeometries(); ++i) {
    const auto &line = result->geometryN(i).as<LineString>();

    // Check for left extension (should reach x=0)
    double minX = std::min(CGAL::to_double(line.startPoint().x()),
                           CGAL::to_double(line.endPoint().x()));
    if (minX < 0.1) {
      foundLeftExtension = true;
    }

    // Check for right extension (should reach x=10)
    double maxX = std::max(CGAL::to_double(line.startPoint().x()),
                           CGAL::to_double(line.endPoint().x()));
    if (maxX > 9.9) {
      foundRightExtension = true;
    }

    // Check for top extension (should reach y=6)
    double maxY = std::max(CGAL::to_double(line.startPoint().y()),
                           CGAL::to_double(line.endPoint().y()));
    if (maxY > 5.9) {
      foundTopExtension = true;
    }
  }

  BOOST_CHECK(foundLeftExtension);
  BOOST_CHECK(foundRightExtension);
  BOOST_CHECK(foundTopExtension);
}

BOOST_AUTO_TEST_CASE(testProjectMedialAxisToEdgesEmptyPolygon)
{
  std::unique_ptr<Geometry>        g(io::readWkt("POLYGON EMPTY"));
  std::unique_ptr<MultiLineString> result(
      algorithm::projectMedialAxisToEdges(*g));

  BOOST_CHECK(result->isEmpty());
  BOOST_CHECK_EQUAL(result->numGeometries(), 0U);
}

BOOST_AUTO_TEST_CASE(testProjectMedialAxisToEdgesInvalidGeometry)
{
  std::unique_ptr<Geometry> g(io::readWkt("LINESTRING (0 0, 1 1, 2 2)"));
  BOOST_CHECK_THROW(algorithm::projectMedialAxisToEdges(*g), Exception);
}

BOOST_AUTO_TEST_SUITE_END()
