/**
 *   SFCGAL
 *
 *   Copyright (C) 2012-2013 Oslandia <infos@oslandia.com>
 *   Copyright (C) 2012-2013 IGN (http://www.ign.fr)
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Library General Public
 *   License as published by the Free Software Foundation; either
 *   version 2 of the License, or (at your option) any later version.
 *
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Library General Public License for more details.

 *   You should have received a copy of the GNU Library General Public
 *   License along with this library; if not, see
 <http://www.gnu.org/licenses/>.
 */
#include "SFCGAL/capi/sfcgal_c.h"
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

#include <boost/test/unit_test.hpp>
#include <memory>

using namespace SFCGAL;

#include "../../../test_config.h"
// always after CGAL
using namespace boost::unit_test;

BOOST_AUTO_TEST_SUITE(SFCGAL_sfcgal_cTest)

bool hasError = false;

auto
on_error(const char * /*msg*/, ...) -> int
{
  hasError = true;
  return 0;
}

BOOST_AUTO_TEST_CASE(testIs3D)
{
  sfcgal_set_error_handlers(printf, on_error);

  // retrieve wkb from geometry via C++ api
  std::unique_ptr<Geometry> const g(io::readWkt("POLYGON ((0 0, 20 0, 20 10, 0 10, 0 0))"));
  // check
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_3d(g.get()));

  std::unique_ptr<Geometry> const g2(
      io::readWkt("POLYGON Z ((0 0 0, 20 0 0, 20 10 0, 0 10 0, 0 0 0))"));
  // check
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_3d(g2.get()));
}

BOOST_AUTO_TEST_CASE(testIsValid)
{
  sfcgal_set_error_handlers(printf, on_error);

  // ========== 2D
  // retrieve wkb from geometry via C++ api
  std::unique_ptr<Geometry> const g(io::readWkt("POLYGON((0 0,10 0,10 0,10 10,0 10,0 0))"));
  char *reason;
  sfcgal_geometry_t *location;
  int result = sfcgal_geometry_is_valid_detail(g.get(), &reason, &location);
  // check
  BOOST_CHECK_EQUAL(true, result);
  BOOST_CHECK_EQUAL(nullptr, reason);
  BOOST_CHECK_EQUAL(nullptr, location);

  std::unique_ptr<Geometry> const g2(io::readWkt("POLYGON((1 2,1 2,1 2,1 2))"));
  result = sfcgal_geometry_is_valid_detail(g2.get(), &reason, &location);
  // check
  BOOST_CHECK_EQUAL(false, result);
  BOOST_CHECK_EQUAL("ring 0 degenerated to a point", reason);
  BOOST_CHECK_EQUAL(nullptr, location);

  sfcgal_free_buffer(reason);
  if (location)
    sfcgal_geometry_delete(location);

  // ========== 3D
  std::unique_ptr<Geometry> const g3(
      io::readWkt("TRIANGLE((1.0 -1.0 -1.0,1.0 1.0 -1.0,1.0 -1.0 1.0,1.0 -1.0 -1.0))"));
  result = sfcgal_geometry_is_valid_detail(g3.get(), &reason, &location);
  // check
  BOOST_CHECK_EQUAL(true, result);
  BOOST_CHECK_EQUAL(nullptr, reason);
  BOOST_CHECK_EQUAL(nullptr, location);

  std::unique_ptr<Geometry> const g4(
      io::readWkt("TRIANGLE((1.0 -1.0 -1.0,1.0 1.0 -1.0,1.0 -1.0 -1.0,1.0 -1.0 -1.0))"));
  result = sfcgal_geometry_is_valid_detail(g4.get(), &reason, &location);
  // check
  BOOST_CHECK_EQUAL(false, result);
  BOOST_CHECK_EQUAL("ring 0 self intersects", reason);
  BOOST_CHECK_EQUAL(nullptr, location);

  sfcgal_free_buffer(reason);
  if (location)
    sfcgal_geometry_delete(location);
}

BOOST_AUTO_TEST_CASE(testIsMeasured)
{
  sfcgal_set_error_handlers(printf, on_error);

  // retrieve wkb from geometry via C++ api
  std::unique_ptr<Geometry> const g(io::readWkt("POLYGON ((0 0, 20 0, 20 10, 0 10, 0 0))"));
  // check
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_measured(g.get()));

  std::unique_ptr<Geometry> const g2(
      io::readWkt("POLYGON Z ((0 0 0, 20 0 0, 20 10 0, 0 10 0, 0 0 0))"));
  // check
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_measured(g2.get()));

  std::unique_ptr<Geometry> const g3(
      io::readWkt("POLYGON M ((0 0 1, 20 0 2, 20 10 3, 0 10 4, 0 0 1))"));
  // check
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_measured(g3.get()));
}

BOOST_AUTO_TEST_CASE(testIsSimple)
{
  sfcgal_set_error_handlers(printf, on_error);

  // retrieve wkb from geometry via C++ api
  std::unique_ptr<Geometry> const g(io::readWkt("LINESTRING (0.0 0.0, 2.0 0.0, 1.0 1.0)"));
  // check
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_simple(g.get()));

  std::unique_ptr<Geometry> const g2(
      io::readWkt("POLYGON Z ((0.0 0.0 0.0, 1.0 0.0 0.0, 1.0 1.0 0.0, "
                  "0.0 1.0 0.0, 0.0 0.0 0.0))"));
  // check
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_simple(g2.get()));

  std::unique_ptr<Geometry> const g3(
      io::readWkt("POLYGON Z ((0.0 0.0 1.0, 1.0 0.0 0.0, 1.0 1.0 0.0, "
                  "0.0 1.0 0.0, 0.0 0.0 1.0))"));
  char *reason;
  // check
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_simple_detail(g3.get(), &reason));
  BOOST_CHECK_EQUAL("Points don't lie in the same plane.", reason);

  sfcgal_free_buffer(reason);
}

/// Coordinate() ;
BOOST_AUTO_TEST_CASE(testErrorOnBadGeometryType)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const l(io::readWkt("LINESTRING (0 0, 0 1)"));
  std::unique_ptr<Geometry>       p(io::readWkt("POINT (0 2)"));
  sfcgal_geometry_t              *gl = l.get();

  hasError = false;
  BOOST_CHECK_EQUAL(2, sfcgal_linestring_num_points(gl)); // should succeed
  BOOST_CHECK(hasError == false);

  hasError = false;
  BOOST_CHECK(sfcgal_triangle_vertex(gl, 0) == nullptr); // should fail
  BOOST_CHECK(hasError == true);

  sfcgal_geometry_t *gp = p.release();
  hasError              = false;
  sfcgal_linestring_add_point(gl, gp); // should succeed
  BOOST_CHECK(hasError == false);

  hasError = false;
  sfcgal_linestring_add_point(gl, gl); // should fail
  BOOST_CHECK(hasError == true);
}

BOOST_AUTO_TEST_CASE(testAsWkb)
{
  sfcgal_set_error_handlers(printf, on_error);

  // retrieve wkb from geometry via C++ api
  std::unique_ptr<Geometry> const geom(
      io::readWkt("POLYGON ((0 0, 20 0, 20 10, 0 10, 0 0))"));

  std::string strGeom;
  strGeom = geom->asWkb();

  // retrieve wkb from C api
  char  *wkbApi;
  size_t wkbLen;
  sfcgal_geometry_as_wkb(geom.get(), &wkbApi, &wkbLen);
  std::string strApi(wkbApi, wkbLen);

  // check
  BOOST_CHECK_EQUAL(strGeom, strApi);
  sfcgal_free_buffer(wkbApi);
}

BOOST_AUTO_TEST_CASE(testStraightSkeletonPolygon)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const g(
      io::readWkt("POLYGON ((0 0, 20 0, 20 10, 0 10, 0 0))"));

  hasError              = false;
  sfcgal_geometry_t *sk = sfcgal_geometry_straight_skeleton(g.get());
  BOOST_CHECK(hasError == false);
  BOOST_CHECK_EQUAL(5, sfcgal_geometry_collection_num_geometries(sk));

  sfcgal_geometry_delete(sk);
}

BOOST_AUTO_TEST_CASE(testStraightSkeletonMultiPolygon)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const g(
      io::readWkt("MULTIPOLYGON (((0 0, 20 0, 20 10, 0 10, 0 0)),((100 0,200 "
                  "0,150 100,100 0)))"));

  hasError              = false;
  sfcgal_geometry_t *sk = sfcgal_geometry_straight_skeleton(g.get());
  BOOST_CHECK(hasError == false);
  BOOST_CHECK_EQUAL(8, sfcgal_geometry_collection_num_geometries(sk));

  sfcgal_geometry_delete(sk);
}

BOOST_AUTO_TEST_CASE(testApproximateMedialAxis)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const g(io::readWkt(
      "POLYGON ((-42 9,-44 9,-42 8,-22 7,-22 21,1 22,-5 13,-5 12,-4 "
      "13,2 23,-23 22,-23 8,-42 9))"));

  hasError              = false;
  sfcgal_geometry_t *sk = sfcgal_geometry_approximate_medial_axis(g.get());
  BOOST_CHECK(hasError == false);
  BOOST_CHECK_EQUAL(11, sfcgal_geometry_collection_num_geometries(sk));
  BOOST_CHECK_EQUAL(71.56, std::round(sfcgal_geometry_length(sk) * 100.0) / 100.0);

  sfcgal_geometry_delete(sk);
}

BOOST_AUTO_TEST_CASE(testCovers)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const g1(
      io::readWkt("MULTIPOLYGON (((0 0, 20 0, 20 10, 0 10, 0 0)),((100 0,200 "
                  "0,150 100,100 0)))"));
  std::unique_ptr<Geometry> const g2(io::readWkt(
      "MULTIPOLYGON (((100 0,200 0,150 100,100 0)), ((0 0, 20 0, 20 "
      "10, 0 10, 0 0)))"));

  BOOST_CHECK(sfcgal_geometry_covers(g1.get(), g2.get()));
}

BOOST_AUTO_TEST_CASE(testLineSubstring)
{
  sfcgal_set_error_handlers(printf, on_error);
  std::unique_ptr<Geometry> const g1(
      io::readWkt("LINESTRING Z (0 0 0, 0 0 10)"));
  std::unique_ptr<Geometry> const g2(
      io::readWkt("LINESTRING Z (0 0 3, 0 0 7)"));
  hasError              = false;
  sfcgal_geometry_t *ls = sfcgal_geometry_line_sub_string(g1.get(), 0.3, 0.7);
  BOOST_CHECK(hasError == false);

  BOOST_CHECK(sfcgal_geometry_covers_3d(ls, g2.get()));

  sfcgal_geometry_delete(ls);
}

BOOST_AUTO_TEST_CASE(testForceRHR)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::string strGeom{"POLYGON ((0 0,0 5,5 5,5 0,0 0),(1 1,2 1,2 2,1 2,1 1))"};

  std::unique_ptr<Geometry> const geom(io::readWkt(strGeom));

  sfcgal_geometry_t *rhr = sfcgal_geometry_force_rhr(geom.get());
  // retrieve wkb from C api
  char  *wkbApi;
  size_t wkbLen;
  sfcgal_geometry_as_text_decim(rhr, 0, &wkbApi, &wkbLen);
  std::string strApi(wkbApi, wkbLen);

  // check
  BOOST_CHECK_EQUAL(strGeom, strApi);
  sfcgal_free_buffer(wkbApi);
  sfcgal_geometry_delete(rhr);
}

BOOST_AUTO_TEST_CASE(testForceLHR)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::string strGeom{"POLYGON ((0 0,0 5,5 5,5 0,0 0),(1 1,2 1,2 2,1 2,1 1))"};
  std::string expectedGeom{
      "POLYGON ((0 0,5 0,5 5,0 5,0 0),(1 1,1 2,2 2,2 1,1 1))"};

  std::unique_ptr<Geometry> const geom(io::readWkt(strGeom));

  sfcgal_geometry_t *lhr = sfcgal_geometry_force_lhr(geom.get());
  // retrieve wkb from C api
  char  *wkbApi;
  size_t wkbLen;
  sfcgal_geometry_as_text_decim(lhr, 0, &wkbApi, &wkbLen);
  std::string strApi(wkbApi, wkbLen);

  // check
  BOOST_CHECK_EQUAL(expectedGeom, strApi);
  sfcgal_free_buffer(wkbApi);
  sfcgal_geometry_delete(lhr);
}

BOOST_AUTO_TEST_CASE(testForceRHR_3D)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::string strGeom{"POLYGON ((0 5 1,0 0 2,5 0 3,5 5 4,0 5 1),(2 1 1,1 1 2,1 "
                      "2 3,2 2 4,2 1 1),(4 3 1,3 3 2,3 4 3,4 4 4,4 3 1))"};
  std::string expectedGeom{
      "POLYGON Z ((0 5 1,5 5 4,5 0 3,0 0 2,0 5 1),(2 1 1,2 2 4,1 2 3,1 1 2,2 1 "
      "1),(4 3 1,4 4 4,3 4 3,3 3 2,4 3 1))"};

  std::unique_ptr<Geometry> const geom(io::readWkt(strGeom));

  sfcgal_geometry_t *rhr = sfcgal_geometry_force_rhr(geom.get());
  // retrieve wkb from C api
  char  *wkbApi;
  size_t wkbLen;
  sfcgal_geometry_as_text_decim(rhr, 0, &wkbApi, &wkbLen);
  std::string strApi(wkbApi, wkbLen);

  // check
  BOOST_CHECK_EQUAL(expectedGeom, strApi);
  sfcgal_free_buffer(wkbApi);
  sfcgal_geometry_delete(rhr);
}

BOOST_AUTO_TEST_CASE(testScaleUniformC)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const g(io::readWkt("POINT (1 2 3)"));

  hasError                  = false;
  sfcgal_geometry_t *scaled = sfcgal_geometry_scale(g.get(), 2.0);
  BOOST_CHECK(hasError == false);

  char  *wkt;
  size_t len;
  sfcgal_geometry_as_text_decim(scaled, 0, &wkt, &len);
  BOOST_CHECK_EQUAL(std::string(wkt), "POINT Z (2 4 6)");

  sfcgal_free_buffer(wkt);
  sfcgal_geometry_delete(scaled);
}

BOOST_AUTO_TEST_CASE(testScaleNonUniformC)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const g(io::readWkt("POINT (1 2 3)"));

  hasError                  = false;
  sfcgal_geometry_t *scaled = sfcgal_geometry_scale_3d(g.get(), 2.0, 3.0, 4.0);
  BOOST_CHECK(hasError == false);

  // FLAKY WKT test
  std::unique_ptr<Geometry> const g1(io::readWkt("POINT Z (2/1 6/1 12/1)"));
  BOOST_CHECK(sfcgal_geometry_covers(g1.get(), scaled));

  sfcgal_geometry_delete(scaled);
}

BOOST_AUTO_TEST_CASE(testScaleAroundCenterC)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const g(io::readWkt("POINT (3 4 5)"));

  hasError = false;
  sfcgal_geometry_t *scaled =
      sfcgal_geometry_scale_3d_around_center(g.get(), 2.0, 2.0, 2.0, 1, 1, 1);
  BOOST_CHECK(hasError == false);

  char  *wkt;
  size_t len;
  sfcgal_geometry_as_text_decim(scaled, 0, &wkt, &len);
  BOOST_CHECK_EQUAL(std::string(wkt), "POINT Z (5 7 9)");

  sfcgal_free_buffer(wkt);
  sfcgal_geometry_delete(scaled);
}

BOOST_AUTO_TEST_CASE(testScaleCubeNonUniformC)
{
  sfcgal_set_error_handlers(printf, on_error);

  // Create a cube with side length 10
  std::string cubeWkt = "POLYHEDRALSURFACE Z ("
                        "((0 0 0, 0 10 0, 10 10 0, 10 0 0, 0 0 0)),"
                        "((0 0 10, 10 0 10, 10 10 10, 0 10 10, 0 0 10)),"
                        "((0 0 0, 10 0 0, 10 0 10, 0 0 10, 0 0 0)),"
                        "((10 0 0, 10 10 0, 10 10 10, 10 0 10, 10 0 0)),"
                        "((0 10 0, 0 10 10, 10 10 10, 10 10 0, 0 10 0)),"
                        "((0 0 0, 0 0 10, 0 10 10, 0 10 0, 0 0 0)))";

  std::unique_ptr<Geometry> const g(io::readWkt(cubeWkt));

  hasError                  = false;
  sfcgal_geometry_t *scaled = sfcgal_geometry_scale_3d(g.get(), 0.5, 1.0, 2.0);
  BOOST_CHECK(hasError == false);

  char  *wkt;
  size_t len;
  sfcgal_geometry_as_text_decim(scaled, 0, &wkt, &len);

  // Check a few key points to ensure the scaling was applied correctly
  std::string scaledWkt(wkt);
  BOOST_CHECK(scaledWkt.find("0 0 0") !=
              std::string::npos); // Origin should remain unchanged
  BOOST_CHECK(scaledWkt.find("5 10 20") !=
              std::string::npos); // (10,10,10) should become (5,10,20)
  BOOST_CHECK(scaledWkt.find("5 0 0") !=
              std::string::npos); // (10,0,0) should become (5,0,0)
  BOOST_CHECK(scaledWkt.find("0 10 20") !=
              std::string::npos); // (0,10,10) should become (0,10,20)

  sfcgal_free_buffer(wkt);
  sfcgal_geometry_delete(scaled);
}

BOOST_AUTO_TEST_CASE(testRotate2D)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const g(io::readWkt("POINT (1 0)"));

  hasError                   = false;
  sfcgal_geometry_t *rotated = sfcgal_geometry_rotate(g.get(), M_PI / 2);
  BOOST_CHECK(hasError == false);

  char  *wkt;
  size_t len;
  sfcgal_geometry_as_text_decim(rotated, 0, &wkt, &len);
  BOOST_CHECK_EQUAL(std::string(wkt), "POINT (0 1)");

  sfcgal_free_buffer(wkt);
  sfcgal_geometry_delete(rotated);
}

BOOST_AUTO_TEST_CASE(testRotate2DAroundPoint)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const g(io::readWkt("POINT (2 0)"));

  hasError = false;
  sfcgal_geometry_t *rotated =
      sfcgal_geometry_rotate_2d(g.get(), M_PI / 2, 1, 0);
  BOOST_CHECK(hasError == false);

  char  *wkt;
  size_t len;
  sfcgal_geometry_as_text_decim(rotated, 0, &wkt, &len);
  BOOST_CHECK_EQUAL(std::string(wkt), "POINT (1 1)");

  sfcgal_free_buffer(wkt);
  sfcgal_geometry_delete(rotated);
}

BOOST_AUTO_TEST_CASE(testRotate3D)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const g(io::readWkt("POINT (1 0 0)"));

  hasError = false;
  sfcgal_geometry_t *rotated =
      sfcgal_geometry_rotate_3d(g.get(), M_PI / 2, 0, 0, 1);
  BOOST_CHECK(hasError == false);

  char  *wkt;
  size_t len;
  sfcgal_geometry_as_text_decim(rotated, 0, &wkt, &len);
  BOOST_CHECK_EQUAL(std::string(wkt), "POINT Z (0 1 0)");

  sfcgal_free_buffer(wkt);
  sfcgal_geometry_delete(rotated);
}

BOOST_AUTO_TEST_CASE(testEnvelope2D)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const g(io::readWkt(
      "MULTIPOLYGON (((0 0, 20 0, 20 10, 0 10, 0 0)), ((25 5, 30 5, 30 15, 25 15, 25 5)))"));

  hasError = false;
  sfcgal_geometry_t *result = sfcgal_geometry_envelope(g.get());
  BOOST_CHECK(hasError == false);

  char *wkt;
  size_t len;
  sfcgal_geometry_as_text_decim(result, 0, &wkt, &len);
  BOOST_CHECK_EQUAL(std::string(wkt), "POLYGON ((0 0,30 0,30 15,0 15,0 0))");

  sfcgal_free_buffer(wkt);
  sfcgal_geometry_delete(result);

  hasError = false;
  result = sfcgal_geometry_envelope_3d(g.get());
  BOOST_CHECK(hasError == false);

  sfcgal_geometry_as_text_decim(result, 0, &wkt, &len);
  BOOST_CHECK_EQUAL(std::string(wkt), "POLYGON ((0 0,30 0,30 15,0 15,0 0))");

  sfcgal_free_buffer(wkt);
  sfcgal_geometry_delete(result);
}

BOOST_AUTO_TEST_CASE(testEnvelope3D)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const g(
      io::readWkt("MULTILINESTRING Z ((0 0 0, 20 0 5, 20 10 10), (25 5 -5, 30 15 20, 25 15 25))"));

  hasError = false;
  sfcgal_geometry_t *result = sfcgal_geometry_envelope(g.get());
  BOOST_CHECK(hasError == false);

  std::string expected = "POLYGON ((0 0,30 0,30 15,0 15,0 0))";
  char *wkt;
  size_t len;
  sfcgal_geometry_as_text_decim(result, 0, &wkt, &len);
  BOOST_CHECK_EQUAL(std::string(wkt), expected);

  sfcgal_free_buffer(wkt);
  sfcgal_geometry_delete(result);

  hasError = false;
  result = sfcgal_geometry_envelope_3d(g.get());
  BOOST_CHECK(hasError == false);

  expected = "POLYHEDRALSURFACE Z (((0 0 -5,0 15 -5,30 15 -5,30 0 -5,0 0 -5)),"
             "((0 0 25,30 0 25,30 15 25,0 15 25,0 0 25)),"
             "((0 0 -5,30 0 -5,30 0 25,0 0 25,0 0 -5)),"
             "((30 15 -5,0 15 -5,0 15 25,30 15 25,30 15 -5)),"
             "((30 0 -5,30 15 -5,30 15 25,30 0 25,30 0 -5)),"
             "((0 0 -5,0 0 25,0 15 25,0 15 -5,0 0 -5)))";
  sfcgal_geometry_as_text_decim(result, 0, &wkt, &len);
  BOOST_CHECK_EQUAL(std::string(wkt), expected);

  sfcgal_free_buffer(wkt);
  sfcgal_geometry_delete(result);
}

BOOST_AUTO_TEST_CASE(testLength2D)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const g(io::readWkt("POLYGON ((0 0,30 0,30 15,0 15,0 0))"));
  hasError = false;
  double result = sfcgal_geometry_length(g.get());
  BOOST_CHECK(hasError == false);
  BOOST_CHECK_EQUAL(0.0, result);

  std::unique_ptr<Geometry> const g2(io::readWkt(
      "MULTIPOLYGON (((0 0, 20 0, 20 10, 0 10, 0 0)), ((25 5, 30 5, 30 15, 25 15, 25 5)))"));
  hasError = false;
  result = sfcgal_geometry_length(g2.get());
  BOOST_CHECK(hasError == false);
  BOOST_CHECK_EQUAL(0.0, result);

  std::unique_ptr<Geometry> const g3(io::readWkt("LINESTRING (0 0, 0 3, 4 3)"));
  hasError = false;
  result = sfcgal_geometry_length(g3.get());
  BOOST_CHECK(hasError == false);
  BOOST_CHECK_EQUAL(7.0, result);

  std::unique_ptr<Geometry> const g4(
      io::readWkt("MULTILINESTRING ((0 0, 0 3, 4 3), (10 0, 10 3, 14 3))"));
  hasError = false;
  result = sfcgal_geometry_length(g4.get());
  BOOST_CHECK(hasError == false);
  BOOST_CHECK_EQUAL(14.0, result);

  std::unique_ptr<Geometry> const g5(io::readWkt("LINESTRING Z (0 0 0, 0 3 10, 4 3 20)"));
  hasError = false;
  result = sfcgal_geometry_length(g5.get());
  BOOST_CHECK(hasError == false);
  BOOST_CHECK_EQUAL(7.0, result);
}

BOOST_AUTO_TEST_CASE(testLength3D)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const g(io::readWkt("POLYGON ((0 0,30 0,30 15,0 15,0 0))"));
  hasError = false;
  double result = sfcgal_geometry_length_3d(g.get());
  BOOST_CHECK(hasError == false);
  BOOST_CHECK_EQUAL(0.0, result);

  std::unique_ptr<Geometry> const g2(io::readWkt(
      "MULTIPOLYGON (((0 0, 20 0, 20 10, 0 10, 0 0)), ((25 5, 30 5, 30 15, 25 15, 25 5)))"));
  hasError = false;
  result = sfcgal_geometry_length_3d(g2.get());
  BOOST_CHECK(hasError == false);
  BOOST_CHECK_EQUAL(0.0, result);

  std::unique_ptr<Geometry> const g3(io::readWkt("LINESTRING (0 0, 0 3, 4 3)"));
  hasError = false;
  result = sfcgal_geometry_length(g3.get());
  BOOST_CHECK(hasError == false);
  BOOST_CHECK_EQUAL(7.0, result);

  std::unique_ptr<Geometry> const g4(
      io::readWkt("MULTILINESTRING ((0 0, 0 3, 4 3), (10 0, 10 3, 14 3))"));
  hasError = false;
  result = sfcgal_geometry_length_3d(g4.get());
  BOOST_CHECK(hasError == false);
  BOOST_CHECK_EQUAL(14.0, result);

  std::unique_ptr<Geometry> const g5(io::readWkt("LINESTRING Z (0 0 0, 0 3 10, 4 3 20)"));
  hasError = false;
  result = sfcgal_geometry_length_3d(g5.get());
  BOOST_CHECK(hasError == false);
  BOOST_CHECK_EQUAL(21.2106, std::round(result * 10000.0) / 10000.0);
}

BOOST_AUTO_TEST_CASE(testCentroid)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const g(io::readWkt(
      "MULTIPOLYGON (((0 0, 20 0, 20 10, 0 10, 0 0)), ((25 5, 30 5, 30 15, 25 15, 25 5)))"));

  hasError = false;
  sfcgal_geometry_t *result = sfcgal_geometry_centroid(g.get());
  BOOST_CHECK(hasError == false);

  char *wkt;
  size_t len;
  sfcgal_geometry_as_text_decim(result, 0, &wkt, &len);
  BOOST_CHECK_EQUAL(std::string(wkt), "POINT (14 6)");

  sfcgal_free_buffer(wkt);
  sfcgal_geometry_delete(result);
}

BOOST_AUTO_TEST_CASE(testCentroid3D)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const g(io::readWkt(
      "MULTIPOLYGON Z (((0 0 0, 20 0 0, 20 10 5, 0 10 5, 0 0 0)), ((25 5 0, 30 5 0, 30 15 -5, 25 15 -5, 25 5 0)))"));

  hasError = false;
  sfcgal_geometry_t *result = sfcgal_geometry_centroid(g.get());
  BOOST_CHECK(hasError == false);

  char *wkt;
  size_t len;
  sfcgal_geometry_as_text_decim(result, 0, &wkt, &len);
  BOOST_CHECK_EQUAL(std::string(wkt), "POINT Z (14 6 2)");

  sfcgal_free_buffer(wkt);
  sfcgal_geometry_delete(result);
}

BOOST_AUTO_TEST_CASE(testRotate3DAroundCenter)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const g(io::readWkt("POINT (2 0 0)"));

  hasError                   = false;
  sfcgal_geometry_t *rotated = sfcgal_geometry_rotate_3d_around_center(
      g.get(), M_PI / 2, 0, 0, 1, 1, 0, 0);
  BOOST_CHECK(hasError == false);

  char  *wkt;
  size_t len;
  sfcgal_geometry_as_text_decim(rotated, 0, &wkt, &len);
  BOOST_CHECK_EQUAL(std::string(wkt), "POINT Z (1 1 0)");

  sfcgal_free_buffer(wkt);
  sfcgal_geometry_delete(rotated);
}

BOOST_AUTO_TEST_CASE(testRotateX)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const g(io::readWkt("POINT (0 1 0)"));

  hasError                   = false;
  sfcgal_geometry_t *rotated = sfcgal_geometry_rotate_x(g.get(), M_PI / 2);
  BOOST_CHECK(hasError == false);

  char  *wkt;
  size_t len;
  sfcgal_geometry_as_text_decim(rotated, 0, &wkt, &len);
  BOOST_CHECK_EQUAL(std::string(wkt), "POINT Z (0 0 1)");

  sfcgal_free_buffer(wkt);
  sfcgal_geometry_delete(rotated);
}

BOOST_AUTO_TEST_CASE(testStraightSkeletonPartitionC)
{

  sfcgal_set_error_handlers(printf, on_error);
  std::unique_ptr<Geometry> const geom(
      io::readWkt("POLYGON ((0 0, 0 2, 1 2, 1 1, 2 1, 2 0, 0 0))"));
  std::string const expectedWKT(
      "MULTIPOLYGON (((0.00 0.00,0.50 0.50,0.50 1.50,0.00 2.00)),((2.00 "
      "0.00,1.50 0.50,0.50 0.50,0.00 0.00)),((2.00 1.00,1.50 0.50,2.00 "
      "0.00)),((1.00 1.00,0.50 0.50,1.50 0.50,2.00 1.00)),((1.00 2.00,0.50 "
      "1.50,0.50 0.50,1.00 1.00)),((0.00 2.00,0.50 1.50,1.00 2.00)))");
  sfcgal_geometry_t *result =
      sfcgal_geometry_straight_skeleton_partition(geom.get(), true);
  // retrieve wkb from C api
  char  *wkbApi;
  size_t wkbLen;
  sfcgal_geometry_as_text_decim(result, 2, &wkbApi, &wkbLen);
  std::string strApi(wkbApi, wkbLen);

  // check
  BOOST_CHECK_EQUAL(expectedWKT, strApi);
  sfcgal_free_buffer(wkbApi);
  sfcgal_geometry_delete(result);
}

BOOST_AUTO_TEST_CASE(testSolidSetExteriorShell)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Solid> const solid = std::make_unique<Solid>();
  BOOST_CHECK(solid->isEmpty());

  std::string const polyhedral1Str =
    "POLYHEDRALSURFACE ("
    "((0 0 0, 0 0 1, 0 1 1, 0 1 0, 0 0 0)),"
    "((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)),"
    "((0 0 0, 1 0 0, 1 0 1, 0 0 1, 0 0 0)),"
    "((1 1 0, 1 1 1, 1 0 1, 1 0 0, 1 1 0)),"
    "((0 1 0, 0 1 1, 1 1 1, 1 1 0, 0 1 0)),"
    "((0 0 1, 1 0 1, 1 1 1, 0 1 1, 0 0 1))"
    ")";

  std::unique_ptr<Geometry> shell1(io::readWkt(polyhedral1Str));
  BOOST_CHECK(!shell1->isEmpty());

  sfcgal_solid_set_exterior_shell(solid.get(), sfcgal_geometry_clone(shell1.get()));

  // check
  BOOST_CHECK(!solid->isEmpty());
  BOOST_CHECK(sfcgal_geometry_covers_3d(sfcgal_solid_shell_n(solid.get(), 0), shell1.get()));
}

BOOST_AUTO_TEST_CASE(testAlphaWrapping3DTest)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::string inputData(SFCGAL_TEST_DIRECTORY);
  inputData += "/data/bunny1000Wkt.txt";
  std::ifstream bunnyFSInput(inputData.c_str());
  BOOST_REQUIRE(bunnyFSInput.good());
  std::ostringstream inputWkt;
  inputWkt << bunnyFSInput.rdbuf();

  std::unique_ptr<Geometry> const geomInput(io::readWkt(inputWkt.str()));
  BOOST_REQUIRE(geomInput->is3D());

  sfcgal_geometry_t *geomAlphaWrapping = sfcgal_geometry_alpha_wrapping_3d(geomInput.get(), 20, 0);

  std::string resultData(SFCGAL_TEST_DIRECTORY);
#if CGAL_VERSION_MAJOR < 6
  resultData += "/data/bunny1000AlphaWrapping20Wkt_cgal5.txt";
#else
  resultData += "/data/bunny1000AlphaWrapping20Wkt_cgal6.txt";
#endif
  std::ifstream bunnyFSResult(resultData.c_str());
  BOOST_REQUIRE(bunnyFSResult.good());
  std::ostringstream resultWkt;
  resultWkt << bunnyFSResult.rdbuf();

  std::unique_ptr<Geometry> alphaWrappingExpectedGeom(io::readWkt(resultWkt.str()));
  BOOST_REQUIRE(alphaWrappingExpectedGeom->is3D());

  BOOST_CHECK(sfcgal_geometry_covers_3d(geomAlphaWrapping, alphaWrappingExpectedGeom.get()));
  sfcgal_geometry_delete(geomAlphaWrapping);
}

BOOST_AUTO_TEST_SUITE_END()
