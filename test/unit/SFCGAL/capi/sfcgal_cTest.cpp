// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/capi/sfcgal_c.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/io/wkt.h"

#include <array>
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

BOOST_AUTO_TEST_CASE(testEmpty)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const emptyGeom(io::readWkt("POLYGON EMPTY"));
  BOOST_CHECK(sfcgal_geometry_is_empty(emptyGeom.get()));
  BOOST_CHECK_EQUAL(sfcgal_geometry_num_geometries(emptyGeom.get()), 0);

  std::unique_ptr<Geometry> const polygonZ(
      io::readWkt("POLYGON Z ((0 0 0, 20 0 0, 20 10 0, 0 10 0, 0 0 0))"));
  BOOST_CHECK(!sfcgal_geometry_is_empty(polygonZ.get()));
  BOOST_CHECK_EQUAL(sfcgal_geometry_num_geometries(polygonZ.get()), 1);
}

BOOST_AUTO_TEST_CASE(testIs3D)
{
  sfcgal_set_error_handlers(printf, on_error);

  // retrieve wkb from geometry via C++ api
  std::unique_ptr<Geometry> const polygon1(
      io::readWkt("POLYGON ((0 0, 20 0, 20 10, 0 10, 0 0))"));
  // check
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_3d(polygon1.get()));

  std::unique_ptr<Geometry> const polygon2(
      io::readWkt("POLYGON Z ((0 0 0, 20 0 0, 20 10 0, 0 10 0, 0 0 0))"));
  // check
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_3d(polygon2.get()));
}

BOOST_AUTO_TEST_CASE(testDimension)
{
  sfcgal_set_error_handlers(printf, on_error);

  // ========== Empty Polygon
  std::unique_ptr<Geometry> const emptyPolygon(io::readWkt("POLYGON EMPTY"));
  BOOST_CHECK(sfcgal_geometry_is_empty(emptyPolygon.get()));
  BOOST_CHECK_EQUAL(sfcgal_geometry_dimension(emptyPolygon.get()), 2);

  // ========== Point
  std::unique_ptr<Geometry> const pointGeom(io::readWkt("POINTZ (4 4 9)"));
  BOOST_CHECK(!sfcgal_geometry_is_empty(pointGeom.get()));
  BOOST_CHECK_EQUAL(sfcgal_geometry_dimension(pointGeom.get()), 0);

  // ========== Line
  std::unique_ptr<Geometry> const lineGeom(
      io::readWkt("LINESTRING Z (-117 33 2, -116 34 4)"));
  BOOST_CHECK(!sfcgal_geometry_is_empty(lineGeom.get()));
  BOOST_CHECK_EQUAL(sfcgal_geometry_dimension(lineGeom.get()), 1);

  // ========== Polygon
  std::unique_ptr<Geometry> const polygonGeom2D(
      io::readWkt("POLYGON ((0 0, 0 3, 3 3, 3 0, 0 0))"));
  BOOST_CHECK(!sfcgal_geometry_is_empty(polygonGeom2D.get()));
  BOOST_CHECK_EQUAL(sfcgal_geometry_dimension(polygonGeom2D.get()), 2);

  std::unique_ptr<Geometry> const polygonGeom3D(
      io::readWkt("POLYGON Z ((0 0 1, 0 3 2, 3 3 3, 3 0 4, 0 0 1))"));
  BOOST_CHECK(!sfcgal_geometry_is_empty(polygonGeom3D.get()));
  BOOST_CHECK_EQUAL(sfcgal_geometry_dimension(polygonGeom3D.get()), 2);

  // ========== PolyhedralSurface
  std::unique_ptr<Geometry> const polyhedralSurface(io::readWkt(
      "POLYHEDRALSURFACE Z (((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)))"));
  BOOST_CHECK(!sfcgal_geometry_is_empty(polyhedralSurface.get()));
  BOOST_CHECK_EQUAL(sfcgal_geometry_dimension(polyhedralSurface.get()), 2);
}

BOOST_AUTO_TEST_CASE(testIsValid)
{
  sfcgal_set_error_handlers(printf, on_error);

  // ========== 2D
  // retrieve wkb from geometry via C++ api
  std::unique_ptr<Geometry> const validPolygon(
      io::readWkt("POLYGON((0 0,10 0,10 0,10 10,0 10,0 0))"));
  char              *reason;
  sfcgal_geometry_t *location;
  int                result =
      sfcgal_geometry_is_valid_detail(validPolygon.get(), &reason, &location);
  // check
  BOOST_CHECK_EQUAL(true, result);
  BOOST_CHECK_EQUAL(nullptr, reason);
  BOOST_CHECK_EQUAL(nullptr, location);

  std::unique_ptr<Geometry> const invalidPolygon(
      io::readWkt("POLYGON((1 2,1 2,1 2,1 2))"));
  result =
      sfcgal_geometry_is_valid_detail(invalidPolygon.get(), &reason, &location);
  // check
  BOOST_CHECK_EQUAL(false, result);
  BOOST_CHECK_EQUAL("ring 0 degenerated to a point", reason);
  BOOST_CHECK_EQUAL(nullptr, location);

  sfcgal_free_buffer(reason);
  if (location != nullptr) {
    sfcgal_geometry_delete(location);
  }

  // ========== 3D
  std::unique_ptr<Geometry> const invalidTriangle1(io::readWkt(
      "TRIANGLE((1.0 -1.0 -1.0,1.0 1.0 -1.0,1.0 -1.0 1.0,1.0 -1.0 -1.0))"));
  result = sfcgal_geometry_is_valid_detail(invalidTriangle1.get(), &reason,
                                           &location);
  // check
  BOOST_CHECK_EQUAL(true, result);
  BOOST_CHECK_EQUAL(nullptr, reason);
  BOOST_CHECK_EQUAL(nullptr, location);

  std::unique_ptr<Geometry> const invalidTriangle2(io::readWkt(
      "TRIANGLE((1.0 -1.0 -1.0,1.0 1.0 -1.0,1.0 -1.0 -1.0,1.0 -1.0 -1.0))"));
  result = sfcgal_geometry_is_valid_detail(invalidTriangle2.get(), &reason,
                                           &location);
  // check
  BOOST_CHECK_EQUAL(false, result);
  BOOST_CHECK_EQUAL("ring 0 self intersects", reason);
  BOOST_CHECK_EQUAL(nullptr, location);

  sfcgal_free_buffer(reason);
  if (location != nullptr) {
    sfcgal_geometry_delete(location);
  }
}

BOOST_AUTO_TEST_CASE(testIsMeasured)
{
  sfcgal_set_error_handlers(printf, on_error);

  // retrieve wkb from geometry via C++ api
  std::unique_ptr<Geometry> const polygon2D(
      io::readWkt("POLYGON ((0 0, 20 0, 20 10, 0 10, 0 0))"));
  // check
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_measured(polygon2D.get()));

  std::unique_ptr<Geometry> const polygon3D(
      io::readWkt("POLYGON Z ((0 0 0, 20 0 0, 20 10 0, 0 10 0, 0 0 0))"));
  // check
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_measured(polygon3D.get()));

  std::unique_ptr<Geometry> const polygonM(
      io::readWkt("POLYGON M ((0 0 1, 20 0 2, 20 10 3, 0 10 4, 0 0 1))"));
  // check
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_measured(polygonM.get()));
}

BOOST_AUTO_TEST_CASE(testDropZM)
{
  sfcgal_set_error_handlers(printf, on_error);

  // retrieve wkb from geometry via C++ api
  std::unique_ptr<Geometry> polygon2D(
      io::readWkt("POLYGON ((0 0, 20 0, 20 10, 0 10, 0 0))"));
  // check
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_3d(polygon2D.get()));
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_measured(polygon2D.get()));
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_drop_z(polygon2D.get()));
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_drop_m(polygon2D.get()));

  // 3D
  std::unique_ptr<Geometry> polygon3D(
      io::readWkt("POLYGON Z ((0 0 2, 20 0 2, 20 10 3, 0 10 2, 0 0 4))"));
  // check
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_3d(polygon3D.get()));
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_measured(polygon3D.get()));
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_drop_m(polygon3D.get()));
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_drop_z(polygon3D.get()));
  char  *wkbApi3D;
  size_t wkb3DLen;
  sfcgal_geometry_as_text_decim(polygon3D.get(), 0, &wkbApi3D, &wkb3DLen);
  std::string strApi3D(wkbApi3D, wkb3DLen);
  sfcgal_free_buffer(wkbApi3D);
  BOOST_CHECK_EQUAL(strApi3D, "POLYGON ((0 0,20 0,20 10,0 10,0 0))");
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_3d(polygon3D.get()));
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_drop_z(polygon3D.get()));

  // M
  std::unique_ptr<Geometry> polygonM(
      io::readWkt("POLYGON M ((0 0 1, 20 0 2, 20 10 3, 0 10 4, 0 0 1))"));
  // check
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_measured(polygonM.get()));
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_3d(polygonM.get()));
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_drop_z(polygonM.get()));
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_drop_m(polygonM.get()));
  char  *wkbApiM;
  size_t wkbMLen;
  sfcgal_geometry_as_text_decim(polygonM.get(), 0, &wkbApiM, &wkbMLen);
  std::string strApiM(wkbApiM, wkbMLen);
  sfcgal_free_buffer(wkbApiM);
  BOOST_CHECK_EQUAL(strApiM, "POLYGON ((0 0,20 0,20 10,0 10,0 0))");
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_3d(polygonM.get()));
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_drop_z(polygonM.get()));

  // ZM
  std::unique_ptr<Geometry> polygonZM(io::readWkt(
      "POLYGON ZM ((0 0 1 2, 20 0 2 2, 20 10 3 2, 0 10 4 2, 0 0 1 2))"));
  // check
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_measured(polygonZM.get()));
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_3d(polygonZM.get()));
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_drop_m(polygonZM.get()));
  char  *wkbApiZM1;
  size_t wkbZM1Len;
  sfcgal_geometry_as_text_decim(polygonZM.get(), 0, &wkbApiZM1, &wkbZM1Len);
  std::string strApiZM1(wkbApiZM1, wkbZM1Len);
  sfcgal_free_buffer(wkbApiZM1);
  BOOST_CHECK_EQUAL(strApiZM1,
                    "POLYGON Z ((0 0 1,20 0 2,20 10 3,0 10 4,0 0 1))");
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_3d(polygonZM.get()));
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_measured(polygonZM.get()));

  BOOST_CHECK_EQUAL(true, sfcgal_geometry_drop_z(polygonZM.get()));
  char  *wkbApiZM2;
  size_t wkbZM2Len;
  sfcgal_geometry_as_text_decim(polygonZM.get(), 0, &wkbApiZM2, &wkbZM2Len);
  std::string strApiZM2(wkbApiZM2, wkbZM2Len);
  sfcgal_free_buffer(wkbApiZM2);
  BOOST_CHECK_EQUAL(strApiZM2, "POLYGON ((0 0,20 0,20 10,0 10,0 0))");
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_3d(polygonZM.get()));
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_measured(polygonZM.get()));
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_drop_z(polygonZM.get()));
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_drop_m(polygonZM.get()));
}

BOOST_AUTO_TEST_CASE(testForceZ)
{
  sfcgal_set_error_handlers(printf, on_error);

  // retrieve wkb from geometry via C++ api
  std::unique_ptr<Geometry> polygon2D(
      io::readWkt("POLYGON ((0 0, 20 0, 20 10, 0 10, 0 0))"));
  // check
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_3d(polygon2D.get()));
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_measured(polygon2D.get()));
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_force_z(polygon2D.get(), 10));
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_3d(polygon2D.get()));
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_measured(polygon2D.get()));
  char  *wkbApi2D;
  size_t wkb2DLen;
  sfcgal_geometry_as_text_decim(polygon2D.get(), 0, &wkbApi2D, &wkb2DLen);
  std::string strApi2D(wkbApi2D, wkb2DLen);
  sfcgal_free_buffer(wkbApi2D);
  BOOST_CHECK_EQUAL(strApi2D,
                    "POLYGON Z ((0 0 10,20 0 10,20 10 10,0 10 10,0 0 10))");

  // 3D
  std::unique_ptr<Geometry> polygon3D(
      io::readWkt("POLYGON Z ((0 0 2, 20 0 2, 20 10 3, 0 10 2, 0 0 4))"));
  // check
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_3d(polygon3D.get()));
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_measured(polygon3D.get()));
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_force_z(polygon3D.get(), 10));
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_3d(polygon3D.get()));
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_measured(polygon3D.get()));
  char  *wkbApi3D;
  size_t wkb3DLen;
  sfcgal_geometry_as_text_decim(polygon3D.get(), 0, &wkbApi3D, &wkb3DLen);
  std::string strApi3D(wkbApi3D, wkb3DLen);
  sfcgal_free_buffer(wkbApi3D);
  BOOST_CHECK_EQUAL(strApi3D,
                    "POLYGON Z ((0 0 2,20 0 2,20 10 3,0 10 2,0 0 4))");

  // M
  std::unique_ptr<Geometry> polygonM(
      io::readWkt("POLYGON M ((0 0 1, 20 0 2, 20 10 3, 0 10 4, 0 0 1))"));
  // check
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_measured(polygonM.get()));
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_3d(polygonM.get()));
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_force_z(polygonM.get(), 2));
  char  *wkbApiM;
  size_t wkbMLen;
  sfcgal_geometry_as_text_decim(polygonM.get(), 0, &wkbApiM, &wkbMLen);
  std::string strApiM(wkbApiM, wkbMLen);
  sfcgal_free_buffer(wkbApiM);
  BOOST_CHECK_EQUAL(
      strApiM, "POLYGON ZM ((0 0 2 1,20 0 2 2,20 10 2 3,0 10 2 4,0 0 2 1))");
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_3d(polygonM.get()));
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_measured(polygonM.get()));

  // ZM
  std::unique_ptr<Geometry> polygonZM(io::readWkt(
      "POLYGON ZM ((0 0 1 2, 20 0 2 2, 20 10 3 2, 0 10 4 2, 0 0 1 2))"));
  // check
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_measured(polygonZM.get()));
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_3d(polygonZM.get()));
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_force_z(polygonM.get(), 0));
  char  *wkbApiZM1;
  size_t wkbZM1Len;
  sfcgal_geometry_as_text_decim(polygonZM.get(), 0, &wkbApiZM1, &wkbZM1Len);
  std::string strApiZM1(wkbApiZM1, wkbZM1Len);
  sfcgal_free_buffer(wkbApiZM1);
  BOOST_CHECK_EQUAL(
      strApiZM1, "POLYGON ZM ((0 0 1 2,20 0 2 2,20 10 3 2,0 10 4 2,0 0 1 2))");
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_3d(polygonZM.get()));
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_measured(polygonZM.get()));
}

BOOST_AUTO_TEST_CASE(testForceM)
{
  sfcgal_set_error_handlers(printf, on_error);

  // retrieve wkb from geometry via C++ api
  std::unique_ptr<Geometry> polygon2D(
      io::readWkt("POLYGON ((0 0, 20 0, 20 10, 0 10, 0 0))"));
  // check
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_3d(polygon2D.get()));
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_measured(polygon2D.get()));
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_force_m(polygon2D.get(), 10));
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_3d(polygon2D.get()));
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_measured(polygon2D.get()));
  char  *wkbApi2D;
  size_t wkb2DLen;
  sfcgal_geometry_as_text_decim(polygon2D.get(), 0, &wkbApi2D, &wkb2DLen);
  std::string strApi2D(wkbApi2D, wkb2DLen);
  sfcgal_free_buffer(wkbApi2D);
  BOOST_CHECK_EQUAL(strApi2D,
                    "POLYGON M ((0 0 10,20 0 10,20 10 10,0 10 10,0 0 10))");

  // 3D
  std::unique_ptr<Geometry> polygon3D(
      io::readWkt("POLYGON Z ((0 0 2, 20 0 2, 20 10 3, 0 10 2, 0 0 4))"));
  // check
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_3d(polygon3D.get()));
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_measured(polygon3D.get()));
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_force_m(polygon3D.get(), 10));
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_3d(polygon3D.get()));
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_measured(polygon3D.get()));
  char  *wkbApi3D;
  size_t wkb3DLen;
  sfcgal_geometry_as_text_decim(polygon3D.get(), 0, &wkbApi3D, &wkb3DLen);
  std::string strApi3D(wkbApi3D, wkb3DLen);
  sfcgal_free_buffer(wkbApi3D);
  BOOST_CHECK_EQUAL(
      strApi3D,
      "POLYGON ZM ((0 0 2 10,20 0 2 10,20 10 3 10,0 10 2 10,0 0 4 10))");

  // M
  std::unique_ptr<Geometry> polygonM(
      io::readWkt("POLYGON M ((0 0 1, 20 0 2, 20 10 3, 0 10 4, 0 0 1))"));
  // check
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_measured(polygonM.get()));
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_3d(polygonM.get()));
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_force_m(polygonM.get(), 2));
  char  *wkbApiM;
  size_t wkbMLen;
  sfcgal_geometry_as_text_decim(polygonM.get(), 0, &wkbApiM, &wkbMLen);
  std::string strApiM(wkbApiM, wkbMLen);
  sfcgal_free_buffer(wkbApiM);
  BOOST_CHECK_EQUAL(strApiM, "POLYGON M ((0 0 1,20 0 2,20 10 3,0 10 4,0 0 1))");
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_3d(polygonM.get()));
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_measured(polygonM.get()));

  // ZM
  std::unique_ptr<Geometry> polygonZM(io::readWkt(
      "POLYGON ZM ((0 0 1 2, 20 0 2 2, 20 10 3 2, 0 10 4 2, 0 0 1 2))"));
  // check
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_measured(polygonZM.get()));
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_3d(polygonZM.get()));
  BOOST_CHECK_EQUAL(false, sfcgal_geometry_force_m(polygonM.get(), 0));
  char  *wkbApiZM1;
  size_t wkbZM1Len;
  sfcgal_geometry_as_text_decim(polygonZM.get(), 0, &wkbApiZM1, &wkbZM1Len);
  std::string strApiZM1(wkbApiZM1, wkbZM1Len);
  sfcgal_free_buffer(wkbApiZM1);
  BOOST_CHECK_EQUAL(
      strApiZM1, "POLYGON ZM ((0 0 1 2,20 0 2 2,20 10 3 2,0 10 4 2,0 0 1 2))");
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_3d(polygonZM.get()));
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_measured(polygonZM.get()));
}

BOOST_AUTO_TEST_CASE(testSwapXY)
{
  sfcgal_set_error_handlers(printf, on_error);

  // retrieve wkb from geometry via C++ api
  std::unique_ptr<Geometry> polygon2D(
      io::readWkt("POLYGON ((0 0, 20 0, 20 10, 0 10, 0 0))"));
  // check
  sfcgal_geometry_swap_xy(polygon2D.get());
  char  *wkbApi2D;
  size_t wkb2DLen;
  sfcgal_geometry_as_text_decim(polygon2D.get(), 0, &wkbApi2D, &wkb2DLen);
  std::string strApi2D(wkbApi2D, wkb2DLen);
  sfcgal_free_buffer(wkbApi2D);
  BOOST_CHECK_EQUAL(strApi2D, "POLYGON ((0 0,0 20,10 20,10 0,0 0))");

  // 3D
  std::unique_ptr<Geometry> polygon3D(
      io::readWkt("POLYGON Z ((0 0 2, 20 0 2, 20 10 3, 0 10 2, 0 0 4))"));
  // check
  sfcgal_geometry_swap_xy(polygon3D.get());
  char  *wkbApi3D;
  size_t wkb3DLen;
  sfcgal_geometry_as_text_decim(polygon3D.get(), 0, &wkbApi3D, &wkb3DLen);
  std::string strApi3D(wkbApi3D, wkb3DLen);
  sfcgal_free_buffer(wkbApi3D);
  BOOST_CHECK_EQUAL(strApi3D,
                    "POLYGON Z ((0 0 2,0 20 2,10 20 3,10 0 2,0 0 4))");

  // M
  std::unique_ptr<Geometry> polygonM(
      io::readWkt("POLYGON M ((0 0 1, 20 0 2, 20 10 3, 0 10 4, 0 0 1))"));
  // check
  sfcgal_geometry_swap_xy(polygonM.get());
  char  *wkbApiM;
  size_t wkbMLen;
  sfcgal_geometry_as_text_decim(polygonM.get(), 0, &wkbApiM, &wkbMLen);
  std::string strApiM(wkbApiM, wkbMLen);
  sfcgal_free_buffer(wkbApiM);
  BOOST_CHECK_EQUAL(strApiM, "POLYGON M ((0 0 1,0 20 2,10 20 3,10 0 4,0 0 1))");

  // ZM
  std::unique_ptr<Geometry> polygonZM(io::readWkt(
      "POLYGON ZM ((0 0 1 2, 20 0 2 2, 20 10 3 2, 0 10 4 2, 0 0 1 2))"));
  // check
  sfcgal_geometry_swap_xy(polygonZM.get());
  char  *wkbApiZM;
  size_t wkbZMLen;
  sfcgal_geometry_as_text_decim(polygonZM.get(), 0, &wkbApiZM, &wkbZMLen);
  std::string strApiZM(wkbApiZM, wkbZMLen);
  sfcgal_free_buffer(wkbApiZM);
  BOOST_CHECK_EQUAL(
      strApiZM, "POLYGON ZM ((0 0 1 2,0 20 2 2,10 20 3 2,10 0 4 2,0 0 1 2))");
}

BOOST_AUTO_TEST_CASE(testType)
{
  // PolyhedralSurface
  std::string const               polySurfaceStr = "POLYHEDRALSURFACE Z ("
                                                   "((0 0 0, 2 0 0, 2 2 1, 0 2 1, 0 0 0)),"
                                                   "((2 0 0, 4 0 0, 4 2 1, 2 2 1, 2 0 0)),"
                                                   "((0 2 1, 2 2 1, 1 3 2, 0 2 1))"
                                                   ")";
  std::unique_ptr<Geometry> const polySurface(io::readWkt(polySurfaceStr));
  char                           *polyType;
  size_t                          polyLen;
  sfcgal_geometry_type(polySurface.get(), &polyType, &polyLen);
  std::string strPolyType(polyType, polyLen);
  sfcgal_free_buffer(polyType);
  BOOST_CHECK_EQUAL(strPolyType, "PolyhedralSurface");

  // LineString
  std::unique_ptr<Geometry> const lineGeom(
      io::readWkt("LINESTRING (0.0 0.0, 2.0 0.0, 1.0 1.0)"));
  char  *lineType;
  size_t lineLen;
  sfcgal_geometry_type(lineGeom.get(), &lineType, &lineLen);
  std::string strLineType(lineType, lineLen);
  sfcgal_free_buffer(lineType);
  BOOST_CHECK_EQUAL(strLineType, "LineString");

  // Point
  std::unique_ptr<Geometry> const pointGeom(
      io::readWkt("POINT Z (0.0 0.0 3.0)"));
  char  *pointType;
  size_t pointLen;
  sfcgal_geometry_type(pointGeom.get(), &pointType, &pointLen);
  std::string strPointType(pointType, pointLen);
  sfcgal_free_buffer(pointType);
  BOOST_CHECK_EQUAL(strPointType, "Point");

  // Polygon
  std::unique_ptr<Geometry> const polygonGeom(
      io::readWkt("POLYGON ((0 0, 0 3, 3 3, 3 0, 0 0))"));
  char  *polygonType;
  size_t polygonLen;
  sfcgal_geometry_type(polygonGeom.get(), &polygonType, &polygonLen);
  std::string strPolygonType(polygonType, polygonLen);
  sfcgal_free_buffer(polygonType);
  BOOST_CHECK_EQUAL(strPolygonType, "Polygon");
}

BOOST_AUTO_TEST_CASE(testIsSimple)
{
  sfcgal_set_error_handlers(printf, on_error);

  // retrieve wkb from geometry via C++ api
  std::unique_ptr<Geometry> const line(
      io::readWkt("LINESTRING (0.0 0.0, 2.0 0.0, 1.0 1.0)"));
  // check
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_simple(line.get()));

  std::unique_ptr<Geometry> const polygonZ1(
      io::readWkt("POLYGON Z ((0.0 0.0 0.0, 1.0 0.0 0.0, 1.0 1.0 0.0, "
                  "0.0 1.0 0.0, 0.0 0.0 0.0))"));
  // check
  BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_simple(polygonZ1.get()));

  std::unique_ptr<Geometry> const polygonZ2(
      io::readWkt("POLYGON Z ((0.0 0.0 1.0, 1.0 0.0 0.0, 1.0 1.0 0.0, "
                  "0.0 1.0 0.0, 0.0 0.0 1.0))"));
  char *reason;
  // check
  BOOST_CHECK_EQUAL(false,
                    sfcgal_geometry_is_simple_detail(polygonZ2.get(), &reason));
  BOOST_CHECK_EQUAL("Points don't lie in the same plane.", reason);

  sfcgal_free_buffer(reason);
}

BOOST_AUTO_TEST_CASE(testIsEqual)
{
  sfcgal_set_error_handlers(printf, on_error);

  {
    // same wkt
    std::unique_ptr<Geometry> const line1(
        io::readWkt("LINESTRING (0.0 0.0, 2.0 0.0, 1.0 1.0)"));
    std::unique_ptr<Geometry> const line2(
        io::readWkt("LINESTRING (0.0 0.0, 2.0 0.0, 1.0 1.0)"));
    // check
    hasError = false;
    BOOST_CHECK_EQUAL(true,
                      sfcgal_geometry_is_equals(line1.get(), line2.get()));
    BOOST_CHECK(hasError == false);
  }

  {
    // 0.1 diff between wkt
    std::unique_ptr<Geometry> const line1(
        io::readWkt("LINESTRING (0.0 0.0, 2.0 0.0, 1.0 1.0)"));
    std::unique_ptr<Geometry> const line2(
        io::readWkt("LINESTRING (0.1 0.1, 2.1 0.1, 1.1 1.1)"));
    // check
    hasError = false;
    BOOST_CHECK_EQUAL(
        true, sfcgal_geometry_is_almost_equals(line1.get(), line2.get(), 0.11));
    BOOST_CHECK(hasError == false);
  }

  {
    // retrieve wkb from geometry via C++ api
    std::unique_ptr<Geometry> const line1(
        io::readWkt("LINESTRING (0.0 0.0, 2.0 0.0, 1.0 1.0)"));
    std::unique_ptr<Geometry> const line2(
        io::readWkt("LINESTRING (0.1 0.1, 2.1 0.1, 1.1 1.1)"));
    // check
    hasError = false;
    BOOST_CHECK_EQUAL(true, sfcgal_geometry_is_almost_equals(
                                line1.get(), line2.get(), 0.100008));
    BOOST_CHECK_EQUAL(false, sfcgal_geometry_is_almost_equals(
                                 line1.get(), line2.get(), 0.099993));
    BOOST_CHECK(hasError == false);
  }
}

/// Coordinate() ;
BOOST_AUTO_TEST_CASE(testErrorOnBadGeometryType)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const line(io::readWkt("LINESTRING (0 0, 0 1)"));
  std::unique_ptr<Geometry>       point(io::readWkt("POINT (0 2)"));
  sfcgal_geometry_t              *lineGeom = line.get();

  hasError = false;
  BOOST_CHECK_EQUAL(2,
                    sfcgal_linestring_num_points(lineGeom)); // should succeed
  BOOST_CHECK(hasError == false);

  hasError = false;
  BOOST_CHECK(sfcgal_triangle_vertex(lineGeom, 0) == nullptr); // should fail
  BOOST_CHECK(hasError == true);

  sfcgal_geometry_t *pointGeom = point.release();
  hasError                     = false;
  sfcgal_linestring_add_point(lineGeom, pointGeom); // should succeed
  BOOST_CHECK(hasError == false);

  hasError = false;
  sfcgal_linestring_add_point(lineGeom, lineGeom); // should fail
  BOOST_CHECK(hasError == true);
}

BOOST_AUTO_TEST_CASE(testGeometryN)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const line(io::readWkt("LINESTRING (0 0, 0 1)"));
  std::unique_ptr<Geometry> const point(io::readWkt("POINT (0 2)"));

  std::string                     geomCollectionStr = "GEOMETRYCOLLECTION ("
                                                      "POINT (1 1),"
                                                      "LINESTRING (0 0, 2 2, 3 3),"
                                                      "POLYGON ((0 0, 4 0, 4 4, 0 4, 0 0)),"
                                                      "MULTIPOINT ((2 2), (3 3), (4 4))"
                                                      ")";
  std::unique_ptr<Geometry> const geomCollection(
      io::readWkt(geomCollectionStr));
  std::unique_ptr<Geometry> const geomCollection0(io::readWkt("POINT (1 1)"));
  std::unique_ptr<Geometry> const geomCollection2(
      io::readWkt("POLYGON ((0 0, 4 0, 4 4, 0 4, 0 0))"));

  std::string                     polySurfaceStr = "POLYHEDRALSURFACE Z ("
                                                   "((0 0 0, 2 0 0, 2 2 1, 0 2 1, 0 0 0)),"
                                                   "((2 0 0, 4 0 0, 4 2 1, 2 2 1, 2 0 0)),"
                                                   "((0 2 1, 2 2 1, 1 3 2, 0 2 1))"
                                                   ")";
  std::unique_ptr<Geometry> const polySurface(io::readWkt(polySurfaceStr));
  std::unique_ptr<Geometry> const polySurface0(
      io::readWkt("POLYGON Z ((0 0 0, 2 0 0, 2 2 1, 0 2 1, 0 0 0))"));
  std::unique_ptr<Geometry> const polySurface2(
      io::readWkt("POLYGON Z ((0 2 1, 2 2 1, 1 3 2, 0 2 1))"));

  std::string                     tinStr = "TIN Z ("
                                           "((0 0 0, 1 0 0, 0.5 1 1, 0 0 0)),"
                                           "((1 0 0, 1.5 1 1, 0.5 1 1, 1 0 0)),"
                                           "((0.5 1 1, 1.5 1 1, 1 2 2, 0.5 1 1))"
                                           ")";
  std::unique_ptr<Geometry> const tin(io::readWkt(tinStr));
  std::unique_ptr<Geometry> const tin0(
      io::readWkt("TRIANGLE Z ((0 0 0, 1 0 0, 0.5 1 1, 0 0 0))"));
  std::unique_ptr<Geometry> const tin2(
      io::readWkt("TRIANGLE Z ((0.5 1 1, 1.5 1 1, 1 2 2, 0.5 1 1))"));

  // GeometryCollection
  hasError = false;
  BOOST_CHECK(!geomCollection->isEmpty());
  BOOST_CHECK(!geomCollection0->isEmpty());
  BOOST_CHECK(!geomCollection2->isEmpty());
  BOOST_CHECK_EQUAL(sfcgal_geometry_num_geometries(geomCollection.get()), 4);
  BOOST_CHECK(sfcgal_geometry_covers(
      sfcgal_geometry_collection_geometry_n(geomCollection.get(), 0),
      geomCollection0.get()));
  BOOST_CHECK(hasError == false);
  hasError = false;
  BOOST_CHECK(sfcgal_geometry_covers(
      sfcgal_geometry_collection_geometry_n(geomCollection.get(), 2),
      geomCollection2.get()));
  BOOST_CHECK(hasError == false);
  sfcgal_geometry_collection_set_geometry_n(geomCollection.get(),
                                            point->clone().release(), 1);
  BOOST_CHECK(sfcgal_geometry_covers(
      sfcgal_geometry_collection_geometry_n(geomCollection.get(), 1),
      point.get()));

  // PolyhedralSurface
  hasError = false;
  BOOST_CHECK(!polySurface->isEmpty());
  BOOST_CHECK(!polySurface0->isEmpty());
  BOOST_CHECK(!polySurface2->isEmpty());
  BOOST_CHECK_EQUAL(sfcgal_geometry_num_geometries(polySurface.get()), 1);
  BOOST_CHECK_EQUAL(sfcgal_polyhedral_surface_num_patches(polySurface.get()),
                    3);
  BOOST_CHECK(sfcgal_geometry_covers_3d(
      sfcgal_polyhedral_surface_patch_n(polySurface.get(), 0),
      polySurface0.get()));
  BOOST_CHECK(hasError == false);
  hasError = false;
  BOOST_CHECK(sfcgal_geometry_covers_3d(
      sfcgal_polyhedral_surface_patch_n(polySurface.get(), 2),
      polySurface2.get()));
  BOOST_CHECK(hasError == false);
  std::unique_ptr<Geometry> const simplePolygon(
      io::readWkt("POLYGON Z ((0 0 0, 2 0 1, 2 2 2, 0 2 1, 0 0 0))"));
  sfcgal_polyhedral_surface_set_patch_n(polySurface.get(),
                                        simplePolygon->clone().release(), 1);
  BOOST_CHECK(sfcgal_geometry_covers_3d(
      sfcgal_polyhedral_surface_patch_n(polySurface.get(), 1),
      simplePolygon.get()));

  // TIN
  hasError = false;
  BOOST_CHECK(!tin->isEmpty());
  BOOST_CHECK(!tin0->isEmpty());
  BOOST_CHECK(!tin2->isEmpty());
  BOOST_CHECK_EQUAL(sfcgal_geometry_num_geometries(tin.get()), 1);
  BOOST_CHECK_EQUAL(sfcgal_triangulated_surface_num_patches(tin.get()), 3);
  BOOST_CHECK(sfcgal_geometry_covers_3d(
      sfcgal_triangulated_surface_patch_n(tin.get(), 0), tin0.get()));
  BOOST_CHECK(hasError == false);
  hasError = false;
  BOOST_CHECK(sfcgal_geometry_covers_3d(
      sfcgal_triangulated_surface_patch_n(tin.get(), 2), tin2.get()));
  BOOST_CHECK(hasError == false);
  std::unique_ptr<Geometry> const simpleTriangle(
      io::readWkt("TRIANGLE Z ((0 0 0, 2 0 1, 2 2 2, 0 0 0))"));
  sfcgal_triangulated_surface_set_patch_n(tin.get(),
                                          simpleTriangle->clone().release(), 1);
  BOOST_CHECK(sfcgal_geometry_covers_3d(
      sfcgal_triangulated_surface_patch_n(tin.get(), 1), simpleTriangle.get()));

  // Line
  hasError = false;
  BOOST_CHECK_EQUAL(sfcgal_geometry_num_geometries(line.get()), 1);
  BOOST_CHECK(hasError == false);

  // Point
  hasError = false;
  BOOST_CHECK_EQUAL(sfcgal_geometry_num_geometries(point.get()), 1);
  BOOST_CHECK(hasError == false);
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

  std::unique_ptr<Geometry> const polygon(
      io::readWkt("POLYGON ((0 0, 20 0, 20 10, 0 10, 0 0))"));

  hasError = false;
  sfcgal_geometry_t *skeleton =
      sfcgal_geometry_straight_skeleton(polygon.get());
  BOOST_CHECK(hasError == false);
  BOOST_CHECK_EQUAL(5, sfcgal_geometry_num_geometries(skeleton));

  sfcgal_geometry_delete(skeleton);
}

BOOST_AUTO_TEST_CASE(testStraightSkeletonMultiPolygon)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const multiPolygon(
      io::readWkt("MULTIPOLYGON (((0 0, 20 0, 20 10, 0 10, 0 0)),((100 0,200 "
                  "0,150 100,100 0)))"));

  hasError = false;
  sfcgal_geometry_t *skeleton =
      sfcgal_geometry_straight_skeleton(multiPolygon.get());
  BOOST_CHECK(hasError == false);
  BOOST_CHECK_EQUAL(8, sfcgal_geometry_num_geometries(skeleton));

  sfcgal_geometry_delete(skeleton);
}

BOOST_AUTO_TEST_CASE(testApproximateMedialAxis)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const polygon(io::readWkt(
      "POLYGON ((-42 9,-44 9,-42 8,-22 7,-22 21,1 22,-5 13,-5 12,-4 "
      "13,2 23,-23 22,-23 8,-42 9))"));

  hasError = false;
  sfcgal_geometry_t *skeleton =
      sfcgal_geometry_approximate_medial_axis(polygon.get());
  BOOST_CHECK(hasError == false);
  BOOST_CHECK_EQUAL(11, sfcgal_geometry_num_geometries(skeleton));
  BOOST_CHECK_EQUAL(
      71.56, std::round(sfcgal_geometry_length(skeleton) * 100.0) / 100.0);

  sfcgal_geometry_delete(skeleton);
}

BOOST_AUTO_TEST_CASE(testCovers)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const multiPolygon1(
      io::readWkt("MULTIPOLYGON (((0 0, 20 0, 20 10, 0 10, 0 0)),((100 0,200 "
                  "0,150 100,100 0)))"));
  std::unique_ptr<Geometry> const multiPolygon2(io::readWkt(
      "MULTIPOLYGON (((100 0,200 0,150 100,100 0)), ((0 0, 20 0, 20 "
      "10, 0 10, 0 0)))"));

  BOOST_CHECK(sfcgal_geometry_covers(multiPolygon1.get(), multiPolygon2.get()));
}

BOOST_AUTO_TEST_CASE(testLineSubstring)
{
  sfcgal_set_error_handlers(printf, on_error);
  std::unique_ptr<Geometry> const line1(
      io::readWkt("LINESTRING Z (0 0 0, 0 0 10)"));
  std::unique_ptr<Geometry> const line2(
      io::readWkt("LINESTRING Z (0 0 3, 0 0 7)"));
  hasError = false;
  sfcgal_geometry_t *subLine =
      sfcgal_geometry_line_sub_string(line1.get(), 0.3, 0.7);
  BOOST_CHECK(hasError == false);

  BOOST_CHECK(sfcgal_geometry_covers_3d(subLine, line2.get()));

  sfcgal_geometry_delete(subLine);
}

BOOST_AUTO_TEST_CASE(testTriangle)
{
  // Create an empty triangle
  sfcgal_geometry_t *emptyTriangle = sfcgal_triangle_create();
  BOOST_CHECK_EQUAL(sfcgal_geometry_type_id(emptyTriangle),
                    SFCGAL_TYPE_TRIANGLE);
  BOOST_CHECK(sfcgal_geometry_is_empty(emptyTriangle));
  sfcgal_geometry_delete(emptyTriangle);

  // Create triangle from points
  sfcgal_geometry_t *pt1 = sfcgal_point_create_from_xyz(0., 0., 0.);
  sfcgal_geometry_t *pt2 = sfcgal_point_create_from_xyz(4., 0., 0.);
  sfcgal_geometry_t *pt3 = sfcgal_point_create_from_xyz(1., 1.5, 1.);
  sfcgal_geometry_t *triangle1 =
      sfcgal_triangle_create_from_points(pt1, pt2, pt3);
  BOOST_CHECK_EQUAL(sfcgal_geometry_type_id(triangle1), SFCGAL_TYPE_TRIANGLE);
  BOOST_CHECK(!sfcgal_geometry_is_empty(triangle1));
  sfcgal_geometry_t const *firstPt = sfcgal_triangle_vertex(triangle1, 0);
  BOOST_CHECK(sfcgal_geometry_covers_3d(pt1, firstPt));
  sfcgal_geometry_t const *secondPt = sfcgal_triangle_vertex(triangle1, 1);
  BOOST_CHECK(sfcgal_geometry_covers_3d(pt2, secondPt));
  sfcgal_geometry_t const *thirdPt = sfcgal_triangle_vertex(triangle1, 2);
  BOOST_CHECK(sfcgal_geometry_covers_3d(pt3, thirdPt));
  sfcgal_geometry_delete(pt1);
  sfcgal_geometry_delete(pt2);
  sfcgal_geometry_delete(pt3);
  sfcgal_geometry_delete(triangle1);

  // triangle from wkt
  sfcgal_set_error_handlers(printf, on_error);
  std::unique_ptr<Geometry> const triangle2(
      io::readWkt("TRIANGLE Z ((0 0 0, 2 0 0, 1.5 1 1, 0 0 0))"));
  BOOST_CHECK_EQUAL(sfcgal_geometry_type_id(triangle2.get()),
                    SFCGAL_TYPE_TRIANGLE);
  sfcgal_geometry_t const *vertex = sfcgal_triangle_vertex(triangle2.get(), 0);
  sfcgal_geometry_t *expectedVertex = sfcgal_point_create_from_xyz(0., 0., 0.);
  BOOST_CHECK(sfcgal_geometry_covers_3d(vertex, expectedVertex));
  sfcgal_geometry_delete(expectedVertex);

  // set vertex
  sfcgal_geometry_t *newVertex = sfcgal_point_create_from_xyz(0.5, 0.5, 0.5);
  sfcgal_triangle_set_vertex(triangle2.get(), 0, newVertex);
  BOOST_CHECK(sfcgal_geometry_covers_3d(
      sfcgal_triangle_vertex(triangle2.get(), 0), newVertex));
  sfcgal_geometry_delete(newVertex);

  sfcgal_triangle_set_vertex_from_xyz(triangle2.get(), 0, 0.4, 0.4, 0.4);
  sfcgal_geometry_t *expectedNewVertex =
      sfcgal_point_create_from_xyz(0.4, 0.4, 0.4);
  BOOST_CHECK(sfcgal_geometry_covers_3d(
      sfcgal_triangle_vertex(triangle2.get(), 0), expectedNewVertex));
  sfcgal_geometry_delete(expectedNewVertex);
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

  std::unique_ptr<Geometry> const point2D(io::readWkt("POINT (1 2 3)"));

  hasError                  = false;
  sfcgal_geometry_t *scaled = sfcgal_geometry_scale(point2D.get(), 2.0);
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

  std::unique_ptr<Geometry> const point2D(io::readWkt("POINT (1 2 3)"));

  hasError = false;
  sfcgal_geometry_t *scaled =
      sfcgal_geometry_scale_3d(point2D.get(), 2.0, 3.0, 4.0);
  BOOST_CHECK(hasError == false);

  // FLAKY WKT test
  std::unique_ptr<Geometry> const pointZ(io::readWkt("POINT Z (2/1 6/1 12/1)"));
  BOOST_CHECK(sfcgal_geometry_covers(pointZ.get(), scaled));

  sfcgal_geometry_delete(scaled);
}

BOOST_AUTO_TEST_CASE(testScaleAroundCenterC)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const point2D(io::readWkt("POINT (3 4 5)"));

  hasError                  = false;
  sfcgal_geometry_t *scaled = sfcgal_geometry_scale_3d_around_center(
      point2D.get(), 2.0, 2.0, 2.0, 1, 1, 1);
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

  std::unique_ptr<Geometry> const cube(io::readWkt(cubeWkt));

  hasError = false;
  sfcgal_geometry_t *scaled =
      sfcgal_geometry_scale_3d(cube.get(), 0.5, 1.0, 2.0);
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

  std::unique_ptr<Geometry> const point2D(io::readWkt("POINT (1 0)"));

  hasError                   = false;
  sfcgal_geometry_t *rotated = sfcgal_geometry_rotate(point2D.get(), M_PI / 2);
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

  std::unique_ptr<Geometry> const point2D(io::readWkt("POINT (2 0)"));

  hasError = false;
  sfcgal_geometry_t *rotated =
      sfcgal_geometry_rotate_2d(point2D.get(), M_PI / 2, 1, 0);
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

  std::unique_ptr<Geometry> const pointZ(io::readWkt("POINT (1 0 0)"));

  hasError = false;
  sfcgal_geometry_t *rotated =
      sfcgal_geometry_rotate_3d(pointZ.get(), M_PI / 2, 0, 0, 1);
  BOOST_CHECK(hasError == false);

  char  *wkt;
  size_t len;
  sfcgal_geometry_as_text_decim(rotated, 0, &wkt, &len);
  BOOST_CHECK_EQUAL(std::string(wkt), "POINT Z (0 1 0)");

  sfcgal_free_buffer(wkt);
  sfcgal_geometry_delete(rotated);
}

BOOST_AUTO_TEST_CASE(testTransform)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const multiPolygon(io::readWkt(
      "MULTIPOLYGON Z (((0 0 0, 20 0 0, 20 10 0, 0 10 0, 0 0 0))," //
      "((25 5 1, 30 5 1, 30 15 1, 25 15 1, 25 5 1)))"));

  // translate
  {
    std::array<float, 16> mat = {1, 0, 0, 0, //
                                 0, 1, 0, 0, //
                                 0, 0, 1, 0, //
                                 1, 2, 3, 1};

    hasError = false;
    sfcgal_geometry_t *result =
        sfcgal_geometry_transform(multiPolygon.get(), mat.data());
    BOOST_CHECK(hasError == false);

    char  *wkt;
    size_t len;
    sfcgal_geometry_as_text_decim(result, 0, &wkt, &len);
    BOOST_CHECK_EQUAL(
        std::string(wkt),
        "MULTIPOLYGON Z (((1 2 3,21 2 3,21 12 3,1 12 3,1 2 3))," //
        "((26 7 4,31 7 4,31 17 4,26 17 4,26 7 4)))");

    sfcgal_free_buffer(wkt);
    sfcgal_geometry_delete(result);
  }

  // scale
  {
    std::array<float, 16> mat = {1, 0, 0, 0, //
                                 0, 2, 0, 0, //
                                 0, 0, 3, 0, //
                                 0, 0, 0, 1};

    hasError = false;
    sfcgal_geometry_t *result =
        sfcgal_geometry_transform(multiPolygon.get(), mat.data());
    BOOST_CHECK(hasError == false);

    char  *wkt;
    size_t len;
    sfcgal_geometry_as_text_decim(result, 0, &wkt, &len);
    BOOST_CHECK_EQUAL(
        std::string(wkt),
        "MULTIPOLYGON Z (((0 0 0,20 0 0,20 20 0,0 20 0,0 0 0))," //
        "((25 10 3,30 10 3,30 30 3,25 30 3,25 10 3)))");

    sfcgal_free_buffer(wkt);
    sfcgal_geometry_delete(result);
  }

  // rotate
  {
    std::array<float, 16> mat = {-1.0e-07, 1.0,      0.0, 0.0, //
                                 -1.0,     -1.0e-07, 0.0, 0.0, //
                                 0.0,      0.0,      1.0, 0.0, //
                                 0.0,      0.0,      0.0, 1.0};

    hasError = false;
    sfcgal_geometry_t *result =
        sfcgal_geometry_transform(multiPolygon.get(), mat.data());
    BOOST_CHECK(hasError == false);

    char  *wkt;
    size_t len;
    sfcgal_geometry_as_text_decim(result, 0, &wkt, &len);
    BOOST_CHECK_EQUAL(
        std::string(wkt),
        "MULTIPOLYGON Z (((0 0 0,0 20 0,-10 20 0,-10 0 0,0 0 0))," //
        "((-5 25 1,-5 30 1,-15 30 1,-15 25 1,-5 25 1)))");

    sfcgal_free_buffer(wkt);
    sfcgal_geometry_delete(result);
  }
}

BOOST_AUTO_TEST_CASE(testEnvelope2D)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const multiPolygon(
      io::readWkt("MULTIPOLYGON (((0 0, 20 0, 20 10, 0 10, 0 0)), ((25 5, 30 "
                  "5, 30 15, 25 15, 25 5)))"));

  hasError                  = false;
  sfcgal_geometry_t *result = sfcgal_geometry_envelope(multiPolygon.get());
  BOOST_CHECK(hasError == false);

  char  *wkt;
  size_t len;
  sfcgal_geometry_as_text_decim(result, 0, &wkt, &len);
  BOOST_CHECK_EQUAL(std::string(wkt), "POLYGON ((0 0,30 0,30 15,0 15,0 0))");

  sfcgal_free_buffer(wkt);
  sfcgal_geometry_delete(result);

  hasError = false;
  result   = sfcgal_geometry_envelope_3d(multiPolygon.get());
  BOOST_CHECK(hasError == false);

  sfcgal_geometry_as_text_decim(result, 0, &wkt, &len);
  BOOST_CHECK_EQUAL(std::string(wkt), "POLYGON ((0 0,30 0,30 15,0 15,0 0))");

  sfcgal_free_buffer(wkt);
  sfcgal_geometry_delete(result);
}

BOOST_AUTO_TEST_CASE(testEnvelope3D)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const multiLine(
      io::readWkt("MULTILINESTRING Z ((0 0 0, 20 0 5, 20 10 10), (25 5 -5, 30 "
                  "15 20, 25 15 25))"));

  hasError                  = false;
  sfcgal_geometry_t *result = sfcgal_geometry_envelope(multiLine.get());
  BOOST_CHECK(hasError == false);

  std::string expected = "POLYGON ((0 0,30 0,30 15,0 15,0 0))";
  char       *wkt;
  size_t      len;
  sfcgal_geometry_as_text_decim(result, 0, &wkt, &len);
  BOOST_CHECK_EQUAL(std::string(wkt), expected);

  sfcgal_free_buffer(wkt);
  sfcgal_geometry_delete(result);

  hasError = false;
  result   = sfcgal_geometry_envelope_3d(multiLine.get());
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

  std::unique_ptr<Geometry> const polygon(
      io::readWkt("POLYGON ((0 0,30 0,30 15,0 15,0 0))"));
  hasError      = false;
  double result = sfcgal_geometry_length(polygon.get());
  BOOST_CHECK(hasError == false);
  BOOST_CHECK_EQUAL(0.0, result);

  std::unique_ptr<Geometry> const multiPolygon(
      io::readWkt("MULTIPOLYGON (((0 0, 20 0, 20 10, 0 10, 0 0)), ((25 5, 30 "
                  "5, 30 15, 25 15, 25 5)))"));
  hasError = false;
  result   = sfcgal_geometry_length(multiPolygon.get());
  BOOST_CHECK(hasError == false);
  BOOST_CHECK_EQUAL(0.0, result);

  std::unique_ptr<Geometry> const line2D(
      io::readWkt("LINESTRING (0 0, 0 3, 4 3)"));
  hasError = false;
  result   = sfcgal_geometry_length(line2D.get());
  BOOST_CHECK(hasError == false);
  BOOST_CHECK_EQUAL(7.0, result);

  std::unique_ptr<Geometry> const multiLine(
      io::readWkt("MULTILINESTRING ((0 0, 0 3, 4 3), (10 0, 10 3, 14 3))"));
  hasError = false;
  result   = sfcgal_geometry_length(multiLine.get());
  BOOST_CHECK(hasError == false);
  BOOST_CHECK_EQUAL(14.0, result);

  std::unique_ptr<Geometry> const lineZ(
      io::readWkt("LINESTRING Z (0 0 0, 0 3 10, 4 3 20)"));
  hasError = false;
  result   = sfcgal_geometry_length(lineZ.get());
  BOOST_CHECK(hasError == false);
  BOOST_CHECK_EQUAL(7.0, result);
}

BOOST_AUTO_TEST_CASE(testLength3D)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const polygon(
      io::readWkt("POLYGON ((0 0,30 0,30 15,0 15,0 0))"));
  hasError      = false;
  double result = sfcgal_geometry_length_3d(polygon.get());
  BOOST_CHECK(hasError == false);
  BOOST_CHECK_EQUAL(0.0, result);

  std::unique_ptr<Geometry> const multiPolygon(
      io::readWkt("MULTIPOLYGON (((0 0, 20 0, 20 10, 0 10, 0 0)), ((25 5, 30 "
                  "5, 30 15, 25 15, 25 5)))"));
  hasError = false;
  result   = sfcgal_geometry_length_3d(multiPolygon.get());
  BOOST_CHECK(hasError == false);
  BOOST_CHECK_EQUAL(0.0, result);

  std::unique_ptr<Geometry> const line2D(
      io::readWkt("LINESTRING (0 0, 0 3, 4 3)"));
  hasError = false;
  result   = sfcgal_geometry_length(line2D.get());
  BOOST_CHECK(hasError == false);
  BOOST_CHECK_EQUAL(7.0, result);

  std::unique_ptr<Geometry> const multiLine(
      io::readWkt("MULTILINESTRING ((0 0, 0 3, 4 3), (10 0, 10 3, 14 3))"));
  hasError = false;
  result   = sfcgal_geometry_length_3d(multiLine.get());
  BOOST_CHECK(hasError == false);
  BOOST_CHECK_EQUAL(14.0, result);

  std::unique_ptr<Geometry> const lineZ(
      io::readWkt("LINESTRING Z (0 0 0, 0 3 10, 4 3 20)"));
  hasError = false;
  result   = sfcgal_geometry_length_3d(lineZ.get());
  BOOST_CHECK(hasError == false);
  BOOST_CHECK_EQUAL(21.2106, std::round(result * 10000.0) / 10000.0);
}

BOOST_AUTO_TEST_CASE(testBoundary)
{
  sfcgal_set_error_handlers(printf, on_error);

  // 2D Polygon
  std::unique_ptr<Geometry> const polygon2D(
      io::readWkt("POLYGON ((0 0,30 0,30 15,0 15,0 0))"));
  hasError = false;
  sfcgal_geometry_t *polygon2DBoundary =
      sfcgal_geometry_boundary(polygon2D.get());
  BOOST_CHECK(hasError == false);

  std::string expectedPolygon2DBoundary =
      "LINESTRING (0 0,30 0,30 15,0 15,0 0)";
  char  *wktPolygon2DBoundary;
  size_t lenPolygon2DBoundary;
  sfcgal_geometry_as_text_decim(polygon2DBoundary, 0, &wktPolygon2DBoundary,
                                &lenPolygon2DBoundary);
  BOOST_CHECK_EQUAL(std::string(wktPolygon2DBoundary),
                    expectedPolygon2DBoundary);

  sfcgal_free_buffer(wktPolygon2DBoundary);
  sfcgal_geometry_delete(polygon2DBoundary);

  // 2D MultiPolygon
  std::unique_ptr<Geometry> const multiPolygon2D(
      io::readWkt("MULTIPOLYGON (((0 0, 20 0, 20 10, 0 10, 0 0)), ((25 5, 30 "
                  "5, 30 15, 25 15, 25 5)))"));
  hasError = false;
  sfcgal_geometry_t *multiPolygon2DBoundary =
      sfcgal_geometry_boundary(multiPolygon2D.get());
  BOOST_CHECK(hasError == false);

  std::string expectedMultiPolygon2DBoundary =
      "MULTILINESTRING ((0 0,20 0),(20 0,20 10),(20 10,0 10),(0 10,0 0),(25 "
      "5,30 5),(30 5,30 15),(30 15,25 15),(25 15,25 5))";
  char  *wktMultiPolygon2DBoundary;
  size_t lenMultiPolygon2DBoundary;
  sfcgal_geometry_as_text_decim(multiPolygon2DBoundary, 0,
                                &wktMultiPolygon2DBoundary,
                                &lenMultiPolygon2DBoundary);
  BOOST_CHECK_EQUAL(std::string(wktMultiPolygon2DBoundary),
                    expectedMultiPolygon2DBoundary);

  sfcgal_free_buffer(wktMultiPolygon2DBoundary);
  sfcgal_geometry_delete(multiPolygon2DBoundary);

  // 3D Polygon
  std::unique_ptr<Geometry> const polygon3D(io::readWkt(
      "POLYGON Z ((0 0,0 10,10 10,10 0,0 0),(1 1 1,1 2 1,2 2 1,2 1 1,1 1 1))"));
  hasError = false;
  sfcgal_geometry_t *polygon3DBoundary =
      sfcgal_geometry_boundary(polygon3D.get());
  BOOST_CHECK(hasError == false);

  std::string expectedPolygon3DBoundary =
      "MULTILINESTRING ((0 0,0 10,10 10,10 0,0 0),(1 1 1,1 2 1,2 2 1,2 1 1,1 1 "
      "1))";
  char  *wktPolygon3DBoundary;
  size_t lenPolygon3DBoundary;
  sfcgal_geometry_as_text_decim(polygon3DBoundary, 0, &wktPolygon3DBoundary,
                                &lenPolygon3DBoundary);
  BOOST_CHECK_EQUAL(std::string(wktPolygon3DBoundary),
                    expectedPolygon3DBoundary);

  sfcgal_free_buffer(wktPolygon3DBoundary);
  sfcgal_geometry_delete(polygon3DBoundary);

  // 2D LineString
  std::unique_ptr<Geometry> const line2D(
      io::readWkt("LINESTRING (0 0, 0 3, 4 3)"));
  hasError                          = false;
  sfcgal_geometry_t *line2DBoundary = sfcgal_geometry_boundary(line2D.get());
  BOOST_CHECK(hasError == false);

  std::string expectedLine2DBoundary = "MULTIPOINT ((0 0),(4 3))";
  char       *wktLine2DBoundary;
  size_t      lenLine2DBoundary;
  sfcgal_geometry_as_text_decim(line2DBoundary, 0, &wktLine2DBoundary,
                                &lenLine2DBoundary);
  BOOST_CHECK_EQUAL(std::string(wktLine2DBoundary), expectedLine2DBoundary);

  sfcgal_free_buffer(wktLine2DBoundary);
  sfcgal_geometry_delete(line2DBoundary);

  // 2D MultiLineString
  std::unique_ptr<Geometry> const multiLine3D(io::readWkt(
      "MULTILINESTRING Z ((0 0 1, 1 1 2, 2 2 3), (3 3 4, 4 4 5, 5 5 6))"));
  hasError = false;
  sfcgal_geometry_t *multiLine3DBoundary =
      sfcgal_geometry_boundary(multiLine3D.get());
  BOOST_CHECK(hasError == false);

  std::string expectedMultiLine3DBoundary =
      "MULTIPOINT Z ((0 0 1),(2 2 3),(3 3 4),(5 5 6))";
  char  *wktMultiLine3DBoundary;
  size_t lenMultiLine3DBoundary;
  sfcgal_geometry_as_text_decim(multiLine3DBoundary, 0, &wktMultiLine3DBoundary,
                                &lenMultiLine3DBoundary);
  BOOST_CHECK_EQUAL(std::string(wktMultiLine3DBoundary),
                    expectedMultiLine3DBoundary);

  sfcgal_free_buffer(wktMultiLine3DBoundary);
  sfcgal_geometry_delete(multiLine3DBoundary);

  // 3D LineString
  std::unique_ptr<Geometry> const line3D(
      io::readWkt("LINESTRING Z (0 0 0, 0 3 10, 4 3 20)"));
  hasError                          = false;
  sfcgal_geometry_t *line3DBoundary = sfcgal_geometry_boundary(line3D.get());
  BOOST_CHECK(hasError == false);

  std::string expectedLine3DBoundary = "MULTIPOINT Z ((0 0 0),(4 3 20))";
  char       *wktLine3DBoundary;
  size_t      lenLine3DBoundary;
  sfcgal_geometry_as_text_decim(line3DBoundary, 0, &wktLine3DBoundary,
                                &lenLine3DBoundary);
  BOOST_CHECK_EQUAL(std::string(wktLine3DBoundary), expectedLine3DBoundary);

  sfcgal_free_buffer(wktLine3DBoundary);
  sfcgal_geometry_delete(line3DBoundary);
}

BOOST_AUTO_TEST_CASE(testCentroid)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const multiPolygon(
      io::readWkt("MULTIPOLYGON (((0 0, 20 0, 20 10, 0 10, 0 0)), ((25 5, 30 "
                  "5, 30 15, 25 15, 25 5)))"));

  hasError                  = false;
  sfcgal_geometry_t *result = sfcgal_geometry_centroid(multiPolygon.get());
  BOOST_CHECK(hasError == false);

  char  *wkt;
  size_t len;
  sfcgal_geometry_as_text_decim(result, 0, &wkt, &len);
  BOOST_CHECK_EQUAL(std::string(wkt), "POINT (14 6)");

  sfcgal_free_buffer(wkt);
  sfcgal_geometry_delete(result);
}

BOOST_AUTO_TEST_CASE(testCentroid3D)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const multiPolygonZ(
      io::readWkt("MULTIPOLYGON Z (((0 0 0, 20 0 0, 20 10 5, 0 10 5, 0 0 0)), "
                  "((25 5 0, 30 5 0, 30 15 -5, 25 15 -5, 25 5 0)))"));

  hasError                  = false;
  sfcgal_geometry_t *result = sfcgal_geometry_centroid(multiPolygonZ.get());
  BOOST_CHECK(hasError == false);

  char  *wkt;
  size_t len;
  sfcgal_geometry_as_text_decim(result, 0, &wkt, &len);
  BOOST_CHECK_EQUAL(std::string(wkt), "POINT Z (14 6 2)");

  sfcgal_free_buffer(wkt);
  sfcgal_geometry_delete(result);
}

BOOST_AUTO_TEST_CASE(testRotate3DAroundCenter)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const pointZ(io::readWkt("POINT (2 0 0)"));

  hasError                   = false;
  sfcgal_geometry_t *rotated = sfcgal_geometry_rotate_3d_around_center(
      pointZ.get(), M_PI / 2, 0, 0, 1, 1, 0, 0);
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

  std::unique_ptr<Geometry> const pointZ(io::readWkt("POINT (0 1 0)"));

  hasError                   = false;
  sfcgal_geometry_t *rotated = sfcgal_geometry_rotate_x(pointZ.get(), M_PI / 2);
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
      "POLYHEDRALSURFACE (((0.00 0.00,0.50 0.50,0.50 1.50,0.00 2.00)),((2.00 "
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

  std::string const polyhedral1Str = "POLYHEDRALSURFACE ("
                                     "((0 0 0, 0 0 1, 0 1 1, 0 1 0, 0 0 0)),"
                                     "((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)),"
                                     "((0 0 0, 1 0 0, 1 0 1, 0 0 1, 0 0 0)),"
                                     "((1 1 0, 1 1 1, 1 0 1, 1 0 0, 1 1 0)),"
                                     "((0 1 0, 0 1 1, 1 1 1, 1 1 0, 0 1 0)),"
                                     "((0 0 1, 1 0 1, 1 1 1, 0 1 1, 0 0 1))"
                                     ")";

  std::unique_ptr<Geometry> shell1(io::readWkt(polyhedral1Str));
  BOOST_CHECK(!shell1->isEmpty());

  sfcgal_solid_set_exterior_shell(solid.get(),
                                  sfcgal_geometry_clone(shell1.get()));

  // check
  BOOST_CHECK(!solid->isEmpty());
  BOOST_CHECK(sfcgal_geometry_covers_3d(sfcgal_solid_shell_n(solid.get(), 0),
                                        shell1.get()));
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

  sfcgal_geometry_t *geomAlphaWrapping =
      sfcgal_geometry_alpha_wrapping_3d(geomInput.get(), 20, 0);

  BOOST_REQUIRE(sfcgal_geometry_is_3d(geomAlphaWrapping));
#if CGAL_VERSION_MAJOR < 6
  // 2304 on Linux
  // 2306 on Mac/FreeBSD maybe other
  BOOST_CHECK_GE(sfcgal_polyhedral_surface_num_patches(geomAlphaWrapping),
                 2304);
#else
  BOOST_CHECK_EQUAL(sfcgal_polyhedral_surface_num_patches(geomAlphaWrapping),
                    2386);
#endif

  sfcgal_geometry_delete(geomAlphaWrapping);
}

BOOST_AUTO_TEST_CASE(testIsClosed)
{
  sfcgal_set_error_handlers(printf, on_error);

  std::unique_ptr<Geometry> const closedLineString(
      io::readWkt("LineString(0 0, 0 1, 1 1, 1 0, 0 0)"));
  BOOST_CHECK(sfcgal_geometry_is_closed(closedLineString.get()));

  std::string wkt = "POLYHEDRALSURFACE Z("
                    "((0 0 0, 2 0 0, 2 2 0, 0 2 0, 0 0 0))," // Base
                    "((0 0 0, 1 1 2, 2 0 0, 0 0 0)),"        // Front
                    "((2 0 0, 1 1 2, 2 2 0, 2 0 0)),"        // Right
                    "((2 2 0, 1 1 2, 0 2 0, 2 2 0)),"        // Back
                    "((0 2 0, 1 1 2, 0 0 0, 0 2 0))"         // Left
                    ")";

  std::unique_ptr<Geometry> const polyhedralSurface(io::readWkt(wkt));
  BOOST_CHECK(sfcgal_geometry_is_closed(polyhedralSurface.get()));
}

BOOST_AUTO_TEST_CASE(testSphereTest)
{
  sfcgal_set_error_handlers(printf, on_error);

  sfcgal_primitive_t *sphere = sfcgal_primitive_create(SFCGAL_TYPE_SPHERE);

  // checks parameter list
  char  *params;
  size_t paramsLen;
  sfcgal_primitive_parameters(sphere, &params, &paramsLen);
  std::string paramsStr(params, paramsLen);
  std::string expectedStr(
      R"([{"name":"center","type":"point3"},{"name":"direction","type":"vector3"},{"name":"num_horizontal","type":"int"},{"name":"num_vertical","type":"int"},{"name":"radius","type":"double"}])");
  BOOST_CHECK(compare_json(paramsStr, expectedStr));
  sfcgal_free_buffer(params);

  // radius parameter
  double radius = sfcgal_primitive_parameter_double(sphere, "radius");
  BOOST_CHECK_CLOSE(radius, 1.0, 1e-6);

  sfcgal_primitive_set_parameter_double(sphere, "radius", 24.2);

  double newRadius = sfcgal_primitive_parameter_double(sphere, "radius");
  BOOST_CHECK_CLOSE(newRadius, 24.2, 1e-6);

  // there is no parameter called min_radius
  sfcgal_primitive_set_parameter_double(sphere, "min_radius", 24.2);
  BOOST_CHECK(hasError);
  hasError = false;

  // center is a point not a double
  sfcgal_primitive_set_parameter_double(sphere, "center", 1.2);
  BOOST_CHECK(hasError);
  hasError = false;

  // radius cannot be negative
  sfcgal_primitive_set_parameter_double(sphere, "radius", -2.2);
  BOOST_CHECK(hasError);
  hasError = false;

  // center parameter
  double *center = sfcgal_primitive_parameter_point(sphere, "center");
  BOOST_CHECK(center != nullptr);
  BOOST_CHECK_CLOSE(center[0], 0.0, 1e-6);
  BOOST_CHECK_CLOSE(center[1], 0.0, 1e-6);
  BOOST_CHECK_CLOSE(center[2], 0.0, 1e-6);
  sfcgal_free_buffer(center);

  std::array<double, 3> expectedCenter{1.0, 2.0, 3.0};
  sfcgal_primitive_set_parameter_point(sphere, "center", expectedCenter.data());

  double *newCenter = sfcgal_primitive_parameter_point(sphere, "center");
  BOOST_CHECK(newCenter != nullptr);
  BOOST_CHECK_CLOSE(newCenter[0], expectedCenter[0], 1e-6);
  BOOST_CHECK_CLOSE(newCenter[1], expectedCenter[1], 1e-6);
  BOOST_CHECK_CLOSE(newCenter[2], expectedCenter[2], 1e-6);
  sfcgal_free_buffer(newCenter);

  // num vertical parameter
  unsigned int numVertical =
      sfcgal_primitive_parameter_int(sphere, "num_vertical");
  BOOST_CHECK_EQUAL(numVertical, 16);

  sfcgal_primitive_set_parameter_int(sphere, "num_vertical", 36);

  double newNumVertical =
      sfcgal_primitive_parameter_int(sphere, "num_vertical");
  BOOST_CHECK_EQUAL(newNumVertical, 36);

  // num horizontal parameter
  unsigned int numHorizontal =
      sfcgal_primitive_parameter_int(sphere, "num_horizontal");
  BOOST_CHECK_EQUAL(numHorizontal, 32);

  sfcgal_primitive_set_parameter_int(sphere, "num_horizontal", 48);

  double newNumHorizontal =
      sfcgal_primitive_parameter_int(sphere, "num_horizontal");
  BOOST_CHECK_EQUAL(newNumHorizontal, 48);

  // direction parameter
  double *direction = sfcgal_primitive_parameter_vector(sphere, "direction");
  BOOST_CHECK(direction != nullptr);
  BOOST_CHECK_CLOSE(direction[0], 0.0, 1e-6);
  BOOST_CHECK_CLOSE(direction[1], 0.0, 1e-6);
  BOOST_CHECK_CLOSE(direction[2], 1.0, 1e-6);
  sfcgal_free_buffer(direction);

  std::array<double, 3> expectedDirection{1.1, 2.1, 3.1};
  sfcgal_primitive_set_parameter_vector(sphere, "direction",
                                        expectedDirection.data());

  double *newDirection = sfcgal_primitive_parameter_vector(sphere, "direction");
  BOOST_CHECK(newDirection != nullptr);
  BOOST_CHECK_CLOSE(newDirection[0], expectedDirection[0], 1e-6);
  BOOST_CHECK_CLOSE(newDirection[1], expectedDirection[1], 1e-6);
  BOOST_CHECK_CLOSE(newDirection[2], expectedDirection[2], 1e-6);
  sfcgal_free_buffer(newDirection);

  // check clone
  sfcgal_primitive_t *sphere2 = sfcgal_primitive_clone(sphere);
  BOOST_CHECK(sfcgal_primitive_is_almost_equals(sphere, sphere2, 0.0));
  sfcgal_primitive_delete(sphere2);

  // check polyhedral conversion
  sfcgal_primitive_set_parameter_int(sphere, "num_vertical", 4);
  sfcgal_primitive_set_parameter_int(sphere, "num_horizontal", 4);
  sfcgal_geometry_t *poly = sfcgal_primitive_as_polyhedral_surface(sphere);
  char              *wkbApi;
  size_t             wkbLen;
  sfcgal_geometry_as_text_decim(poly, 0, &wkbApi, &wkbLen);
  std::string strApi(wkbApi, wkbLen);
  sfcgal_geometry_delete(poly);

  BOOST_CHECK_EQUAL(
      "POLYHEDRALSURFACE Z (((8 15 22,0 0 27,-9 19 17,8 15 22)),((8 15 22,21 3 "
      "17,0 0 27,8 15 22)),((8 15 22,12 23 6,21 3 17,8 15 22)),((8 15 22,-9 19 "
      "17,12 23 6,8 15 22)),((-9 19 17,0 0 27,-8 -15 18,-9 19 17)),((-9 19 "
      "17,-8 -15 18,-20 13 3,-9 19 17)),((0 0 27,21 3 17,22 -9 3,0 0 27)),((0 "
      "0 27,22 -9 3,-8 -15 18,0 0 27)),((21 3 17,12 23 6,10 19 -12,21 3 "
      "17)),((21 3 17,10 19 -12,22 -9 3,21 3 17)),((12 23 6,-9 19 17,-20 13 "
      "3,12 23 6)),((12 23 6,-20 13 3,10 19 -12,12 23 6)),((-20 13 3,-8 -15 "
      "18,-10 -19 0,-20 13 3)),((-20 13 3,-10 -19 0,-19 0 -11,-20 13 3)),((-8 "
      "-15 18,22 -9 3,11 -15 -11,-8 -15 18)),((-8 -15 18,11 -15 -11,-10 -19 "
      "0,-8 -15 18)),((22 -9 3,10 19 -12,2 5 -21,22 -9 3)),((22 -9 3,2 5 "
      "-21,11 -15 -11,22 -9 3)),((10 19 -12,-20 13 3,-19 0 -11,10 19 "
      "-12)),((10 19 -12,-19 0 -11,2 5 -21,10 19 -12)),((-6 -11 -16,-19 0 "
      "-11,-10 -19 0,-6 -11 -16)),((-6 -11 -16,-10 -19 0,11 -15 -11,-6 -11 "
      "-16)),((-6 -11 -16,11 -15 -11,2 5 -21,-6 -11 -16)),((-6 -11 -16,2 5 "
      "-21,-19 0 -11,-6 -11 -16)))",
      strApi);
  sfcgal_free_buffer(wkbApi);

  sfcgal_primitive_set_parameter_double(sphere, "radius", 2.0);
  // check volume
  BOOST_CHECK_CLOSE(sfcgal_primitive_volume(sphere, false), 33.51, 0.01);
  BOOST_CHECK_CLOSE(sfcgal_primitive_volume(sphere, true), 33.51, 0.01);

  // check area
  BOOST_CHECK_CLOSE(sfcgal_primitive_area(sphere, false), 50.265, 0.01);
  BOOST_CHECK_CLOSE(sfcgal_primitive_area(sphere, true), 50.265, 0.01);

  sfcgal_primitive_delete(sphere);
}

BOOST_AUTO_TEST_CASE(testCylinderTest)
{
  sfcgal_set_error_handlers(printf, on_error);

  sfcgal_primitive_t *cylinder = sfcgal_primitive_create(SFCGAL_TYPE_CYLINDER);

  // base_center parameter
  double *baseCenter =
      sfcgal_primitive_parameter_point(cylinder, "base_center");
  BOOST_CHECK(baseCenter != nullptr);
  BOOST_CHECK_CLOSE(baseCenter[0], 0.0, 1e-6);
  BOOST_CHECK_CLOSE(baseCenter[1], 0.0, 1e-6);
  BOOST_CHECK_CLOSE(baseCenter[2], 0.0, 1e-6);
  sfcgal_free_buffer(baseCenter);

  std::array<double, 3> expectedBaseCenter{1.0, 2.0, 3.0};
  sfcgal_primitive_set_parameter_point(cylinder, "base_center",
                                       expectedBaseCenter.data());

  double *newBaseCenter =
      sfcgal_primitive_parameter_point(cylinder, "base_center");
  BOOST_CHECK(newBaseCenter != nullptr);
  BOOST_CHECK_CLOSE(newBaseCenter[0], expectedBaseCenter[0], 1e-6);
  BOOST_CHECK_CLOSE(newBaseCenter[1], expectedBaseCenter[1], 1e-6);
  BOOST_CHECK_CLOSE(newBaseCenter[2], expectedBaseCenter[2], 1e-6);
  sfcgal_free_buffer(newBaseCenter);

  // axis parameter
  double *axis = sfcgal_primitive_parameter_vector(cylinder, "axis");
  BOOST_CHECK(axis != nullptr);
  BOOST_CHECK_CLOSE(axis[0], 0.0, 1e-6);
  BOOST_CHECK_CLOSE(axis[1], 0.0, 1e-6);
  BOOST_CHECK_CLOSE(axis[2], 1.0, 1e-6);
  sfcgal_free_buffer(axis);

  std::array<double, 3> expectedAxis{1.1, 2.1, 3.1};
  sfcgal_primitive_set_parameter_vector(cylinder, "axis", expectedAxis.data());

  double *newAxis = sfcgal_primitive_parameter_vector(cylinder, "axis");
  BOOST_CHECK(newAxis != nullptr);
  BOOST_CHECK_CLOSE(newAxis[0], expectedAxis[0], 1e-6);
  BOOST_CHECK_CLOSE(newAxis[1], expectedAxis[1], 1e-6);
  BOOST_CHECK_CLOSE(newAxis[2], expectedAxis[2], 1e-6);
  sfcgal_free_buffer(newAxis);

  // radius parameter
  double radius = sfcgal_primitive_parameter_double(cylinder, "radius");
  BOOST_CHECK_CLOSE(radius, 1.0, 1e-6);

  sfcgal_primitive_set_parameter_double(cylinder, "radius", 24.2);

  double newRadius = sfcgal_primitive_parameter_double(cylinder, "radius");
  BOOST_CHECK_CLOSE(newRadius, 24.2, 1e-6);

  // there is no parameter called min_radius
  sfcgal_primitive_set_parameter_double(cylinder, "min_radius", 24.2);
  BOOST_CHECK(hasError);
  hasError = false;

  // base_center is a point not a double
  sfcgal_primitive_set_parameter_double(cylinder, "base_center", 1.2);
  BOOST_CHECK(hasError);
  hasError = false;

  // radius cannot be negative
  sfcgal_primitive_set_parameter_double(cylinder, "radius", -3.7);
  BOOST_CHECK(hasError);
  hasError = false;

  // height parameter
  double height = sfcgal_primitive_parameter_double(cylinder, "height");
  BOOST_CHECK_CLOSE(height, 1.0, 1e-6);

  sfcgal_primitive_set_parameter_double(cylinder, "height", 24.2);

  double newHeight = sfcgal_primitive_parameter_double(cylinder, "height");
  BOOST_CHECK_CLOSE(newHeight, 24.2, 1e-6);

  // num vertical parameter
  unsigned int numRadial =
      sfcgal_primitive_parameter_int(cylinder, "num_radial");
  BOOST_CHECK_EQUAL(numRadial, 32);

  sfcgal_primitive_set_parameter_int(cylinder, "num_radial", 36);

  double newNumRadial = sfcgal_primitive_parameter_int(cylinder, "num_radial");
  BOOST_CHECK_EQUAL(newNumRadial, 36);

  // check clone
  sfcgal_primitive_t *cylinder2 = sfcgal_primitive_clone(cylinder);
  BOOST_CHECK(sfcgal_primitive_is_almost_equals(cylinder, cylinder2, 0.0));
  sfcgal_primitive_delete(cylinder2);

  // polyhedral surface generation
  sfcgal_geometry_t *polySurface =
      sfcgal_primitive_as_polyhedral_surface(cylinder);

  BOOST_REQUIRE(sfcgal_geometry_is_valid(polySurface));

  char  *wktApi;
  size_t wkbLen;
  sfcgal_geometry_as_text_decim(polySurface, 2, &wktApi, &wkbLen);
  const std::string strApi(wktApi, wkbLen);
  sfcgal_free_buffer(wktApi);

  std::string expectedWkt(SFCGAL_TEST_DIRECTORY);
  expectedWkt += "/data/cylinder_expected.wkt";
  std::ifstream efs(expectedWkt.c_str());
  BOOST_REQUIRE(efs.good());
  std::getline(efs, expectedWkt);

  BOOST_CHECK_EQUAL(strApi, expectedWkt);

  sfcgal_geometry_delete(polySurface);

  sfcgal_primitive_set_parameter_double(cylinder, "radius", 2.0);
  sfcgal_primitive_set_parameter_double(cylinder, "height", 10.0);
  // check volume
  BOOST_CHECK_CLOSE(sfcgal_primitive_volume(cylinder, false), 125.66, 0.01);
  BOOST_CHECK_CLOSE(sfcgal_primitive_volume(cylinder, true), 125.66, 0.01);

  // check area
  BOOST_CHECK_CLOSE(sfcgal_primitive_area(cylinder, false), 150.79, 0.01);
  BOOST_CHECK_CLOSE(sfcgal_primitive_area(cylinder, true), 150.79, 0.01);

  // checks parameter list
  char  *params;
  size_t paramsLen;
  sfcgal_primitive_parameters(cylinder, &params, &paramsLen);
  std::string paramsStr(params, paramsLen);
  std::string expectedStr(
      R"([{"name":"num_radial","type":"int"},{"name":"height","type":"double"},{"name":"radius","type":"double"},{"name":"axis","type":"vector3"},{"name":"base_center","type":"point3"}])");
  BOOST_CHECK(compare_json(paramsStr, expectedStr));
  sfcgal_free_buffer(params);

  // checks get parameter value generic
  char  *param;
  size_t paramLen;
  sfcgal_primitive_parameter(cylinder, "num_radial", &param, &paramLen);
  std::string paramStr = std::string(param, paramLen);
  std::string expectedParamStr(
      R"({"name":"num_radial","type":"int","value":36})");
  BOOST_CHECK(compare_json(paramStr, expectedParamStr));

  sfcgal_free_buffer(param);

  sfcgal_primitive_parameter(cylinder, "axis", &param, &paramLen);
  paramStr = std::string(param, paramLen);
  expectedParamStr =
      std::string(R"({"name":"axis","type":"vector3","value":[1.1,2.1,3.1]})");
  BOOST_CHECK(compare_json(paramStr, expectedParamStr));
  sfcgal_free_buffer(param);

  // checks set parameter value generic
  sfcgal_primitive_set_parameter(
      cylinder, "num_radial",
      R"({"name":"num_radial","type":"int","value":250})");
  unsigned int numRadial0 =
      sfcgal_primitive_parameter_int(cylinder, "num_radial");
  BOOST_CHECK_EQUAL(numRadial0, 250);

  sfcgal_primitive_set_parameter(
      cylinder, "axis",
      R"({"name":"axis","type":"vector3","value":[2.0,3.0,4.0]})");
  newAxis = sfcgal_primitive_parameter_vector(cylinder, "axis");
  BOOST_CHECK(newAxis != nullptr);
  BOOST_CHECK(!hasError);
  BOOST_CHECK_CLOSE(newAxis[0], 2.0, 1e-6);
  BOOST_CHECK_CLOSE(newAxis[1], 3.0, 1e-6);
  BOOST_CHECK_CLOSE(newAxis[2], 4.0, 1e-6);
  sfcgal_free_buffer(newAxis);

  sfcgal_primitive_delete(cylinder);
}

BOOST_AUTO_TEST_CASE(testTorusTest)
{
  sfcgal_set_error_handlers(printf, on_error);

  sfcgal_primitive_t *torus = sfcgal_primitive_create(SFCGAL_TYPE_TORUS);

  // checks parameter list
  char  *params;
  size_t paramsLen;
  sfcgal_primitive_parameters(torus, &params, &paramsLen);
  std::string paramsStr(params, paramsLen);
  std::string expectedStr(
      R"([{"name":"tube_num_radial","type":"int"},{"name":"main_num_radial","type":"int"},{"name":"tube_radius","type":"double"},{"name":"main_radius","type":"double"}])");
  BOOST_CHECK(compare_json(paramsStr, expectedStr));
  sfcgal_free_buffer(params);

  // main_radius parameter
  double mainRadius = sfcgal_primitive_parameter_double(torus, "main_radius");
  BOOST_CHECK_CLOSE(mainRadius, 10.0, 1e-6);

  sfcgal_primitive_set_parameter_double(torus, "main_radius", 11.2);
  double newMainRadius =
      sfcgal_primitive_parameter_double(torus, "main_radius");
  BOOST_CHECK_CLOSE(newMainRadius, 11.2, 1e-6);

  // there is no parameter called min_radius
  sfcgal_primitive_set_parameter_double(torus, "min_radius", 24.2);
  BOOST_CHECK(hasError);
  hasError = false;

  // main_num_radial is an unsigned int not a double
  sfcgal_primitive_set_parameter_double(torus, "main_num_radial", 1.2);
  BOOST_CHECK(hasError);
  hasError = false;

  // tube_radius parameter
  double tubeRadius = sfcgal_primitive_parameter_double(torus, "tube_radius");
  BOOST_CHECK_CLOSE(tubeRadius, 2.0, 1e-6);

  sfcgal_primitive_set_parameter_double(torus, "tube_radius", 5.2);
  double newTubeRadius =
      sfcgal_primitive_parameter_double(torus, "tube_radius");
  BOOST_CHECK_CLOSE(newTubeRadius, 5.2, 1e-6);

  // main num radial parameter
  unsigned int mainNumRadial =
      sfcgal_primitive_parameter_int(torus, "main_num_radial");
  BOOST_CHECK_EQUAL(mainNumRadial, 32);
  sfcgal_primitive_set_parameter_int(torus, "main_num_radial", 36);

  double newMainNumRadial =
      sfcgal_primitive_parameter_int(torus, "main_num_radial");
  BOOST_CHECK_EQUAL(newMainNumRadial, 36);

  // tube num radial parameter
  unsigned int tubeNumRadial =
      sfcgal_primitive_parameter_int(torus, "tube_num_radial");
  BOOST_CHECK_EQUAL(tubeNumRadial, 16);
  sfcgal_primitive_set_parameter_int(torus, "tube_num_radial", 18);

  double newTubeNumRadial =
      sfcgal_primitive_parameter_int(torus, "tube_num_radial");
  BOOST_CHECK_EQUAL(newTubeNumRadial, 18);

  // check clone
  sfcgal_primitive_t *torus2 = sfcgal_primitive_clone(torus);
  BOOST_CHECK(sfcgal_primitive_is_almost_equals(torus, torus2, 0.0));
  sfcgal_primitive_delete(torus2);

  // check polyhedral conversion
  sfcgal_primitive_set_parameter_int(torus, "main_num_radial", 4);
  sfcgal_primitive_set_parameter_int(torus, "tube_num_radial", 4);
  sfcgal_geometry_t *poly = sfcgal_primitive_as_polyhedral_surface(torus);
  char              *wkbApi;
  size_t             wkbLen;
  sfcgal_geometry_as_text_decim(poly, 0, &wkbApi, &wkbLen);
  std::string strApi(wkbApi, wkbLen);
  sfcgal_geometry_delete(poly);

  BOOST_CHECK_EQUAL(
      "POLYHEDRALSURFACE Z (((16 0 0,0 16 0,0 11 5,11 0 5,16 0 0)),((11 0 5,0 "
      "11 5,0 6 0,6 0 0,11 0 5)),((6 0 0,0 6 0,0 11 -5,11 0 -5,6 0 0)),((11 0 "
      "-5,0 11 -5,0 16 0,16 0 0,11 0 -5)),((0 16 0,-16 0 0,-11 0 5,0 11 5,0 16 "
      "0)),((0 11 5,-11 0 5,-6 0 0,0 6 0,0 11 5)),((0 6 0,-6 0 0,-11 0 -5,0 11 "
      "-5,0 6 0)),((0 11 -5,-11 0 -5,-16 0 0,0 16 0,0 11 -5)),((-16 0 0,0 -16 "
      "0,0 -11 5,-11 0 5,-16 0 0)),((-11 0 5,0 -11 5,0 -6 0,-6 0 0,-11 0 "
      "5)),((-6 0 0,0 -6 0,0 -11 -5,-11 0 -5,-6 0 0)),((-11 0 -5,0 -11 -5,0 "
      "-16 0,-16 0 0,-11 0 -5)),((0 -16 0,16 0 0,11 0 5,0 -11 5,0 -16 0)),((0 "
      "-11 5,11 0 5,6 0 0,0 -6 0,0 -11 5)),((0 -6 0,6 0 0,11 0 -5,0 -11 -5,0 "
      "-6 0)),((0 -11 -5,11 0 -5,16 0 0,0 -16 0,0 -11 -5)))",
      strApi);
  sfcgal_free_buffer(wkbApi);

  sfcgal_primitive_set_parameter_double(torus, "main_radius", 2.0);
  sfcgal_primitive_set_parameter_double(torus, "tube_radius", 1.0);
  // check volume
  BOOST_CHECK_CLOSE(sfcgal_primitive_volume(torus, false), 221.07, 0.01);
  BOOST_CHECK_CLOSE(sfcgal_primitive_volume(torus, true), 221.07, 0.01);

  // check area
  BOOST_CHECK_CLOSE(sfcgal_primitive_area(torus, false), 442.15, 0.01);
  BOOST_CHECK_CLOSE(sfcgal_primitive_area(torus, true), 442.15, 0.01);

  sfcgal_primitive_delete(torus);
}

BOOST_AUTO_TEST_CASE(testBoxTest)
{
  sfcgal_set_error_handlers(printf, on_error);

  sfcgal_primitive_t *box = sfcgal_primitive_create(SFCGAL_TYPE_BOX);

  // checks parameter list
  char  *params;
  size_t paramsLen;
  sfcgal_primitive_parameters(box, &params, &paramsLen);
  std::string paramsStr(params, paramsLen);
  std::string expectedParamsStr(
      R"([{"name":"z_extent","type":"double"},{"name":"y_extent","type":"double"},{"name":"x_extent","type":"double"}])");
  BOOST_CHECK(compare_json(paramsStr, expectedParamsStr));
  sfcgal_free_buffer(params);

  // x_extent parameter
  double xExtent = sfcgal_primitive_parameter_double(box, "x_extent");
  BOOST_CHECK_CLOSE(xExtent, 1.0, 1e-6);

  sfcgal_primitive_set_parameter_double(box, "x_extent", 11.2);
  double newXExtent = sfcgal_primitive_parameter_double(box, "x_extent");
  BOOST_CHECK_CLOSE(newXExtent, 11.2, 1e-6);

  // y_extent parameter
  double yExtent = sfcgal_primitive_parameter_double(box, "y_extent");
  BOOST_CHECK_CLOSE(yExtent, 1.0, 1e-6);

  sfcgal_primitive_set_parameter_double(box, "y_extent", 5.2);
  double newYExtent = sfcgal_primitive_parameter_double(box, "y_extent");
  BOOST_CHECK_CLOSE(newYExtent, 5.2, 1e-6);

  // z_extent parameter
  double zExtent = sfcgal_primitive_parameter_double(box, "z_extent");
  BOOST_CHECK_CLOSE(zExtent, 1.0, 1e-6);

  sfcgal_primitive_set_parameter_double(box, "z_extent", 1.2);
  double newZExtent = sfcgal_primitive_parameter_double(box, "z_extent");
  BOOST_CHECK_CLOSE(newZExtent, 1.2, 1e-6);

  // there is no parameter called radius
  sfcgal_primitive_set_parameter_double(box, "radius", 24.2);
  BOOST_CHECK(hasError);
  hasError = false;

  // x_extent is a double nont an unsigned int
  sfcgal_primitive_set_parameter_int(box, "x_extent", 5);
  BOOST_CHECK(hasError);
  hasError = false;

  // check clone
  sfcgal_primitive_t *box2 = sfcgal_primitive_clone(box);
  BOOST_CHECK(sfcgal_primitive_is_almost_equals(box, box2, 0.0));
  sfcgal_primitive_delete(box2);

  // check polyhedral conversion
  sfcgal_geometry_t *poly = sfcgal_primitive_as_polyhedral_surface(box);
  char              *wkbApi;
  size_t             wkbLen;
  sfcgal_geometry_as_text_decim(poly, 0, &wkbApi, &wkbLen);
  std::string strApi(wkbApi, wkbLen);
  sfcgal_geometry_delete(poly);

  BOOST_CHECK_EQUAL(
      "POLYHEDRALSURFACE Z (((0 0 0,0 5 0,11 5 0,11 0 0,0 0 0)),((0 0 1,11 0 "
      "1,11 5 1,0 5 1,0 0 1)),((0 0 0,11 0 0,11 0 1,0 0 1,0 0 0)),((0 5 0,0 5 "
      "1,11 5 1,11 5 0,0 5 0)),((11 0 0,11 5 0,11 5 1,11 0 1,11 0 0)),((0 0 "
      "0,0 0 1,0 5 1,0 5 0,0 0 0)))",
      strApi);
  sfcgal_free_buffer(wkbApi);

  sfcgal_primitive_set_parameter_double(box, "x_extent", 10.0);
  sfcgal_primitive_set_parameter_double(box, "y_extent", 5.0);
  sfcgal_primitive_set_parameter_double(box, "z_extent", 2.0);
  // check volume
  BOOST_CHECK_CLOSE(sfcgal_primitive_volume(box, false), 100.0, 0.01);
  BOOST_CHECK_CLOSE(sfcgal_primitive_volume(box, true), 100.0, 0.01);

  // check area
  BOOST_CHECK_CLOSE(sfcgal_primitive_area(box, false), 160.0, 0.01);
  BOOST_CHECK_CLOSE(sfcgal_primitive_area(box, true), 160.0, 0.01);

  sfcgal_primitive_delete(box);
}

BOOST_AUTO_TEST_CASE(testCubeTest)
{
  sfcgal_set_error_handlers(printf, on_error);

  sfcgal_primitive_t *cube = sfcgal_primitive_create(SFCGAL_TYPE_CUBE);

  // checks parameter list
  char  *params;
  size_t paramsLen;
  sfcgal_primitive_parameters(cube, &params, &paramsLen);
  std::string paramsStr(params, paramsLen);
  std::string expectedParamsStr(R"([{"name":"size","type":"double"}])");
  BOOST_CHECK(compare_json(paramsStr, expectedParamsStr));
  sfcgal_free_buffer(params);

  // size parameter
  double size = sfcgal_primitive_parameter_double(cube, "size");
  BOOST_CHECK_CLOSE(size, 1.0, 1e-6);

  sfcgal_primitive_set_parameter_double(cube, "size", 11.2);
  double newXExtent = sfcgal_primitive_parameter_double(cube, "size");
  BOOST_CHECK_CLOSE(newXExtent, 11.2, 1e-6);

  // there is no parameter called radius
  sfcgal_primitive_set_parameter_double(cube, "radius", 24.2);
  BOOST_CHECK(hasError);
  hasError = false;

  // size is a double nont an unsigned int
  sfcgal_primitive_set_parameter_int(cube, "size", 3);
  BOOST_CHECK(hasError);
  hasError = false;

  // size cannot be negative
  sfcgal_primitive_set_parameter_double(cube, "size", -2.2);
  BOOST_CHECK(hasError);
  hasError = false;

  // check clone
  sfcgal_primitive_t *cube2 = sfcgal_primitive_clone(cube);
  BOOST_CHECK(sfcgal_primitive_is_almost_equals(cube, cube2, 0.0));
  sfcgal_primitive_delete(cube2);

  // check polyhedral conversion
  sfcgal_primitive_set_parameter_double(cube, "size", 2.0);
  sfcgal_geometry_t *poly = sfcgal_primitive_as_polyhedral_surface(cube);
  char              *wkbApi;
  size_t             wkbLen;
  sfcgal_geometry_as_text_decim(poly, 0, &wkbApi, &wkbLen);
  std::string strApi(wkbApi, wkbLen);
  sfcgal_geometry_delete(poly);

  BOOST_CHECK_EQUAL("POLYHEDRALSURFACE Z (((0 0 0,0 2 0,2 2 0,2 0 0,0 0 0)),"
                    "((0 0 2,2 0 2,2 2 2,0 2 2,0 0 2)),"
                    "((0 0 0,2 0 0,2 0 2,0 0 2,0 0 0)),"
                    "((0 2 0,0 2 2,2 2 2,2 2 0,0 2 0)),"
                    "((2 0 0,2 2 0,2 2 2,2 0 2,2 0 0)),"
                    "((0 0 0,0 0 2,0 2 2,0 2 0,0 0 0)))",
                    strApi);
  sfcgal_free_buffer(wkbApi);

  // check volume
  BOOST_CHECK_CLOSE(sfcgal_primitive_volume(cube, false), 8.0, 0.01);
  BOOST_CHECK_CLOSE(sfcgal_primitive_volume(cube, true), 8.0, 0.01);

  // check area
  BOOST_CHECK_CLOSE(sfcgal_primitive_area(cube, false), 24.0, 0.01);
  BOOST_CHECK_CLOSE(sfcgal_primitive_area(cube, true), 24.0, 0.01);

  sfcgal_primitive_delete(cube);
}

BOOST_AUTO_TEST_CASE(testConeTest)
{
  sfcgal_set_error_handlers(printf, on_error);

  sfcgal_primitive_t *cone = sfcgal_primitive_create(SFCGAL_TYPE_CONE);

  // checks parameter list
  char  *params;
  size_t paramsLen;
  sfcgal_primitive_parameters(cone, &params, &paramsLen);
  std::string paramsStr(params, paramsLen);
  std::string expectedParamsStr(
      R"([{"name":"num_radial","type":"int"},{"name":"height","type":"double"},{"name":"top_radius","type":"double"},{"name":"bottom_radius","type":"double"}])");
  BOOST_CHECK(compare_json(paramsStr, expectedParamsStr));
  sfcgal_free_buffer(params);

  // bottom_radius parameter
  double bottomRadius =
      sfcgal_primitive_parameter_double(cone, "bottom_radius");
  BOOST_CHECK_CLOSE(bottomRadius, 1.0, 1e-6);

  sfcgal_primitive_set_parameter_double(cone, "bottom_radius", 14.31);
  double newBottomRadius =
      sfcgal_primitive_parameter_double(cone, "bottom_radius");
  BOOST_CHECK_CLOSE(newBottomRadius, 14.31, 1e-6);

  // there is no parameter called min_radius
  sfcgal_primitive_set_parameter_double(cone, "min_radius", 24.2);
  BOOST_CHECK(hasError);
  hasError = false;

  // num_radial is an unsigned int not a double
  sfcgal_primitive_set_parameter_double(cone, "num_radial", 1.2);
  BOOST_CHECK(hasError);
  hasError = false;

  // top_radius parameter
  double topRadius = sfcgal_primitive_parameter_double(cone, "top_radius");
  BOOST_CHECK_CLOSE(topRadius, 0.0, 1e-6);

  sfcgal_primitive_set_parameter_double(cone, "top_radius", 0.3);
  double newTubeRadius = sfcgal_primitive_parameter_double(cone, "top_radius");
  BOOST_CHECK_CLOSE(newTubeRadius, 0.3, 1e-6);

  // top_radius cannot be negative
  sfcgal_primitive_set_parameter_double(cone, "top_radius", -2.3);
  BOOST_CHECK(hasError);
  hasError = false;

  // num radial parameter
  unsigned int numRadial = sfcgal_primitive_parameter_int(cone, "num_radial");
  BOOST_CHECK_EQUAL(numRadial, 32);
  sfcgal_primitive_set_parameter_int(cone, "num_radial", 36);

  double newNumRadial = sfcgal_primitive_parameter_int(cone, "num_radial");
  BOOST_CHECK_EQUAL(newNumRadial, 36);

  // check clone
  sfcgal_primitive_t *cone2 = sfcgal_primitive_clone(cone);
  BOOST_CHECK(sfcgal_primitive_is_almost_equals(cone, cone2, 0.0));
  sfcgal_primitive_delete(cone2);

  // check polyhedral conversion
  sfcgal_primitive_set_parameter_int(cone, "num_radial", 4);
  sfcgal_geometry_t *poly = sfcgal_primitive_as_polyhedral_surface(cone);
  char              *wkbApi;
  size_t             wkbLen;
  sfcgal_geometry_as_text_decim(poly, 0, &wkbApi, &wkbLen);
  std::string strApi(wkbApi, wkbLen);
  sfcgal_geometry_delete(poly);

  // check
  BOOST_CHECK_EQUAL(
      "POLYHEDRALSURFACE Z (((14 0 0,0 -14 0,-14 0 0,0 14 0,14 0 0)),((0 0 1,0 "
      "0 1,0 0 1,0 0 1,0 0 1)),((14 0 0,0 0 1,0 0 1,0 -14 0,14 0 0)),((0 -14 "
      "0,0 0 1,0 0 1,-14 0 0,0 -14 0)),((-14 0 0,0 0 1,0 0 1,0 14 0,-14 0 "
      "0)),((0 14 0,0 0 1,0 0 1,14 0 0,0 14 0)))",
      strApi);
  sfcgal_free_buffer(wkbApi);

  sfcgal_primitive_set_parameter_double(cone, "bottom_radius", 2.0);
  sfcgal_primitive_set_parameter_double(cone, "top_radius", 0.01);
  sfcgal_primitive_set_parameter_double(cone, "height", 10.0);
  // check volume
  BOOST_CHECK_CLOSE(sfcgal_primitive_volume(cone, false), 42.098, 0.01);
  BOOST_CHECK_CLOSE(sfcgal_primitive_volume(cone, true), 26.80, 0.01);

  // check area
  BOOST_CHECK_CLOSE(sfcgal_primitive_area(cone, false), 76.95, 0.01);
  BOOST_CHECK_CLOSE(sfcgal_primitive_area(cone, true), 65.41, 0.01);

  sfcgal_primitive_delete(cone);
}

BOOST_AUTO_TEST_SUITE_END()
