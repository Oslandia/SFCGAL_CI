// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/capi/sfcgal_c.h"
#include "SFCGAL/Solid.h"
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
                                            point->clone(), 1);
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
                                        simplePolygon->clone(), 1);
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
  sfcgal_triangulated_surface_set_patch_n(tin.get(), simpleTriangle->clone(),
                                          1);
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

BOOST_AUTO_TEST_SUITE_END()
