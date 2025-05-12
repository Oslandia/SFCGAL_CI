// Copyright (c) 2025-2025, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/Geometry.h"
#include "SFCGAL/Kernel.h"

#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/io/wkt.h"

#include <boost/test/unit_test.hpp>
#include <memory>
#include <string>
using namespace boost::unit_test;

using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_PolyhedralSurfaceTest)

BOOST_AUTO_TEST_CASE(setPatchNTest)
{
  std::unique_ptr<Geometry> emptyGeom(io::readWkt("POLYHEDRALSURFACE EMPTY"));
  BOOST_CHECK(emptyGeom->is<PolyhedralSurface>());
  BOOST_CHECK(emptyGeom->isEmpty());
  BOOST_CHECK_EQUAL(emptyGeom->numGeometries(), 0);
  BOOST_CHECK_EQUAL(emptyGeom->as<PolyhedralSurface>().numPatches(), 0);
  BOOST_CHECK_EQUAL(emptyGeom->geometryN(0).asText(), emptyGeom->asText());

  std::string const polyhedralStr = "POLYHEDRALSURFACE Z ("
                                    "((0 0 0, 10 0 0, 10 10 0, 0 10 0, 0 0 0)),"
                                    "((0 0 0, 10 0 0, 5 0 5, 0 0 0)),"
                                    "((0 0 0, 0 10 0, 5 5 5, 0 0 0))"
                                    ")";

  std::unique_ptr<Geometry> geom(io::readWkt(polyhedralStr));
  BOOST_CHECK(!geom->isEmpty());
  BOOST_CHECK_EQUAL(geom->numGeometries(), 1);
  BOOST_CHECK_EQUAL(geom->as<PolyhedralSurface>().numPatches(), 3);
  BOOST_CHECK_EQUAL(geom->as<PolyhedralSurface>().patchN(0).asText(0),
                    "POLYGON Z ((0 0 0,10 0 0,10 10 0,0 10 0,0 0 0))");
  BOOST_CHECK_EQUAL(geom->as<PolyhedralSurface>().patchN(1).asText(0),
                    "POLYGON Z ((0 0 0,10 0 0,5 0 5,0 0 0))");
  BOOST_CHECK_EQUAL(geom->as<PolyhedralSurface>().patchN(2).asText(0),
                    "POLYGON Z ((0 0 0,0 10 0,5 5 5,0 0 0))");

  // set new Polygon at index 1 from a Geometry object
  std::string const newPolygonStr =
      "POLYGON Z ((0 0 0, 0 5 0, 5 5 5, 5 0 5, 0 0 0))";
  std::unique_ptr<Geometry> newGeom(io::readWkt(newPolygonStr));
  geom->as<PolyhedralSurface>().setPatchN(newGeom->clone(), 1);

  BOOST_CHECK_EQUAL(geom->numGeometries(), 1);
  BOOST_CHECK_EQUAL(geom->as<PolyhedralSurface>().numPatches(), 3);
  BOOST_CHECK_EQUAL(geom->as<PolyhedralSurface>().patchN(0).asText(0),
                    "POLYGON Z ((0 0 0,10 0 0,10 10 0,0 10 0,0 0 0))");
  BOOST_CHECK_EQUAL(geom->as<PolyhedralSurface>().patchN(1).asText(),
                    newGeom->asText());
  BOOST_CHECK_EQUAL(geom->as<PolyhedralSurface>().patchN(2).asText(0),
                    "POLYGON Z ((0 0 0,0 10 0,5 5 5,0 0 0))");

  // set New Polygon at index 2 from a Polygon
  std::string const newPolygonStr2 =
      "POLYGON Z ((0 0 0, 0 15 0, 15 15 15, 15 0 15, 0 0 0))";
  std::unique_ptr<Geometry> newGeom2(io::readWkt(newPolygonStr2));
  Polygon *newPolygon2 = dynamic_cast<Polygon *>(newGeom2.get());
  BOOST_CHECK(newPolygon2);
  geom->as<PolyhedralSurface>().setPatchN(newPolygon2->clone(), 2);

  BOOST_CHECK_EQUAL(geom->numGeometries(), 1);
  BOOST_CHECK_EQUAL(geom->as<PolyhedralSurface>().numPatches(), 3);
  BOOST_CHECK_EQUAL(geom->as<PolyhedralSurface>().patchN(0).asText(0),
                    "POLYGON Z ((0 0 0,10 0 0,10 10 0,0 10 0,0 0 0))");
  BOOST_CHECK_EQUAL(geom->as<PolyhedralSurface>().patchN(1).asText(),
                    newGeom->asText());
  BOOST_CHECK_EQUAL(geom->as<PolyhedralSurface>().patchN(2).asText(),
                    newGeom2->asText());
}

BOOST_AUTO_TEST_CASE(dropZMTest)
{
  std::unique_ptr<Geometry> emptyGeom(io::readWkt("POLYHEDRALSURFACE EMPTY"));
  BOOST_CHECK(!emptyGeom->is3D());
  BOOST_CHECK(!emptyGeom->isMeasured());
  BOOST_CHECK(!emptyGeom->dropZ());
  BOOST_CHECK(!emptyGeom->dropM());

  std::string const polyhedral3DStr =
      "POLYHEDRALSURFACE Z ("
      "((0 0 0, 10 0 0, 10 10 0, 0 10 0, 0 0 0)),"
      "((0 0 0, 10 0 0, 5 0 5, 0 0 0)),"
      "((0 0 0, 0 10 0, 5 5 5, 0 0 0))"
      ")";

  std::unique_ptr<Geometry> geom3D(io::readWkt(polyhedral3DStr));
  BOOST_CHECK(geom3D->is3D());
  BOOST_CHECK(!geom3D->isMeasured());
  BOOST_CHECK(!geom3D->dropM());
  BOOST_CHECK(geom3D->dropZ());
  BOOST_CHECK_EQUAL(geom3D->asText(0),
                    "POLYHEDRALSURFACE (((0 0,10 0,10 10,0 10,0 0)),"
                    "((0 0,10 0,5 0,0 0)),((0 0,0 10,5 5,0 0)))");
  BOOST_CHECK(!geom3D->is3D());
  BOOST_CHECK(!geom3D->dropZ());
  BOOST_CHECK(!geom3D->isMeasured());
  BOOST_CHECK(!geom3D->dropM());

  std::string const polyhedralMStr =
      "POLYHEDRALSURFACE M (((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)))";
  std::unique_ptr<Geometry> geomM(io::readWkt(polyhedralMStr));
  BOOST_CHECK(!geomM->is3D());
  BOOST_CHECK(geomM->isMeasured());
  BOOST_CHECK(!geomM->dropZ());
  BOOST_CHECK(geomM->dropM());
  BOOST_CHECK_EQUAL(geomM->asText(0),
                    "POLYHEDRALSURFACE (((0 0,0 1,1 1,1 0,0 0)))");
  BOOST_CHECK(!geomM->is3D());
  BOOST_CHECK(!geomM->isMeasured());
  BOOST_CHECK(!geomM->dropZ());
  BOOST_CHECK(!geomM->dropM());

  std::string const polyhedralZMStr =
      "POLYHEDRALSURFACE ZM (((0 0 0 1,0 1 0 2,1 1 0 3,1 0 0 4,0 0 0 1)))";
  std::unique_ptr<Geometry> geomZM(io::readWkt(polyhedralZMStr));
  BOOST_CHECK(geomZM->is3D());
  BOOST_CHECK(geomZM->isMeasured());
  BOOST_CHECK(geomZM->dropM());
  BOOST_CHECK(geomZM->is3D());
  BOOST_CHECK(!geomZM->isMeasured());
  BOOST_CHECK_EQUAL(geomZM->asText(0),
                    "POLYHEDRALSURFACE Z (((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)))");
  BOOST_CHECK(!geomZM->dropM());
  BOOST_CHECK(geomZM->dropZ());
  BOOST_CHECK_EQUAL(geomZM->asText(0),
                    "POLYHEDRALSURFACE (((0 0,0 1,1 1,1 0,0 0)))");
  BOOST_CHECK(!geomZM->isMeasured());
  BOOST_CHECK(!geomZM->is3D());
  BOOST_CHECK(!geomZM->dropZ());
  BOOST_CHECK(!geomZM->dropM());
}

BOOST_AUTO_TEST_CASE(swapXYTest)
{
  std::unique_ptr<Geometry> emptyGeom(io::readWkt("POLYHEDRALSURFACE EMPTY"));
  BOOST_CHECK(emptyGeom->isEmpty());
  emptyGeom->swapXY();
  BOOST_CHECK(emptyGeom->isEmpty());

  std::string const polyhedral3DStr =
      "POLYHEDRALSURFACE Z ("
      "((0 0 0, 10 0 0, 10 10 0, 0 10 0, 0 0 0)),"
      "((0 0 0, 10 0 0, 5 0 5, 0 0 0)),"
      "((0 0 0, 0 10 0, 5 5 5, 0 0 0))"
      ")";
  std::unique_ptr<Geometry> geom3D(io::readWkt(polyhedral3DStr));
  geom3D->swapXY();
  BOOST_CHECK_EQUAL(geom3D->asText(0), "POLYHEDRALSURFACE Z "
                                       "(((0 0 0,0 10 0,10 10 0,10 0 0,0 0 0)),"
                                       "((0 0 0,0 10 0,0 5 5,0 0 0)),"
                                       "((0 0 0,10 0 0,5 5 5,0 0 0)))");

  std::string const polyhedralMStr =
      "POLYHEDRALSURFACE M (((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)))";
  std::unique_ptr<Geometry> geomM(io::readWkt(polyhedralMStr));
  geomM->swapXY();
  BOOST_CHECK_EQUAL(geomM->asText(0),
                    "POLYHEDRALSURFACE M (((0 0 0,1 0 0,1 1 0,0 1 0,0 0 0)))");

  std::string const polyhedralZMStr =
      "POLYHEDRALSURFACE ZM (((0 0 0 1,0 1 0 2,1 1 0 3,1 0 0 4,0 0 0 1)))";
  std::unique_ptr<Geometry> geomZM(io::readWkt(polyhedralZMStr));
  geomZM->swapXY();
  BOOST_CHECK_EQUAL(
      geomZM->asText(0),
      "POLYHEDRALSURFACE ZM (((0 0 0 1,1 0 0 2,1 1 0 3,0 1 0 4,0 0 0 1)))");
}

BOOST_AUTO_TEST_CASE(getCoordinateType)
{
  BOOST_CHECK_EQUAL(io::readWkt("POLYHEDRALSURFACE (((0 0, 1 0, 0 1, 0 0)))")
                        ->getCoordinateType(),
                    CoordinateType::COORDINATE_XY);
  BOOST_CHECK_EQUAL(
      io::readWkt("POLYHEDRALSURFACE Z(((0 0 1, 1 0 1, 0 1 1, 0 0 1)))")
          ->getCoordinateType(),
      CoordinateType::COORDINATE_XYZ);
  BOOST_CHECK_EQUAL(
      io::readWkt("POLYHEDRALSURFACE M(((0 0 1, 1 0 1, 0 1 1, 0 0 1)))")
          ->getCoordinateType(),
      CoordinateType::COORDINATE_XYM);
  BOOST_CHECK_EQUAL(
      io::readWkt(
          "POLYHEDRALSURFACE ZM(((0 0 1 2, 1 0 1 2, 0 1 1 2, 0 0 1 2)))")
          ->getCoordinateType(),
      CoordinateType::COORDINATE_XYZM);
}

BOOST_AUTO_TEST_SUITE_END()
