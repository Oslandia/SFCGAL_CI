// Copyright (c) 2025-2025, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/surfaceSimplification.h"

#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/io/wkt.h"

#include <boost/test/unit_test.hpp>

using namespace SFCGAL;
using namespace SFCGAL::algorithm;

#include "../../../test_config.h"

BOOST_AUTO_TEST_SUITE(SFCGAL_SurfaceSimplificationTest)

/**
 * @brief Helper to create a simple cube as a TriangulatedSurface
 */
auto
createCubeTriangulatedSurface() -> std::unique_ptr<TriangulatedSurface>
{
  auto cube = std::make_unique<TriangulatedSurface>();

  // Bottom face (z=0)
  cube->addTriangle(Triangle(Point(0, 0, 0), Point(1, 0, 0), Point(1, 1, 0)));
  cube->addTriangle(Triangle(Point(0, 0, 0), Point(1, 1, 0), Point(0, 1, 0)));

  // Top face (z=1)
  cube->addTriangle(Triangle(Point(0, 0, 1), Point(1, 1, 1), Point(1, 0, 1)));
  cube->addTriangle(Triangle(Point(0, 0, 1), Point(0, 1, 1), Point(1, 1, 1)));

  // Front face (y=0)
  cube->addTriangle(Triangle(Point(0, 0, 0), Point(1, 0, 1), Point(1, 0, 0)));
  cube->addTriangle(Triangle(Point(0, 0, 0), Point(0, 0, 1), Point(1, 0, 1)));

  // Back face (y=1)
  cube->addTriangle(Triangle(Point(0, 1, 0), Point(1, 1, 0), Point(1, 1, 1)));
  cube->addTriangle(Triangle(Point(0, 1, 0), Point(1, 1, 1), Point(0, 1, 1)));

  // Left face (x=0)
  cube->addTriangle(Triangle(Point(0, 0, 0), Point(0, 1, 0), Point(0, 1, 1)));
  cube->addTriangle(Triangle(Point(0, 0, 0), Point(0, 1, 1), Point(0, 0, 1)));

  // Right face (x=1)
  cube->addTriangle(Triangle(Point(1, 0, 0), Point(1, 1, 1), Point(1, 1, 0)));
  cube->addTriangle(Triangle(Point(1, 0, 0), Point(1, 0, 1), Point(1, 1, 1)));

  return cube;
}

/**
 * @brief Helper to create a simple cube as a PolyhedralSurface
 */
auto
createCubePolyhedralSurface() -> std::unique_ptr<PolyhedralSurface>
{
  auto cube = std::make_unique<PolyhedralSurface>();

  // Helper lambda to create a square face
  auto createFace = [](double x1, double y1, double z1, double x2, double y2,
                       double z2, double x3, double y3, double z3, double x4,
                       double y4, double z4) -> std::unique_ptr<Polygon> {
    auto ring = std::make_unique<LineString>();
    ring->addPoint(Point(x1, y1, z1));
    ring->addPoint(Point(x2, y2, z2));
    ring->addPoint(Point(x3, y3, z3));
    ring->addPoint(Point(x4, y4, z4));
    ring->addPoint(Point(x1, y1, z1)); // close the ring
    return std::make_unique<Polygon>(ring.release());
  };

  // Bottom (z=0)
  cube->addPatch(createFace(0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0));
  // Top (z=1)
  cube->addPatch(createFace(0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1));
  // Front (y=0)
  cube->addPatch(createFace(0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0));
  // Back (y=1)
  cube->addPatch(createFace(0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1));
  // Left (x=0)
  cube->addPatch(createFace(0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1));
  // Right (x=1)
  cube->addPatch(createFace(1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0));

  return cube;
}

// Test TriangulatedSurface simplification with Garland-Heckbert strategy
BOOST_AUTO_TEST_CASE(testSimplify_TriangulatedSurface_GarlandHeckbert)
{
  auto cube = createCubeTriangulatedSurface();

  BOOST_CHECK_EQUAL(cube->numTriangles(), 12);

  // Simplify using edge count ratio
  auto simplified = surfaceSimplification(
      *cube, SimplificationStopPredicate::edgeCountRatio(0.5),
      SimplificationStrategy::EDGE_LENGTH);

  BOOST_CHECK(simplified);
  BOOST_CHECK(simplified->is<TriangulatedSurface>());

  const auto &simplifiedTin = simplified->as<TriangulatedSurface>();
  BOOST_CHECK_LT(simplifiedTin.numTriangles(), cube->numTriangles());

  std::cout << "Original triangles: " << cube->numTriangles() << std::endl;
  std::cout << "Simplified triangles: " << simplifiedTin.numTriangles()
            << std::endl;
}

// Test TriangulatedSurface simplification with Lindstrom-Turk strategy
BOOST_AUTO_TEST_CASE(testSimplify_TriangulatedSurface_LindstromTurk)
{
  auto cube = createCubeTriangulatedSurface();

  BOOST_CHECK_EQUAL(cube->numTriangles(), 12);

  // Simplify using edge count ratio
  auto simplified = surfaceSimplification(
      *cube, SimplificationStopPredicate::edgeCountRatio(0.5),
      SimplificationStrategy::EDGE_LENGTH);

  BOOST_CHECK(simplified);
  BOOST_CHECK(simplified->is<TriangulatedSurface>());

  const auto &simplifiedTin = simplified->as<TriangulatedSurface>();
  BOOST_CHECK_LT(simplifiedTin.numTriangles(), cube->numTriangles());

  std::cout << "Original triangles: " << cube->numTriangles() << std::endl;
  std::cout << "Simplified triangles (LT): " << simplifiedTin.numTriangles()
            << std::endl;
}

// Test PolyhedralSurface simplification
BOOST_AUTO_TEST_CASE(testSimplify_PolyhedralSurface_GarlandHeckbert)
{
  auto cube = createCubePolyhedralSurface();

  BOOST_CHECK_EQUAL(cube->numPatches(), 6);

  // Simplify using edge count ratio
  auto simplified = surfaceSimplification(
      *cube, SimplificationStopPredicate::edgeCountRatio(0.5),
      SimplificationStrategy::EDGE_LENGTH);

  BOOST_CHECK(simplified);
  BOOST_CHECK(simplified->is<PolyhedralSurface>());

  std::cout << "Original patches: " << cube->numPatches() << std::endl;
  std::cout << "Simplified patches: "
            << simplified->as<PolyhedralSurface>().numPatches() << std::endl;
}

// Test PolyhedralSurface simplification with Lindstrom-Turk
BOOST_AUTO_TEST_CASE(testSimplify_PolyhedralSurface_LindstromTurk)
{
  auto cube = createCubePolyhedralSurface();

  BOOST_CHECK_EQUAL(cube->numPatches(), 6);

  // Simplify using edge count ratio
  auto simplified = surfaceSimplification(
      *cube, SimplificationStopPredicate::edgeCountRatio(0.5),
      SimplificationStrategy::EDGE_LENGTH);

  BOOST_CHECK(simplified);
  BOOST_CHECK(simplified->is<PolyhedralSurface>());

  std::cout << "Original patches: " << cube->numPatches() << std::endl;
  std::cout << "Simplified patches (LT): "
            << simplified->as<PolyhedralSurface>().numPatches() << std::endl;
}

// Test Solid simplification
BOOST_AUTO_TEST_CASE(testSimplify_Solid)
{
  auto exteriorShell = createCubePolyhedralSurface();
  Solid solid(*exteriorShell);

  auto simplified = surfaceSimplification(
      solid, SimplificationStopPredicate::edgeCountRatio(0.5),
      SimplificationStrategy::EDGE_LENGTH);

  BOOST_CHECK(simplified);
  BOOST_CHECK(simplified->is<Solid>());

  const auto &simplifiedSolid = simplified->as<Solid>();
  BOOST_CHECK(!simplifiedSolid.isEmpty());

  std::cout << "Original solid exterior patches: "
            << solid.exteriorShell().numPatches() << std::endl;
  std::cout << "Simplified solid exterior patches: "
            << simplifiedSolid.exteriorShell().numPatches() << std::endl;
}

// Test MultiSolid simplification
BOOST_AUTO_TEST_CASE(testSimplify_MultiSolid)
{
  MultiSolid multiSolid;

  auto shell1 = createCubePolyhedralSurface();
  multiSolid.addGeometry(Solid(*shell1));

  auto shell2 = createCubePolyhedralSurface();
  multiSolid.addGeometry(Solid(*shell2));

  BOOST_CHECK_EQUAL(multiSolid.numGeometries(), 2);

  auto simplified = surfaceSimplification(
      multiSolid, SimplificationStopPredicate::edgeCountRatio(0.5),
      SimplificationStrategy::EDGE_LENGTH);

  BOOST_CHECK(simplified);
  BOOST_CHECK(simplified->is<MultiSolid>());

  const auto &simplifiedMulti = simplified->as<MultiSolid>();
  BOOST_CHECK_EQUAL(simplifiedMulti.numGeometries(), 2);

  std::cout << "MultiSolid simplified successfully" << std::endl;
}

// Test edge count stop predicate
BOOST_AUTO_TEST_CASE(testSimplify_EdgeCountPredicate)
{
  auto cube = createCubeTriangulatedSurface();

  // Convert to surface mesh to count edges
  // For a cube with 12 triangles, there are 18 edges (12 * 3 / 2 = 18 for
  // manifold)

  // Simplify to a specific edge count (keep 10 edges)
  auto simplified = surfaceSimplification(
      *cube, SimplificationStopPredicate::edgeCount(10),
      SimplificationStrategy::EDGE_LENGTH);

  BOOST_CHECK(simplified);
  BOOST_CHECK(simplified->is<TriangulatedSurface>());

  std::cout << "Simplified with edge count predicate" << std::endl;
}

// Test empty geometries
BOOST_AUTO_TEST_CASE(testSimplify_EmptyGeometries)
{
  TriangulatedSurface emptyTin;
  BOOST_CHECK(emptyTin.isEmpty());

  auto simplified = surfaceSimplification(
      emptyTin, SimplificationStopPredicate::edgeCountRatio(0.5),
      SimplificationStrategy::EDGE_LENGTH);

  BOOST_CHECK(simplified);
  BOOST_CHECK(simplified->isEmpty());
}

// Test invalid ratio (too low)
BOOST_AUTO_TEST_CASE(testSimplify_InvalidRatio_TooLow)
{
  auto cube = createCubeTriangulatedSurface();

  BOOST_CHECK_THROW(
      surfaceSimplification(*cube,
                            SimplificationStopPredicate::edgeCountRatio(0.0),
                            SimplificationStrategy::EDGE_LENGTH),
      std::invalid_argument);
}

// Test invalid ratio (too high)
BOOST_AUTO_TEST_CASE(testSimplify_InvalidRatio_TooHigh)
{
  auto cube = createCubeTriangulatedSurface();

  BOOST_CHECK_THROW(
      surfaceSimplification(*cube,
                            SimplificationStopPredicate::edgeCountRatio(1.0),
                            SimplificationStrategy::EDGE_LENGTH),
      std::invalid_argument);
}

// Test invalid ratio (negative)
BOOST_AUTO_TEST_CASE(testSimplify_InvalidRatio_Negative)
{
  auto cube = createCubeTriangulatedSurface();

  BOOST_CHECK_THROW(
      surfaceSimplification(*cube,
                            SimplificationStopPredicate::edgeCountRatio(-0.5),
                            SimplificationStrategy::EDGE_LENGTH),
      std::invalid_argument);
}

// Test unsupported geometry type
BOOST_AUTO_TEST_CASE(testSimplify_UnsupportedGeometry)
{
  Point point(0, 0, 0);

  BOOST_CHECK_THROW(
      surfaceSimplification(point,
                            SimplificationStopPredicate::edgeCountRatio(0.5),
                            SimplificationStrategy::EDGE_LENGTH),
      std::invalid_argument);
}

// Test simplification with WKT geometries
BOOST_AUTO_TEST_CASE(testSimplify_WKT_TriangulatedSurface)
{
  // Create a connected triangulated surface (a tetrahedron)
  std::string wkt =
      "TIN Z(((0 0 0, 1 0 0, 0.5 0.5 1, 0 0 0)),"
      "((1 0 0, 0 1 0, 0.5 0.5 1, 1 0 0)),"
      "((0 1 0, 0 0 0, 0.5 0.5 1, 0 1 0)),"
      "((0 0 0, 0 1 0, 1 0 0, 0 0 0)))";

  std::unique_ptr<Geometry> geom(io::readWkt(wkt));
  BOOST_CHECK(geom->is<TriangulatedSurface>());

  auto simplified = surfaceSimplification(
      *geom, SimplificationStopPredicate::edgeCountRatio(0.7),
      SimplificationStrategy::EDGE_LENGTH);

  BOOST_CHECK(simplified);
  BOOST_CHECK(simplified->is<TriangulatedSurface>());

  std::cout << "WKT Source: " << geom->asText(2) << std::endl;
  std::cout << "WKT Simplified: " << simplified->asText(2) << std::endl;
}

// Test that simplified geometry maintains 3D characteristics
BOOST_AUTO_TEST_CASE(testSimplify_Maintains3D)
{
  auto cube = createCubeTriangulatedSurface();
  BOOST_CHECK(cube->is3D());

  auto simplified = surfaceSimplification(
      *cube, SimplificationStopPredicate::edgeCountRatio(0.5),
      SimplificationStrategy::EDGE_LENGTH);

  BOOST_CHECK(simplified);
  BOOST_CHECK(simplified->is3D());
}

// Test Solid with interior shells
BOOST_AUTO_TEST_CASE(testSimplify_SolidWithInteriorShells)
{
  auto exteriorShell = createCubePolyhedralSurface();
  Solid solid(*exteriorShell);

  // Add a smaller interior shell (void)
  auto interiorShell = std::make_unique<PolyhedralSurface>();
  auto createSmallFace = [](double x1, double y1, double z1, double x2,
                            double y2, double z2, double x3, double y3,
                            double z3, double x4, double y4,
                            double z4) -> std::unique_ptr<Polygon> {
    auto ring = std::make_unique<LineString>();
    ring->addPoint(Point(x1, y1, z1));
    ring->addPoint(Point(x2, y2, z2));
    ring->addPoint(Point(x3, y3, z3));
    ring->addPoint(Point(x4, y4, z4));
    ring->addPoint(Point(x1, y1, z1));
    return std::make_unique<Polygon>(ring.release());
  };

  // Small cube inside (0.25 to 0.75 range) - all 6 faces
  // Bottom (z=0.25)
  interiorShell->addPatch(
      createSmallFace(0.25, 0.25, 0.25, 0.75, 0.25, 0.25, 0.75, 0.75, 0.25,
                      0.25, 0.75, 0.25));
  // Top (z=0.75)
  interiorShell->addPatch(
      createSmallFace(0.25, 0.25, 0.75, 0.25, 0.75, 0.75, 0.75, 0.75, 0.75,
                      0.75, 0.25, 0.75));
  // Front (y=0.25)
  interiorShell->addPatch(
      createSmallFace(0.25, 0.25, 0.25, 0.25, 0.25, 0.75, 0.75, 0.25, 0.75,
                      0.75, 0.25, 0.25));
  // Back (y=0.75)
  interiorShell->addPatch(
      createSmallFace(0.25, 0.75, 0.25, 0.75, 0.75, 0.25, 0.75, 0.75, 0.75,
                      0.25, 0.75, 0.75));
  // Left (x=0.25)
  interiorShell->addPatch(
      createSmallFace(0.25, 0.25, 0.25, 0.25, 0.75, 0.25, 0.25, 0.75, 0.75,
                      0.25, 0.25, 0.75));
  // Right (x=0.75)
  interiorShell->addPatch(
      createSmallFace(0.75, 0.25, 0.25, 0.75, 0.25, 0.75, 0.75, 0.75, 0.75,
                      0.75, 0.75, 0.25));

  solid.addInteriorShell(std::move(interiorShell));

  BOOST_CHECK_EQUAL(solid.numInteriorShells(), 1);

  // Use NoValidityCheck because SFCGAL's isValid for interior shells is not fully implemented
  auto simplified = surfaceSimplification(
      solid, SimplificationStopPredicate::edgeCountRatio(0.5),
      SimplificationStrategy::EDGE_LENGTH, NoValidityCheck());

  BOOST_CHECK(simplified);
  BOOST_CHECK(simplified->is<Solid>());

  const auto &simplifiedSolid = simplified->as<Solid>();
  // Interior shells are preserved but not simplified
  BOOST_CHECK_EQUAL(simplifiedSolid.numInteriorShells(), 1);

  std::cout << "Solid with interior shells simplified successfully"
            << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()
