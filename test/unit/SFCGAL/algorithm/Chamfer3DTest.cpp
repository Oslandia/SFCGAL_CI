// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/chamfer3D.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/io/wkt.h"
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <memory>

using namespace SFCGAL;

BOOST_AUTO_TEST_SUITE(SFCGAL_algorithm_Chamfer3DTest)

// Helper function to create a unit cube solid
auto
createUnitCube() -> std::unique_ptr<Solid>
{
  // Create a cube from (0,0,0) to (1,1,1)
  auto ps = std::make_unique<PolyhedralSurface>();

  // Bottom face (z=0)
  {
    auto *ring = new LineString();
    ring->addPoint(Point(0, 0, 0));
    ring->addPoint(Point(1, 0, 0));
    ring->addPoint(Point(1, 1, 0));
    ring->addPoint(Point(0, 1, 0));
    ring->addPoint(Point(0, 0, 0));
    ps->addPatch(Polygon(ring));
  }

  // Top face (z=1)
  {
    auto *ring = new LineString();
    ring->addPoint(Point(0, 0, 1));
    ring->addPoint(Point(0, 1, 1));
    ring->addPoint(Point(1, 1, 1));
    ring->addPoint(Point(1, 0, 1));
    ring->addPoint(Point(0, 0, 1));
    ps->addPatch(Polygon(ring));
  }

  // Front face (y=0)
  {
    auto *ring = new LineString();
    ring->addPoint(Point(0, 0, 0));
    ring->addPoint(Point(0, 0, 1));
    ring->addPoint(Point(1, 0, 1));
    ring->addPoint(Point(1, 0, 0));
    ring->addPoint(Point(0, 0, 0));
    ps->addPatch(Polygon(ring));
  }

  // Back face (y=1)
  {
    auto *ring = new LineString();
    ring->addPoint(Point(0, 1, 0));
    ring->addPoint(Point(1, 1, 0));
    ring->addPoint(Point(1, 1, 1));
    ring->addPoint(Point(0, 1, 1));
    ring->addPoint(Point(0, 1, 0));
    ps->addPatch(Polygon(ring));
  }

  // Left face (x=0)
  {
    auto *ring = new LineString();
    ring->addPoint(Point(0, 0, 0));
    ring->addPoint(Point(0, 1, 0));
    ring->addPoint(Point(0, 1, 1));
    ring->addPoint(Point(0, 0, 1));
    ring->addPoint(Point(0, 0, 0));
    ps->addPatch(Polygon(ring));
  }

  // Right face (x=1)
  {
    auto *ring = new LineString();
    ring->addPoint(Point(1, 0, 0));
    ring->addPoint(Point(1, 0, 1));
    ring->addPoint(Point(1, 1, 1));
    ring->addPoint(Point(1, 1, 0));
    ring->addPoint(Point(1, 0, 0));
    ps->addPatch(Polygon(ring));
  }

  return std::make_unique<Solid>(*ps);
}

BOOST_AUTO_TEST_CASE(testChamfer3D_InvalidInput_2D)
{
  // 2D geometry should throw
  std::unique_ptr<Geometry> polygon(io::readWkt("POLYGON((0 0, 1 0, 1 1, 0 1, 0 0))"));

  std::vector<algorithm::EdgeIdentifier> edges;
  edges.emplace_back(Kernel::Point_3(0, 0, 0), Kernel::Point_3(1, 0, 0));
  auto params = algorithm::ChamferParameters::symmetric(0.1);

  BOOST_CHECK_THROW(
      algorithm::chamfer3D(*polygon, edges, params),
      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(testChamfer3D_InvalidInput_Empty)
{
  PolyhedralSurface empty;

  std::vector<algorithm::EdgeIdentifier> edges;
  edges.emplace_back(Kernel::Point_3(0, 0, 0), Kernel::Point_3(1, 0, 0));
  auto params = algorithm::ChamferParameters::symmetric(0.1);

  BOOST_CHECK_THROW(
      algorithm::chamfer3D(empty, edges, params),
      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(testChamfer3D_InvalidDistance)
{
  auto cube = createUnitCube();

  std::vector<algorithm::EdgeIdentifier> edges;
  edges.emplace_back(Kernel::Point_3(0, 0, 0), Kernel::Point_3(1, 0, 0));

  auto params_zero = algorithm::ChamferParameters::symmetric(0.0);
  BOOST_CHECK_THROW(
      algorithm::chamfer3D(*cube, edges, params_zero),
      std::invalid_argument);

  auto params_neg = algorithm::ChamferParameters::symmetric(-0.1);
  BOOST_CHECK_THROW(
      algorithm::chamfer3D(*cube, edges, params_neg),
      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(testChamfer3D_SingleEdge)
{
  auto cube = createUnitCube();

  // Chamfer one edge
  std::vector<algorithm::EdgeIdentifier> edges;
  edges.emplace_back(Kernel::Point_3(0, 0, 0), Kernel::Point_3(1, 0, 0));

  auto params = algorithm::ChamferParameters::symmetric(0.1);

  std::unique_ptr<Geometry> result = algorithm::chamfer3D(*cube, edges, params);

  BOOST_CHECK(result != nullptr);
  BOOST_CHECK(!result->isEmpty());
  BOOST_CHECK(result->is<Solid>() || result->is<PolyhedralSurface>());
}

BOOST_AUTO_TEST_CASE(testChamfer3D_FindSharpEdges)
{
  auto cube = createUnitCube();

  algorithm::Chamfer3D chamfer(*cube);

  // All edges of a cube should be "sharp" with 90-degree threshold
  // A cube has 12 edges
  auto sharpEdges = chamfer.findSharpEdges(M_PI / 4); // 45 degrees

  BOOST_CHECK_EQUAL(sharpEdges.size(), 12U);
}

BOOST_AUTO_TEST_CASE(testChamfer3D_ComputeDihedralAngle)
{
  auto cube = createUnitCube();

  algorithm::Chamfer3D chamfer(*cube);

  // Test an edge of the cube - should have 90 degree dihedral angle
  algorithm::EdgeIdentifier edge(Kernel::Point_3(0, 0, 0), Kernel::Point_3(1, 0, 0));

  double angle = chamfer.computeDihedralAngle(edge);

  // Cube edges have 90 degree (PI/2) dihedral angle
  // Allow some tolerance for floating point comparison
  BOOST_CHECK_CLOSE(angle, M_PI / 2, 1.0); // 1% tolerance
}

BOOST_AUTO_TEST_CASE(testChamfer3D_EdgeSelector_Explicit)
{
  auto cube = createUnitCube();

  // Select only one edge to chamfer
  std::vector<algorithm::EdgeIdentifier> edges;
  edges.emplace_back(Kernel::Point_3(0, 0, 0), Kernel::Point_3(1, 0, 0));

  auto params = algorithm::ChamferParameters::symmetric(0.1);

  std::unique_ptr<Geometry> result = algorithm::chamfer3D(*cube, edges, params);

  BOOST_CHECK(result != nullptr);
}

BOOST_AUTO_TEST_CASE(testChamfer3D_Asymmetric)
{
  auto cube = createUnitCube();

  // Chamfer one edge with asymmetric parameters
  std::vector<algorithm::EdgeIdentifier> edges;
  edges.emplace_back(Kernel::Point_3(0, 0, 0), Kernel::Point_3(1, 0, 0));

  auto selector = algorithm::EdgeSelector::explicit_(edges);
  auto params = algorithm::ChamferParameters::asymmetric(0.1, 0.2);

  algorithm::Chamfer3D chamfer(*cube);
  std::unique_ptr<Geometry> result = chamfer.chamferEdges(selector, params);

  BOOST_CHECK(result != nullptr);
}

BOOST_AUTO_TEST_CASE(testChamfer3D_VertexChamfer)
{
  auto cube = createUnitCube();

  // Chamfer one vertex
  std::vector<Kernel::Point_3> vertices;
  vertices.emplace_back(0, 0, 0);

  auto params = algorithm::VertexChamferParameters::create(0.1);

  std::unique_ptr<Geometry> result = algorithm::chamfer3DVertices(*cube, vertices, params);

  BOOST_CHECK(result != nullptr);
}

BOOST_AUTO_TEST_CASE(testChamfer3D_PolyhedralSurface)
{
  // Create a simple polyhedral surface (open box without top)
  auto ps = std::make_unique<PolyhedralSurface>();

  // Bottom face
  {
    auto *ring = new LineString();
    ring->addPoint(Point(0, 0, 0));
    ring->addPoint(Point(1, 0, 0));
    ring->addPoint(Point(1, 1, 0));
    ring->addPoint(Point(0, 1, 0));
    ring->addPoint(Point(0, 0, 0));
    ps->addPatch(Polygon(ring));
  }

  // Front face
  {
    auto *ring = new LineString();
    ring->addPoint(Point(0, 0, 0));
    ring->addPoint(Point(0, 0, 1));
    ring->addPoint(Point(1, 0, 1));
    ring->addPoint(Point(1, 0, 0));
    ring->addPoint(Point(0, 0, 0));
    ps->addPatch(Polygon(ring));
  }

  // Try to chamfer the shared edge - this edge is internal to both faces
  std::vector<algorithm::EdgeIdentifier> edges;
  edges.emplace_back(Kernel::Point_3(0, 0, 0), Kernel::Point_3(1, 0, 0));

  auto params = algorithm::ChamferParameters::symmetric(0.1);
  std::unique_ptr<Geometry> result = algorithm::chamfer3D(*ps, edges, params);

  BOOST_CHECK(result != nullptr);
  BOOST_CHECK(!result->isEmpty());
}

BOOST_AUTO_TEST_CASE(testChamfer3D_EdgeIdentifier_Equality)
{
  // Test edge identifier equality (order-independent)
  algorithm::EdgeIdentifier e1(Kernel::Point_3(0, 0, 0), Kernel::Point_3(1, 0, 0));
  algorithm::EdgeIdentifier e2(Kernel::Point_3(1, 0, 0), Kernel::Point_3(0, 0, 0));
  algorithm::EdgeIdentifier e3(Kernel::Point_3(0, 0, 0), Kernel::Point_3(0, 1, 0));

  BOOST_CHECK(e1 == e2); // Same edge, different order
  BOOST_CHECK(!(e1 == e3)); // Different edge
}

BOOST_AUTO_TEST_CASE(testChamfer3D_ChamferParameters)
{
  auto sym = algorithm::ChamferParameters::symmetric(0.5);
  BOOST_CHECK(sym.type == algorithm::ChamferParameters::Type::SYMMETRIC);
  BOOST_CHECK_CLOSE(sym.distance1, 0.5, 0.001);
  BOOST_CHECK_CLOSE(sym.distance2, 0.5, 0.001);

  auto asym = algorithm::ChamferParameters::asymmetric(0.3, 0.7);
  BOOST_CHECK(asym.type == algorithm::ChamferParameters::Type::ASYMMETRIC);
  BOOST_CHECK_CLOSE(asym.distance1, 0.3, 0.001);
  BOOST_CHECK_CLOSE(asym.distance2, 0.7, 0.001);
}

BOOST_AUTO_TEST_CASE(testChamfer3D_EdgeSelector_Modes)
{
  auto sel1 = algorithm::EdgeSelector::all();
  BOOST_CHECK(sel1.mode == algorithm::EdgeSelector::Mode::ALL);

  auto sel2 = algorithm::EdgeSelector::byAngle(M_PI / 6);
  BOOST_CHECK(sel2.mode == algorithm::EdgeSelector::Mode::BY_ANGLE);
  BOOST_CHECK_CLOSE(sel2.angleThreshold, M_PI / 6, 0.001);

  std::vector<algorithm::EdgeIdentifier> edges;
  edges.emplace_back(Kernel::Point_3(0, 0, 0), Kernel::Point_3(1, 0, 0));
  auto sel3 = algorithm::EdgeSelector::explicit_(edges);
  BOOST_CHECK(sel3.mode == algorithm::EdgeSelector::Mode::EXPLICIT);
  BOOST_CHECK_EQUAL(sel3.edges.size(), 1U);
}

BOOST_AUTO_TEST_CASE(testChamfer3D_ContinuousEdges)
{
  auto cube = createUnitCube();

  // Chamfer two adjacent edges (sharing vertex at (0,0,0))
  std::vector<algorithm::EdgeIdentifier> edges;
  edges.emplace_back(Kernel::Point_3(0, 0, 0), Kernel::Point_3(1, 0, 0)); // bottom front
  edges.emplace_back(Kernel::Point_3(0, 0, 0), Kernel::Point_3(0, 1, 0)); // bottom left

  auto params = algorithm::ChamferParameters::symmetric(0.1);

  std::unique_ptr<Geometry> result = algorithm::chamfer3D(*cube, edges, params);

  BOOST_CHECK(result != nullptr);
  BOOST_CHECK(!result->isEmpty());
}

BOOST_AUTO_TEST_SUITE_END()
