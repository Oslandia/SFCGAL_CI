// Copyright (c) 2025-2025, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/surfaceSimplification.h"

#include "SFCGAL/Kernel.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Surface.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/algorithm/isValid.h"

#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

#include <map>
#include <stdexcept>

namespace SFCGAL::algorithm {

namespace SMS = CGAL::Surface_mesh_simplification;

namespace detail {

/**
 * @brief Apply edge collapse simplification using Edge Length cost and Midpoint placement
 */
template <typename StopPredicate>
auto
simplifyMesh(Surface_mesh_3 &mesh, const StopPredicate &stop) -> size_t
{
  using Cost = SMS::Edge_length_cost<Surface_mesh_3>;
  using Placement = SMS::Midpoint_placement<Surface_mesh_3>;

  return SMS::edge_collapse(mesh, stop,
                            CGAL::parameters::get_cost(Cost())
                                .get_placement(Placement()));
}

/**
 * @brief Simplify a surface mesh
 */
auto
simplifySurfaceMesh(Surface_mesh_3                     &mesh,
                    const SimplificationStopPredicate &stopPredicate) -> size_t
{
  if (stopPredicate.type == SimplificationStopPredicate::Type::EDGE_COUNT) {
    // Stop after collapsing a specific number of edges
    auto const initialEdgeCount =
        static_cast<size_t>(mesh.number_of_edges());
    auto const targetEdgeCount = static_cast<size_t>(stopPredicate.value);

    if (targetEdgeCount >= initialEdgeCount) {
      return 0; // No simplification needed
    }

    size_t const edgesToRemove = initialEdgeCount - targetEdgeCount;
    SMS::Edge_count_stop_predicate<Surface_mesh_3> stop(edgesToRemove);

    return simplifyMesh(mesh, stop);

  } else {
    // Stop when edge count ratio is reached
    double const ratio = stopPredicate.value;

    if (ratio <= 0.0 || ratio >= 1.0) {
      throw std::invalid_argument(
          "Edge count ratio must be in the range (0.0, 1.0)");
    }

    SMS::Edge_count_ratio_stop_predicate<Surface_mesh_3> stop(ratio);

    return simplifyMesh(mesh, stop);
  }
}

/**
 * @brief Simplify a TriangulatedSurface
 */
auto
simplifyTriangulatedSurface(const TriangulatedSurface &surface,
                            const SimplificationStopPredicate &stopPredicate)
    -> std::unique_ptr<TriangulatedSurface>
{
  if (surface.isEmpty()) {
    return std::make_unique<TriangulatedSurface>(surface);
  }

  // Convert to Surface_mesh
  Surface_mesh_3 mesh = surface.toSurfaceMesh();

  // Simplify
  simplifySurfaceMesh(mesh, stopPredicate);

  // Convert back to PolyhedralSurface first, then to TriangulatedSurface
  PolyhedralSurface polyhedral(mesh);
  return std::make_unique<TriangulatedSurface>(
      polyhedral.toTriangulatedSurface());
}

/**
 * @brief Simplify a PolyhedralSurface
 */
auto
simplifyPolyhedralSurface(const PolyhedralSurface &surface,
                          const SimplificationStopPredicate &stopPredicate)
    -> std::unique_ptr<PolyhedralSurface>
{
  if (surface.isEmpty()) {
    return std::make_unique<PolyhedralSurface>(surface);
  }

  // Convert to Surface_mesh
  Surface_mesh_3 mesh = surface.toSurfaceMesh();

  // Simplify
  simplifySurfaceMesh(mesh, stopPredicate);

  // Convert back to PolyhedralSurface
  return std::make_unique<PolyhedralSurface>(mesh);
}

/**
 * @brief Simplify a Solid (only exterior shell)
 */
auto
simplifySolid(const Solid &solid, const SimplificationStopPredicate &stopPredicate)
    -> std::unique_ptr<Solid>
{
  if (solid.isEmpty()) {
    return std::make_unique<Solid>(solid);
  }

  // Simplify the exterior shell
  auto simplifiedExterior =
      simplifyPolyhedralSurface(solid.exteriorShell(), stopPredicate);

  // Create new solid with simplified exterior shell
  auto result = std::make_unique<Solid>(*simplifiedExterior);

  // Copy interior shells (not simplified)
  for (size_t i = 0; i < solid.numInteriorShells(); ++i) {
    result->addInteriorShell(solid.interiorShellN(i).clone());
  }

  return result;
}

/**
 * @brief Simplify a MultiSolid
 */
auto
simplifyMultiSolid(const MultiSolid &multiSolid,
                   const SimplificationStopPredicate &stopPredicate)
    -> std::unique_ptr<MultiSolid>
{
  auto result = std::make_unique<MultiSolid>();

  for (size_t i = 0; i < multiSolid.numGeometries(); ++i) {
    const Solid &solid = multiSolid.geometryN(i).as<Solid>();
    result->addGeometry(simplifySolid(solid, stopPredicate));
  }

  return result;
}

} // namespace detail

/**
 * @brief Simplify a surface mesh using CGAL edge collapse algorithm
 * @param geometry The input geometry to simplify (must be a surface or solid)
 * @param stopPredicate When to stop the simplification process
 * @param strategy The cost and placement strategy to use
 * @return A simplified copy of the input geometry
 * @pre The input geometry must be valid and non-empty
 * @pre For EDGE_COUNT_RATIO, the ratio must be in the range (0.0, 1.0)
 * @pre The geometry must be 3-dimensional
 */
auto
surfaceSimplification(const Geometry                 &geometry,
                      const SimplificationStopPredicate &stopPredicate,
                      SimplificationStrategy strategy)
    -> std::unique_ptr<Geometry>
{
  // Validate geometry
  switch (geometry.geometryTypeId()) {
  case TYPE_TRIANGULATEDSURFACE:
  case TYPE_POLYHEDRALSURFACE:
  case TYPE_SOLID:
  case TYPE_MULTISOLID:
    SFCGAL_ASSERT_GEOMETRY_VALIDITY_3D(geometry);
    break;
  default:
    break;
  }

  std::unique_ptr<Geometry> result(
      surfaceSimplification(geometry, stopPredicate, SimplificationStrategy::EDGE_LENGTH, NoValidityCheck()));
  propagateValidityFlag(*result, true);
  return result;
}

auto
surfaceSimplification(const Geometry                 &geometry,
                      const SimplificationStopPredicate &stopPredicate,
                      SimplificationStrategy strategy,
                      NoValidityCheck /*unused*/) -> std::unique_ptr<Geometry>
{
  if (geometry.isEmpty()) {
    return geometry.clone();
  }

  // Validate ratio for EDGE_COUNT_RATIO
  if (stopPredicate.type == SimplificationStopPredicate::Type::EDGE_COUNT_RATIO) {
    if (stopPredicate.value <= 0.0 || stopPredicate.value >= 1.0) {
      throw std::invalid_argument(
          "Edge count ratio must be in the range (0.0, 1.0)");
    }
  }

  switch (geometry.geometryTypeId()) {
  case TYPE_TRIANGULATEDSURFACE:
    return detail::simplifyTriangulatedSurface(
        geometry.as<TriangulatedSurface>(), stopPredicate);

  case TYPE_POLYHEDRALSURFACE:
    return detail::simplifyPolyhedralSurface(geometry.as<PolyhedralSurface>(),
                                             stopPredicate);

  case TYPE_SOLID:
    return detail::simplifySolid(geometry.as<Solid>(), stopPredicate);

  case TYPE_MULTISOLID:
    return detail::simplifyMultiSolid(geometry.as<MultiSolid>(), stopPredicate);

  default:
    throw std::invalid_argument(
        "surfaceSimplification only supports TriangulatedSurface, "
        "PolyhedralSurface, Solid, and MultiSolid geometries");
  }
}

} // namespace SFCGAL::algorithm
