// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/orientation.h"

#include "SFCGAL/algorithm/ConsistentOrientationBuilder.h"

#include "SFCGAL/detail/graph/GeometryGraph.h"
#include "SFCGAL/detail/graph/GeometryGraphBuilder.h"
#include "SFCGAL/detail/graph/algorithm/isHalfEdge.h"

#include "SFCGAL/detail/graph/algorithm/orientation.h"

namespace SFCGAL::algorithm {

auto
makeValidOrientation(CGAL::Polygon_2<Kernel> &polygon) -> void
{
  if (polygon.orientation() != CGAL::COUNTERCLOCKWISE) {
    polygon.reverse_orientation();
  }
}

auto
makeValidOrientation(CGAL::Polygon_with_holes_2<Kernel> &polygon) -> void
{
  if (polygon.outer_boundary().orientation() != CGAL::COUNTERCLOCKWISE) {
    polygon.outer_boundary().reverse_orientation();
  }

  for (auto it = polygon.holes_begin(); it != polygon.holes_end(); ++it) {
    if (it->orientation() != CGAL::CLOCKWISE) {
      it->reverse_orientation();
    }
  }
}

auto
makeValidOrientation(Polygon &polygon) -> void
{
  for (size_t i = 0; i < polygon.numRings(); i++) {
    LineString &ring = polygon.ringN(i);

    if (i == 0) {
      if (ring.toPolygon_2().orientation() != CGAL::COUNTERCLOCKWISE) {
        ring.reverse();
      }
    } else {
      if (ring.toPolygon_2().orientation() != CGAL::CLOCKWISE) {
        ring.reverse();
      }
    }
  }
}

/// @private
auto
hasConsistentOrientation3D(const TriangulatedSurface &g) -> bool
{
  using namespace graph;

  if (g.isEmpty()) {
    return true;
  }

  GeometryGraph        graph;
  GeometryGraphBuilder graphBuilder(graph);
  graphBuilder.addTriangulatedSurface(g);
  return graph::algorithm::isHalfEdge(graph);
}

/// @private
auto
hasConsistentOrientation3D(const PolyhedralSurface &g) -> bool
{
  using namespace graph;

  if (g.isEmpty()) {
    return true;
  }

  GeometryGraph        graph;
  GeometryGraphBuilder graphBuilder(graph);
  graphBuilder.addPolyhedralSurface(g);
  return graph::algorithm::isHalfEdge(graph);
}

auto
makeConsistentOrientation3D(TriangulatedSurface &g) -> void
{
  ConsistentOrientationBuilder builder;
  builder.addTriangulatedSurface(g);
  g = builder.buildTriangulatedSurface();
}

/// @private
auto
isCounterClockWiseOriented(const LineString &lineString) -> bool
{
  // Use CGAL's built-in orientation() which properly handles lazy exact numbers
  // The manual Newell's formula builds a huge DAG with lazy exact arithmetic
  // causing stack overflow on Windows with polygons having 100,000+ points
  // (see PostGIS issue #3651)
  return lineString.toPolygon_2(false).orientation() == CGAL::COUNTERCLOCKWISE;
}

/// @private
auto
isCounterClockWiseOriented(const Triangle &triangle) -> bool
{
  if (triangle.isEmpty()) {
    return false;
  }
  // Compute the 'z' part of the cross product

  return (triangle.vertex(2).x() - triangle.vertex(1).x()) *
                 (triangle.vertex(0).y() - triangle.vertex(1).y()) -
             (triangle.vertex(2).y() - triangle.vertex(1).y()) *
                 (triangle.vertex(0).x() - triangle.vertex(1).x()) >
         0;
}

/// @private
auto
isCounterClockWiseOriented(const Polygon &polygon) -> bool
{
  if (polygon.isEmpty()) {
    return false;
  }
  return isCounterClockWiseOriented(polygon.exteriorRing());
}

} // namespace SFCGAL::algorithm
