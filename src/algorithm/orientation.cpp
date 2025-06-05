// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/orientation.h"

#include "SFCGAL/algorithm/ConsistentOrientationBuilder.h"

#include "SFCGAL/detail/graph/GeometryGraph.h"
#include "SFCGAL/detail/graph/GeometryGraphBuilder.h"
#include "SFCGAL/detail/graph/algorithm/isHalfEdge.h"

#include "SFCGAL/detail/graph/algorithm/orientation.h"

namespace SFCGAL::algorithm {

void
makeValidOrientation(CGAL::Polygon_2<Kernel> &polygon)
{
  if (polygon.orientation() != CGAL::COUNTERCLOCKWISE) {
    polygon.reverse_orientation();
  }
}

void
makeValidOrientation(CGAL::Polygon_with_holes_2<Kernel> &polygon)
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

void
makeValidOrientation(Polygon &polygon)
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

void
makeConsistentOrientation3D(TriangulatedSurface &g)
{
  ConsistentOrientationBuilder builder;
  builder.addTriangulatedSurface(g);
  g = builder.buildTriangulatedSurface();
}

auto
isCounterClockWiseOriented(const LineString &ls) -> bool
{
  // Compute the 'z' part of the Newell's formula
  // and test against 0
  Kernel::FT z = 0;

  for (size_t i = 0; i < ls.numSegments(); ++i) {
    const Point &pi = ls.pointN(i);
    const Point &pj = ls.pointN(i + 1);
    z += (pi.x() - pj.x()) * (pi.y() + pj.y());
  }

  return z > 0;
}

auto
isCounterClockWiseOriented(const Triangle &tri) -> bool
{
  // Compute the 'z' part of the cross product

  return (tri.vertex(2).x() - tri.vertex(1).x()) *
                 (tri.vertex(0).y() - tri.vertex(1).y()) -
             (tri.vertex(2).y() - tri.vertex(1).y()) *
                 (tri.vertex(0).x() - tri.vertex(1).x()) >
         0;
}

auto
isCounterClockWiseOriented(const Polygon &poly) -> bool
{
  return isCounterClockWiseOriented(poly.exteriorRing());
}

} // namespace SFCGAL::algorithm
