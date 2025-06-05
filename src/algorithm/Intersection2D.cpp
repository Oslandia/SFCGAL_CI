// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/intersections.h>

#include "SFCGAL/Exception.h"
#include "SFCGAL/algorithm/intersection.h"
#include "SFCGAL/algorithm/intersects.h"
#include "SFCGAL/detail/GeometrySet.h"
#include "SFCGAL/detail/triangulate/triangulateInGeometrySet.h"

using namespace SFCGAL::detail;

namespace SFCGAL::algorithm {

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection

// local function : get the number of intersection points between rings of a
// polygon
auto
numIntersectionPoints(const CGAL::Polygon_with_holes_2<Kernel> &poly) -> int
{
  int  numIntersectionPoints = 0;
  auto hit                   = poly.holes_begin();

  for (int i = 0; i == 0 || hit != poly.holes_end(); ++i) {
    GeometrySet<2> ringI;

    if (i == 0) {
      ringI.addSegments(poly.outer_boundary().edges_begin(),
                        poly.outer_boundary().edges_end());
    } else {
      ringI.addSegments(hit->edges_begin(), hit->edges_end());
      hit++;
    }

    for (auto hjt = hit; hjt != poly.holes_end(); ++hjt) {
      GeometrySet<2> ringJ;
      GeometrySet<2> inter;
      ringJ.addSegments(hjt->edges_begin(), hjt->edges_end());

      algorithm::intersection(ringI, ringJ, inter);
      numIntersectionPoints += inter.points().size();
    }
  }

  return numIntersectionPoints;
}

/// @} end of private section

// ----------------------------------------------------------------------------------
// -- public interface
// ----------------------------------------------------------------------------------
/// @publicsection

/// must be called with pa's dimension larger than pb's
void
intersection(const PrimitiveHandle<2> &pa, const PrimitiveHandle<2> &pb,
             GeometrySet<2> &output, dim_t<2> /*unused*/)
{
  // everything vs a point
  if (pb.handle.which() == PrimitivePoint) {
    if (algorithm::intersects(pa, pb)) {
      output.addPrimitive(*pb.as<CGAL::Point_2<Kernel>>());
    }
  } else if (pa.handle.which() == PrimitiveSurface &&
             pb.handle.which() == PrimitiveSurface) {
    const auto *poly1 = pa.as<CGAL::Polygon_with_holes_2<Kernel>>();
    const auto *poly2 = pb.as<CGAL::Polygon_with_holes_2<Kernel>>();

    // shortcut for triangles
    if (poly1->holes_begin() == poly1->holes_end() &&
        poly1->outer_boundary().size() == 3 &&
        poly2->holes_begin() == poly2->holes_end() &&
        poly2->outer_boundary().size() == 3) {
      auto vit1 = poly1->outer_boundary().vertices_begin();
      CGAL::Point_2<Kernel> const    pa1 = *vit1++;
      CGAL::Point_2<Kernel> const    pa2 = *vit1++;
      CGAL::Point_2<Kernel> const    pa3 = *vit1;
      CGAL::Triangle_2<Kernel> const tri1(pa1, pa2, pa3);

      auto vit2 = poly2->outer_boundary().vertices_begin();
      CGAL::Point_2<Kernel> const    pb1 = *vit2++;
      CGAL::Point_2<Kernel> const    pb2 = *vit2++;
      CGAL::Point_2<Kernel> const    pb3 = *vit2;
      CGAL::Triangle_2<Kernel> const tri2(pb1, pb2, pb3);

      CGAL::Object const interObj = CGAL::intersection(tri1, tri2);
      output.addPrimitive(interObj, /* pointsAsRing */ true);
      return;
    }

    // CGAL::intersection does not work when the intersection is a point or a
    // segment We have to call intersection on boundaries first

    GeometrySet<2> gpoly1;
    GeometrySet<2> gpoly2;
    gpoly1.addBoundary(*poly1);
    gpoly2.addBoundary(*poly2);

    algorithm::intersection(gpoly1, gpoly2, output);

    // CGAL::intersection does not work when rings intersect themselves
    // However, touching by a single point is valid for OGC
    // FIXME: not implemented yet
    if (numIntersectionPoints(*poly1) > 0) {
      BOOST_THROW_EXCEPTION(NotImplementedException(
          "Intersection does not support polygon with connected rings"));
    }

    if (numIntersectionPoints(*poly2) > 0) {
      BOOST_THROW_EXCEPTION(NotImplementedException(
          "Intersection does not support polygon with connected rings"));
    }

    // now call on polygon's interiors
    CGAL::intersection(*poly1, *poly2, std::back_inserter(output.surfaces()));
  } else if (pa.handle.which() == PrimitiveSegment &&
             pb.handle.which() == PrimitiveSegment) {
    const auto        *seg1     = pa.as<CGAL::Segment_2<Kernel>>();
    const auto        *seg2     = pb.as<CGAL::Segment_2<Kernel>>();
    CGAL::Object const interObj = CGAL::intersection(*seg1, *seg2);
    output.addPrimitive(interObj);
  } else if (pa.handle.which() == PrimitiveSurface &&
             pb.handle.which() == PrimitiveSegment) {
    const auto *poly = pa.as<CGAL::Polygon_with_holes_2<Kernel>>();
    const auto *seg  = pb.as<CGAL::Segment_2<Kernel>>();

    // shortcut for triangle
    if (poly->holes_begin() == poly->holes_end() &&
        poly->outer_boundary().size() == 3) {
      // no holes and 3 vertices => it is a triangle
      auto                        vit = poly->outer_boundary().vertices_begin();
      CGAL::Point_2<Kernel> const p1(*vit++);
      CGAL::Point_2<Kernel> const p2(*vit++);
      CGAL::Point_2<Kernel> const p3(*vit++);
      CGAL::Triangle_2<Kernel> const tri(p1, p2, p3);
      CGAL::Object const             interObj = CGAL::intersection(tri, *seg);
      output.addPrimitive(interObj);
      return;
    }

    // if it s a regulat polygon, triangulate it and recurse call
    GeometrySet<2> triangles;
    GeometrySet<2> g;
    triangulate::triangulate(*poly, triangles);
    g.addPrimitive(pb);

    // recurse call
    algorithm::intersection(triangles, g, output);
  }
}
} // namespace SFCGAL::algorithm
