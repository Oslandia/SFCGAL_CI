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

#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/intersections.h>

#include <SFCGAL/Exception.h>
#include <SFCGAL/algorithm/intersection.h>
#include <SFCGAL/algorithm/intersects.h>
#include <SFCGAL/detail/GeometrySet.h>
#include <SFCGAL/detail/triangulate/triangulateInGeometrySet.h>

using namespace SFCGAL::detail;

namespace SFCGAL {
namespace algorithm {

// local function : get the number of intersection points between rings of a
// polygon
int
numIntersectionPoints(const CGAL::Polygon_with_holes_2<Kernel> &poly)
{
  int numIntersectionPoints = 0;
  CGAL::Polygon_with_holes_2<Kernel>::Hole_const_iterator hit =
      poly.holes_begin();

  for (int i = 0; i == 0 || hit != poly.holes_end(); ++i) {
    GeometrySet<2> ringI;

    if (i == 0) {
      ringI.addSegments(poly.outer_boundary().edges_begin(),
                        poly.outer_boundary().edges_end());
    } else {
      ringI.addSegments(hit->edges_begin(), hit->edges_end());
      hit++;
    }

    for (CGAL::Polygon_with_holes_2<Kernel>::Hole_const_iterator hjt = hit;
         hjt != poly.holes_end(); ++hjt) {
      GeometrySet<2> ringJ, inter;
      ringJ.addSegments(hjt->edges_begin(), hjt->edges_end());

      algorithm::intersection(ringI, ringJ, inter);
      numIntersectionPoints += inter.points().size();
    }
  }

  return numIntersectionPoints;
}

//
// must be called with pa's dimension larger than pb's
void
intersection(const PrimitiveHandle<2> &pa, const PrimitiveHandle<2> &pb,
             GeometrySet<2> &output, dim_t<2>)
{
  // everything vs a point
  if (pb.handle.which() == PrimitivePoint) {
    if (algorithm::intersects(pa, pb)) {
      output.addPrimitive(*pb.as<CGAL::Point_2<Kernel>>());
    }
  } else if (pa.handle.which() == PrimitiveSurface &&
             pb.handle.which() == PrimitiveSurface) {
    const CGAL::Polygon_with_holes_2<Kernel> *poly1 =
        pa.as<CGAL::Polygon_with_holes_2<Kernel>>();
    const CGAL::Polygon_with_holes_2<Kernel> *poly2 =
        pb.as<CGAL::Polygon_with_holes_2<Kernel>>();

    // shortcut for triangles
    if (poly1->holes_begin() == poly1->holes_end() &&
        poly1->outer_boundary().size() == 3 &&
        poly2->holes_begin() == poly2->holes_end() &&
        poly2->outer_boundary().size() == 3) {
      CGAL::Polygon_2<Kernel>::Vertex_iterator vit1 =
          poly1->outer_boundary().vertices_begin();
      CGAL::Point_2<Kernel>    pa1 = *vit1++;
      CGAL::Point_2<Kernel>    pa2 = *vit1++;
      CGAL::Point_2<Kernel>    pa3 = *vit1;
      CGAL::Triangle_2<Kernel> tri1(pa1, pa2, pa3);

      CGAL::Polygon_2<Kernel>::Vertex_iterator vit2 =
          poly2->outer_boundary().vertices_begin();
      CGAL::Point_2<Kernel>    pb1 = *vit2++;
      CGAL::Point_2<Kernel>    pb2 = *vit2++;
      CGAL::Point_2<Kernel>    pb3 = *vit2;
      CGAL::Triangle_2<Kernel> tri2(pb1, pb2, pb3);

      CGAL::Object interObj = CGAL::intersection(tri1, tri2);
      output.addPrimitive(interObj, /* pointsAsRing */ true);
      return;
    }

    // CGAL::intersection does not work when the intersection is a point or a
    // segment We have to call intersection on boundaries first

    GeometrySet<2> gpoly1, gpoly2;
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
    const CGAL::Segment_2<Kernel> *seg1     = pa.as<CGAL::Segment_2<Kernel>>();
    const CGAL::Segment_2<Kernel> *seg2     = pb.as<CGAL::Segment_2<Kernel>>();
    CGAL::Object                   interObj = CGAL::intersection(*seg1, *seg2);
    output.addPrimitive(interObj);
  } else if (pa.handle.which() == PrimitiveSurface &&
             pb.handle.which() == PrimitiveSegment) {
    const CGAL::Polygon_with_holes_2<Kernel> *poly =
        pa.as<CGAL::Polygon_with_holes_2<Kernel>>();
    const CGAL::Segment_2<Kernel> *seg = pb.as<CGAL::Segment_2<Kernel>>();

    // shortcut for triangle
    if (poly->holes_begin() == poly->holes_end() &&
        poly->outer_boundary().size() == 3) {
      // no holes and 3 vertices => it is a triangle
      CGAL::Polygon_2<Kernel>::Vertex_iterator vit =
          poly->outer_boundary().vertices_begin();
      CGAL::Point_2<Kernel>    p1(*vit++);
      CGAL::Point_2<Kernel>    p2(*vit++);
      CGAL::Point_2<Kernel>    p3(*vit++);
      CGAL::Triangle_2<Kernel> tri(p1, p2, p3);
      CGAL::Object             interObj = CGAL::intersection(tri, *seg);
      output.addPrimitive(interObj);
      return;
    }

    // if it s a regulat polygon, triangulate it and recurse call
    GeometrySet<2> triangles, g;
    triangulate::triangulate(*poly, triangles);
    g.addPrimitive(pb);

    // recurse call
    algorithm::intersection(triangles, g, output);
  }
}
} // namespace algorithm
} // namespace SFCGAL
