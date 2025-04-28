// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <map>
#include <memory>
#include <sstream>

#include "SFCGAL/Envelope.h"
#include "SFCGAL/Exception.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/algorithm/connection.h"
#include "SFCGAL/algorithm/covers.h"
#include "SFCGAL/algorithm/intersection.h"
#include "SFCGAL/algorithm/intersects.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/detail/GeometrySet.h"
#include "SFCGAL/detail/triangulate/triangulateInGeometrySet.h"

#include <CGAL/box_intersection_d.h>

#include <CGAL/Side_of_triangle_mesh.h>

using namespace SFCGAL::detail;

namespace SFCGAL::algorithm {

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection

/// Type of pa must be of larger dimension than type of pb
auto
_intersects(const PrimitiveHandle<2> &pa, const PrimitiveHandle<2> &pb) -> bool
{
  //
  // Point vs. Point
  //

  if (pa.handle.which() == PrimitivePoint &&
      pb.handle.which() == PrimitivePoint) {
    return *boost::get<const CGAL::Point_2<Kernel> *>(pa.handle) ==
           *boost::get<const CGAL::Point_2<Kernel> *>(pb.handle);
  }

  //
  // Segment vs. Point
  //

  if (pa.handle.which() == PrimitiveSegment &&
      pb.handle.which() == PrimitivePoint) {
    const auto *seg = pa.as<CGAL::Segment_2<Kernel>>();
    const auto *pt  = pb.as<CGAL::Point_2<Kernel>>();
    return seg->has_on(*pt);
  }

  //
  // Segment vs. Segment
  //

  if (pa.handle.which() == PrimitiveSegment &&
      pb.handle.which() == PrimitiveSegment) {
    const auto *seg1 = pa.as<CGAL::Segment_2<Kernel>>();
    const auto *seg2 = pb.as<CGAL::Segment_2<Kernel>>();
    return CGAL::do_intersect(*seg1, *seg2);
  }

  //
  // Polygon vs. Point
  //

  if (pa.handle.which() == PrimitiveSurface &&
      pb.handle.which() == PrimitivePoint) {
    // Polygon versus Point
    const auto *poly = pa.as<CGAL::Polygon_with_holes_2<Kernel>>();
    const auto *pt   = pb.as<CGAL::Point_2<Kernel>>();

    int const b1 = poly->outer_boundary().bounded_side(*pt);

    if (b1 == CGAL::ON_BOUNDARY) {
      return true;
    }

    if (b1 == CGAL::ON_BOUNDED_SIDE) {
      CGAL::Polygon_with_holes_2<Kernel>::Hole_const_iterator it;

      for (it = poly->holes_begin(); it != poly->holes_end(); ++it) {
        int const b = it->bounded_side(*pt);

        if (b == CGAL::ON_BOUNDED_SIDE) {
          return false;
        }
      }
    } else {
      return false;
    }

    return true;
  }

  //
  // Polygon vs. Segment
  //

  if (pa.handle.which() == PrimitiveSurface &&
      pb.handle.which() == PrimitiveSegment) {
    const auto *poly = pa.as<CGAL::Polygon_with_holes_2<Kernel>>();
    const auto *seg  = pb.as<CGAL::Segment_2<Kernel>>();

    // 1. if the segment intersects a boundary of the polygon, returns true
    // 2. else, if one of the point of the segment intersects the polygon,
    // returns true

    GeometrySet<2> segment;
    segment.addSegments(seg, seg + 1);

    // 1.
    GeometrySet<2> boundary;
    boundary.addSegments(poly->outer_boundary().edges_begin(),
                         poly->outer_boundary().edges_end());

    // recurse call
    if (intersects(boundary, segment)) {
      return true;
    }

    for (auto it = poly->holes_begin(); it != poly->holes_end(); ++it) {
      GeometrySet<2> hole;
      hole.addSegments(it->edges_begin(), it->edges_end());

      // recurse call
      if (intersects(hole, segment)) {
        return true;
      }
    }

    // 2. call the polygon, point version
    CGAL::Point_2<Kernel> const pt = seg->source();
    PrimitiveHandle<2> const    ppoly(poly);
    PrimitiveHandle<2> const    ppt(&pt);
    return intersects(ppoly, ppt);
  }

  //
  // Polygon vs. Polygon
  //

  if (pa.handle.which() == PrimitiveSurface &&
      pb.handle.which() == PrimitiveSurface) {
    const auto *poly1 = pa.as<CGAL::Polygon_with_holes_2<Kernel>>();
    const auto *poly2 = pb.as<CGAL::Polygon_with_holes_2<Kernel>>();

    // 1. if rings intersects, returns true
    // 2. else, if poly1 is inside poly2 or poly1 inside poly2 (but not in
    // holes), returns true

    GeometrySet<2> rings1;
    GeometrySet<2> rings2;
    rings1.addSegments(poly1->outer_boundary().edges_begin(),
                       poly1->outer_boundary().edges_end());

    for (auto it = poly1->holes_begin(); it != poly1->holes_end(); ++it) {
      rings1.addSegments(it->edges_begin(), it->edges_end());
    }

    rings2.addSegments(poly2->outer_boundary().edges_begin(),
                       poly2->outer_boundary().edges_end());

    for (auto it = poly2->holes_begin(); it != poly2->holes_end(); ++it) {
      rings2.addSegments(it->edges_begin(), it->edges_end());
    }

    // 1.
    if (intersects(rings1, rings2)) {
      return true;
    }

    // 2.
    CGAL::Bbox_2 box1;
    CGAL::Bbox_2 box2;
    box1 = poly1->bbox();
    box2 = poly2->bbox();
    Envelope const e1(box1.xmin(), box1.xmax(), box1.ymin(), box1.ymax());
    Envelope const e2(box2.xmin(), box2.xmax(), box2.ymin(), box2.ymax());

    // if pa is inside pb
    if (Envelope::contains(e2, e1)) {
      // is pa inside one of pb's holes ?
      CGAL::Point_2<Kernel> const pt =
          *poly1->outer_boundary().vertices_begin();

      for (auto it = poly2->holes_begin(); it != poly2->holes_end(); ++it) {
        CGAL::Bounded_side const b2 = it->bounded_side(pt);

        if (b2 == CGAL::ON_BOUNDED_SIDE) {
          return false;
        }
      }

      return true;
    }

    // if pb is inside pa
    if (Envelope::contains(e1, e2)) {
      // is pa inside one of pb's holes ?
      CGAL::Point_2<Kernel> const pt =
          *poly2->outer_boundary().vertices_begin();

      for (auto it = poly1->holes_begin(); it != poly1->holes_end(); ++it) {
        CGAL::Bounded_side const b2 = it->bounded_side(pt);

        if (b2 == CGAL::ON_BOUNDED_SIDE) {
          return false;
        }
      }

      return true;
    }

    return false;
  }

  return false;
}

/// intersects of a volume with any other type
struct intersects_volume_x : public boost::static_visitor<bool> {
  const MarkedPolyhedron *polyhedron;

  intersects_volume_x(const MarkedPolyhedron *vol) : polyhedron(vol) {}

  template <class T>
  auto
  operator()(const T *geometry) const -> bool
  {
    // intersection between a solid and a geometry
    // 1. either one of the geometry' point lies inside the solid
    // 2. or the geometry intersects one of the surfaces

    // 1.

    if (polyhedron->is_closed()) {
      // this test is needed only if its a volume
      // if the polyhedron is not closed, this is not a volume, actually

      CGAL::Side_of_triangle_mesh<MarkedPolyhedron, Kernel> const is_in_poly(
          *polyhedron);

      GeometrySet<3> points;
      points.collectPoints(geometry);

      for (const auto &pit : points.points()) {
        if (is_in_poly(pit.primitive()) != CGAL::ON_UNBOUNDED_SIDE) {
          return true;
        }
      }
    }

    // 2.

    // triangulate the polyhedron_3
    GeometrySet<3> g;
    g.addPrimitive(*geometry);

    GeometrySet<3> triangles;
    triangulate::triangulate(*polyhedron, triangles);

    return intersects(g, triangles);
  }
};

/// Type of pa must be of larger dimension than type of pb
auto
_intersects(const PrimitiveHandle<3> &pa, const PrimitiveHandle<3> &pb) -> bool
{
  if (pa.handle.which() == PrimitivePoint &&
      pb.handle.which() == PrimitivePoint) {
    return *boost::get<const CGAL::Point_3<Kernel> *>(pa.handle) ==
           *boost::get<const CGAL::Point_3<Kernel> *>(pb.handle);
  }
  if (pa.handle.which() == PrimitiveSegment &&
      pb.handle.which() == PrimitivePoint) {
    const auto *seg = pa.as<CGAL::Segment_3<Kernel>>();
    const auto *pt  = pb.as<CGAL::Point_3<Kernel>>();
    return seg->has_on(*pt);
  }
  if (pa.handle.which() == PrimitiveSegment &&
      pb.handle.which() == PrimitiveSegment) {
    const auto *sega = pa.as<CGAL::Segment_3<Kernel>>();
    const auto *segb = pb.as<CGAL::Segment_3<Kernel>>();
    return CGAL::do_intersect(*sega, *segb);
  }

  if (pa.handle.which() == PrimitiveVolume) {
    intersects_volume_x visitor(pa.as<MarkedPolyhedron>());
    return boost::apply_visitor(visitor, pb.handle);
  }

  if (pa.handle.which() == PrimitiveSurface &&
      pb.handle.which() == PrimitivePoint) {
    const auto *tri = pa.as<CGAL::Triangle_3<Kernel>>();
    const auto *pt  = pb.as<CGAL::Point_3<Kernel>>();
    return tri->has_on(*pt);
  }

  if (pa.handle.which() == PrimitiveSurface &&
      pb.handle.which() == PrimitiveSegment) {
    const auto *tri = pa.as<CGAL::Triangle_3<Kernel>>();
    const auto *seg = pb.as<CGAL::Segment_3<Kernel>>();
    return CGAL::do_intersect(*tri, *seg);
  }

  if (pa.handle.which() == PrimitiveSurface &&
      pb.handle.which() == PrimitiveSurface) {
    const auto *tri1 = pa.as<CGAL::Triangle_3<Kernel>>();
    const auto *tri2 = pb.as<CGAL::Triangle_3<Kernel>>();
    return CGAL::do_intersect(*tri1, *tri2);
  }

  return false;
}

//
// We deal here with symmetric call
template <int Dim>
auto
dispatch_intersects_sym(const PrimitiveHandle<Dim> &pa,
                        const PrimitiveHandle<Dim> &pb) -> bool
{
  // assume types are ordered by dimension within the boost::variant
  if (pa.handle.which() >= pb.handle.which()) {
    return _intersects(pa, pb);
  }
  return _intersects(pb, pa);
}

/// @} end of private section

// ----------------------------------------------------------------------------------
// -- public interface
// ----------------------------------------------------------------------------------
/// @publicsection

template <int Dim>
auto
intersects(const PrimitiveHandle<Dim> &pa, const PrimitiveHandle<Dim> &pb)
    -> bool
{
  return dispatch_intersects_sym(pa, pb);
}

struct found_an_intersection {};

template <int Dim>
struct intersects_cb {
  void
  operator()(const typename PrimitiveBox<Dim>::Type &a,
             const typename PrimitiveBox<Dim>::Type &b)
  {
    if (dispatch_intersects_sym(*a.handle(), *b.handle())) {
      throw found_an_intersection();
    }
  }
};

template <int Dim>
auto
intersects(const GeometrySet<Dim> &a, const GeometrySet<Dim> &b) -> bool
{
  typename SFCGAL::detail::HandleCollection<Dim>::Type ahandles;
  typename SFCGAL::detail::HandleCollection<Dim>::Type bhandles;
  typename SFCGAL::detail::BoxCollection<Dim>::Type    aboxes;
  typename SFCGAL::detail::BoxCollection<Dim>::Type    bboxes;
  a.computeBoundingBoxes(ahandles, aboxes);
  b.computeBoundingBoxes(bhandles, bboxes);

  try {
    intersects_cb<Dim> const cb;
    CGAL::box_intersection_d(aboxes.begin(), aboxes.end(), bboxes.begin(),
                             bboxes.end(), cb);
  } catch (found_an_intersection &e) {
    return true;
  }

  return false;
}

/// @private
template bool
intersects<2>(const GeometrySet<2> &a, const GeometrySet<2> &b);
/// @private
template bool
intersects<3>(const GeometrySet<3> &a, const GeometrySet<3> &b);

/// @private
template bool
intersects<2>(const PrimitiveHandle<2> &a, const PrimitiveHandle<2> &b);
/// @private
template bool
intersects<3>(const PrimitiveHandle<3> &a, const PrimitiveHandle<3> &b);

auto
intersects(const Geometry &ga, const Geometry &gb) -> bool
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(ga);
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(gb);

  GeometrySet<2> const gsa(ga);
  GeometrySet<2> const gsb(gb);

  return intersects(gsa, gsb);
}

auto
intersects3D(const Geometry &ga, const Geometry &gb) -> bool
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_3D(ga);
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_3D(gb);

  GeometrySet<3> const gsa(ga);
  GeometrySet<3> const gsb(gb);

  return intersects(gsa, gsb);
}

auto
intersects(const Geometry &ga, const Geometry &gb, NoValidityCheck /*unused*/)
    -> bool
{
  GeometrySet<2> const gsa(ga);
  GeometrySet<2> const gsb(gb);

  return intersects(gsa, gsb);
}

auto
intersects3D(const Geometry &ga, const Geometry &gb, NoValidityCheck /*unused*/)
    -> bool
{
  GeometrySet<3> const gsa(ga);
  GeometrySet<3> const gsb(gb);

  return intersects(gsa, gsb);
}

/// @private
template <int Dim>
auto
selfIntersectsImpl(const LineString &line) -> bool
{

  if (line.numSegments() < 2) {
    return false; // one segment cannot intersect
  }

  // note: zero length segments are a pain, to avoid algorithm complexity
  // we start by filtering them out
  const size_t numPoints = line.numPoints();
  LineString   l;

  for (size_t i = 0; i != numPoints; ++i) {
    if (i == 0 || l.endPoint() != line.pointN(i)) {
      l.addPoint(line.pointN(i));
    }
  }

  const size_t numSegments = l.numSegments();

  // test any two pairs of segments
  for (size_t i = 0; i != numSegments; ++i) {
    // first line segment is point i and i+1
    for (size_t j = i + 1; j < numSegments; ++j) {
      /** @todo find a way to avoid ugly copy/paste here, toPoint_d< Dim > can
       * be used, but I dont know what to do with Kernel::Segment_Dim and
       * Kernel::Point_Dim
       */
      std::unique_ptr<Geometry> inter; // null if no intersection

      if (Dim == 2) {
        const CGAL::Segment_2<Kernel> s1(l.pointN(i).toPoint_2(),
                                         l.pointN(i + 1).toPoint_2());
        const CGAL::Segment_2<Kernel> s2(l.pointN(j).toPoint_2(),
                                         l.pointN(j + 1).toPoint_2());
        const CGAL::Object            out = CGAL::intersection(s1, s2);

        if (out.is<Kernel::Point_2>()) {
          inter =
              std::make_unique<Point>(CGAL::object_cast<Kernel::Point_2>(out));
        } else if (out.is<Kernel::Segment_2>()) {
          const Kernel::Segment_2 &s =
              CGAL::object_cast<Kernel::Segment_2>(out);
          inter = std::make_unique<LineString>(s.point(0), s.point(1));
        }
      } else {
        const CGAL::Segment_3<Kernel> s1(l.pointN(i).toPoint_3(),
                                         l.pointN(i + 1).toPoint_3());
        const CGAL::Segment_3<Kernel> s2(l.pointN(j).toPoint_3(),
                                         l.pointN(j + 1).toPoint_3());
        const CGAL::Object            out = CGAL::intersection(s1, s2);

        if (out.is<Kernel::Point_3>()) {
          inter =
              std::make_unique<Point>(CGAL::object_cast<Kernel::Point_3>(out));
        } else if (out.is<Kernel::Segment_3>()) {
          const Kernel::Segment_3 &s =
              CGAL::object_cast<Kernel::Segment_3>(out);
          inter = std::make_unique<LineString>(s.point(0), s.point(1));
        }
      }

      if (inter.get() && inter->is<LineString>()) {
        return true; // segments overlap
      }
      if (inter.get() && inter->is<Point>() &&
          !(i + 1 == j) // one contact point between consecutive segments is ok
          && !((i == 0) && (j + 1 == numSegments) &&
               inter->as<Point>() == l.startPoint() &&
               inter->as<Point>() == l.endPoint())) {
        return true; // contact point that is not a contact between startPoint
                     // and endPoint
      }
    }
  }

  return false;
}

auto
selfIntersects(const LineString &l) -> bool
{
  return selfIntersectsImpl<2>(l);
}
auto
selfIntersects3D(const LineString &l) -> bool
{
  return selfIntersectsImpl<3>(l);
}

/// @private
template <int Dim>
auto
selfIntersectsImpl(const PolyhedralSurface &s, const SurfaceGraph &graph)
    -> bool
{
  size_t const numPatchs = s.numPatchs();

  for (size_t pi = 0; pi != numPatchs; ++pi) {
    for (size_t pj = pi + 1; pj < numPatchs; ++pj) {
      std::unique_ptr<Geometry> inter =
          Dim == 3 ? intersection3D(s.patchN(pi), s.patchN(pj))
                   : intersection(s.patchN(pi), s.patchN(pj));

      if (!inter->isEmpty()) {
        // two cases:
        // - neighbors can have a line as intersection
        // - non neighbors can only have a point or a set of points
        using Iterator = SurfaceGraph::FaceGraph::adjacency_iterator;
        std::pair<Iterator, Iterator> const neighbors =
            boost::adjacent_vertices(pi, graph.faceGraph());

        if (neighbors.second !=
            std::find(neighbors.first, neighbors.second, pj)) {
          // neighbor
          // std::cerr << pi << " " << pj << " neighbor\n";
          if (!inter->is<LineString>()) {
            return true;
          }
        } else {
          // not a neighbor
          // std::cerr << pi << " " << pj << " not neighbor\n";
          if (inter->dimension() != 0) {
            return true;
          }
        }
      }
    }
  }

  return false;
}

auto
selfIntersects(const PolyhedralSurface &s, const SurfaceGraph &g) -> bool
{
  return selfIntersectsImpl<2>(s, g);
}

auto
selfIntersects3D(const PolyhedralSurface &s, const SurfaceGraph &g) -> bool
{
  return selfIntersectsImpl<3>(s, g);
}

/// @private
template <int Dim>
auto
selfIntersectsImpl(const TriangulatedSurface &tin, const SurfaceGraph &graph)
    -> bool
{
  size_t const numPatchs = tin.numPatchs();

  for (size_t ti = 0; ti != numPatchs; ++ti) {
    for (size_t tj = ti + 1; tj < numPatchs; ++tj) {
      std::unique_ptr<Geometry> inter =
          Dim == 3 ? intersection3D(tin.patchN(ti), tin.patchN(tj))
                   : intersection(tin.patchN(ti), tin.patchN(tj));

      if (!inter->isEmpty()) {
        // two cases:
        // - neighbors can have a line as intersection
        // - non neighbors can only have a point or a set of points
        using Iterator = SurfaceGraph::FaceGraph::adjacency_iterator;
        std::pair<Iterator, Iterator> const neighbors =
            boost::adjacent_vertices(ti, graph.faceGraph());

        if (neighbors.second !=
            std::find(neighbors.first, neighbors.second, tj)) {
          // neighbor
          // std::cerr << ti << " " << tj << " neighbor\n";
          if (!inter->is<LineString>()) {
            return true;
          }
        } else {
          // not a neighbor
          // std::cerr << ti << " " << tj << " not neighbor\n";
          if (inter->dimension() != 0) {
            return true;
          }
        }
      }
    }
  }

  return false;
}

auto
selfIntersects(const TriangulatedSurface &tin, const SurfaceGraph &g) -> bool
{
  return selfIntersectsImpl<2>(tin, g);
}

auto
selfIntersects3D(const TriangulatedSurface &tin, const SurfaceGraph &g) -> bool
{
  return selfIntersectsImpl<3>(tin, g);
}

} // namespace SFCGAL::algorithm
