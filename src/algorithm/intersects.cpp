// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
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
_intersects(const PrimitiveHandle<2> &handle1,
            const PrimitiveHandle<2> &handle2) -> bool
{
  //
  // Point vs. Point
  //

  if (handle1.handle.which() == PrimitivePoint &&
      handle2.handle.which() == PrimitivePoint) {
    return *boost::get<const CGAL::Point_2<Kernel> *>(handle1.handle) ==
           *boost::get<const CGAL::Point_2<Kernel> *>(handle2.handle);
  }

  //
  // Segment vs. Point
  //

  if (handle1.handle.which() == PrimitiveSegment &&
      handle2.handle.which() == PrimitivePoint) {
    const auto *seg = handle1.as<CGAL::Segment_2<Kernel>>();
    const auto *pt  = handle2.as<CGAL::Point_2<Kernel>>();
    return seg->has_on(*pt);
  }

  //
  // Segment vs. Segment
  //

  if (handle1.handle.which() == PrimitiveSegment &&
      handle2.handle.which() == PrimitiveSegment) {
    const auto *seg1 = handle1.as<CGAL::Segment_2<Kernel>>();
    const auto *seg2 = handle2.as<CGAL::Segment_2<Kernel>>();
    return CGAL::do_intersect(*seg1, *seg2);
  }

  //
  // Polygon vs. Point
  //

  if (handle1.handle.which() == PrimitiveSurface &&
      handle2.handle.which() == PrimitivePoint) {
    // Polygon versus Point
    const auto *poly = handle1.as<CGAL::Polygon_with_holes_2<Kernel>>();
    const auto *pt   = handle2.as<CGAL::Point_2<Kernel>>();

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

  if (handle1.handle.which() == PrimitiveSurface &&
      handle2.handle.which() == PrimitiveSegment) {
    const auto *poly = handle1.as<CGAL::Polygon_with_holes_2<Kernel>>();
    const auto *seg  = handle2.as<CGAL::Segment_2<Kernel>>();

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

  if (handle1.handle.which() == PrimitiveSurface &&
      handle2.handle.which() == PrimitiveSurface) {
    const auto *poly1 = handle1.as<CGAL::Polygon_with_holes_2<Kernel>>();
    const auto *poly2 = handle2.as<CGAL::Polygon_with_holes_2<Kernel>>();

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
_intersects(const PrimitiveHandle<3> &handle1,
            const PrimitiveHandle<3> &handle2) -> bool
{
  if (handle1.handle.which() == PrimitivePoint &&
      handle2.handle.which() == PrimitivePoint) {
    return *boost::get<const CGAL::Point_3<Kernel> *>(handle1.handle) ==
           *boost::get<const CGAL::Point_3<Kernel> *>(handle2.handle);
  }
  if (handle1.handle.which() == PrimitiveSegment &&
      handle2.handle.which() == PrimitivePoint) {
    const auto *seg = handle1.as<CGAL::Segment_3<Kernel>>();
    const auto *pt  = handle2.as<CGAL::Point_3<Kernel>>();
    return seg->has_on(*pt);
  }
  if (handle1.handle.which() == PrimitiveSegment &&
      handle2.handle.which() == PrimitiveSegment) {
    const auto *sega = handle1.as<CGAL::Segment_3<Kernel>>();
    const auto *segb = handle2.as<CGAL::Segment_3<Kernel>>();
    return CGAL::do_intersect(*sega, *segb);
  }

  if (handle1.handle.which() == PrimitiveVolume) {
    intersects_volume_x visitor(handle1.as<MarkedPolyhedron>());
    return boost::apply_visitor(visitor, handle2.handle);
  }

  if (handle1.handle.which() == PrimitiveSurface &&
      handle2.handle.which() == PrimitivePoint) {
    const auto *tri = handle1.as<CGAL::Triangle_3<Kernel>>();
    const auto *pt  = handle2.as<CGAL::Point_3<Kernel>>();
    return tri->has_on(*pt);
  }

  if (handle1.handle.which() == PrimitiveSurface &&
      handle2.handle.which() == PrimitiveSegment) {
    const auto *tri = handle1.as<CGAL::Triangle_3<Kernel>>();
    const auto *seg = handle2.as<CGAL::Segment_3<Kernel>>();
    return CGAL::do_intersect(*tri, *seg);
  }

  if (handle1.handle.which() == PrimitiveSurface &&
      handle2.handle.which() == PrimitiveSurface) {
    const auto *tri1 = handle1.as<CGAL::Triangle_3<Kernel>>();
    const auto *tri2 = handle2.as<CGAL::Triangle_3<Kernel>>();
    return CGAL::do_intersect(*tri1, *tri2);
  }

  return false;
}

//
// We deal here with symmetric call
template <int Dim>
auto
dispatch_intersects_sym(const PrimitiveHandle<Dim> &handle1,
                        const PrimitiveHandle<Dim> &handle2) -> bool
{
  // assume types are ordered by dimension within the boost::variant
  if (handle1.handle.which() >= handle2.handle.which()) {
    return _intersects(handle1, handle2);
  }
  return _intersects(handle2, handle1);
}

struct found_an_intersection {};

template <int Dim>
struct intersects_callback {
  void
  operator()(const typename PrimitiveBox<Dim>::Type &box1,
             const typename PrimitiveBox<Dim>::Type &box2)
  {
    if (dispatch_intersects_sym(*box1.handle(), *box2.handle())) {
      throw found_an_intersection();
    }
  }
};

/// @private
template <int Dim>
auto
intersects(const PrimitiveHandle<Dim> &handle1,
           const PrimitiveHandle<Dim> &handle2) -> bool
{
  return dispatch_intersects_sym(handle1, handle2);
}

/// @private
template <int Dim>
auto
intersects(const GeometrySet<Dim> &geometrySet1,
           const GeometrySet<Dim> &geometrySet2) -> bool
{
  typename SFCGAL::detail::HandleCollection<Dim>::Type ahandles;
  typename SFCGAL::detail::HandleCollection<Dim>::Type bhandles;
  typename SFCGAL::detail::BoxCollection<Dim>::Type    aboxes;
  typename SFCGAL::detail::BoxCollection<Dim>::Type    bboxes;
  geometrySet1.computeBoundingBoxes(ahandles, aboxes);
  geometrySet2.computeBoundingBoxes(bhandles, bboxes);

  try {
    intersects_callback<Dim> const callback;
    CGAL::box_intersection_d(aboxes.begin(), aboxes.end(), bboxes.begin(),
                             bboxes.end(), callback);
  } catch (found_an_intersection &e) {
    return true;
  }

  return false;
}

/// @private
template bool
intersects<2>(const GeometrySet<2> &geometrySet1,
              const GeometrySet<2> &geometrySet2);
/// @private
template bool
intersects<3>(const GeometrySet<3> &geometrySet1,
              const GeometrySet<3> &geometrySet2);

/// @private
template bool
intersects<2>(const PrimitiveHandle<2> &handle1,
              const PrimitiveHandle<2> &handle2);
/// @private
template bool
intersects<3>(const PrimitiveHandle<3> &handle1,
              const PrimitiveHandle<3> &handle2);

/// @private
auto
intersects(const Geometry &geometry1, const Geometry &geometry2) -> bool
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(geometry1);
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(geometry2);

  GeometrySet<2> const gsa(geometry1);
  GeometrySet<2> const gsb(geometry2);

  return intersects(gsa, gsb);
}

/// @private
auto
intersects3D(const Geometry &geometry1, const Geometry &geometry2) -> bool
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_3D(geometry1);
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_3D(geometry2);

  GeometrySet<3> const gsa(geometry1);
  GeometrySet<3> const gsb(geometry2);

  return intersects(gsa, gsb);
}

/// @private
auto
intersects(const Geometry &geometry1, const Geometry &geometry2,
           NoValidityCheck /*unused*/) -> bool
{
  GeometrySet<2> const gsa(geometry1);
  GeometrySet<2> const gsb(geometry2);

  return intersects(gsa, gsb);
}

/// @private
auto
intersects3D(const Geometry &geometry1, const Geometry &geometry2,
             NoValidityCheck /*unused*/) -> bool
{
  GeometrySet<3> const gsa(geometry1);
  GeometrySet<3> const gsb(geometry2);

  return intersects(gsa, gsb);
}

/// @private
template <int Dim>
auto
selfIntersectsImpl(const LineString &lineString) -> bool
{

  if (lineString.numSegments() < 2) {
    return false; // one segment cannot intersect
  }

  // note: zero length segments are a pain, to avoid algorithm complexity
  // we start by filtering them out
  const size_t numPoints = lineString.numPoints();
  LineString   l;

  for (size_t i = 0; i != numPoints; ++i) {
    if (i == 0 || l.endPoint() != lineString.pointN(i)) {
      l.addPoint(lineString.pointN(i));
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

/// @private
auto
selfIntersects(const LineString &lineString) -> bool
{
  return selfIntersectsImpl<2>(lineString);
}
/// @private
auto
selfIntersects3D(const LineString &lineString) -> bool
{
  return selfIntersectsImpl<3>(lineString);
}

/// @private
template <int Dim>
auto
selfIntersectsImpl(const PolyhedralSurface &surface, const SurfaceGraph &graph)
    -> bool
{
  size_t const numPatches = surface.numPatches();

  for (size_t pi = 0; pi != numPatches; ++pi) {
    GeometrySet<Dim> geomSetI(surface.patchN(pi));
    for (size_t pj = pi + 1; pj < numPatches; ++pj) {
      GeometrySet<Dim> geomSetJ(surface.patchN(pj));
      bool             hasIntersection = intersects<Dim>(geomSetI, geomSetJ);
      if (hasIntersection) {
        GeometrySet<Dim> geomSetInter;
        // call intersection() version without recompose:
        intersection<Dim>(geomSetI, geomSetJ, geomSetInter);

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

          // apply recompose code from intersection(const Geometry, const
          // Geometry, NoValidityCheck):
          GeometrySet<Dim> filtered;
          geomSetInter.filterCovered(filtered);
          std::unique_ptr<Geometry> interRecomposed = filtered.recompose();
          if (!interRecomposed->is<LineString>()) {
            return true;
          }
        } else {
          // not a neighbor
          // std::cerr << pi << " " << pj << " not neighbor\n";
          if (geomSetInter.dimension() != 0) {
            return true;
          }
        }
      }
    }
  }

  return false;
}

/// @private
auto
selfIntersects(const PolyhedralSurface &surface, const SurfaceGraph &graph)
    -> bool
{
  return selfIntersectsImpl<2>(surface, graph);
}

/// @private
auto
selfIntersects3D(const PolyhedralSurface &surface, const SurfaceGraph &graph)
    -> bool
{
  return selfIntersectsImpl<3>(surface, graph);
}

/// @private
template <int Dim>
auto
selfIntersectsImpl(const TriangulatedSurface &triangulatedSurface,
                   const SurfaceGraph        &graph) -> bool
{
  size_t const numPatches = triangulatedSurface.numPatches();

  for (size_t ti = 0; ti != numPatches; ++ti) {
    GeometrySet<Dim> geomSetI(triangulatedSurface.patchN(ti));
    for (size_t tj = ti + 1; tj < numPatches; ++tj) {
      GeometrySet<Dim> geomSetJ(triangulatedSurface.patchN(tj));
      bool             hasIntersection = intersects<Dim>(geomSetI, geomSetJ);

      if (hasIntersection) {
        GeometrySet<Dim> geomSetInter;
        // call intersection() version without recompose:
        intersection<Dim>(geomSetI, geomSetJ, geomSetInter);

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

          // apply recompose code from intersection(const Geometry, const
          // Geometry, NoValidityCheck):
          GeometrySet<Dim> filtered;
          geomSetInter.filterCovered(filtered);
          std::unique_ptr<Geometry> interRecomposed = filtered.recompose();
          if (!interRecomposed->is<LineString>()) {
            return true;
          }
        } else {
          // not a neighbor
          // std::cerr << ti << " " << tj << " not neighbor\n";
          if (geomSetInter.dimension() != 0) {
            return true;
          }
        }
      }
    }
  }

  return false;
}

/// @private
auto
selfIntersects(const TriangulatedSurface &triangulatedSurface,
               const SurfaceGraph        &graph) -> bool
{
  return selfIntersectsImpl<2>(triangulatedSurface, graph);
}

/// @private
auto
selfIntersects3D(const TriangulatedSurface &triangulatedSurface,
                 const SurfaceGraph        &graph) -> bool
{
  return selfIntersectsImpl<3>(triangulatedSurface, graph);
}

} // namespace SFCGAL::algorithm
