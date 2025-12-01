// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_DIFFERENCEPRIMITIVES_H_
#define SFCGAL_ALGORITHM_DIFFERENCEPRIMITIVES_H_

#include "SFCGAL/Exception.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/detail/GeometrySet.h"
#include "SFCGAL/triangulate/triangulatePolygon.h"

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/box_intersection_d.h>

namespace SFCGAL::algorithm {

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection
using NoVolume         = detail::NoVolume;
using MarkedPolyhedron = detail::MarkedPolyhedron;

/**
 * @brief Compute intersection of two CGAL triangles
 * @param triangle1 First triangle
 * @param triangle2 Second triangle
 * @return CGAL::Object containing the intersection result
 */
auto
intersection(const CGAL::Triangle_3<Kernel> &triangle1,
             const CGAL::Triangle_3<Kernel> &triangle2) -> CGAL::Object;

inline auto
do_intersect(const Point_2 &point, const Polygon_with_holes_2 &polygon) -> bool
{
  // point intersects if it's inside the ext ring and outside all holes

  if (CGAL::bounded_side_2(polygon.outer_boundary().vertices_begin(),
                           polygon.outer_boundary().vertices_end(), point,
                           Kernel()) == CGAL::ON_UNBOUNDED_SIDE) {
    return false;
  }

  for (auto hit = polygon.holes_begin(); hit != polygon.holes_end(); ++hit) {
    if (CGAL::bounded_side_2(hit->vertices_begin(), hit->vertices_end(), point,
                             Kernel()) != CGAL::ON_UNBOUNDED_SIDE) {
      return false;
    }
  }

  return true;
}

template <typename PointOutputIteratorType>
auto
difference(const Point_2 &a, const Point_2 &b, PointOutputIteratorType out)
    -> PointOutputIteratorType
{
  if (a != b) {
    *out++ = a;
  }

  return out;
}

template <typename PointOutputIteratorType>
auto
difference(const Point_2 &a, const Segment_2 &b, PointOutputIteratorType out)
    -> PointOutputIteratorType
{
  if (!CGAL::do_intersect(a, b)) {
    *out++ = a;
  }

  return out;
}

template <typename PointOutputIteratorType>
auto
difference(const Point_2 &a, const Polygon_with_holes_2 &b,
           PointOutputIteratorType out) -> PointOutputIteratorType
{
  if (!do_intersect(a, b)) {
    *out++ = a;
  }

  return out;
}

template <typename PointOutputIteratorType>
auto
difference(const Point_2 & /*unused*/, const NoVolume & /*unused*/,
           PointOutputIteratorType out) -> PointOutputIteratorType
{
  BOOST_ASSERT(false);
  return out;
}

template <typename SegmentOutputIteratorType>
auto
difference(const Segment_2 & /*unused*/, const NoVolume & /*unused*/,
           SegmentOutputIteratorType out) -> SegmentOutputIteratorType
{
  BOOST_ASSERT(false);
  return out;
}

template <typename SurfaceOutputIteratorType>
auto
difference(const Polygon_with_holes_2 & /*unused*/, const NoVolume & /*unused*/,
           SurfaceOutputIteratorType out) -> SurfaceOutputIteratorType
{
  BOOST_ASSERT(false);
  return out;
}

template <typename PointOutputIteratorType>
auto
difference(const Point_3 &a, const Point_3 &b, PointOutputIteratorType out)
    -> PointOutputIteratorType
{
  if (a != b) {
    *out++ = a;
  }

  return out;
}

template <typename PointOutputIteratorType>
auto
difference(const Point_3 &a, const Segment_3 &b, PointOutputIteratorType out)
    -> PointOutputIteratorType
{
  if (!b.has_on(a)) {
    *out++ = a;
  }

  return out;
}

template <typename PointOutputIteratorType>
auto
difference(const Point_3 &a, const Triangle_3 &b, PointOutputIteratorType out)
    -> PointOutputIteratorType
{
  if (!b.has_on(a)) {
    *out++ = a;
  }

  return out;
}

template <typename PointOutputIteratorType>
auto
difference(const Point_3 &a, const MarkedPolyhedron &b,
           PointOutputIteratorType out) -> PointOutputIteratorType
{
  CGAL::Side_of_triangle_mesh<MarkedPolyhedron, Kernel> is_in_poly(b);

  if (CGAL::ON_UNBOUNDED_SIDE == is_in_poly(a)) {
    *out++ = a;
  }

  return out;
}

template <typename SegmentType, typename SegmentOrSurfaceType,
          typename SegmentOutputIteratorType>
auto
difference(const SegmentType &segA, const SegmentOrSurfaceType &segB,
           SegmentOutputIteratorType out) -> SegmentOutputIteratorType
{
  CGAL::Object       inter    = CGAL::intersection(segA, segB);
  const SegmentType *interSeg = CGAL::object_cast<SegmentType>(&inter);

  if (interSeg !=
      nullptr) { // there maybe zero, one or two segments as a result
    if (CGAL::squared_distance(segA.source(), interSeg->source()) <
        CGAL::squared_distance(segA.source(), interSeg->target())) {
      if (segA.source() != interSeg->source()) {
        *out++ = SegmentType(segA.source(), interSeg->source());
      }

      if (interSeg->target() != segA.target()) {
        *out++ = SegmentType(interSeg->target(), segA.target());
      }
    } else {
      if (segA.source() != interSeg->target()) {
        *out++ = SegmentType(segA.source(), interSeg->target());
      }

      if (interSeg->source() != segA.target()) {
        *out++ = SegmentType(interSeg->source(), segA.target());
      }
    }
  } else { // intersection is a point or nothing, leave a unchanged
    *out++ = segA;
  }

  return out;
}

template <typename PointType>
struct Nearer {
  explicit Nearer(const PointType &reference) : _ref(reference) {}
  auto
  operator()(const PointType &lhs, const PointType &rhs) const -> bool
  {
    return CGAL::squared_distance(_ref, lhs) <
           CGAL::squared_distance(_ref, rhs);
  }

private:
  const PointType _ref;
};

// NOLINTBEGIN(readability-function-cognitive-complexity)
template <typename SegmentOutputIteratorType>
auto
difference(const Segment_2 &segment, const Polygon_with_holes_2 &polygon,
           SegmentOutputIteratorType out) -> SegmentOutputIteratorType
{
  // we could triangulate the polygon and substract each triangle
  //
  // we could also cut the line by polygon contours and test if the middle of
  // the segment is inside but if the segment lies on the contour it's a special
  // case we first substract the contours to take care of this special case, we
  // obtain a vector of segments, for each segment of this vector, we subdivide
  // it with the intersection points with the rings once done, we check, for
  // each subdivision that has distinct end-points if the middle is in or out.

  std::vector<Segment_2> result(1, segment);
  std::vector<Polygon_2> rings(1, polygon.outer_boundary());
  rings.insert(rings.end(), polygon.holes_begin(), polygon.holes_end());

  for (auto &ring : rings) {
    for (auto target = ring.vertices_begin(); target != ring.vertices_end();
         ++target) {
      const Segment_2 segmentCut(target == ring.vertices_begin()
                                     ? *(ring.vertices_end() - 1)
                                     : *(target - 1),
                                 *target);
      std::vector<Segment_2> tmp;

      for (const auto &seg : result) {
        difference(seg, segmentCut, std::back_inserter(tmp));
      }

      tmp.swap(result);
    }
  }

  for (const auto &seg : result) {
    std::vector<Point_2> points(1, seg.source());

    for (auto &ring : rings) {
      for (auto target = ring.vertices_begin(); target != ring.vertices_end();
           ++target) {
        Segment_2 segmentCut(target == ring.vertices_begin()
                                 ? *(ring.vertices_end() - 1)
                                 : *(target - 1),
                             *target);
        CGAL::Object inter = CGAL::intersection(seg, segmentCut);
        const auto  *point = CGAL::object_cast<Point_2>(&inter);

        if (point != nullptr) {
          points.push_back(*point);
        }
      }
    }

    points.push_back(seg.target());
    // order point according to the distance from source
    const Nearer<Point_2> nearer(seg.source());
    std::sort(points.begin() + 1, points.end() - 1, nearer);

    // append segments that has length and wich midpoint is outside polygon to
    // result
    for (auto pit = points.begin(); pit != points.end() - 1; ++pit) {
      auto qit = pit + 1;

      if (*pit != *qit && !do_intersect(CGAL::midpoint(*pit, *qit), polygon)) {
        *out++ = Segment_2(*pit, *qit);
      }
    }
  }

  return out;
}
// NOLINTEND(readability-function-cognitive-complexity)

// assuming two disjoint (except at a point) polygon, test if the first a hole
// of the second
inline auto
isHoleOf(const Polygon_2 &hole, const Polygon_2 &poly) -> bool
{
  return CGAL::bounded_side_2(poly.vertices_begin(), poly.vertices_end(),
                              *hole.vertices_begin(),
                              Kernel()) == CGAL::ON_BOUNDED_SIDE ||
         CGAL::bounded_side_2(poly.vertices_begin(), poly.vertices_end(),
                              *(hole.vertices_begin() + 1),
                              Kernel()) == CGAL::ON_BOUNDED_SIDE;
}

// NOLINTBEGIN(readability-function-cognitive-complexity)
template <typename PolygonOutputIteratorType>
auto
fix_cgal_valid_polygon(const Polygon_with_holes_2 &poly,
                       PolygonOutputIteratorType   out)
    -> PolygonOutputIteratorType
{
  const Polygon_2 &outer = poly.outer_boundary();

  // std::cerr << "in fix outer " << outer << "\n";
  if (!outer.is_simple()) {
    // the holes are simple, so we need to find the intersection points, then
    // split the outer ring into simple components and put holes in the right
    // one note that the hole may touch the outer boundary at a point, so if the
    // tested hole point falls on the boundary, we test the next

    std::vector<Polygon_2>            boundaries;
    std::vector<std::vector<Point_2>> stack(1);

    for (auto vit = outer.vertices_begin(); vit != outer.vertices_end();
         ++vit) {
      if (!stack.back().empty() && stack.back()[0] == *vit) { // closing ring
        boundaries.emplace_back(stack.back().begin(), stack.back().end());
        stack.pop_back();
      } else if (std::find(vit + 1, outer.vertices_end(), *vit) !=
                 outer.vertices_end()) { // split point
        stack.back().push_back(*vit);
        stack.resize(stack.size() + 1);
        stack.back().push_back(*vit);
      } else {
        stack.back().push_back(*vit);
      }
    }

    if (!stack.empty()) {
      boundaries.emplace_back(stack.back().begin(), stack.back().end());
    }

    // std::cerr << "in fix boundaries " << boundaries.size() << "\n";

    std::vector<Polygon_2> holes(poly.holes_begin(), poly.holes_end());

    // one of the boundaries may be a hole
    std::vector<Polygon_2> clockwiseBoundaries;
    std::vector<Polygon_2> counterClockwiseBoundaries;

    for (const auto &boundary : boundaries) {
      if (boundary.orientation() == CGAL::CLOCKWISE) {
        clockwiseBoundaries.push_back(boundary);
      } else {
        counterClockwiseBoundaries.push_back(boundary);
      }
    }

    // if we have holes, check the orientation of the first hole to see
    // what is a hole orientation
    // if we don't have holes, we test if the first ccw is a hole of any
    // of the cw, if not, then the other are holes
    bool holesAreCCW = false;

    if (clockwiseBoundaries.empty()) {
      holesAreCCW = false;
    } else if (counterClockwiseBoundaries.empty()) {
      holesAreCCW = true;
    } else if (!holes.empty()) {
      holesAreCCW = holes[0].orientation() != CGAL::CLOCKWISE;
    } else {
      for (const auto &cwBoundary : clockwiseBoundaries) {
        if (isHoleOf(counterClockwiseBoundaries[0], cwBoundary)) {
          holesAreCCW = true;
          break;
        }
      }
    }

    if (holesAreCCW) {
      holes.insert(holes.end(), counterClockwiseBoundaries.begin(),
                   counterClockwiseBoundaries.end());
      boundaries.swap(clockwiseBoundaries);
    } else {
      holes.insert(holes.end(), clockwiseBoundaries.begin(),
                   clockwiseBoundaries.end());
      boundaries.swap(counterClockwiseBoundaries);
    }

    std::vector<std::vector<Polygon_2>> sortedHoles(
        boundaries.size()); // 1/1 with boundaries

    for (const auto &hole : holes) {
      for (size_t idx = 0; idx < boundaries.size(); ++idx) {
        if (isHoleOf(hole, boundaries[idx])) {
          sortedHoles[idx].push_back(hole);
        }
      }
    }

    for (unsigned idx = 0; idx < boundaries.size(); idx++) {
      *out++ = Polygon_with_holes_2(boundaries[idx], sortedHoles[idx].begin(),
                                    sortedHoles[idx].end());
    }
  } else {
    *out++ = poly;
  }

  return out;
}
// NOLINTEND(readability-function-cognitive-complexity)

inline auto
fix_sfs_valid_polygon(const Polygon_with_holes_2 &poly) -> Polygon_with_holes_2
{
  CGAL::Gps_segment_traits_2<Kernel> traits;

  if (are_holes_and_boundary_pairwise_disjoint(poly, traits)) {
    return poly;
  }

  // a polygon is valid for sfs and invalid for CGAL when two rings intersect
  // on a point that is not a ring vertex
  // we add this vertex to fix the polygon
  // for each ring segment
  //    for every other ring point
  //        add point to segment

  // put all rings in a vector to avoid distinction between outer and holes
  std::vector<Polygon_2> rings(1, poly.outer_boundary());
  rings.insert(rings.end(), poly.holes_begin(), poly.holes_end());

  std::vector<Polygon_2> out;

  for (size_t ringIdx = 0; ringIdx < rings.size(); ++ringIdx) {
    out.emplace_back();

    for (auto target = rings[ringIdx].vertices_begin();
         target != rings[ringIdx].vertices_end(); ++target) {
      Segment_2 segment(target == rings[ringIdx].vertices_begin()
                            ? *(rings[ringIdx].vertices_end() - 1)
                            : *(target - 1),
                        *target);

      // for every other ring
      for (size_t otherIdx = 0; otherIdx < rings.size(); ++otherIdx) {
        if (ringIdx == otherIdx) {
          continue;
        }

        for (auto vertex = rings[otherIdx].vertices_begin();
             vertex != rings[otherIdx].vertices_end(); ++vertex) {
          if (CGAL::do_intersect(*vertex, segment)) {
            out.back().push_back(*vertex);
          }
        }
      }

      out.back().push_back(*target);
    }
  }

  return {out[0], out.begin() + 1, out.end()};
}

template <typename OutputIteratorType>
auto
difference(const Triangle_3 &triangleP, const Triangle_3 &triangleQ,
           OutputIteratorType out) -> OutputIteratorType
{

  const Plane_3 plane = triangleP.supporting_plane();

  if (plane != triangleQ.supporting_plane() ||
      !CGAL::do_intersect(triangleP, triangleQ)) {
    *out++ = triangleP;
  } else {
    // project on plane
    // difference between polygons
    // triangulate the result

    Polygon_with_holes_2 pProj;
    Polygon_with_holes_2 qProj;

    for (int idx = 0; idx < 3; idx++) {
      pProj.outer_boundary().push_back(plane.to_2d(triangleP.vertex(idx)));
      qProj.outer_boundary().push_back(plane.to_2d(triangleQ.vertex(idx)));
    }

    std::vector<Polygon_with_holes_2> res;
    difference(pProj, qProj, std::back_inserter(res));

    for (const auto &pwh : res) {
      const Polygon       poly(pwh);
      TriangulatedSurface triangSurface;
      triangulate::triangulatePolygon3D(poly, triangSurface);

      for (auto &triangleElem : triangSurface) {
        *out++ = Triangle_3(plane.to_3d(triangleElem.vertex(0).toPoint_2()),
                            plane.to_3d(triangleElem.vertex(1).toPoint_2()),
                            plane.to_3d(triangleElem.vertex(2).toPoint_2()));
      }
    }
  }

  return out;
}

template <typename VolumeOutputIteratorType>
auto
difference(const MarkedPolyhedron &polyA, const MarkedPolyhedron &polyB,
           VolumeOutputIteratorType out) -> VolumeOutputIteratorType
{
  MarkedPolyhedron meshP = polyA;
  MarkedPolyhedron meshQ = polyB;
  if (CGAL::Polygon_mesh_processing::corefine_and_compute_difference(
          meshP, meshQ, meshP)) {
    if (std::next(vertices(meshP).first) != vertices(meshP).second) {
      *out++ = meshP;
    }
  }
  return out;
}

using FaceBboxBase = CGAL::Box_intersection_d::Box_with_handle_d<
    double, 3, MarkedPolyhedron::Halfedge_around_facet_const_circulator>;

struct FaceBbox : FaceBboxBase {
  /**
   * @brief Axis-aligned bounding box for a polygon face.
   *
   * Computes and stores the 3D bounding box that encompasses all vertices
   * of a face represented by a halfedge circulator.
   */
  struct Bbox : CGAL::Bbox_3 {
    /**
     * @brief Constructor that computes the bounding box of a face.
     * @param handle Halfedge circulator pointing to the face vertices.
     */
    Bbox(MarkedPolyhedron::Halfedge_around_facet_const_circulator handle)
        : CGAL::Bbox_3(handle->vertex()->point().bbox())
    {
      const MarkedPolyhedron::Halfedge_around_facet_const_circulator end =
          handle;

      do {
        // @note with CGAL 4.5 you would write simply
        // *this += (++handle)->vertex()->point().bbox();
        this->CGAL::Bbox_3::operator=(*this +
                                      (++handle)->vertex()->point().bbox());
      } while (handle != end);
    }
  };

  FaceBbox(const MarkedPolyhedron::Facet &facet)
      : FaceBboxBase(Bbox(facet.facet_begin()), facet.facet_begin())
  {
  }
};

struct FaceSegmentCollide {
  using CollisionVector =
      std::vector<MarkedPolyhedron::Halfedge_around_facet_const_circulator>;
  explicit FaceSegmentCollide(CollisionVector &list) : _list(list) {}
  auto
  operator()(const FaceBboxBase & /*unused*/, const FaceBboxBase &face) -> void
  {
    _list.push_back(face.handle());
  }

private:
  CollisionVector &_list;
};

template <typename TriangleOutputIteratorType>
auto
collidingTriangles(const FaceSegmentCollide::CollisionVector &collisions,
                   TriangleOutputIteratorType out) -> TriangleOutputIteratorType
{
  for (const auto &collision : collisions) {
    auto               circulator = collision;
    std::vector<Point> points(1, circulator->vertex()->point());

    do {
      points.emplace_back((++circulator)->vertex()->point());
    } while (circulator != collision);

    if (points.size() == 3) {
      *out++ = Triangle_3(points[0].toPoint_3(), points[1].toPoint_3(),
                          points[2].toPoint_3());
    } else {
      const Polygon       poly(points);
      TriangulatedSurface triangSurface;
      triangulate::triangulatePolygon3D(poly, triangSurface);

      for (auto &triangleElem : triangSurface) {
        *out++ = Triangle_3(triangleElem.vertex(0).toPoint_3(),
                            triangleElem.vertex(1).toPoint_3(),
                            triangleElem.vertex(2).toPoint_3());
      }
    }
  }

  return out;
}

// NOLINTBEGIN(readability-function-cognitive-complexity)
template <typename SegmentOutputIteratorType>
auto
difference(const Segment_3 &segment, const MarkedPolyhedron &polyhedron,
           SegmentOutputIteratorType out) -> SegmentOutputIteratorType
{
  // this is a bit of a pain
  // the algo should follow the same lines as the Segment_2 -
  // Polygon_with_holes_2 namely, remove the pieces of the segment were it
  // touches facets, then compute the intersections with facets to cut the
  // segments and create segments for output were the middle point is inside
  //
  // to speed thing up we put facets in AABB-Tree

  std::vector<FaceBbox>     bboxes(polyhedron.facets_begin(),
                                   polyhedron.facets_end());
  std::vector<FaceBboxBase> bbox(
      1, FaceBboxBase(segment.bbox(),
                      polyhedron.facets_begin()->facet_begin())); // nevermind
                                                                  // the facet
                                                                  // handle,
                                                                  // it's not
                                                                  // used anyway
  FaceSegmentCollide::CollisionVector collisions;
  FaceSegmentCollide                  cb(collisions);
  CGAL::box_intersection_d(bbox.begin(), bbox.end(), bboxes.begin(),
                           bboxes.end(), cb);

  if (collisions.empty()) {
    // completely in or out, we just test one point
    CGAL::Side_of_triangle_mesh<MarkedPolyhedron, Kernel> is_in_poly(
        polyhedron);

    if (CGAL::ON_UNBOUNDED_SIDE == is_in_poly(segment.source())) {
      *out++ = segment;
    }
  } else {
    std::vector<Triangle_3> triangles;
    collidingTriangles(collisions, std::back_inserter(triangles));

    // first step, substract faces
    std::vector<Segment_3> res1(1, segment);

    for (const auto &tri : triangles) {
      std::vector<Segment_3> tmp;

      for (const auto &seg : res1) {
        difference(seg, tri, std::back_inserter(tmp));
      }

      res1.swap(tmp);
    }

    // second step, for each segment, add intersection points and test each
    // middle point to know if it's in or out
    for (const auto &seg : res1) {
      std::vector<Point_3> points(1, seg.source());

      for (const auto &tri : triangles) {
        CGAL::Object inter = CGAL::intersection(seg, tri);
        const auto  *point = CGAL::object_cast<Point_3>(&inter);

        if (point != nullptr) {
          points.push_back(*point);
        }
      }

      points.push_back(seg.target());
      // order point according to the distance from source

      const Nearer<Point_3> nearer(seg.source());
      std::sort(points.begin() + 1, points.end() - 1, nearer);

      CGAL::Side_of_triangle_mesh<MarkedPolyhedron, Kernel> is_in_poly(
          polyhedron);

      // append segments that has length and wich midpoint is outside polyhedron
      // to result
      for (auto pit = points.begin(); pit != points.end() - 1; ++pit) {
        auto qit = pit + 1;

        if (*pit != *qit &&
            CGAL::ON_UNBOUNDED_SIDE == is_in_poly(CGAL::midpoint(*pit, *qit))) {
          *out++ = Segment_3(*pit, *qit);
        }
      }
    }
  }

  return out;
}
// NOLINTEND(readability-function-cognitive-complexity)

// @TODO put that in a proper header
auto
_intersection_solid_triangle(const MarkedPolyhedron &pa, const Triangle_3 &tri,
                             detail::GeometrySet<3> &output) -> void;

template <typename TriangleOutputIteratorType>
auto
difference(const Triangle_3 &triangle, const MarkedPolyhedron &polyhedron,
           TriangleOutputIteratorType out) -> TriangleOutputIteratorType
{
  std::vector<Triangle_3> inter;
  // call _intersection_solid_triangle
  detail::GeometrySet<3> interSet;
  _intersection_solid_triangle(polyhedron, triangle, interSet);

  for (const auto &surface : interSet.surfaces()) {
    inter.push_back(surface.primitive());
  }

  std::vector<Triangle_3> res(1, triangle);

  // GOTCHA for intersection points (volume touching triangle) , need to
  // retriangulate
  for (const auto &pointElem : interSet.points()) {
    std::vector<Triangle_3> tmp;

    for (const auto &tri : res) {
      const Point_3 &point = pointElem.primitive();

      for (int side = 0; side < 3; side++) {
        if (point != tri.vertex(side) && point != tri.vertex((side + 1) % 3) &&
            Segment_3(tri.vertex(side), tri.vertex((side + 1) % 3))
                .has_on(point)) {
          tmp.emplace_back(tri.vertex(side), point, tri.vertex((side + 2) % 3));
          tmp.emplace_back(point, tri.vertex((side + 1) % 3),
                           tri.vertex((side + 2) % 3));
        }
      }
    }

    tmp.swap(res);
  }

  for (const auto &interTri : inter) {
    std::vector<Triangle_3> tmp;

    for (const auto &tri : res) {
      difference(tri, interTri, std::back_inserter(tmp));
    }

    tmp.swap(res);
  }

  for (const auto &tri : res) {
    *out++ = tri;
  }

  return out;

  /*
      std::vector< FaceBbox > bboxes(polyhedron.facets_begin(),
     polyhedron.facets_end() ); std::vector< FaceBboxBase > bbox( 1,
     FaceBboxBase(triangle.bbox(),polyhedron.facets_begin()->facet_begin()) );
     // nevermind the facet handle, it's not used anyway
      FaceSegmentCollide::CollisionVector collisions;
      FaceSegmentCollide cb(collisions);
      CGAL::box_intersection_d( bbox.begin(), bbox.end(),
                                bboxes.begin(), bboxes.end(),
                                cb );

      if ( !collisions.size() ){
          // completely in or out, we just test one point
          CGAL::Point_inside_polyhedron_3<MarkedPolyhedron, Kernel> is_in_poly(
     polyhedron ); if ( CGAL::ON_UNBOUNDED_SIDE == is_in_poly(
     triangle.vertex(0) ) ) *out++ = triangle;
      }
      else {
          // now we first transform bboxes colliding faces into triangles
          // then we test for intersection and store resulting segments in a
     vector
          // we also store resulting polygons as segments
          //
          // we need to convert the resulting segments to a multipolygon of sort
          //
          // finally we triangulate the result and substract those triangles
          //
          std::vector< Triangle_3 > interTriangles;
          collidingTriangles( collisions, std::back_inserter( interTriangles )
     );

          std::vector< Segment_3 > intersectionCountours;

          BOOST_THROW_EXCEPTION(NotImplementedException("Triangle_3 - Volume is
     not implemented") );
      }
      return out;
  */
}

template <typename PolygonOutputIteratorType>
auto
difference(const Polygon_with_holes_2 &a, const Polygon_with_holes_2 &b,
           PolygonOutputIteratorType out) -> PolygonOutputIteratorType
{
  CGAL::Gps_segment_traits_2<Kernel> traits;

  std::vector<Polygon_with_holes_2> temp;
  CGAL::difference(are_holes_and_boundary_pairwise_disjoint(a, traits)
                       ? a
                       : fix_sfs_valid_polygon(a),
                   are_holes_and_boundary_pairwise_disjoint(b, traits)
                       ? b
                       : fix_sfs_valid_polygon(b),
                   std::back_inserter(temp));

  // polygon outer rings from difference can self intersect at points
  // therefore we need to split the generated polygons so that they are valid
  // for SFS
  for (const auto &poly : temp) {
    out = fix_cgal_valid_polygon(poly, out);
  }

  return out;
}

/// @} end of private section

} // namespace SFCGAL::algorithm
#endif
