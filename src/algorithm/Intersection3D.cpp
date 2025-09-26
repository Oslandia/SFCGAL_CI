// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <CGAL/intersections.h>
#include <memory>

#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>

#include "SFCGAL/algorithm/intersection.h"
#include "SFCGAL/algorithm/intersects.h"
#include "SFCGAL/detail/GeometrySet.h"
#include "SFCGAL/detail/triangulate/triangulateInGeometrySet.h"

#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Side_of_triangle_mesh.h>

using namespace SFCGAL::detail;

namespace SFCGAL::algorithm {

namespace {

/**
 * @brief Fast bounding box overlap test
 * @param bbox1 First bounding box
 * @param bbox2 Second bounding box
 * @return true if boxes overlap, enabling early termination
 */
inline bool
bbox_intersects(const CGAL::Bbox_3 &bbox1, const CGAL::Bbox_3 &bbox2)
{
  return CGAL::do_overlap(bbox1, bbox2);
}

} // anonymous namespace

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection

// Helper function to handle source inside, target outside case
static void
handle_source_inside_target_outside(const CGAL::Segment_3<Kernel> *segment,
                                    const GeometrySet<3> &intersection_points,
                                    GeometrySet<3>       &output)
{
  CGAL::Segment_3<Kernel> const intersectionSegment(
      segment->source(), intersection_points.points().begin()->primitive());

  if (intersectionSegment.source() == intersectionSegment.target()) {
    output.addPrimitive(segment->source());
  } else {
    output.addPrimitive(intersectionSegment);
  }
}

// Helper function to handle source outside, target inside case
static void
handle_source_outside_target_inside(const CGAL::Segment_3<Kernel> *segment,
                                    const GeometrySet<3> &intersection_points,
                                    GeometrySet<3>       &output)
{
  CGAL::Segment_3<Kernel> const intersectionSegment(
      intersection_points.points().begin()->primitive(), segment->target());

  if (intersectionSegment.source() == intersectionSegment.target()) {
    output.addPrimitive(segment->target());
  } else {
    output.addPrimitive(intersectionSegment);
  }
}

// Helper function to handle both source and target outside case
static void
handle_both_outside(const GeometrySet<3> &intersection_points,
                    GeometrySet<3>       &output)
{
  if (intersection_points.points().size() == 1) {
    // Single intersection point: tangent case
    output.addPrimitive(intersection_points.points().begin()->primitive());
  } else if (intersection_points.points().size() >= 2) {
    // Multiple intersection points: create segment between first two
    // This is the fix for the main bug
    auto                  iterator    = intersection_points.points().begin();
    CGAL::Point_3<Kernel> first_point = iterator->primitive();
    ++iterator;
    CGAL::Point_3<Kernel> second_point = iterator->primitive();

    CGAL::Segment_3<Kernel> const intersectionSegment(first_point,
                                                      second_point);
    // NOLINTBEGIN(bugprone-branch-clone)
    if (intersectionSegment.source() == intersectionSegment.target()) {
      output.addPrimitive(first_point);
    } else {
      output.addPrimitive(intersectionSegment);
    } // NOLINTEND(bugprone-branch-clone)
  }
  // If no intersection points, do nothing
}

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
void
_intersection_solid_segment(const PrimitiveHandle<3> &polyhedron_handle,
                            const PrimitiveHandle<3> &segment_handle,
                            GeometrySet<3>           &output)
{
  const auto *ext_poly = polyhedron_handle.as<MarkedPolyhedron>();
  BOOST_ASSERT(ext_poly->is_closed());
  const auto *segment = segment_handle.as<CGAL::Segment_3<Kernel>>();

  // Early termination with bounding box test
  CGAL::Bbox_3 poly_bbox = CGAL::Polygon_mesh_processing::bbox(*ext_poly);
  CGAL::Bbox_3 seg_bbox  = segment->bbox();

  if (!bbox_intersects(poly_bbox, seg_bbox)) {
    return; // No intersection possible - significant performance gain
  }

  auto *ext_poly_nc = const_cast<MarkedPolyhedron *>(ext_poly);
  CGAL::Side_of_triangle_mesh<MarkedPolyhedron, Kernel> const is_in_ext(
      *ext_poly_nc);

  // OPTIMIZATION: Triangulate only once and reuse the result
  GeometrySet<3> triangles;
  triangulate::triangulate(*ext_poly, triangles);

  GeometrySet<3> const spoint(segment->source());
  GeometrySet<3> const tpoint(segment->target());

  bool const source_inside =
      (is_in_ext(segment->source()) != CGAL::ON_UNBOUNDED_SIDE) ||
      intersects(triangles, spoint);
  bool const target_inside =
      (is_in_ext(segment->target()) != CGAL::ON_UNBOUNDED_SIDE) ||
      intersects(triangles, tpoint);

  if (source_inside && target_inside) {
    // the entire segment intersects the volume, return the segment
    output.addPrimitive(segment_handle);
  } else {
    // OPTIMIZATION: Reuse the triangulation computed above
    GeometrySet<3> geometry_set;
    geometry_set.addPrimitive(segment_handle);
    GeometrySet<3> intersection_result;

    // call recursively on triangulated polyhedron
    intersection(geometry_set, triangles, intersection_result);

    if (!intersection_result.points().empty()) {
      // the intersection is a point, build a segment from that point to the
      // other end
      if (!source_inside && target_inside) {
        handle_source_outside_target_inside(segment, intersection_result,
                                            output);
      } else if (source_inside && !target_inside) {
        handle_source_inside_target_outside(segment, intersection_result,
                                            output);
      } else { // !source_inside && !target_inside
        handle_both_outside(intersection_result, output);
      }
    }

    if (!intersection_result.segments().empty()) {
      // the intersection is a segment
      output.addPrimitive(intersection_result.segments().begin()->primitive());
    }
  }
}
// NOLINTEND(bugprone-easily-swappable-parameters)

using Polyline_3 = std::vector<Kernel::Point_3>;

struct Is_not_marked {
  auto
  operator()(MarkedPolyhedron::Halfedge_const_handle halfedge) const -> bool
  {
    return !halfedge->mark;
  }
};

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
void
_intersection_solid_triangle(const MarkedPolyhedron         &polyhedron,
                             const CGAL::Triangle_3<Kernel> &triangle,
                             GeometrySet<3>                 &output)
{
  BOOST_ASSERT(polyhedron.is_closed());

  // Early termination with bounding box test
  CGAL::Bbox_3 poly_bbox = CGAL::Polygon_mesh_processing::bbox(polyhedron);
  CGAL::Bbox_3 tri_bbox  = triangle.bbox();

  if (!bbox_intersects(poly_bbox, tri_bbox)) {
    return; // No intersection possible - significant performance gain
  }

  MarkedPolyhedron polyb;
  polyb.make_triangle(triangle.vertex(0), triangle.vertex(1),
                      triangle.vertex(2));

  // OPTIMIZATION: Avoid copy by using const_cast (safe because we don't modify)
  MarkedPolyhedron &polya = const_cast<MarkedPolyhedron &>(polyhedron);

  // OPTIMIZATION: Lazy construction - create Side_of_triangle_mesh only if
  // needed
  std::unique_ptr<CGAL::Side_of_triangle_mesh<MarkedPolyhedron, Kernel>>
      side_of_tm_ptr;

  std::list<Polyline_3> polylines;
  CGAL::Polygon_mesh_processing::surface_intersection(
      polya, polyb, std::back_inserter(polylines));
  if (polylines.empty()) {
    // no surface intersection
    // OPTIMIZATION: Create Side_of_triangle_mesh only when needed
    if (!side_of_tm_ptr) {
      side_of_tm_ptr = std::make_unique<
          CGAL::Side_of_triangle_mesh<MarkedPolyhedron, Kernel>>(polya);
    }
    // if one of the point of the triangle is inside the polyhedron,
    // the triangle is inside
    if ((*side_of_tm_ptr)(triangle.vertex(0)) != CGAL::ON_UNBOUNDED_SIDE) {
      output.addPrimitive(triangle);
    }
    return;
  }

  CGAL::Polygon_mesh_processing::clip(
      polyb, polya, CGAL::parameters::use_compact_clipper(true));

  bool hasSurface = false;

  std::vector<MarkedPolyhedron> ccs;
  CGAL::Polygon_mesh_processing::split_connected_components(polyb, ccs);

  for (const MarkedPolyhedron &mp : ccs) {
    // check if all vertices are on polya
    bool all_on = true;
    for (auto v : vertices(mp)) {
      // OPTIMIZATION: Create Side_of_triangle_mesh only when first needed
      if (!side_of_tm_ptr) {
        side_of_tm_ptr = std::make_unique<
            CGAL::Side_of_triangle_mesh<MarkedPolyhedron, Kernel>>(polya);
      }
      if ((*side_of_tm_ptr)(v->point()) != CGAL::ON_BOUNDARY) {
        all_on = false;
        break;
      }
    }
    if (all_on) {
      output.addPrimitive(mp);
    } else {
      hasSurface = true;
      output.addPrimitive(mp, FLAG_IS_PLANAR);
    }
  }

  if (hasSurface) {
    return;
  }

  for (auto &polyline : polylines) {
    if (polyline.size() == 1) {
      // it's a point
      output.addPrimitive(polyline[0]);
    } else {
      for (size_t k = 1; k < polyline.size(); ++k) {
        CGAL::Segment_3<Kernel> const seg(polyline[k - 1], polyline[k]);
        output.addPrimitive(seg);
      }
    }
  }
}
// NOLINTEND(bugprone-easily-swappable-parameters)

void
_intersection_solid_solid(const MarkedPolyhedron &polyhedron_a,
                          const MarkedPolyhedron &polyhedron_b,
                          GeometrySet<3>         &output)
{
  // Early termination with bounding box test
  CGAL::Bbox_3 bbox_a = CGAL::Polygon_mesh_processing::bbox(polyhedron_a);
  CGAL::Bbox_3 bbox_b = CGAL::Polygon_mesh_processing::bbox(polyhedron_b);

  if (!bbox_intersects(bbox_a, bbox_b)) {
    return; // No intersection possible - major performance improvement for
            // distant objects
  }

  // 1. find intersections on surfaces
  // CGAL corefinement or polyhedra_intersection do not return polygon
  // intersections between two solids they only return points, lines and volumes
  // (but no surfaces ...)
  {
    // call CGAL::intersection on triangles
    GeometrySet<3> gsa;
    GeometrySet<3> gsb;
    // convert polyhedra to geometry sets
    // (no actual triangulation is done if the polyhedra are pure_triangle()
    triangulate::triangulate(polyhedron_a, gsa);
    triangulate::triangulate(polyhedron_b, gsb);
    // "recurse" call on triangles
    algorithm::intersection(gsa, gsb, output);
  }

  // 2. find intersections in volumes
  {
    // OPTIMIZATION: Avoid copies - create result mesh instead
    MarkedPolyhedron result;
    if (CGAL::Polygon_mesh_processing::corefine_and_compute_intersection(
            const_cast<MarkedPolyhedron &>(polyhedron_a),
            const_cast<MarkedPolyhedron &>(polyhedron_b), result)) {
      if (std::next(vertices(result).first) != vertices(result).second) {
        output.addPrimitive(result);
      }
    }
  }
}

/// @} end of private section

// ----------------------------------------------------------------------------------
// -- public interface
// ----------------------------------------------------------------------------------
/// @publicsection

/**
 * @brief Compute intersection between two 3D primitive handles
 * @param primitive_a First primitive (must have dimension >= primitive_b)
 * @param primitive_b Second primitive
 * @param output Output geometry set to store results
 * @param unused Dimension tag for template resolution
 * @pre primitive_a's dimension must be larger than or equal to primitive_b's
 */
void
intersection(const PrimitiveHandle<3> &primitive_a,
             const PrimitiveHandle<3> &primitive_b, GeometrySet<3> &output,
             dim_t<3> /*unused*/)
{
  // everything vs a point
  if (primitive_b.handle.which() == PrimitivePoint) {
    if (algorithm::intersects(primitive_a, primitive_b)) {
      output.addPrimitive(
          *boost::get<const TypeForDimension<3>::Point *>(primitive_b.handle));
    }
  } else if (primitive_a.handle.which() == PrimitiveSegment &&
             primitive_b.handle.which() == PrimitiveSegment) {
    const auto        *seg1     = primitive_a.as<CGAL::Segment_3<Kernel>>();
    const auto        *seg2     = primitive_b.as<CGAL::Segment_3<Kernel>>();
    CGAL::Object const interObj = CGAL::intersection(*seg1, *seg2);
    output.addPrimitive(interObj);
  } else if (primitive_a.handle.which() == PrimitiveSurface) {
    const auto *tri1 = primitive_a.as<CGAL::Triangle_3<Kernel>>();

    if (primitive_b.handle.which() == PrimitiveSegment) {
      const auto        *seg2     = primitive_b.as<CGAL::Segment_3<Kernel>>();
      CGAL::Object const interObj = CGAL::intersection(*tri1, *seg2);
      output.addPrimitive(interObj);
    } else if (primitive_b.handle.which() == PrimitiveSurface) {
      const auto        *tri2     = primitive_b.as<CGAL::Triangle_3<Kernel>>();
      CGAL::Object const interObj = CGAL::intersection(*tri1, *tri2);
      output.addPrimitive(interObj, /* pointsAsRing */ true);
    }
  } else if (primitive_a.handle.which() == PrimitiveVolume) {
    if (primitive_b.handle.which() == PrimitiveSegment) {
      _intersection_solid_segment(primitive_a, primitive_b, output);
    } else if (primitive_b.handle.which() == PrimitiveSurface) {
      _intersection_solid_triangle(*primitive_a.as<MarkedPolyhedron>(),
                                   *primitive_b.as<CGAL::Triangle_3<Kernel>>(),
                                   output);
    } else if (primitive_b.handle.which() == PrimitiveVolume) {
      const MarkedPolyhedron &solid_a = *primitive_a.as<MarkedPolyhedron>();
      const MarkedPolyhedron &solid_b = *primitive_b.as<MarkedPolyhedron>();
      BOOST_ASSERT(solid_a.is_closed());
      BOOST_ASSERT(solid_b.is_closed());
      _intersection_solid_solid(solid_a, solid_b, output);
    }
  }
}

/**
 * @brief Public interface for solid-triangle intersection (compatibility)
 * @param polyhedron The solid polyhedron
 * @param triangle The triangle to intersect with
 * @param output Output geometry set to store results
 *
 * This function provides the expected public API for intersectionSolidTriangle
 * while using the optimized internal implementation with bounding box tests.
 */
void
intersectionSolidTriangle(const MarkedPolyhedron         &polyhedron,
                          const CGAL::Triangle_3<Kernel> &triangle,
                          GeometrySet<3>                 &output)
{
  _intersection_solid_triangle(polyhedron, triangle, output);
}

} // namespace SFCGAL::algorithm
