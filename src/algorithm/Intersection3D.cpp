// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <CGAL/intersections.h>

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

void
_intersection_solid_segment(const PrimitiveHandle<3> &pa,
                            const PrimitiveHandle<3> &pb,
                            GeometrySet<3>           &output)
{
  // typedef CGAL::Polyhedral_mesh_domain_3<MarkedPolyhedron, Kernel>
  // Mesh_domain;

  const auto *ext_poly = pa.as<MarkedPolyhedron>();
  BOOST_ASSERT(ext_poly->is_closed());
  const auto *segment = pb.as<CGAL::Segment_3<Kernel>>();

  auto *ext_poly_nc = const_cast<MarkedPolyhedron *>(ext_poly);
  CGAL::Side_of_triangle_mesh<MarkedPolyhedron, Kernel> const is_in_ext(
      *ext_poly_nc);

  GeometrySet<3>       triangles;
  GeometrySet<3> const spoint(segment->source());
  GeometrySet<3> const tpoint(segment->target());
  triangulate::triangulate(*ext_poly, triangles);

  bool const source_inside =
      (is_in_ext(segment->source()) != CGAL::ON_UNBOUNDED_SIDE) ||
      intersects(triangles, spoint);
  bool const target_inside =
      (is_in_ext(segment->target()) != CGAL::ON_UNBOUNDED_SIDE) ||
      intersects(triangles, tpoint);

  if (source_inside && target_inside) {
    // the entire segment intersects the volume, return the segment
    output.addPrimitive(pb);
  } else {
    GeometrySet<3> triangles;
    GeometrySet<3> g;
    triangulate::triangulate(*ext_poly, triangles);
    g.addPrimitive(pb);
    GeometrySet<3> inter;

    // call recursively on triangulated polyhedron
    intersection(g, triangles, inter);

    if (!inter.points().empty()) {
      // the intersection is a point, build a segment from that point to the
      // other end
      if (!source_inside && target_inside) {
        CGAL::Segment_3<Kernel> const interSeg(
            inter.points().begin()->primitive(), segment->target());

        if (interSeg.source() == interSeg.target()) {
          output.addPrimitive(segment->target());
        } else {
          output.addPrimitive(interSeg);
        }
      } else if (source_inside && !target_inside) {
        CGAL::Segment_3<Kernel> const interSeg(
            segment->source(), inter.points().begin()->primitive());

        if (interSeg.source() == interSeg.target()) {
          output.addPrimitive(segment->source());
        } else {
          output.addPrimitive(interSeg);
        }
      } else { // !source_inside && !target_inside => intersection on a point
        output.addPrimitive(inter.points().begin()->primitive());
      }
    }

    if (!inter.segments().empty()) {
      // the intersection is a segment
      output.addPrimitive(inter.segments().begin()->primitive());
    }
  }
}

using Polyline_3 = std::vector<Kernel::Point_3>;

struct Is_not_marked {
  auto
  operator()(MarkedPolyhedron::Halfedge_const_handle h) const -> bool
  {
    return !h->mark;
  }
};

void
_intersection_solid_triangle(const MarkedPolyhedron         &pa,
                             const CGAL::Triangle_3<Kernel> &tri,
                             GeometrySet<3>                 &output)
{
  BOOST_ASSERT(pa.is_closed());

  MarkedPolyhedron polyb;
  polyb.make_triangle(tri.vertex(0), tri.vertex(1), tri.vertex(2));

  MarkedPolyhedron                                            polya = pa;
  CGAL::Side_of_triangle_mesh<MarkedPolyhedron, Kernel> const side_of_tm(polya);

  std::list<Polyline_3> polylines;
  CGAL::Polygon_mesh_processing::surface_intersection(
      polya, polyb, std::back_inserter(polylines));
  if (polylines.empty()) {
    // no surface intersection
    // if one of the point of the triangle is inside the polyhedron,
    // the triangle is inside
    if (side_of_tm(tri.vertex(0)) != CGAL::ON_UNBOUNDED_SIDE) {
      output.addPrimitive(tri);
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
      if (side_of_tm(v->point()) != CGAL::ON_BOUNDARY) {
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

void
_intersection_solid_solid(const MarkedPolyhedron &pa,
                          const MarkedPolyhedron &pb, GeometrySet<3> &output)
{
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
    triangulate::triangulate(pa, gsa);
    triangulate::triangulate(pb, gsb);
    // "recurse" call on triangles
    algorithm::intersection(gsa, gsb, output);
  }

  // 2. find intersections in volumes
  {
    MarkedPolyhedron polya = pa;
    MarkedPolyhedron polyb = pb;
    if (CGAL::Polygon_mesh_processing::corefine_and_compute_intersection(
            polya, polyb, polya)) {
      if (std::next(vertices(polya).first) != vertices(polya).second) {
        output.addPrimitive(polya);
      }
    }
  }
}

//
// must be called with pa's dimension larger than pb's
void
intersection(const PrimitiveHandle<3> &pa, const PrimitiveHandle<3> &pb,
             GeometrySet<3> &output, dim_t<3> /*unused*/)
{
  // everything vs a point
  if (pb.handle.which() == PrimitivePoint) {
    if (algorithm::intersects(pa, pb)) {
      output.addPrimitive(
          *boost::get<const TypeForDimension<3>::Point *>(pb.handle));
    }
  } else if (pa.handle.which() == PrimitiveSegment &&
             pb.handle.which() == PrimitiveSegment) {
    const auto        *seg1     = pa.as<CGAL::Segment_3<Kernel>>();
    const auto        *seg2     = pb.as<CGAL::Segment_3<Kernel>>();
    CGAL::Object const interObj = CGAL::intersection(*seg1, *seg2);
    output.addPrimitive(interObj);
  } else if (pa.handle.which() == PrimitiveSurface) {
    const auto *tri1 = pa.as<CGAL::Triangle_3<Kernel>>();

    if (pb.handle.which() == PrimitiveSegment) {
      const auto        *seg2     = pb.as<CGAL::Segment_3<Kernel>>();
      CGAL::Object const interObj = CGAL::intersection(*tri1, *seg2);
      output.addPrimitive(interObj);
    } else if (pb.handle.which() == PrimitiveSurface) {
      const auto        *tri2     = pb.as<CGAL::Triangle_3<Kernel>>();
      CGAL::Object const interObj = CGAL::intersection(*tri1, *tri2);
      output.addPrimitive(interObj, /* pointsAsRing */ true);
    }
  } else if (pa.handle.which() == PrimitiveVolume) {
    if (pb.handle.which() == PrimitiveSegment) {
      _intersection_solid_segment(pa, pb, output);
    } else if (pb.handle.which() == PrimitiveSurface) {
      _intersection_solid_triangle(*pa.as<MarkedPolyhedron>(),
                                   *pb.as<CGAL::Triangle_3<Kernel>>(), output);
    } else if (pb.handle.which() == PrimitiveVolume) {
      const MarkedPolyhedron &sa = *pa.as<MarkedPolyhedron>();
      const MarkedPolyhedron &sb = *pb.as<MarkedPolyhedron>();
      BOOST_ASSERT(sa.is_closed());
      BOOST_ASSERT(sb.is_closed());
      _intersection_solid_solid(sa, sb, output);
    }
  }
}

} // namespace SFCGAL::algorithm
