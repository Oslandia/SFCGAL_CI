// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#include <SFCGAL/Polygon.h>
#include <SFCGAL/TriangulatedSurface.h>
#include <SFCGAL/detail/triangulate/ConstraintDelaunayTriangulation.h>
#include <SFCGAL/detail/triangulate/triangulateInGeometrySet.h>
#include <SFCGAL/triangulate/triangulatePolygon.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

using namespace SFCGAL::detail;

namespace SFCGAL {
namespace triangulate {

/**
 * vertex information with original coordinates
 */
struct VertexInfo2 {
  VertexInfo2() : original(NULL) {}

  const CGAL::Point_3<Kernel> *original;
};

typedef CGAL::Triangulation_vertex_base_with_info_2<VertexInfo2, Kernel>
    triangulation_vertex_base;
typedef CGAL::Triangulation_data_structure_2<triangulation_vertex_base>
    triangulation_data_structure;

typedef CGAL::Delaunay_triangulation_2<Kernel, triangulation_data_structure>
    Triangulation;

///
/// input polyhedron must have its planes computed
///
void
triangulate(const MarkedPolyhedron &polyhedron, GeometrySet<3> &geometry)
{
  if (polyhedron.is_pure_triangle()) {
    for (MarkedPolyhedron::Facet_const_iterator fit = polyhedron.facets_begin();
         fit != polyhedron.facets_end(); ++fit) {
      MarkedPolyhedron::Facet::Halfedge_around_facet_const_circulator pit;

      pit                              = fit->facet_begin();
      const CGAL::Point_3<Kernel> &pt1 = pit->vertex()->point();
      ++pit;
      const CGAL::Point_3<Kernel> &pt2 = pit->vertex()->point();
      ++pit;
      const CGAL::Point_3<Kernel> &pt3 = pit->vertex()->point();

      CGAL::Triangle_3<Kernel> tri(pt1, pt2, pt3);
      geometry.addPrimitive(tri);
    }

    return;
  }

  Triangulation triangulation;

  for (MarkedPolyhedron::Facet_const_iterator fit = polyhedron.facets_begin();
       fit != polyhedron.facets_end(); ++fit) {
    MarkedPolyhedron::Facet::Halfedge_around_facet_const_circulator pit;

    triangulation.clear();

    CGAL::Plane_3<Kernel> plane = fit->plane();
    BOOST_ASSERT(!plane.is_degenerate());

    pit = fit->facet_begin();

    do {
      const CGAL::Point_3<Kernel> &pt3 = pit->vertex()->point();
      CGAL::Point_2<Kernel>        pt2 = plane.to_2d(pt3);

      Triangulation::Vertex_handle vh = triangulation.insert(pt2);
      vh->info().original             = &pt3;

      pit++;
    } while (pit != fit->facet_begin());

    for (Triangulation::Finite_faces_iterator tit =
             triangulation.finite_faces_begin();
         tit != triangulation.finite_faces_end(); ++tit) {
      const CGAL::Point_3<Kernel> *a = tit->vertex(0)->info().original;
      const CGAL::Point_3<Kernel> *b = tit->vertex(1)->info().original;
      const CGAL::Point_3<Kernel> *c = tit->vertex(2)->info().original;

      CGAL::Triangle_3<Kernel> tri(*a, *b, *c);
      geometry.addPrimitive(tri);
    }
  }
}

///
///
///
void
triangulate(const CGAL::Polygon_with_holes_2<Kernel> &polygon,
            GeometrySet<2>                           &output)
{
  // FIXME: lots of copies here, can we avoid this ?
  Polygon             poly(polygon);
  TriangulatedSurface surf;
  triangulatePolygon3D(poly, surf);

  for (size_t i = 0; i < surf.numTriangles(); ++i) {
    output.addGeometry(surf.triangleN(i));
  }
}

} // namespace triangulate
} // namespace SFCGAL
