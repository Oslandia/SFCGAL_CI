// Copyright (c) 2023, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/triangulate/triangulatePolygonalDomain.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#include <algorithm>
#include <iostream>
#include <string>

namespace SFCGAL {
namespace triangulate {

struct FaceInfo2 {
  FaceInfo2() {}
  int nesting_level;
  bool
  in_domain()
  {
    return nesting_level % 2 == 1;
  }
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel      K;
typedef CGAL::Triangulation_vertex_base_2<K>                     Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, K>  Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<K, Fbb>      Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>             TDS;
typedef CGAL::Exact_predicates_tag                               Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> CDT;
typedef CDT::Point                                               Point;
typedef CDT::Segment                                             Segment;
typedef CGAL::Polygon_2<K>                                       Polygon_2;
typedef CGAL::Point_2<K>                                         Point_2;

void
mark_domains(CDT &ct, CDT::Face_handle start, int index,
             std::list<CDT::Edge> &border, const std::vector<Segment> &rings)
{
  if (start->info().nesting_level != -1) {
    return;
  }
  std::list<CDT::Face_handle> queue;
  queue.push_back(start);
  while (!queue.empty()) {
    CDT::Face_handle fh = queue.front();
    queue.pop_front();
    if (fh->info().nesting_level == -1) {
      fh->info().nesting_level = index;
      for (int i = 0; i < 3; i++) {
        CDT::Edge        e(fh, i);
        CDT::Face_handle n = fh->neighbor(i);
        if (n->info().nesting_level == -1) {
          if (ct.is_constrained(e) && std::find(rings.begin(), rings.end(),
                                                ct.segment(e)) != rings.end()) {
            border.push_back(e);
          } else {
            queue.push_back(n);
          }
        }
      }
    }
  }
}

// explore set of facets connected with non constrained edges,
// and attribute to each such set a nesting level.
// We start from facets incident to the infinite vertex, with a nesting
// level of 0. Then we recursively consider the non-explored facets incident
// to constrained edges bounding the former set and increase the nesting level
// by 1. Facets in the domain are those with an odd nesting level.
void
mark_domains(CDT &cdt, const std::vector<Segment> &rings)
{
  for (CDT::All_faces_iterator it = cdt.all_faces_begin();
       it != cdt.all_faces_end(); ++it) {
    it->info().nesting_level = -1;
  }
  std::list<CDT::Edge> border;
  mark_domains(cdt, cdt.infinite_face(), 0, border, rings);
  while (!border.empty()) {
    CDT::Edge e = border.front();
    border.pop_front();
    CDT::Face_handle n = e.first->neighbor(e.second);
    if (n->info().nesting_level == -1) {
      mark_domains(cdt, n, e.first->info().nesting_level + 1, border, rings);
    }
  }
}

static Polygon_2
polygon_from_ring(const LineString &ring)
{
  Polygon_2 poly;
  for (size_t i = 0; i < ring.numPoints(); i++)
    poly.push_back(Point_2(CGAL::to_double(ring.pointN(i).x()),
                           CGAL::to_double(ring.pointN(i).y())));
  return std::move(poly);
}

auto
constrainedPolygonalDomain(const Polygon         &poly,
                           const MultiLineString &multiline,
                           const MultiPoint &multipoint) -> TriangulatedSurface
{
  TriangulatedSurface result;
  CDT                 cdt;

  for (auto &r : poly) {
    Polygon_2 poly2(polygon_from_ring(r));
    cdt.insert_constraint(poly2.vertices_begin(), poly2.vertices_end(), true);
  }

  // insert line constrains
  for (auto &l : multiline) {
    const LineString line = l.as<LineString>();
    for (size_t i = 1; i < line.numPoints(); ++i) {
      const Point_2 pa(CGAL::to_double(line.pointN(i - 1).x()),
                       CGAL::to_double(line.pointN(i - 1).y()));
      const Point_2 pb(CGAL::to_double(line.pointN(i).x()),
                       CGAL::to_double(line.pointN(i).y()));
      cdt.insert_constraint(pa, pb);
    }
  }

  // insert points
  for (size_t i = 0; i < multipoint.numGeometries(); i++) {
    const Point_2 point(CGAL::to_double(multipoint.pointN(i).x()),
                        CGAL::to_double(multipoint.pointN(i).y()));
    cdt.push_back(point);
  }

  std::vector<Segment> border;
  for (CDT::Constrained_edges_iterator e = cdt.constrained_edges_begin();
       e != cdt.constrained_edges_end(); e++) {
    const Segment s(cdt.segment(*e));
    border.push_back(s);
    border.push_back(s.opposite());
  }
  mark_domains(cdt, border);

  int i = 0;
  for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
       fit != cdt.finite_faces_end(); ++fit) {
    std::cout << "faces: " << i;
    if (fit->info().in_domain()) {
      std::cout << "inserted.";
      const SFCGAL::Point t0(fit->vertex(0 % 3)->point().x(),
                             fit->vertex(0 % 3)->point().y());
      const SFCGAL::Point t1(fit->vertex(1 % 3)->point().x(),
                             fit->vertex(1 % 3)->point().y());
      const SFCGAL::Point t2(fit->vertex(2 % 3)->point().x(),
                             fit->vertex(2 % 3)->point().y());
      result.addTriangle(Triangle(t0, t1, t2));
    }
    std::cout << "\n";
  }
  return result;
}

} // namespace triangulate
} // namespace SFCGAL
