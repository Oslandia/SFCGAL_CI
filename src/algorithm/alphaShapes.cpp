/**
 * Part of the My Project, under the MIT license.
 * See LICENSE.txt for license information.
 *
 * SPDX-License-Identifier: MIT
 */

/**
 * \file        : alphaShapes
 * \author      : Loïc Bartoletti (loic dot bartoletti @ oslandia dot com)
 * \created     : novembre 2021
 */

#include <SFCGAL/algorithm/alphaShapes.h>

#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/MultiPolygon.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Triangle.h>
#include <SFCGAL/TriangulatedSurface.h>

#include <SFCGAL/detail/GetPointsVisitor.h>

#include <SFCGAL/Exception.h>
#include <SFCGAL/Kernel.h>
#include <boost/format.hpp>

#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/algorithm.h>

#include <CGAL/Arr_non_caching_segment_basic_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <vector>

namespace SFCGAL {
namespace algorithm {

typedef CGAL::Alpha_shape_vertex_base_2<Kernel> Vb;
typedef CGAL::Alpha_shape_face_base_2<Kernel> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<Kernel, Tds> Triangulation_2;
typedef CGAL::Point_2<Kernel> Point_2;
typedef CGAL::Segment_2<Kernel> Segment_2;
typedef CGAL::Alpha_shape_2<Triangulation_2> Alpha_shape_2;

typedef Alpha_shape_2::Alpha_shape_edges_iterator Alpha_shape_edges_iterator;

typedef CGAL::Arr_non_caching_segment_basic_traits_2<Kernel> Traits_2;
// typedef CGAL::Arr_segment_traits_2<Kernel> Traits_2;
typedef CGAL::Arrangement_2<Traits_2> Arrangement;

template <class OutputIterator>
void alpha_edges(const Alpha_shape_2 &A, OutputIterator out) {
  Alpha_shape_edges_iterator it = A.alpha_shape_edges_begin(),
                             end = A.alpha_shape_edges_end();
  for (; it != end; ++it) {
    if (A.classify(*it) == 2)
      *out++ = A.segment(*it);
  }
}

static double computeAlpha(const Geometry &g, Alpha_shape_2 &alphaShape,
                           double alpha = 0, size_t nb_components = 1) {
  using CGAL::object_cast;
  double result = -1.0;

  if (g.isEmpty()) {
    return result;
  }

  SFCGAL::detail::GetPointsVisitor getPointVisitor;
  const_cast<Geometry &>(g).accept(getPointVisitor);

  // collect points

  if (getPointVisitor.points.size() == 0) {
    return result;
  }

  std::vector<Point_2> points;

  for (size_t i = 0; i < getPointVisitor.points.size(); i++) {
    points.push_back(getPointVisitor.points[i]->toPoint_2());
  }

  std::vector<Segment_2> segments;
  alphaShape.make_alpha_shape(points.begin(), points.end());
  alphaShape.set_alpha(Kernel::FT(alpha));
  alpha_edges(alphaShape, std::back_inserter(segments));

  if (segments.size() == 0)
    return result;

  result = CGAL::to_double(*alphaShape.find_optimal_alpha(nb_components));

  return result;
}

static std::unique_ptr<Geometry> alpha_to_geometry(const Alpha_shape_2 &A,
                                                   bool allow_holes) {
  std::vector<Segment_2> segments;
  alpha_edges(A, std::back_inserter(segments));

  Arrangement arr;

  CGAL::insert_non_intersecting_curves(arr, segments.begin(), segments.end());
  Polygon *poly = new Polygon;
  for (auto f = arr.faces_begin(); f != arr.faces_end(); f++) {
    LineString *ring = new LineString;
    for (auto h = f->holes_begin(); h != f->holes_end(); h++) {
      auto he = *h;
      do {
        ring->addPoint(he->source()->point());
      } while (++he != *h);
    }

    if (ring->numPoints() > 3) {
      ring->addPoint(ring->startPoint());
      if (f->is_unbounded()) {
        poly->setExteriorRing(ring);
      } else if (allow_holes) {
        poly->addInteriorRing(ring);
      }
    }
  }

  return std::unique_ptr<Geometry>(poly);
}

std::unique_ptr<Geometry> optimal_alpha_shapes(const Geometry &g,
                                               bool allow_holes,
                                               size_t nb_components) {
  Alpha_shape_2 A;
  const double optimalAlpha = computeAlpha(g, A, 10000, nb_components);
  if (optimalAlpha < 0) {
    return std::unique_ptr<Geometry>(new GeometryCollection());
  }

  A.set_alpha(optimalAlpha);

  return alpha_to_geometry(A, allow_holes);
}

std::unique_ptr<Geometry> alphaShapes(const Geometry &g, double alpha,
                                      bool allow_holes) {
  using CGAL::object_cast;
  Alpha_shape_2 A;
  const double optimalAlpha = computeAlpha(g, A, alpha);
  if (optimalAlpha < 0) {
    return std::unique_ptr<Geometry>(new GeometryCollection());
  }

  return alpha_to_geometry(A, allow_holes);
}

} // namespace algorithm
} // namespace SFCGAL
