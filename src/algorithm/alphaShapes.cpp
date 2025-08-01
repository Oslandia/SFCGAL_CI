// Copyright (c) 2022-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/alphaShapes.h"

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"

#include "SFCGAL/detail/GetPointsVisitor.h"

#include "SFCGAL/Exception.h"
#include "SFCGAL/Kernel.h"
#include <boost/format.hpp>

#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/algorithm.h>

#include <CGAL/Arr_non_caching_segment_basic_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <vector>

namespace SFCGAL::algorithm {

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection

using Vb              = CGAL::Alpha_shape_vertex_base_2<Kernel>;
using Fb              = CGAL::Alpha_shape_face_base_2<Kernel>;
using Tds             = CGAL::Triangulation_data_structure_2<Vb, Fb>;
using Triangulation_2 = CGAL::Delaunay_triangulation_2<Kernel, Tds>;
using Alpha_shape_2   = CGAL::Alpha_shape_2<Triangulation_2>;

using Alpha_shape_edges_iterator = Alpha_shape_2::Alpha_shape_edges_iterator;

using Traits_2    = CGAL::Arr_non_caching_segment_basic_traits_2<Kernel>;
using Arrangement = CGAL::Arrangement_2<Traits_2>;

template <class OutputIterator>
void
alpha_edges(const Alpha_shape_2 &A, OutputIterator out)
{
  auto it  = A.alpha_shape_edges_begin();
  auto end = A.alpha_shape_edges_end();
  for (; it != end; ++it) {
    if (A.classify(*it) == 2) {
      *out++ = A.segment(*it);
    }
  }
}

static auto
computeAlpha(const Geometry &g, Alpha_shape_2 &alphaShape, double alpha = 0,
             size_t nb_components = 1) -> double
{
  using CGAL::object_cast;
  double result = -1.0;

  if (g.isEmpty()) {
    return result;
  }

  SFCGAL::detail::GetPointsVisitor getPointVisitor;
  const_cast<Geometry &>(g).accept(getPointVisitor);

  // collect points
  if (getPointVisitor.points.size() < 4) {
    return result;
  }

  std::vector<Point_2> points;

  points.reserve(getPointVisitor.points.size());
  for (auto &point : getPointVisitor.points) {
    points.push_back(point->toPoint_2());
  }

  std::vector<Segment_2> segments;
  alphaShape.make_alpha_shape(points.begin(), points.end());
  alphaShape.set_alpha(Kernel::FT(alpha));
  alpha_edges(alphaShape, std::back_inserter(segments));

  // Ensure we compare the iterator from find_optimal_alpha(nb_components)
  // against alpha_end() before dereferencing to avoid potential crash when
  // no valid alpha is found.
  auto it_alpha = alphaShape.find_optimal_alpha(nb_components);
  if (it_alpha != alphaShape.alpha_end()) {
    result = CGAL::to_double(*it_alpha);
  } else {
    std::ostringstream oss;
    oss << "Can't find alpha value.";
    throw std::runtime_error(oss.str());
  }

  return result;
}

static auto
alpha_to_geometry(const Alpha_shape_2 &A, bool allow_holes)
    -> std::unique_ptr<Geometry>
{
  std::vector<Segment_2> segments;
  alpha_edges(A, std::back_inserter(segments));

  Arrangement arr;

  CGAL::insert_non_intersecting_curves(arr, segments.begin(), segments.end());
  auto poly{std::make_unique<Polygon>()};
  for (auto f = arr.faces_begin(); f != arr.faces_end(); f++) {
    auto ring{std::make_unique<LineString>()};
    for (auto h = f->holes_begin(); h != f->holes_end(); h++) {
      auto he = *h;
      do {
        ring->addPoint(he->source()->point());
      } while (++he != *h);
    }

    if (ring->numPoints() > 3) {
      ring->addPoint(ring->startPoint());
      if (f->is_unbounded()) {
        poly->setExteriorRing(ring.release());
      } else if (allow_holes) {
        poly->addInteriorRing(ring.release());
      }
    }
  }

  std::unique_ptr<Geometry> result = std::move(poly);

  return result;
}

/// @} end of private section

// ----------------------------------------------------------------------------------
// -- public interface
// ----------------------------------------------------------------------------------
/// @publicsection

auto
optimal_alpha_shapes(const Geometry &g, bool allow_holes, size_t nb_components)
    -> std::unique_ptr<Geometry>
{
  Alpha_shape_2 A;
  const double  optimalAlpha{computeAlpha(g, A, 10000, nb_components)};
  if (optimalAlpha < 0) {
    return std::unique_ptr<Geometry>(new GeometryCollection());
  }

  A.set_alpha(optimalAlpha);

  return alpha_to_geometry(A, allow_holes);
}

auto
alphaShapes(const Geometry &g, double alpha, bool allow_holes)
    -> std::unique_ptr<Geometry>
{
  using CGAL::object_cast;
  Alpha_shape_2 A;
  const double  optimalAlpha{computeAlpha(g, A, alpha)};
  if (optimalAlpha < 0) {
    return std::unique_ptr<Geometry>(new GeometryCollection());
  }

  return alpha_to_geometry(A, allow_holes);
}

} // namespace SFCGAL::algorithm
