// Copyright (c) 2024-2024, SFCGAL Contributors and Oslandia
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/alphaShapes3D.h"
#include "SFCGAL/detail/GetPointsVisitor.h"
#include "SFCGAL/detail/transform/ForceOrderPoints.h"

#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

namespace SFCGAL::algorithm {

using Vb3  = CGAL::Alpha_shape_vertex_base_3<Kernel>;
using Fb3  = CGAL::Alpha_shape_cell_base_3<Kernel>;
using Tds3 = CGAL::Triangulation_data_structure_3<Vb3, Fb3>;
using Delaunay3 =
    CGAL::Delaunay_triangulation_3<Kernel, Tds3, CGAL::Fast_location>;
using Alpha_shape_3 = CGAL::Alpha_shape_3<Delaunay3>;
using Point_3       = Kernel::Point_3;

//! Reads the geometry and compute alphaShape
auto
buildAlphaShape(const Geometry &geom, AlphaShape3DMode mode)
    -> std::unique_ptr<Alpha_shape_3>
{
  if (geom.isEmpty()) {
    return nullptr;
  }

  // Collect points from geometry
  SFCGAL::detail::GetPointsVisitor getPointVisitor;
  const_cast<Geometry &>(geom).accept(getPointVisitor);

  // Need at least 4 points for 3D alpha shapes
  if (getPointVisitor.points.size() < 4) {
    return nullptr;
  }

  // Create points vector
  std::vector<Point_3> points;
  points.reserve(getPointVisitor.points.size());
  for (const auto &point : getPointVisitor.points) {
    points.push_back(point->toPoint_3());
  }

  // Create and configure alpha shape
  auto alphaShape =
      std::make_unique<Alpha_shape_3>(points.begin(), points.end());
  alphaShape->set_mode(mode == AlphaShape3DMode::REGULARIZED
                           ? Alpha_shape_3::REGULARIZED
                           : Alpha_shape_3::GENERAL);

  // Set optimal alpha value for one solid component
  auto optimalAlphaValue = alphaShape->find_optimal_alpha(1);
  alphaShape->set_alpha(*optimalAlphaValue);

  return alphaShape;
}

//! Converts alphaShape to a PolyhedralSurface
auto
alphaShape3D_to_polyhedralSurface(const Alpha_shape_3 &alphaShape)
    -> std::unique_ptr<PolyhedralSurface>
{
  auto resultSurface = std::make_unique<PolyhedralSurface>();

  // Iterate through all facets of the alpha shape
  for (auto facetIterator = alphaShape.alpha_shape_facets_begin();
       facetIterator != alphaShape.alpha_shape_facets_end(); ++facetIterator) {

    if (alphaShape.classify(*facetIterator) == Alpha_shape_3::REGULAR) {
      // Get facet vertices
      auto cell = facetIterator->first;
      int  idx  = facetIterator->second;

      const Point_3 pt1 = cell->vertex((idx + 1) & 3)->point();
      const Point_3 pt2 = cell->vertex((idx + 2) & 3)->point();
      const Point_3 pt3 = cell->vertex((idx + 3) & 3)->point();

      // Create triangle polygon
      auto poly = std::make_unique<Polygon>();
      auto ring = std::make_unique<LineString>();

      // create a closed ring
      ring->addPoint(Point(pt1));
      ring->addPoint(Point(pt2));
      ring->addPoint(Point(pt3));
      ring->addPoint(Point(pt1));

      poly->setExteriorRing(ring.release());

      // SFCGAL::transform::ForceOrderPoints force(false);
      // poly->accept(force);

      resultSurface->addPolygon(poly.release());
    }
  }

  // SFCGAL::transform::ForceOrderPoints force(true);
  // resultSurface->accept(force);

  return resultSurface;
}

auto
alphaShapes3D(const Geometry &geom, AlphaShape3DMode mode)
    -> std::unique_ptr<PolyhedralSurface>
{
  std::unique_ptr<Alpha_shape_3> alphaShape = buildAlphaShape(geom, mode);
  if (!alphaShape) {
    return std::make_unique<PolyhedralSurface>();
  }
  return alphaShape3D_to_polyhedralSurface(*alphaShape);
}

} // namespace SFCGAL::algorithm
