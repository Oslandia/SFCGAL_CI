// Copyright (c) 2024-2024, SFCGAL Contributors and Oslandia
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/alphaShapes3D.h"
#include "SFCGAL/detail/GetPointsVisitor.h"

#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>

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

  // Iterate through all facets of the alpha shape
  std::vector<Point_3>                    points;
  std::vector<std::array<std::size_t, 3>> polygons;
  for (auto facetIterator = alphaShape.alpha_shape_facets_begin();
       facetIterator != alphaShape.alpha_shape_facets_end(); ++facetIterator) {

    if (alphaShape.classify(*facetIterator) == Alpha_shape_3::REGULAR) {
      // Get facet vertices
      auto cell     = facetIterator->first;
      int  facetIdx = facetIterator->second;

      const size_t currentIdx = points.size();
      polygons.push_back({currentIdx, currentIdx + 1, currentIdx + 2});
      points.push_back(cell->vertex((facetIdx + 1) & 3)->point());
      points.push_back(cell->vertex((facetIdx + 2) & 3)->point());
      points.push_back(cell->vertex((facetIdx + 3) & 3)->point());

      // // Create triangle polygon
      // auto poly = std::make_unique<Polygon>();
      // auto ring = std::make_unique<LineString>();

      // // create a closed ring
      // ring->addPoint(Point(pt1));
      // ring->addPoint(Point(pt2));
      // ring->addPoint(Point(pt3));
      // ring->addPoint(Point(pt1));

      // poly->setExteriorRing(ring.release());

      // resultSurface->addPolygon(poly.release());
    }
  }

  CGAL::Polygon_mesh_processing::repair_polygon_soup(points, polygons);
  CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);
  CGAL::Surface_mesh<Point_3> surfaceMesh;
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons,
                                                              surfaceMesh);

  return std::make_unique<PolyhedralSurface>(PolyhedralSurface(surfaceMesh));
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
