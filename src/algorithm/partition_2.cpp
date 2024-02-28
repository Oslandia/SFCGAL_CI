// Copyright (c) 2022-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/partition_2.h"

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/Polygon.h"

#include "SFCGAL/detail/GetPointsVisitor.h"

#include "SFCGAL/Exception.h"
#include "SFCGAL/Kernel.h"
#include <boost/format.hpp>

#include <CGAL/algorithm.h>
#include <CGAL/partition_2.h>

#include <iterator>
#include <list>
#include <vector>

namespace SFCGAL::algorithm {

using Traits     = CGAL::Partition_traits_2<Kernel>;
using TPoint_2   = Traits::Point_2;
using TPolygon_2 = Traits::Polygon_2;

static auto
toTPolygon_2(const Polygon &poly) -> CGAL::Partition_traits_2<Kernel>::Polygon_2
{
  if (poly.isEmpty()) {
    return {};
  }

  LineString::const_iterator pend = poly.exteriorRing().end();
  // skip the last point
  pend--;

  std::list<TPoint_2> points;
  Point               lastP;

  for (LineString::const_iterator pit = poly.exteriorRing().begin();
       pit != pend; ++pit) {
    if (pit == poly.exteriorRing().begin()) {
      lastP = *pit;
      TPoint_2 const point2(pit->coordinate().x(), pit->coordinate().y());
      points.push_back(point2);
    }

    if (lastP != *pit) {
      TPoint_2 const point2(pit->coordinate().x(), pit->coordinate().y());
      points.push_back(point2);
    }

    lastP = *pit;
  }

  TPolygon_2 result(points.begin(), points.end());

  return result;
}

static auto
polygons_to_geometry(const std::list<TPolygon_2> &polys)
    -> std::unique_ptr<Geometry>
{

  auto *geoms = new GeometryCollection;

  std::list<TPolygon_2>::const_iterator poly_it;
  for (poly_it = polys.begin(); poly_it != polys.end(); ++poly_it) {
    auto *poly = new Polygon;
    for (const TPoint_2 &p : poly_it->container()) {
      poly->exteriorRing().addPoint(p);
    }
    poly->exteriorRing().addPoint(*(poly_it->vertices_begin()));
    geoms->addGeometry(poly);
  }
  return std::unique_ptr<Geometry>(geoms);
}

auto
partition_2(const Geometry &g, PartitionAlgorithm alg)
    -> std::unique_ptr<Geometry>
{
  using CGAL::object_cast;

  if (g.isEmpty()) {
    return std::unique_ptr<Geometry>(new GeometryCollection());
  }

  if (g.geometryTypeId() != TYPE_POLYGON) {
    return std::unique_ptr<Geometry>(new GeometryCollection());
  }

  const TPolygon_2 poly{toTPolygon_2(g.as<Polygon>())};

  std::list<TPolygon_2> partition_polys;
  switch (alg) {
  case y_monotone:
    y_monotone_partition_2(poly.vertices_begin(), poly.vertices_end(),
                           std::back_inserter(partition_polys));
    break;
  case optimal_convex:
    optimal_convex_partition_2(poly.vertices_begin(), poly.vertices_end(),
                               std::back_inserter(partition_polys));
    break;
  case greene_approx_convex:
    greene_approx_convex_partition_2(poly.vertices_begin(), poly.vertices_end(),
                                     std::back_inserter(partition_polys));
    break;
  case approx_convex:
    approx_convex_partition_2(poly.vertices_begin(), poly.vertices_end(),
                              std::back_inserter(partition_polys));
    break;
  default:
    break;
  }
  return polygons_to_geometry(partition_polys);
}

} // namespace SFCGAL::algorithm
