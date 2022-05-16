// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/algorithm/polygonPartitioning.h>
#include <SFCGAL/algorithm/isValid.h>

#include <SFCGAL/detail/GetPointsVisitor.h>

#include <SFCGAL/Exception.h>
#include <SFCGAL/Kernel.h>
#include <SFCGAL/MultiPolygon.h>
#include <SFCGAL/Polygon.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>

namespace SFCGAL {
namespace algorithm {

using Traits       = CGAL::Partition_traits_2<Kernel>;
using Polygon_2    = Traits::Polygon_2;
using Point_2      = Traits::Point_2;
using Polygon_list = std::list<Polygon_2>;

static auto
detail_greene_approx_convex_partition(const Geometry &g) -> std::unique_ptr<Geometry>
{
  using CGAL::object_cast;

  if (g.isEmpty()) {
    return std::unique_ptr<Geometry>(g.clone());
  }

  SFCGAL::detail::GetPointsVisitor getPointVisitor;
  const_cast<Geometry &>(g).accept(getPointVisitor);

  // collect points

  if (getPointVisitor.points.size() == 0) {
    return std::unique_ptr<Geometry>(new GeometryCollection());
  }

  Polygon_2    polygon;
  Polygon_list partition_polys;

  for (auto &point : getPointVisitor.points) {
    polygon.push_back(point->toPoint_2());
  }

  // remove last point to avoid an auto-intersection:
  // pre-condition is_simple
  polygon.container().pop_back();
  // reverse orientation since the algorithm requires the inverse of the OGC one:
  // pre-condition: orientation
  polygon.reverse_orientation();

  CGAL::greene_approx_convex_partition_2(polygon.vertices_begin(),
                                         polygon.vertices_end(),
                                         std::back_inserter(partition_polys));

  auto *partitions = new MultiPolygon;
  for (auto &partition_poly : partition_polys) {
    auto *poly = new Polygon;

    for (std::list<Point_2>::const_iterator it = partition_poly.begin();
         it != partition_poly.end(); ++it) {
      poly->exteriorRing().addPoint(*it);
    }

    // add back the first point to close the ring
    poly->exteriorRing().addPoint(*partition_poly.begin());
    partitions->addGeometry(poly);
  }

  return std::unique_ptr<Geometry>(partitions);
}

auto
greene_approx_convex_partition(const Geometry &g, NoValidityCheck) -> std::unique_ptr<Geometry>
{
  switch (g.geometryTypeId()) {
  case TYPE_POINT:
  case TYPE_LINESTRING:
  case TYPE_TRIANGLE:
  case TYPE_MULTIPOINT:
  case TYPE_MULTILINESTRING:
  case TYPE_POLYHEDRALSURFACE: 
  case TYPE_SOLID: 
  case TYPE_MULTISOLID:
  case TYPE_GEOMETRYCOLLECTION: 
  case TYPE_TRIANGULATEDSURFACE: 
    return std::unique_ptr<Geometry>(g.clone());

  case TYPE_POLYGON:
  {
    return detail_greene_approx_convex_partition(g);
  }

  case TYPE_MULTIPOLYGON:{
    std::unique_ptr<GeometryCollection> result(new GeometryCollection);

    for (size_t i = 0; i < g.numGeometries(); ++i) {
      result->addGeometry(detail_greene_approx_convex_partition(g.geometryN(i)).release());
    }

    return std::unique_ptr<Geometry>(result.release());
  }
  default:
    break;
  }
  return std::unique_ptr<Geometry>(g.clone());
}

auto
greene_approx_convex_partition(const Geometry &g) -> std::unique_ptr<Geometry>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY(g);

  return greene_approx_convex_partition(g, NoValidityCheck());
}

} // namespace algorithm
} // namespace SFCGAL
